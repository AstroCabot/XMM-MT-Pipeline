import csv
import json
import os
import shutil
import subprocess
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.time import Time, TimeDelta

import psf_operator

INSTRUMENT = {"pn": "EPN", "mos1": "EMOS1", "mos2": "EMOS2"}
STACK_GRID_KEYS = (
    "CTYPE1",
    "CTYPE2",
    "CRPIX1",
    "CRPIX2",
    "CRVAL1",
    "CRVAL2",
    "CDELT1",
    "CDELT2",
)


def read_tsv(path):
    with Path(path).open(newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def event_mjd(seconds, header):
    ref = float(
        header["MJDREF"]
        if "MJDREF" in header
        else header["MJDREFI"] + header.get("MJDREFF", 0)
    )
    scale = str(header.get("TIMESYS", "TT")).strip().lower()
    return (
        Time(ref, format="mjd", scale=scale) + TimeDelta(seconds, format="sec")
    ).utc.mjd


def event_gti(event):
    with fits.open(event, memmap=True) as hdul:
        header = hdul["EVENTS"].header.copy()
        ontime = float(header["ONTIME"])
        candidates = []
        for hdu in hdul:
            if (
                hdu.data is None
                or not getattr(hdu.data, "names", None)
                or not {"START", "STOP"} <= set(hdu.data.names)
            ):
                continue
            spans = np.column_stack((hdu.data["START"], hdu.data["STOP"])).astype(float)
            duration = float(np.sum(spans[:, 1] - spans[:, 0]))
            candidates.append((abs(duration - ontime), spans))
    if not candidates:
        raise RuntimeError(f"no event GTI: {event}")
    difference, spans = min(candidates, key=lambda value: value[0])
    if (
        not len(spans)
        or np.any(spans[:, 1] <= spans[:, 0])
        or np.any(spans[1:, 0] < spans[:-1, 1])
        or difference > 0.01
    ):
        raise RuntimeError(f"invalid event GTI: {event}")
    return spans, header


def gti_samples(frame, event, step_s):
    spans, header = event_gti(event)
    if step_s <= 0 or (
        frame.get("ontime") is not None
        and abs(float(header["ONTIME"]) - float(frame["ontime"])) > 0.051
    ):
        raise RuntimeError(f"invalid GTI sampling for {frame['frame']}")
    times, weights = [], []
    for start, stop in spans:
        edges = np.linspace(
            start, stop, max(1, int(np.ceil((stop - start) / step_s))) + 1
        )
        times.extend((edges[:-1] + edges[1:]) / 2)
        weights.extend(np.diff(edges))
    if not weights or max(weights) > step_s + 1e-9:
        raise RuntimeError(f"invalid GTI sampling for {frame['frame']}")
    return np.asarray(times), np.asarray(weights), header


def comet_track(frame, event, ephemeris, step_s=5.0):
    times, weights, header = gti_samples(frame, event, step_s)
    mjd_sample = event_mjd(times, header)
    mjd = ephemeris[:, 0]
    if mjd_sample.min() < mjd[0] or mjd_sample.max() > mjd[-1]:
        raise RuntimeError("Task 01 ephemeris does not cover a retained frame")
    ra = (
        np.rad2deg(np.interp(mjd_sample, mjd, np.unwrap(np.deg2rad(ephemeris[:, 1]))))
        % 360
    )
    dec = np.interp(mjd_sample, mjd, ephemeris[:, 2])
    return ra, dec, weights


def mean_position(ra, dec, weights):
    mean_ra = np.rad2deg(np.average(np.unwrap(np.deg2rad(ra)), weights=weights)) % 360
    return float(mean_ra), float(np.average(dec, weights=weights))


def run_psfgen(output, image, detector, ra, dec, energy_ev, work, sas_env):
    output = Path(output)
    output.unlink(missing_ok=True)
    pfiles = Path(work) / ".pfiles"
    pfiles.mkdir(parents=True, exist_ok=True)
    env = os.environ.copy()
    env.update(sas_env)
    system = env.get("PFILES", "").split(";")[-1]
    env["PFILES"] = f"{pfiles};{system}"
    argv = [
        "psfgen",
        f"instrument={INSTRUMENT[detector]}",
        f"energy={energy_ev:.8g}",
        "level=ELLBETA",
        "coordtype=EQPOS",
        f"x={ra:.9f}",
        f"y={dec:.9f}",
        f"image={image}",
        "xsize=1024",
        "ysize=1024",
        f"output={output}",
    ]
    result = subprocess.run(argv, cwd=work, env=env, text=True, capture_output=True)
    if result.returncode:
        raise RuntimeError(result.stderr or result.stdout)
    return output


def exposure_maps(reduction, frames, soft_band_ev):
    reduction = Path(reduction)
    low, high = map(int, soft_band_ev)
    rows = read_tsv(reduction / "05_background/background.tsv")
    keys = [(row["frame"], row["band"]) for row in rows]
    expected = {(frame, band) for frame in frames for band in ("soft", "hard")}
    if len(keys) != len(expected) or set(keys) != expected:
        raise RuntimeError("Task 01 background table is incomplete or duplicated")
    index = dict(zip(keys, rows))
    maps = {}
    for frame in frames:
        detector, exposure, _ = frame.split("-", 2)
        row = index[frame, "soft"]
        path = (
            reduction
            / "05_background"
            / frame
            / "soft"
            / (f"{detector}{exposure.upper()}-expimsky-{low}-{high}.fits")
        )
        if (
            row["det"] != frame_detector(frame)
            or Path(row["exposure_map"]).resolve() != path.resolve()
        ):
            raise RuntimeError(f"Task 01 exposure-map mismatch for {frame}")
        if not path.is_file() or not path.stat().st_size:
            raise RuntimeError(f"missing Task 01 exposure map for {frame}")
        data = np.asarray(fits.getdata(path), float)
        if (
            data.ndim != 2
            or not np.isfinite(data).all()
            or np.any(data < 0)
            or not np.any(data > 0)
        ):
            raise RuntimeError(f"invalid Task 01 exposure map for {frame}")
        maps[frame] = path
    return maps


def frame_detector(frame):
    detector = frame.split("-", 1)[0]
    if detector not in INSTRUMENT:
        raise RuntimeError(f"invalid EPIC frame name: {frame}")
    return detector


def event_path(reduction, frame, row):
    path = Path(reduction) / "04_sources/events" / f"{frame}.fits"
    if Path(row["event"]).resolve() != path.resolve() or not path.is_file():
        raise RuntimeError(f"Task 01 event-path mismatch for {frame}")
    return path


def detector_exposure_stacks(reduction, soft_band_ev):
    stacks, grid = {}, None
    for detector in INSTRUMENT:
        path = Path(reduction) / "07_comove/stack" / (f"{detector}-exposure-soft.fits")
        data, header = fits.getdata(path, header=True)
        data = np.asarray(data, float)
        current = (data.shape, tuple(header.get(key) for key in STACK_GRID_KEYS))
        if (
            header.get("DETECTOR") != detector
            or header.get("BAND") != "soft"
            or [header.get("E_MIN"), header.get("E_MAX")] != list(soft_band_ev)
            or not header.get("COMOVING", False)
            or data.ndim != 2
            or np.any(~np.isfinite(data))
            or not np.any(data > 0)
        ):
            raise RuntimeError(f"invalid Task 01 exposure stack for {detector}")
        if grid is None:
            grid = current
        elif current != grid:
            raise RuntimeError("Task 01 exposure stacks use different comet grids")
        stacks[detector] = (path, data, header)
    return stacks


def frame_exposure_weights(
    reduction,
    frames,
    events,
    maps,
    ephemeris,
    soft_band_ev,
    annulus_edges_arcsec,
    step_s,
):
    stacks = detector_exposure_stacks(reduction, soft_band_ev)
    rows, centers = [], []
    rebuilt = {detector: np.zeros_like(item[1]) for detector, item in stacks.items()}
    for frame in frames:
        name, detector = frame["frame"], frame_detector(frame["frame"])
        event = event_path(reduction, name, events[name])
        ra, dec, sample_weights = comet_track(frame, event, ephemeris, step_s)
        centers.append(mean_position(ra, dec, sample_weights))
        exposure, header = fits.getdata(maps[name], header=True)
        _, _, output_header = stacks[detector]
        plane = psf_operator.comoving_exposure(
            exposure,
            header,
            ra,
            dec,
            sample_weights,
            output_header,
            rebuilt[detector].shape,
        )
        rebuilt[detector] += plane
        rows.append(
            psf_operator.annulus_exposure(plane, output_header, annulus_edges_arcsec)
        )
    closure = []
    for detector, (_, expected_plane, header) in stacks.items():
        expected = psf_operator.annulus_exposure(
            expected_plane, header, annulus_edges_arcsec
        )
        actual = psf_operator.annulus_exposure(
            rebuilt[detector], header, annulus_edges_arcsec
        )
        closure.append(
            psf_operator.fractional_closure(actual, expected, f"Task 01 {detector}")
        )
    return np.asarray(rows), np.asarray(centers), np.asarray(closure)


def radial_density(image, header):
    image = np.asarray(image, float).copy()
    image[~np.isfinite(image) | (image < 0)] = 0
    total = image.sum()
    if total <= 0:
        raise RuntimeError("empty ELLBETA image")
    image /= total
    rows, columns = np.indices(image.shape)
    center = np.asarray(image.shape, float) / 2
    scale = np.abs([header.get("CDELT2", 0), header.get("CDELT1", 0)]) * 3600
    crpix = [header.get("CRPIX2"), header.get("CRPIX1")]
    units = [str(header.get(key, "")).strip().lower() for key in ("CUNIT2", "CUNIT1")]
    if (
        units != ["deg", "deg"]
        or not np.allclose(scale, 1, rtol=0, atol=1e-9)
        or not np.allclose(crpix, center, rtol=0, atol=1e-9)
    ):
        raise RuntimeError("unexpected ELLBETA grid")
    # SAS places the ELLBETA peak at zero-based N/2.
    peaks = np.argwhere(image == image.max())
    if len(peaks) != 1 or not np.array_equal(peaks[0], center):
        raise RuntimeError("off-center ELLBETA image")
    radius = np.hypot(columns - center[1], rows - center[0])
    index = np.floor(radius).astype(int)
    count = min(image.shape) // 2
    probability = np.bincount(index.ravel(), image.ravel(), minlength=count)[:count]
    edges = np.arange(count + 1, dtype=float)
    density = probability / (np.pi * (edges[1:] ** 2 - edges[:-1] ** 2))
    return 0.5 * (edges[:-1] + edges[1:]), density


def generate(reduction, spectrum, work, soft_band_ev, annulus_edges_arcsec):
    reduction = Path(reduction)
    frames = read_tsv(reduction / "03_events/frames.tsv")
    event_rows = read_tsv(reduction / "04_sources/events.tsv")
    events = {row["frame"]: row for row in event_rows}
    settings = json.loads((reduction.parent / "settings.json").read_text())
    sas_env = {
        "SAS_CCF": str(reduction / "01_init/ccf.cif"),
        "SAS_CCFPATH": settings["ccf"],
        "SAS_ODF": str(reduction / "01_init/odf" / Path(settings["summary"]).name),
    }
    frame_names = [row["frame"] for row in frames]
    expected_frames = tuple(settings["expected_frames"])
    maps = exposure_maps(reduction, set(expected_frames), soft_band_ev)
    expected = tuple(spectrum["detector_K_ct_s_per_erg_cm2_s"])
    if (
        tuple(frame_names) != expected_frames
        or len(event_rows) != len(expected_frames)
        or set(events) != set(expected_frames)
        or set(expected) != set(INSTRUMENT)
    ):
        raise RuntimeError(
            f"PSF construction requires the {len(expected_frames)} retained frames"
        )
    eph_rows = read_tsv(reduction / "07_comove/ephemeris.tsv")
    ephemeris = np.array(
        [
            [float(row[name]) for name in ("mjd", "ra_deg", "dec_deg")]
            for row in eph_rows
        ]
    )
    if (
        len(eph_rows) < 2
        or np.any(~np.isfinite(ephemeris))
        or np.any(np.diff(ephemeris[:, 0]) <= 0)
    ):
        raise RuntimeError("Task 01 ephemeris is invalid")
    for frame in frames:
        detector = frame_detector(frame["frame"])
        if frame["det"] != detector or events[frame["frame"]]["det"] != detector:
            raise RuntimeError(f"Task 01 detector mismatch for {frame['frame']}")
    raw_weights, centers, closure = frame_exposure_weights(
        reduction,
        frames,
        events,
        maps,
        ephemeris,
        soft_band_ev,
        annulus_edges_arcsec,
        float(settings["comove"]["map_step_s"]),
    )
    factors = spectrum["detector_K_ct_s_per_erg_cm2_s"]
    scales = np.array([factors[frame["det"]] / factors["mos2"] for frame in frames])
    weights = raw_weights * scales[:, None]
    work = Path(work)
    work.mkdir(parents=True, exist_ok=True)
    densities, names, psf_radius = [], [], None
    for index, frame in enumerate(frames):
        name = frame["frame"]
        detector = frame_detector(name)
        exposure_path = maps[name]
        ra, dec = centers[index]
        path = run_psfgen(
            work / f"{name}.fits",
            exposure_path,
            detector,
            ra,
            dec,
            spectrum["mean_photon_energy_eV"],
            work,
            sas_env,
        )
        image, header = fits.getdata(path, header=True)
        radius, density = radial_density(image, header)
        path.unlink()
        if psf_radius is None:
            psf_radius = radius
        elif not np.array_equal(radius, psf_radius):
            raise RuntimeError("ELLBETA images have different radial grids")
        densities.append(density)
        names.append(name)
    shutil.rmtree(work / ".pfiles")
    if np.any(~np.isfinite(weights)) or np.any(weights.sum(axis=0) <= 0):
        raise RuntimeError("empty per-annulus PSF exposure weights")
    return (psf_radius, np.asarray(densities), weights, np.asarray(names), closure)
