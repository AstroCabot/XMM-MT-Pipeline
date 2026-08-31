import json
import re
import subprocess
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.time import Time, TimeDelta
from astropy.wcs import WCS
from scipy.ndimage import map_coordinates
from scipy.signal import fftconvolve

from solar_wind import read_tsv, write_tsv

HERE = Path(__file__).resolve().parent
REPRO = HERE.parent
REDUCTION = REPRO / "analysis_01_reduction"
OUTPUT = HERE / "outputs"
WORK = HERE / "work"
PARAMS = json.loads((REPRO / "parameters.json").read_text())
SETTINGS = json.loads((REDUCTION / "settings.json").read_text())

XMM_MJDREF = 50814.0
BAND_EV = (200, 1000)
PATTERN = 0
MAP_STEP_S = 5.0
PIXEL_BIN = 80
EXPECTED = tuple(f"pn-s003-p{i:02d}" for i in (1, 2, 3, 4, 5, 8, 10, 11))


def valid_fits(path):
    if not path.exists() or not path.stat().st_size:
        return False
    try:
        with fits.open(path):
            return True
    except OSError:
        return False


def frame_fits(frame):
    root = WORK / frame
    return (
        WORK / "exposure" / f"{frame}.fits",
        root / "image.fits",
        root / "exposure.fits",
    )


def remove_frame_fits(frame):
    paths = frame_fits(frame)
    root = WORK.resolve()
    if any(root not in path.resolve().parents for path in paths):
        raise RuntimeError(f"refusing to remove a file outside {WORK}")
    for path in paths:
        path.unlink(missing_ok=True)


def pn_frames():
    frame_rows = [
        row
        for row in read_tsv(REDUCTION / "outputs/03_events/frames.tsv")
        if row["det"] == "pn"
    ]
    masked_rows = [
        row
        for row in read_tsv(REDUCTION / "outputs/04_sources/events.tsv")
        if row["det"] == "pn"
    ]
    frame_names = [row["frame"] for row in frame_rows]
    masked_names = [row["frame"] for row in masked_rows]
    if (
        len(frame_names) != len(EXPECTED)
        or set(frame_names) != set(EXPECTED)
        or len(masked_names) != len(EXPECTED)
        or set(masked_names) != set(EXPECTED)
    ):
        raise RuntimeError("Task 01 PN frame set changed")
    frames = dict(zip(frame_names, frame_rows))
    masked = dict(zip(masked_names, masked_rows))
    return [{**frames[name], **masked[name]} for name in EXPECTED]


def event_gti(path):
    with fits.open(path) as hdul:
        ontime = float(hdul["EVENTS"].header["ONTIME"])
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
    difference, spans = min(candidates, key=lambda value: value[0])
    if (
        not len(spans)
        or np.any(spans[:, 1] <= spans[:, 0])
        or np.any(spans[1:, 0] < spans[:-1, 1])
        or difference > 0.01
    ):
        raise RuntimeError(f"invalid event GTI: {path}")
    return spans, ontime


def gti_samples(row):
    """Sample the filtered event GTI at MAP_STEP_S with duration weights."""
    spans, ontime = event_gti(Path(row["event"]))
    if abs(ontime - float(row["ontime"])) > 0.051:
        raise RuntimeError(f'ONTIME changed for {row["frame"]}')
    times, weights = [], []
    for start, stop in spans:
        count = max(1, int(np.ceil((stop - start) / MAP_STEP_S)))
        edges = np.linspace(start, stop, count + 1)
        times.extend((edges[:-1] + edges[1:]) / 2)
        weights.extend(np.diff(edges))
    return np.asarray(times), np.asarray(weights)


def event_wcs(header):
    wcs = WCS(naxis=2)
    wcs.wcs.crval = [header["REFXCRVL"], header["REFYCRVL"]]
    wcs.wcs.crpix = [header["REFXCRPX"], header["REFYCRPX"]]
    wcs.wcs.cdelt = [header["REFXCDLT"], header["REFYCDLT"]]
    wcs.wcs.ctype = [header["REFXCTYP"], header["REFYCTYP"]]
    return wcs


def sky_wcs(header):
    """Normalize RADECSYS before constructing the celestial WCS."""
    if "RADECSYS" in header:
        header["RADESYS"] = header["RADECSYS"]
        del header["RADECSYS"]
    return WCS(header, fix=False).celestial


def event_mjd(seconds, header):
    ref = float(
        header.get(
            "MJDREF", header.get("MJDREFI", XMM_MJDREF) + header.get("MJDREFF", 0)
        )
    )
    return (
        Time(ref, format="mjd", scale=str(header.get("TIMESYS", "TT")).lower())
        + TimeDelta(seconds, format="sec")
    ).utc.mjd


def ephemeris():
    path = REDUCTION / "outputs/07_comove/ephemeris.tsv"
    if not path.exists():
        raise FileNotFoundError("Task 01 comoving ephemeris is not complete")
    rows = read_tsv(path)
    value = {
        key: np.array([float(row[key]) for row in rows])
        for key in ("mjd", "ra_deg", "dec_deg")
    }
    if (
        len(rows) < 2
        or np.any(~np.isfinite(np.column_stack(list(value.values()))))
        or np.any(np.diff(value["mjd"]) <= 0)
    ):
        raise RuntimeError("Task 01 ephemeris is invalid")
    return value


def comet_position(times, eph):
    times = np.asarray(times, float)
    if len(times) and (times.min() < eph["mjd"][0] or times.max() > eph["mjd"][-1]):
        raise RuntimeError("Task 01 ephemeris does not cover the selected events")
    ra = (
        np.degrees(np.interp(times, eph["mjd"], np.unwrap(np.radians(eph["ra_deg"]))))
        % 360
    )
    return ra, np.interp(times, eph["mjd"], eph["dec_deg"])


def sources():
    rows = read_tsv(REDUCTION / "outputs/04_sources/sources.tsv")
    values = [
        (float(row["ra"]), float(row["dec"]), float(row["radius_arcsec"]))
        for row in rows
    ]
    if (
        len(values) != 15
        or len({item[:2] for item in values}) != 15
        or not np.isfinite(values).all()
        or any(item[2] != PARAMS["source_mask_radius_arcsec"] for item in values)
    ):
        raise RuntimeError("Task 01 source catalog is invalid")
    return values


def separation(ra1, dec1, ra2, dec2):
    dra = np.radians(ra1 - ra2)
    d1, d2 = np.radians(dec1), np.radians(dec2)
    value = np.sin((d1 - d2) / 2) ** 2 + np.cos(d1) * np.cos(d2) * np.sin(dra / 2) ** 2
    return np.degrees(2 * np.arcsin(np.sqrt(np.clip(value, 0, 1)))) * 3600


def counts(path, eph):
    with fits.open(path, memmap=True) as hdul:
        data, header = hdul["EVENTS"].data, hdul["EVENTS"].header
        use = (
            (data["PI"] >= BAND_EV[0])
            & (data["PI"] <= BAND_EV[1])
            & (data["PATTERN"] == PATTERN)
            & (data["FLAG"] == 0)
        )
        ra, dec = event_wcs(header).wcs_pix2world(data["X"][use], data["Y"][use], 1)
        radius = separation(
            ra, dec, *comet_position(event_mjd(data["TIME"][use], header), eph)
        )
    inner, outer = PARAMS["lightcurve_background_annulus_arcsec"]
    return (
        int(np.sum(radius < PARAMS["science_aperture_radius_arcsec"])),
        int(np.sum((radius >= inner) & (radius < outer))),
    )


def sas_environment():
    command = (
        f'export HEADAS="{SETTINGS["headas"]}" && source "$HEADAS/headas-init.sh" && '
        f'source "{SETTINGS["sas"]}" && '
        f'export LD_LIBRARY_PATH="{SETTINGS["python_lib"]}:${{LD_LIBRARY_PATH:-}}" && '
        "env -0"
    )
    result = subprocess.run(["bash", "-lc", command], capture_output=True)
    if result.returncode:
        raise RuntimeError(result.stderr.decode(errors="replace"))
    env = dict(
        item.split("=", 1) for item in result.stdout.decode().split("\0") if "=" in item
    )
    env.update(
        SAS_CCFPATH=SETTINGS["ccf"], SAS_CCF=str(REDUCTION / "outputs/01_init/ccf.cif")
    )
    summary = REDUCTION / "outputs/01_init/odf" / Path(SETTINGS["summary"]).name
    if not summary.is_file():
        raise RuntimeError("Task 01 ODF summary is missing")
    env["SAS_ODF"] = str(summary)
    env.update(
        TMPDIR=str(WORK / "tmp"),
        XDG_CACHE_HOME=str(WORK / "cache"),
        XDG_CONFIG_HOME=str(WORK / "config"),
    )
    for path in (WORK / "tmp", WORK / "cache", WORK / "config", WORK / "pfiles"):
        path.mkdir(parents=True, exist_ok=True)
    env["PFILES"] = f'{WORK / "pfiles"};{env.get("PFILES", "").split(";")[-1]}'
    return env


def sas(tool, cwd, env, **parameters):
    argv = [tool] + [
        f"{key}={'yes' if value is True else 'no' if value is False else value}"
        for key, value in parameters.items()
    ]
    result = subprocess.run(argv, cwd=cwd, env=env, text=True, capture_output=True)
    text = result.stdout + result.stderr
    if result.returncode or re.search(
        r"^\*\* .*:\s+(?:fatal\s+)?error\b", text, flags=re.I | re.M
    ):
        raise RuntimeError(
            f"{tool} failed\n{result.stdout[-3000:]}\n{result.stderr[-3000:]}"
        )


def exposure_map(row, env):
    """Unvignetted exposure seconds including QE, filter, gaps, bad pixels, and source holes."""
    destination, image, temporary = frame_fits(row["frame"])
    remove_frame_fits(row["frame"])
    root = image.parent
    root.mkdir(parents=True, exist_ok=True)
    event = Path(row["event"])
    sas(
        "evselect",
        root,
        env,
        table=f"{event}:EVENTS",
        withimageset=True,
        imageset=image,
        imagebinning="binSize",
        xcolumn="X",
        ycolumn="Y",
        ximagebinsize=PIXEL_BIN,
        yimagebinsize=PIXEL_BIN,
        squarepixels=True,
        expression=f"(PI>={BAND_EV[0]})&&(PI<={BAND_EV[1]})"
        f"&&(PATTERN=={PATTERN})&&(FLAG==0)",
    )
    sas(
        "eexpmap",
        root,
        env,
        imageset=image,
        eventset=f"{event}:EVENTS",
        attitudeset=REDUCTION / "outputs/02_reprocess/atthk.fits",
        expimageset=temporary,
        withdetcoords=False,
        withvignetting=False,
        usefastpixelization=False,
        attrebin=4,
        pimin=BAND_EV[0],
        pimax=BAND_EV[1],
    )
    with fits.open(temporary) as hdul:
        hdu = next(h for h in hdul if h.data is not None and h.data.ndim == 2)
        data, header = np.asarray(hdu.data, float), hdu.header.copy()
    data[~np.isfinite(data) | (data < 0)] = 0
    y, x = np.indices(data.shape)
    ra, dec = sky_wcs(header).pixel_to_world_values(x, y)
    for source_ra, source_dec, radius_arcsec in sources():
        data[separation(ra, dec, source_ra, source_dec) < radius_arcsec] = 0
    if not np.any(data > 0):
        raise RuntimeError(f'empty exposure map for {row["frame"]}')
    header.update(
        VIGNET=False,
        SRCMASK=True,
        PI_MIN=BAND_EV[0],
        PI_MAX=BAND_EV[1],
        PATTERN=PATTERN,
    )
    partial = destination.with_suffix(".partial")
    destination.parent.mkdir(parents=True, exist_ok=True)
    partial.unlink(missing_ok=True)
    fits.PrimaryHDU(data.astype("f4"), header).writeto(partial)
    if not valid_fits(partial):
        raise RuntimeError(f"invalid exposure map: {partial}")
    partial.replace(destination)
    return destination


def disk(radius):
    half = int(np.ceil(radius))
    y, x = np.ogrid[-half : half + 1, -half : half + 1]
    return ((x * x + y * y) <= radius * radius).astype(float)


def annulus(inner, outer):
    kernel = disk(outer)
    edge = (kernel.shape[0] - disk(inner).shape[0]) // 2
    kernel[edge:-edge, edge:-edge] -= disk(inner)
    return kernel


def region_exposure(path, times, weights, eph):
    with fits.open(path) as hdul:
        hdu = next(h for h in hdul if h.data is not None and h.data.ndim == 2)
        data, wcs = np.asarray(hdu.data, float), sky_wcs(hdu.header.copy())
    scale = np.sqrt(abs(np.linalg.det(wcs.pixel_scale_matrix))) * 3600
    pixel_area = (scale / 60) ** 2
    inner, outer = PARAMS["lightcurve_background_annulus_arcsec"]
    kernels = (
        disk(PARAMS["science_aperture_radius_arcsec"] / scale),
        annulus(inner / scale, outer / scale),
    )
    x, y = wcs.world_to_pixel_values(*comet_position(times, eph))
    values = [
        np.average(
            map_coordinates(
                fftconvolve(data, kernel, mode="same") * pixel_area,
                [y, x],
                order=1,
                mode="constant",
            ),
            weights=weights,
        )
        for kernel in kernels
    ]
    if min(values) <= 0:
        raise RuntimeError(f"non-positive region exposure for {path}")
    return values


def rate(source, background, source_oot, background_oot, source_exp, background_exp):
    factor = PARAMS["pn_oot_fraction"]
    area = np.pi * (PARAMS["science_aperture_radius_arcsec"] / 60) ** 2
    net = (source - factor * source_oot, background - factor * background_oot)
    variance = (
        source + factor**2 * source_oot,
        background + factor**2 * background_oot,
    )
    value = area * (net[0] / source_exp - net[1] / background_exp)
    error = area * np.sqrt(
        variance[0] / source_exp**2 + variance[1] / background_exp**2
    )
    return value, error


def build(frames):
    eph, env, rows = ephemeris(), sas_environment(), []
    for frame in frames:
        try:
            times, weights = gti_samples(frame)
            science = counts(Path(frame["event"]), eph)
            oot = counts(Path(frame["oot"]), eph)
            times_mjd = event_mjd(times, fits.getheader(frame["event"], "EVENTS"))
            exposure = region_exposure(
                exposure_map(frame, env), times_mjd, weights, eph
            )
            value, error = rate(*science, *oot, *exposure)
            row = {
                "frame": frame["frame"],
                "midpoint_utc": Time(
                    np.average(times_mjd, weights=weights), format="mjd"
                ).isot,
                "gti_s": f"{weights.sum():.3f}",
                "source_counts": science[0],
                "background_counts": science[1],
                "source_oot": oot[0],
                "background_oot": oot[1],
                "source_exposure_s_arcmin2": f"{exposure[0]:.6f}",
                "background_exposure_s_arcmin2": f"{exposure[1]:.6f}",
                "rate_ct_s": f"{value:.8f}",
                "error_ct_s": f"{error:.8f}",
            }
        finally:
            remove_frame_fits(frame["frame"])
        rows.append(row)
    write_tsv(OUTPUT / "lightcurve.tsv", tuple(rows[0]), rows)
    return rows
