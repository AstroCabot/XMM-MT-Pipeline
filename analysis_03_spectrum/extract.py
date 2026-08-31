import json
import os
import re
import shutil
import subprocess
import time
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.time import Time, TimeDelta
from astropy.wcs import WCS
from astropy.wcs.utils import proj_plane_pixel_area

import spectrum

PHYSICAL_PIXEL_ARCSEC = 0.05
ROOT = Path(__file__).resolve().parent.parent
WORK = ROOT / "work" / "analysis_03_spectrum"
REDUCTION = ROOT / "analysis_01_reduction" / "outputs"
SETTINGS = json.loads((ROOT / "analysis_01_reduction/settings.json").read_text())
SAS_ENV = {
    "SAS_CCF": str(REDUCTION / "01_init/ccf.cif"),
    "SAS_CCFPATH": SETTINGS["ccf"],
    "SAS_ODF": str(REDUCTION / "01_init/odf" / Path(SETTINGS["summary"]).name),
}


def run_tool(tool, cwd, **parameters):
    argv = [tool] + [f"{key}={value}" for key, value in parameters.items()]
    cwd = Path(cwd)
    log_dir = WORK / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)
    stamp = f"{tool}-{time.time_ns()}"
    log = log_dir / f"{stamp}.log"
    runtime = WORK / "runtime" / stamp
    local = runtime / "pfiles"
    local.mkdir(parents=True)
    env = os.environ.copy()
    env.update(SAS_ENV)
    env["PFILES"] = f"{local};{env.get('PFILES', '').split(';')[-1]}"
    scratch = runtime / "tmp"
    scratch.mkdir()
    env["TMPDIR"] = str(scratch)
    try:
        with log.open("wb") as stream:
            result = subprocess.run(
                argv, cwd=cwd, env=env, stdout=stream, stderr=subprocess.STDOUT
            )
    finally:
        shutil.rmtree(runtime)
    text = log.read_text(errors="replace")
    fatal = any(
        line.startswith("** ") and re.search(r":\s+(?:fatal\s+)?error\b", line, re.I)
        for line in text.splitlines()
    )
    if result.returncode or fatal:
        raise RuntimeError(f"{tool} failed; see {log}")


def event_wcs(path):
    header = fits.getheader(path, "EVENTS")
    wcs = WCS(naxis=2)
    wcs.wcs.crval = [header["REFXCRVL"], header["REFYCRVL"]]
    wcs.wcs.crpix = [header["REFXCRPX"], header["REFYCRPX"]]
    wcs.wcs.cdelt = [header["REFXCDLT"], header["REFYCDLT"]]
    wcs.wcs.ctype = [header["REFXCTYP"], header["REFYCTYP"]]
    return wcs


def event_gti(path):
    with fits.open(path, memmap=True) as hdul:
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
    if not candidates:
        raise RuntimeError(f"no event GTI: {path}")
    difference, spans = min(candidates, key=lambda value: value[0])
    if (
        not len(spans)
        or np.any(spans[:, 1] <= spans[:, 0])
        or np.any(spans[1:, 0] < spans[:-1, 1])
        or difference > 0.01
    ):
        raise RuntimeError(f"invalid event GTI: {path}")
    return spans


def center(frame, ephemeris, step_s=5.0):
    spans = event_gti(frame["event"])
    times, weights = [], []
    for start, stop in spans:
        edges = np.linspace(
            start, stop, max(1, int(np.ceil((stop - start) / step_s))) + 1
        )
        times.extend((edges[:-1] + edges[1:]) / 2)
        weights.extend(np.diff(edges))
    if not weights or max(weights) > step_s + 1e-9:
        raise RuntimeError(f"invalid GTI sampling for {frame['frame']}")
    header = fits.getheader(frame["event"], "EVENTS")
    ref = float(
        header.get(
            "MJDREF", header.get("MJDREFI", 50814.0) + header.get("MJDREFF", 0.0)
        )
    )
    scale = str(header.get("TIMESYS", "TT")).strip().lower()
    mjd = np.asarray(
        (Time(ref, format="mjd", scale=scale) + TimeDelta(times, format="sec")).utc.mjd
    )
    if mjd.min() < ephemeris["mjd"].min() or mjd.max() > ephemeris["mjd"].max():
        raise RuntimeError("ephemeris does not cover a retained GTI")
    ra = np.interp(mjd, ephemeris["mjd"], np.unwrap(np.deg2rad(ephemeris["ra_deg"])))
    values = {
        "ra_deg": np.rad2deg(np.average(ra, weights=weights)) % 360,
        "dec_deg": np.average(
            np.interp(mjd, ephemeris["mjd"], ephemeris["dec_deg"]), weights=weights
        ),
        "r_au": np.average(
            np.interp(mjd, ephemeris["mjd"], ephemeris["r_au"]), weights=weights
        ),
        "distance_au": np.average(
            np.interp(mjd, ephemeris["mjd"], ephemeris["delta_au"]), weights=weights
        ),
        "duration_s": sum(weights),
    }
    return {key: float(value) for key, value in values.items()}


def separation(ra1, dec1, ra2, dec2):
    return float(
        np.hypot((ra1 - ra2) * np.cos(np.deg2rad((dec1 + dec2) / 2)), dec1 - dec2)
        * 3600
    )


def circle_overlap(radius, hole_radius, distance):
    if distance >= radius + hole_radius:
        return 0.0
    if distance <= abs(radius - hole_radius):
        return np.pi * min(radius, hole_radius) ** 2
    a = np.arccos((distance**2 + radius**2 - hole_radius**2) / (2 * distance * radius))
    b = np.arccos(
        (distance**2 + hole_radius**2 - radius**2) / (2 * distance * hole_radius)
    )
    cross = (
        np.sqrt(
            (-distance + radius + hole_radius)
            * (distance + radius - hole_radius)
            * (distance - radius + hole_radius)
            * (distance + radius + hole_radius)
        )
        / 2
    )
    return radius**2 * a + hole_radius**2 * b - cross


def solid_angle(region, position, sources):
    inner, outer = region
    area = np.pi * (outer**2 - inner**2)
    for source in sources:
        distance = separation(
            position["ra_deg"], position["dec_deg"], source["ra"], source["dec"]
        )
        area -= circle_overlap(outer, source["radius"], distance)
        area += circle_overlap(inner, source["radius"], distance)
    return area / 3600


def selection(event, detector, region, position, sources, chips):
    wcs = event_wcs(event)
    x, y = wcs.wcs_world2pix(position["ra_deg"], position["dec_deg"], 1)
    inner, outer = (value / PHYSICAL_PIXEL_ARCSEC for value in region)
    aperture = (
        f"(X,Y) IN circle({x:.4f},{y:.4f},{outer:.4f})"
        if inner == 0
        else f"(X,Y) IN annulus({x:.4f},{y:.4f},{inner:.4f},{outer:.4f})"
    )
    holes = []
    for source in sources:
        distance = separation(
            position["ra_deg"], position["dec_deg"], source["ra"], source["dec"]
        )
        if (
            distance - source["radius"] >= region[1]
            or distance + source["radius"] <= region[0]
        ):
            continue
        sx, sy = wcs.wcs_world2pix(source["ra"], source["dec"], 1)
        holes.append(
            f"!((X,Y) IN circle({sx:.4f},{sy:.4f},"
            f"{source['radius'] / PHYSICAL_PIXEL_ARCSEC:.4f}))"
        )
    pattern = "PATTERN==0" if detector == "pn" else "PATTERN<=12"
    states = chips.split()
    expected = 4 if detector == "pn" else 7
    if (
        len(states) != expected
        or set(states) - {"T", "F"}
        or (detector == "pn" and "F" in states)
    ):
        raise RuntimeError(f"invalid chip selection for {detector}: {chips}")
    excluded = (
        ()
        if detector == "pn"
        else tuple(
            f"CCDNR!={ccd}" for ccd, state in enumerate(states, 1) if state == "F"
        )
    )
    return "&&".join(("FLAG==0", pattern, *excluded, aperture, *holes))


def map_region(path, region, position):
    with fits.open(path, memmap=True) as hdul:
        hdu = next(h for h in hdul if h.data is not None and h.data.ndim == 2)
        data, header = np.asarray(hdu.data, float), hdu.header.copy()
    if "RADECSYS" in header:
        if "RADESYS" not in header:
            header["RADESYS"] = header["RADECSYS"]
        del header["RADECSYS"]
    wcs = WCS(header, fix=False).celestial
    y, x = np.indices(data.shape)
    ra, dec = wcs.pixel_to_world_values(x, y)
    radius = (
        np.hypot(
            (ra - position["ra_deg"]) * np.cos(np.deg2rad(position["dec_deg"])),
            dec - position["dec_deg"],
        )
        * 3600
    )
    use = (radius >= region[0]) & (radius < region[1]) & np.isfinite(data)
    pixel_area = abs(float(proj_plane_pixel_area(wcs))) * 3600
    return data, use, pixel_area


def map_fraction(path, region, position):
    data, use, _ = map_region(path, region, position)
    total = float(np.nansum(data))
    return float(np.nansum(data[use]) / total) if total > 0 else 0.0


def mean_exposure(path, region, position, omega):
    data, use, pixel_area = map_region(path, region, position)
    return float(np.nansum(data[use]) * pixel_area / omega)


def link(path, target):
    if path.is_symlink() and path.resolve() == Path(target).resolve():
        return
    path.unlink(missing_ok=True)
    path.symlink_to(target)


def one(frame, region_name, region, sources, channel_max, oot_scale, fit_band, output):
    work = output / "frames" / f"{frame['frame']}-{region_name}"
    work.mkdir(parents=True, exist_ok=True)
    event, oot = work / "events.fits", work / "oot.fits"
    link(event, frame["event"])
    if frame["oot"]:
        link(oot, frame["oot"])
    expression = selection(
        event,
        frame["det"],
        region,
        frame["position"],
        sources,
        frame["chips"],
    )
    source, rmf, arf = work / "source.pi", work / "source.rmf", work / "source.arf"
    args = dict(
        table="events.fits:EVENTS",
        withspectrumset="yes",
        spectrumset=source.name,
        energycolumn="PI",
        spectralbinsize=5,
        specchannelmin=0,
        specchannelmax=channel_max,
        withspecranges="yes",
        expression=expression,
    )
    run_tool("evselect", work, **args)
    run_tool("backscale", work, spectrumset=source.name, badpixlocation=event.name)
    first_backscal = spectrum.read(source).backscal
    run_tool(
        "arfgen",
        work,
        spectrumset=source.name,
        arfset=arf.name,
        extendedsource="yes",
        detmaptype="flat",
        withbadpixcorr="yes",
        detxbins=400,
        detybins=400,
        badpixlocation=event.name,
        filterdss="yes",
        setbackscale="yes",
        withrmfset="no",
        modelootcorr="no",
    )
    run_tool(
        "rmfgen",
        work,
        spectrumset=source.name,
        rmfset=rmf.name,
        detmaptype="flat",
        format="var",
        threshold=1e-6,
    )
    if not np.isclose(first_backscal, spectrum.read(source).backscal, rtol=0.01):
        raise RuntimeError(f"BACKSCAL disagreement for {frame['frame']}/{region_name}")
    if frame["oot"]:
        run_tool(
            "evselect",
            work,
            **{**args, "table": "oot.fits:EVENTS", "spectrumset": "oot.pi"},
        )
    fractions = (
        map_fraction(frame["qpb_map"], region, frame["position"]),
        map_fraction(frame["sp_map"], region, frame["position"]),
    )
    background = work / "background.pi"
    spectrum.background(
        background,
        source,
        frame["qpb_pha"],
        frame["sp_pha"],
        fractions,
        work / "oot.pi" if frame["oot"] else None,
        oot_scale,
    )
    pha = spectrum.read(source)
    omega = solid_angle(region, frame["position"], sources)
    mapped_exposure = mean_exposure(
        frame["exposure_map"], region, frame["position"], omega
    )
    throughput = mapped_exposure / pha.exposure
    if not 0 < throughput <= 1.01:
        raise RuntimeError(
            f"invalid exposure-map throughput for {frame['frame']}/{region_name}"
        )
    low_ev, high_ev = (1000 * float(value) for value in fit_band)
    return {
        "frame": frame["frame"],
        "det": frame["det"],
        "region": region_name,
        "source_pha": source,
        "background_pha": background,
        "rmf": rmf,
        "arf": arf,
        "exposure_s": pha.exposure,
        "omega_arcmin2": omega,
        "mean_exposure_s": mapped_exposure,
        "response_throughput": throughput,
        "qpb_fraction": fractions[0],
        "sp_fraction": fractions[1],
        "source_counts": spectrum.band_counts(source, low_ev, high_ev),
        "background_counts": spectrum.band_counts(background, low_ev, high_ev),
    }
