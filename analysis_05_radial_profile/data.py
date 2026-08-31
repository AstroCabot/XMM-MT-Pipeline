import json
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS

from analysis_02_pps_luminosity import run as pps
from analysis_02_pps_luminosity.selection import band_mask

HERE = Path(__file__).resolve().parent
REPRO = HERE.parent
IMAGE = REPRO / "analysis_04_morphology/outputs/image.fits"
PARAMS = json.loads((REPRO / "parameters.json").read_text())
SOFT_BAND = tuple(map(int, PARAMS["soft_band_ev"]))
HARD_BAND = tuple(map(int, PARAMS["hard_band_ev"]))
PPS_GRID_KEYS = (
    "CTYPE1",
    "CTYPE2",
    "CRPIX1",
    "CRPIX2",
    "CRVAL1",
    "CRVAL2",
    "CDELT1",
    "CDELT2",
)


def energy_masks(pi):
    return (
        band_mask(pi, *SOFT_BAND, include_high=True),
        band_mask(pi, *HARD_BAND, include_high=True),
    )


def radii(header, shape):
    rows, columns = np.indices(shape)
    wcs = WCS(header)
    cx, cy = wcs.wcs.crpix - 1
    scale = abs(float(wcs.wcs.cdelt[0])) * 60
    dx, dy = (cx - columns) * scale, (rows - cy) * scale
    return np.hypot(dx, dy), np.arctan2(dx, dy)


def epic_pixels():
    with fits.open(IMAGE, memmap=True) as hdul:
        planes = {
            name: np.asarray(hdul[name].data, float)
            for name in ("EVENTS", "QPB", "SP", "OOT", "EXPOSURE")
        }
        header = hdul["EVENTS"].header.copy()
    if (
        header.get("DETECTOR") != "EPIC"
        or header.get("BAND") != "soft"
        or [header.get("E_MIN"), header.get("E_MAX")] != list(SOFT_BAND)
        or not header.get("COMOVING", False)
    ):
        raise RuntimeError("Task 04 image metadata differs from the shared soft band")
    radius, angle = radii(header, planes["EVENTS"].shape)
    planes["radius"] = radius
    planes["sector"] = np.floor((angle + np.pi) / (np.pi / 6)).astype(int) % 12
    planes["pixel_area_arcmin2"] = abs(float(np.prod(WCS(header).wcs.cdelt[:2]))) * 3600
    return planes


def aggregate_epic(planes, edges, excluded_sector=None):
    radius = planes["radius"]
    ring = np.searchsorted(edges, radius, side="right") - 1
    use = (ring >= 0) & (ring < len(edges) - 1) & (planes["EXPOSURE"] > 0)
    if excluded_sector is not None:
        use &= planes["sector"] != excluded_sector
    index = ring[use]
    count = lambda values: np.bincount(index, values[use], minlength=len(edges) - 1)
    events = count(planes["EVENTS"])
    qpb, sp, oot = (count(planes[name]) for name in ("QPB", "SP", "OOT"))
    weight = planes["EXPOSURE"][use] * planes["pixel_area_arcmin2"]
    geometric = np.bincount(
        ring[(ring >= 0) & (ring < len(edges) - 1)], minlength=len(edges) - 1
    )
    covered = np.bincount(index, minlength=len(edges) - 1)
    return {
        "counts": events,
        "qpb": qpb,
        "sp": sp,
        "oot": oot,
        "fixed": qpb + sp + PARAMS["pn_oot_fraction"] * oot,
        "variance": events + qpb + sp + PARAMS["pn_oot_fraction"] ** 2 * oot,
        "area": np.bincount(index, weight, minlength=len(edges) - 1),
        "coverage": covered / np.maximum(geometric, 1),
        "radius": radius[use],
        "weight": weight,
        "ring": index,
    }


def _gti_mask(data, gtis):
    keep = np.zeros(len(data), bool)
    for ccd, intervals in enumerate(gtis, 1):
        rows = np.flatnonzero(data["CCDNR"] == ccd)
        if not len(rows):
            continue
        times = np.asarray(data["TIME"][rows], float)
        accepted = np.zeros(len(rows), bool)
        for start, stop in intervals:
            accepted |= (times >= start) & (times <= stop)
        keep[rows[accepted]] = True
    return keep


def pps_maps():
    maps, grid = [], None
    for name in pps.EXPMAPS:
        image, header = fits.getdata(pps.MT_PPS / name, header=True)
        image = np.asarray(image, float)
        if "RADESYS" not in header and "RADECSYS" in header:
            header["RADESYS"] = header.pop("RADECSYS")
        wcs = WCS(header, fix=False)
        center = wcs.all_world2pix([pps.REFERENCE], 0)[0]
        linear = tuple(np.asarray(wcs.pixel_scale_matrix, float).ravel())
        current = (image.shape, tuple(header.get(key) for key in PPS_GRID_KEYS), linear)
        expected_linear = np.diag([-pps.MAP_PIXEL / 3600, pps.MAP_PIXEL / 3600])
        if np.max(np.abs(center - pps.MAP_CENTER)) > 1e-4 or not np.allclose(
            wcs.pixel_scale_matrix, expected_linear, rtol=0, atol=1e-14
        ):
            raise RuntimeError(f"unexpected PPS exposure grid: {name}")
        if grid is None:
            grid = current
        elif current != grid:
            raise RuntimeError("PPS exposure maps are on different grids")
        maps.append(image)
    return maps


def pps_profile(edges):
    event_path = pps.MT_PPS / pps.EVENT
    gtis = pps.read_gtis(pps.MT_PPS / pps.GTI_IMAGE)
    sources, ephemeris, _ = pps.motion_inputs()
    with fits.open(event_path, memmap=True) as hdul:
        data, header = hdul["EVENTS"].data, hdul["EVENTS"].header.copy()
        quality = (
            (data["RAWY"] > 12)
            & (data["PATTERN"] <= 4)
            & (data["FLAG"] == 0)
            & _gti_mask(data, gtis)
        )
        pi = np.asarray(data["PI"], int)
        soft_energy, hard_energy = energy_masks(pi)
        rows = np.flatnonzero(quality & (soft_energy | hard_energy))
        time = np.asarray(data["TIME"][rows], float)
        mjd = pps.event_mjd(time, header)
        keep = pps.source_keep(
            data["X"][rows],
            data["Y"][rows],
            mjd,
            sources,
            ephemeris,
            PARAMS["source_mask_radius_arcsec"],
        )
        radius_event = (
            np.hypot(
                data["X"][rows] - pps.EVENT_CENTER[0],
                data["Y"][rows] - pps.EVENT_CENTER[1],
            )
            * pps.EVENT_PIXEL
            / 60
        )
        selected_pi = pi[rows]
    radius_event, selected_pi = radius_event[keep], selected_pi[keep]
    event_ring = np.searchsorted(edges, radius_event, side="right") - 1
    within = (event_ring >= 0) & (event_ring < len(edges) - 1)
    soft_energy, hard_energy = energy_masks(selected_pi)
    soft, hard = within & soft_energy, within & hard_energy

    maps = pps_maps()
    coverage = pps.exposure_coverage(
        gtis,
        header,
        sources,
        ephemeris,
        maps[0].shape,
        PARAMS["source_mask_radius_arcsec"],
    )
    yy, xx = np.indices(maps[0].shape)
    radius_map = (
        np.hypot(xx - pps.MAP_CENTER[0], yy - pps.MAP_CENTER[1]) * pps.MAP_PIXEL / 60
    )
    map_area = (pps.MAP_PIXEL / 60) ** 2
    inside = radius_map < edges[-1]
    rates = []
    for index, ((lo, hi), image) in enumerate(zip(pps.BANDS, maps, strict=True)):
        selected = band_mask(
            selected_pi, lo, hi, include_high=index == len(pps.BANDS) - 1
        )
        count = np.count_nonzero((radius_event < edges[-1]) & selected)
        rates.append(count / (np.sum(image[inside] * coverage[inside]) * map_area))
    weights = np.asarray(rates) / np.sum(rates)
    exposure = (
        sum(weight * image for weight, image in zip(weights, maps, strict=True))
        * coverage
    )
    map_ring = np.searchsorted(edges, radius_map, side="right") - 1
    map_use = (map_ring >= 0) & (map_ring < len(edges) - 1) & (exposure > 0)
    index, weight = map_ring[map_use], exposure[map_use] * map_area
    return {
        "soft": np.bincount(event_ring[soft], minlength=len(edges) - 1),
        "hard": np.bincount(event_ring[hard], minlength=len(edges) - 1),
        "area": np.bincount(index, weight, minlength=len(edges) - 1),
        "radius": radius_map[map_use],
        "weight": weight,
        "ring": index,
        "exposure_weights": weights,
    }
