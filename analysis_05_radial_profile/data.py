import json
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS

HERE = Path(__file__).resolve().parent
REPRO = HERE.parent
IMAGE = REPRO / "analysis_04_morphology/outputs/image.fits"
PARAMS = json.loads((REPRO / "parameters.json").read_text())
SOFT_BAND = tuple(map(int, PARAMS["soft_band_ev"]))


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
        "pixel_area_arcmin2": planes["pixel_area_arcmin2"],
        "weight": weight,
        "ring": index,
    }
