import csv
import json
import os
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPRO = HERE.parent
RUNTIME = REPRO / "work" / "runtime" / "analysis_04_morphology"
for name in ("tmp", "cache", "config", "matplotlib"):
    (RUNTIME / name).mkdir(parents=True, exist_ok=True)
os.environ.update(
    TMPDIR=str(RUNTIME / "tmp"),
    XDG_CACHE_HOME=str(RUNTIME / "cache"),
    XDG_CONFIG_HOME=str(RUNTIME / "config"),
    MPLCONFIGDIR=str(RUNTIME / "matplotlib"),
)

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from scipy.ndimage import label, median_filter

from figure import CONTOURS, draw
from smooth import SETTINGS, adaptive

STACK = REPRO / "analysis_01_reduction" / "outputs" / "07_comove" / "stack"
SPECTRUM = REPRO / "analysis_03_spectrum" / "outputs" / "result.json"
OUTPUT = HERE / "outputs"
FINAL_PRODUCTS = (
    "image.fits",
    "slice.tsv",
    "result.json",
    "figure1.png",
    "figure1.pdf",
)
PARAMS = json.loads((REPRO / "parameters.json").read_text())
DETECTORS = ("pn", "mos1", "mos2")
SUN_PA_DEG = 113.510
MOTION_PA_DEG = 289.90867
UTC_MID = "2025-12-03T18:59:16.860"
KM_PER_ARCMIN = 81920.50
CONTOUR_SEED_RADIUS_ARCMIN = 0.75
GRID_KEYS = ("CRPIX1", "CRPIX2", "CRVAL1", "CRVAL2", "CDELT1", "CDELT2")


def load_image(path):
    data, header = fits.getdata(path, header=True)
    return np.asarray(data, float), header


def combine(products, conversion):
    events = sum(products[d]["events"] for d in DETECTORS)
    qpb = sum(products[d]["qpb"] for d in DETECTORS)
    sp = sum(products[d]["sp"] for d in DETECTORS)
    oot = products["pn"]["oot"]
    reference = conversion["mos2"]
    exposure = sum(
        products[d]["exposure"] * conversion[d] / reference for d in DETECTORS
    )
    scale = PARAMS["pn_oot_fraction"]
    net = events - qpb - sp - scale * oot
    variance = events + qpb + sp + scale**2 * oot
    return events, qpb, sp, oot, exposure, net, variance


def coordinates(shape, header):
    rows, cols = np.indices(shape)
    wcs = WCS(header)
    cx, cy = wcs.wcs.crpix - 1
    scale = abs(float(wcs.wcs.cdelt[0])) * 60
    return (cx - cols) * scale, (rows - cy) * scale


def display_mask(exposure, pn_exposure):
    local = median_filter(pn_exposure, size=PARAMS["pn_gap_box_pixels"])
    gap = (pn_exposure < PARAMS["pn_gap_local_ratio"] * local) & (
        local > PARAMS["display_exposure_floor_fraction"] * pn_exposure.max()
    )
    return (
        exposure >= PARAMS["display_exposure_floor_fraction"] * exposure.max()
    ) & ~gap


def source_extent(detected, dx, dy):
    components, _ = label(detected, structure=np.ones((3, 3), int))
    radius = np.hypot(dx, dy)
    seeds = np.unique(components[detected & (radius <= CONTOUR_SEED_RADIUS_ARCMIN)])
    seeds = seeds[seeds > 0]
    if not len(seeds):
        raise RuntimeError("lowest contour does not reach the comet")
    connected = np.isin(components, seeds)
    return float(radius[connected].max())


def slice_rows(
    net, variance, exposure, smooth, mask, dx, dy, pixel_area, max_radius=None
):
    # Position angle runs north through east, so the sunward unit vector
    # is (sin, cos) in (east, north).
    angle = np.deg2rad(SUN_PA_DEG)
    along = dx * np.sin(angle) + dy * np.cos(angle)
    across = dx * np.cos(angle) - dy * np.sin(angle)
    width = float(PARAMS["profile_bin_width_arcmin"])
    radius = PARAMS["profile_fit_arcmin"][1] if max_radius is None else max_radius
    bins = round(float(radius) / width)
    edges = np.arange(-bins, bins + 1) * width
    rows = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        whole = (abs(across) <= 0.5) & (along >= lo) & (along < hi)
        use = whole & mask
        area_exposure = exposure[use].sum() * pixel_area
        sb = net[use].sum() / area_exposure if area_exposure else np.nan
        error = (
            np.sqrt(variance[use].sum()) / area_exposure if area_exposure else np.nan
        )
        smoothed = (
            np.average(smooth[use], weights=exposure[use]) if use.any() else np.nan
        )
        rows.append(
            [
                0.5 * (lo + hi),
                0.5 * (lo + hi) * KM_PER_ARCMIN,
                sb,
                error,
                use.sum() / max(whole.sum(), 1),
                smoothed,
            ]
        )
    reflected = np.array([r[-1] for r in rows])[::-1]
    return [row + [reflected[i]] for i, row in enumerate(rows)]


def centroid_sunward(net, exposure, mask, dx, dy, pixel_area, radius=2):
    symmetric = mask & mask[::-1, ::-1]
    use = symmetric & (np.hypot(dx, dy) <= radius) & (exposure > 0)
    weight = net[use] / (exposure[use] * pixel_area)
    total = weight.sum()
    east = np.sum(weight * dx[use]) / total
    north = np.sum(weight * dy[use]) / total
    angle = np.deg2rad(SUN_PA_DEG)
    return 60 * (east * np.sin(angle) + north * np.cos(angle))


def write_products(arrays, header, profile, conversion, extent):
    hdus = [fits.PrimaryHDU()]
    for name, data, unit in arrays:
        h = header.copy()
        for key in (
            "SIMPLE",
            "XTENSION",
            "BITPIX",
            "NAXIS",
            "NAXIS1",
            "NAXIS2",
            "PCOUNT",
            "GCOUNT",
            "EXTEND",
        ):
            h.remove(key, ignore_missing=True)
        h["BUNIT"] = unit
        dtype = {"EVENTS": np.int32, "DISPLAY_MASK": np.uint8}.get(name, np.float32)
        hdus.append(fits.ImageHDU(np.asarray(data, dtype), h, name=name))
    fits.HDUList(hdus).writeto(OUTPUT / "image.fits", overwrite=True)
    columns = (
        "distance_arcmin",
        "distance_km",
        "surface_brightness",
        "error",
        "coverage",
        "smooth",
        "reflected_smooth",
    )
    with (OUTPUT / "slice.tsv").open("w", newline="") as stream:
        writer = csv.writer(stream, delimiter="\t")
        writer.writerow(columns)
        writer.writerows(profile)
    result = {
        "band_eV": PARAMS["soft_band_ev"],
        "midpoint_utc": UTC_MID,
        "sun_position_angle_deg": SUN_PA_DEG,
        "motion_position_angle_deg": MOTION_PA_DEG,
        "km_per_arcmin": KM_PER_ARCMIN,
        "detector_exposure_weights": {
            d: conversion[d] / conversion["mos2"] for d in DETECTORS
        },
        "display_exposure_floor_fraction": PARAMS["display_exposure_floor_fraction"],
        "contour_seed_radius_arcmin": CONTOUR_SEED_RADIUS_ARCMIN,
        "contour_surface_brightness": list(CONTOURS),
        "lowest_contour_extent_arcmin": extent,
    }
    (OUTPUT / "result.json").write_text(json.dumps(result, indent=2) + "\n")


def load_stack():
    """Load all soft-band comet-frame planes on a shared grid."""
    planes = [
        (d, k) for d in DETECTORS for k in ("events", "qpb", "sp", "exposure")
    ] + [("pn", "oot")]
    products, header, grid = {}, None, None
    for detector, kind in planes:
        data, current = load_image(STACK / f"{detector}-{kind}-soft.fits")
        if (
            current.get("DETECTOR") != detector
            or current.get("BAND") != "soft"
            or (current.get("E_MIN"), current.get("E_MAX"))
            != tuple(PARAMS["soft_band_ev"])
            or not current.get("COMOVING", False)
        ):
            raise RuntimeError(f"invalid stack metadata for {detector}-{kind}")
        shape = (data.shape, tuple(current[key] for key in GRID_KEYS))
        if grid is None:
            header, grid = current, shape
        elif shape != grid:
            raise RuntimeError(f"{detector}-{kind}-soft.fits is on a different grid")
        products.setdefault(detector, {})[kind] = data
    pixel = np.abs([float(header["CDELT1"]), float(header["CDELT2"])]) * 3600
    if not np.allclose(pixel, SETTINGS["comove"]["pixel_arcsec"], rtol=0, atol=1e-9):
        raise RuntimeError("comet stack pixel scale differs from Task 01 settings")
    header["DETECTOR"] = "EPIC"
    return products, header


def main():
    OUTPUT.mkdir(parents=True, exist_ok=True)
    for name in FINAL_PRODUCTS:
        (OUTPUT / name).unlink(missing_ok=True)
    if not (STACK.parent / ".done").is_file():
        raise RuntimeError("Task 01 comet stack is incomplete")
    spectral = json.loads(SPECTRUM.read_text())
    expected_band = np.asarray(PARAMS["soft_band_ev"], float) / 1000
    if not np.allclose(spectral["fit_band_keV"], expected_band, rtol=0, atol=1e-12):
        raise RuntimeError("Task 03 fit band differs from the shared soft band")
    conversion = spectral["detector_K_ct_s_per_erg_cm2_s"]
    products, header = load_stack()
    events, qpb, sp, oot, exposure, net, variance = combine(products, conversion)
    if not np.allclose(events, np.rint(events), atol=1e-5):
        raise RuntimeError("event stack is not integral")
    coverage = exposure >= PARAMS["display_exposure_floor_fraction"] * exposure.max()
    mask = display_mask(exposure, products["pn"]["exposure"])
    smooth = adaptive(events, qpb, sp, oot, exposure, header, mask)
    display_smooth = adaptive(events, qpb, sp, oot, exposure, header)
    dx, dy = coordinates(events.shape, header)
    pixel_area = abs(float(header["CDELT1"] * header["CDELT2"])) * 3600
    profile = slice_rows(net, variance, exposure, smooth, mask, dx, dy, pixel_area)
    figure_profile = slice_rows(
        net,
        variance,
        exposure,
        display_smooth,
        coverage,
        dx,
        dy,
        pixel_area,
        max_radius=20,
    )
    centroid = centroid_sunward(net, exposure, mask, dx, dy, pixel_area)
    arrays = (
        ("EVENTS", np.rint(events), "count"),
        ("QPB", qpb, "count"),
        ("SP", sp, "count"),
        ("OOT", oot, "count"),
        ("EXPOSURE", exposure, "s"),
        ("DISPLAY_MASK", mask, "1"),
        ("NET", net, "count"),
        ("VARIANCE", variance, "count2"),
        ("SMOOTH", smooth, "count/s/arcmin2"),
    )
    detected = mask & np.isfinite(smooth) & (smooth >= CONTOURS[0])
    extent = source_extent(detected, dx, dy)
    write_products(arrays, header, profile, conversion, extent)
    draw(
        display_smooth,
        mask,
        dx,
        dy,
        figure_profile,
        SUN_PA_DEG,
        MOTION_PA_DEG,
        PARAMS["science_aperture_radius_arcsec"] / 60,
        OUTPUT,
        KM_PER_ARCMIN,
        centroid,
        coverage,
    )


if __name__ == "__main__":
    main()
