import sys

import numpy as np

from _common import DATA, OUT, REPRO, read_json, restore_env, stage_dir, write_json

sys.path.insert(0, str(REPRO / "analysis_04_morphology"))

import run as run4

restore_env()

SMOOTH_COUNTS = 50
KERNEL_RADII_PX = (1, 2, 3, 4, 6, 8, 11, 16, 22, 32, 45, 64, 90)


def adaptive_surrogate(events, qpb, sp, oot, exposure, header, mask=None):
    from scipy.signal import fftconvolve

    scale = run4.PARAMS["pn_oot_fraction"]
    net = events - qpb - sp - scale * oot
    pixel_area = abs(float(header["CDELT1"] * header["CDELT2"])) * 3600
    valid = exposure > 0 if mask is None else (np.asarray(mask, bool) & (exposure > 0))
    counts = np.where(valid, events, 0.0)
    net = np.where(valid, net, 0.0)
    live = np.where(valid, exposure, 0.0)
    output = np.full(events.shape, np.nan)
    remaining = valid.copy()
    for radius in KERNEL_RADII_PX:
        yy, xx = np.ogrid[-radius : radius + 1, -radius : radius + 1]
        disk = (xx * xx + yy * yy <= radius * radius).astype(float)
        total = fftconvolve(counts, disk, "same")
        signal = fftconvolve(net, disk, "same")
        area = fftconvolve(live, disk, "same") * pixel_area
        done = remaining & (
            (total >= SMOOTH_COUNTS - 0.5) | (radius == KERNEL_RADII_PX[-1])
        )
        done &= area > 0
        output[done] = signal[done] / area[done]
        remaining &= ~done
        if not remaining.any():
            break
    if mask is not None:
        output[~np.asarray(mask, bool)] = np.nan
    return output


def main():
    output = stage_dir(4)
    run4.STACK = DATA / "stack04" / "stack"
    run4.SPECTRUM = OUT / "03" / "result.json"
    run4.OUTPUT = output
    run4.DETECTORS = ("pn",)
    run4.adaptive = adaptive_surrogate
    run4.main()

    conversion = read_json(OUT / "03" / "result.json")["detector_K_ct_s_per_erg_cm2_s"]
    result = read_json(output / "result.json")
    result["detector_exposure_weights"] = {
        name: conversion[name] / conversion["mos2"] for name in ("pn", "mos1", "mos2")
    }
    write_json(output / "result.json", result)
    print(
        f"stage 04: lowest-contour extent {result['lowest_contour_extent_arcmin']:.2f}"
        " arcmin, Figure 1 written"
    )


if __name__ == "__main__":
    main()
