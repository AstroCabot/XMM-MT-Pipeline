#!/usr/bin/env python3
"""Diagnose counts-vs-exposure alignment in PN per-detector images.

Reads the PN counts and exposure maps, measures:
1. FITS header WCS comparison (CRPIX, CRVAL, CDELT)
2. Cross-correlation offset between counts FOV and exposure FOV
3. RA_NOM/DEC_NOM vs RA_PNT/DEC_PNT from event file headers
4. Edge profile asymmetry (counts vs exposure boundary)

Produces a diagnostic PNG with overlay and offset measurements.
"""
import sys
import json
import math
import numpy as np
from pathlib import Path
from astropy.io import fits
from astropy.coordinates import SkyCoord
import astropy.units as u
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.ndimage import shift as ndshift


def read_image(path):
    with fits.open(path) as hdul:
        for ext in hdul:
            if ext.data is not None and ext.data.ndim == 2:
                return ext.data.astype(np.float64), ext.header
    raise RuntimeError(f"No 2D image in {path}")


def fov_mask(data, threshold_frac=0.01):
    """Boolean mask of the illuminated FOV."""
    finite = np.isfinite(data) & (data > 0)
    if not np.any(finite):
        return finite
    thresh = np.nanpercentile(data[finite], threshold_frac * 100)
    return data > thresh


def centroid(data, mask):
    """Weighted centroid of data within mask."""
    yy, xx = np.mgrid[:data.shape[0], :data.shape[1]]
    w = np.where(mask, data, 0.0)
    total = w.sum()
    if total == 0:
        return np.nan, np.nan
    return (xx * w).sum() / total, (yy * w).sum() / total


def cross_correlate_fov(counts, expo):
    """Cross-correlate the FOV boundaries to find pixel offset."""
    # Use binary FOV masks — this isolates the geometric offset
    # from the brightness distribution
    c_fov = fov_mask(counts, 0.01).astype(float)
    e_fov = fov_mask(expo, 0.01).astype(float)

    # FFT cross-correlation
    from numpy.fft import fft2, ifft2
    f1 = fft2(c_fov)
    f2 = fft2(e_fov)
    cc = np.real(ifft2(f1 * np.conj(f2)))
    # Peak of cross-correlation
    peak = np.unravel_index(np.argmax(cc), cc.shape)
    dy = peak[0] if peak[0] < cc.shape[0] // 2 else peak[0] - cc.shape[0]
    dx = peak[1] if peak[1] < cc.shape[1] // 2 else peak[1] - cc.shape[1]
    return dx, dy


def edge_profiles(data, mask):
    """Sum along rows and columns to get edge profiles."""
    masked = np.where(mask, data, 0.0)
    row_sum = masked.sum(axis=1)
    col_sum = masked.sum(axis=0)
    return row_sum, col_sum


def main():
    base = Path("/mnt/c/Users/cabot/OneDrive/Documents/Research/Cambridge/X3I/reduction_v4/output")
    counts_path = base / "images" / "soft" / "PN_counts.fits"
    exp_path = base / "images" / "soft" / "PN_exp.fits"
    evt_path = base / "comet" / "PN_comet.fits"
    grid_path = base / "images" / "grid.json"
    outdir = base / "qc"

    print("=" * 60)
    print("COUNTS vs EXPOSURE ALIGNMENT DIAGNOSTIC")
    print("=" * 60)

    # --- 1. Read images and compare WCS headers ---
    counts, c_hdr = read_image(str(counts_path))
    expo, e_hdr = read_image(str(exp_path))
    print(f"\nImage shapes: counts={counts.shape}, expo={expo.shape}")

    wcs_keys = ["CRPIX1", "CRPIX2", "CRVAL1", "CRVAL2", "CDELT1", "CDELT2",
                "REFXCRPX", "REFYCRPX", "REFXCRVL", "REFYCRVL"]
    print("\nWCS header comparison:")
    print(f"  {'key':<12} {'counts':<20} {'exposure':<20} {'diff':<15}")
    for key in wcs_keys:
        cv = c_hdr.get(key)
        ev = e_hdr.get(key)
        if cv is not None or ev is not None:
            if isinstance(cv, (int, float)) and isinstance(ev, (int, float)):
                diff = cv - ev
                print(f"  {key:<12} {cv:<20.8f} {ev:<20.8f} {diff:<15.8f}")
            else:
                print(f"  {key:<12} {str(cv):<20} {str(ev):<20}")

    # Also check RA_NOM, DEC_NOM in the image headers
    for key in ["RA_NOM", "DEC_NOM", "RA_PNT", "DEC_PNT"]:
        cv = c_hdr.get(key)
        ev = e_hdr.get(key)
        if cv is not None or ev is not None:
            if isinstance(cv, (int, float)) and isinstance(ev, (int, float)):
                diff = (cv - ev) * 3600
                print(f"  {key:<12} {cv:<20.8f} {ev:<20.8f} {diff:<15.4f} arcsec")
            else:
                print(f"  {key:<12} {str(cv):<20} {str(ev):<20}")

    # --- 2. Event file header keywords ---
    print("\nEvent file header (PN_comet.fits):")
    with fits.open(str(evt_path)) as hdul:
        for ext in hdul:
            if ext.name == "EVENTS":
                h = ext.header
                ra_nom = h.get("RA_NOM", None)
                dec_nom = h.get("DEC_NOM", None)
                ra_pnt = h.get("RA_PNT", None)
                dec_pnt = h.get("DEC_PNT", None)
                print(f"  RA_NOM  = {ra_nom}")
                print(f"  DEC_NOM = {dec_nom}")
                print(f"  RA_PNT  = {ra_pnt}")
                print(f"  DEC_PNT = {dec_pnt}")
                if all(v is not None for v in [ra_nom, dec_nom, ra_pnt, dec_pnt]):
                    nom = SkyCoord(ra=ra_nom*u.deg, dec=dec_nom*u.deg)
                    pnt = SkyCoord(ra=ra_pnt*u.deg, dec=dec_pnt*u.deg)
                    print(f"  NOM-PNT sep = {nom.separation(pnt).to_value(u.arcsec):.2f} arcsec")
                    dlon, dlat = nom.spherical_offsets_to(pnt)
                    print(f"  NOM→PNT offset: dRA*cos(dec)={dlon.to_value(u.arcsec):.2f}\", "
                          f"dDec={dlat.to_value(u.arcsec):.2f}\"")
                break

    # --- 3. FOV centroid comparison ---
    c_mask = fov_mask(counts, 0.5)
    e_mask = fov_mask(expo, 0.5)
    cx, cy = centroid(counts, c_mask)
    ex, ey = centroid(expo, e_mask)
    print(f"\nFOV centroids:")
    print(f"  Counts:   ({cx:.1f}, {cy:.1f})")
    print(f"  Exposure: ({ex:.1f}, {ey:.1f})")
    print(f"  Offset:   dx={cx-ex:.2f} pix, dy={cy-ey:.2f} pix")

    # Get pixel scale
    pixel_arcsec = 2.0  # default
    if grid_path.exists():
        grid = json.loads(grid_path.read_text())
        scale = float(grid.get("scale_arcsec_per_phys", 1.0))
        bin_phys = float(grid.get("bin_phys", 1.0))
        pixel_arcsec = scale * bin_phys
        print(f"  Pixel scale: {pixel_arcsec:.3f} arcsec/pixel")
    dx_arcsec = (cx - ex) * pixel_arcsec
    dy_arcsec = (cy - ey) * pixel_arcsec
    print(f"  Offset:   dx={dx_arcsec:.1f}\", dy={dy_arcsec:.1f}\"")

    # --- 4. Cross-correlation of FOV masks ---
    dx_cc, dy_cc = cross_correlate_fov(counts, expo)
    print(f"\nCross-correlation (FOV mask):")
    print(f"  dx={dx_cc} pix, dy={dy_cc} pix")
    print(f"  dx={dx_cc * pixel_arcsec:.1f}\", dy={dy_cc * pixel_arcsec:.1f}\"")

    # --- 5. Edge profile comparison ---
    # Project FOV onto rows (Y axis) — this should show the offset most clearly
    c_row = (counts > 0).astype(float).sum(axis=1)
    e_row = (expo > 0).astype(float).sum(axis=1)
    # Find the FOV edges (where the row sum drops below 50% of max)
    c_half = c_row.max() * 0.5
    e_half = e_row.max() * 0.5
    c_top = np.max(np.where(c_row > c_half)) if np.any(c_row > c_half) else 0
    c_bot = np.min(np.where(c_row > c_half)) if np.any(c_row > c_half) else 0
    e_top = np.max(np.where(e_row > e_half)) if np.any(e_row > e_half) else 0
    e_bot = np.min(np.where(e_row > e_half)) if np.any(e_row > e_half) else 0
    print(f"\nFOV Y-extent (row indices at 50% illumination):")
    print(f"  Counts: bottom={c_bot}, top={c_top}, height={c_top-c_bot}")
    print(f"  Exposure: bottom={e_bot}, top={e_top}, height={e_top-e_bot}")
    print(f"  Top edge shift: {c_top - e_top} pix = {(c_top-e_top)*pixel_arcsec:.1f}\"")
    print(f"  Bottom edge shift: {c_bot - e_bot} pix = {(c_bot-e_bot)*pixel_arcsec:.1f}\"")
    print(f"  Center shift: {((c_top+c_bot)/2 - (e_top+e_bot)/2):.1f} pix "
          f"= {((c_top+c_bot)/2 - (e_top+e_bot)/2)*pixel_arcsec:.1f}\"")

    # Same for X
    c_col = (counts > 0).astype(float).sum(axis=0)
    e_col = (expo > 0).astype(float).sum(axis=0)
    c_halfx = c_col.max() * 0.5
    e_halfx = e_col.max() * 0.5
    c_right = np.max(np.where(c_col > c_halfx)) if np.any(c_col > c_halfx) else 0
    c_left = np.min(np.where(c_col > c_halfx)) if np.any(c_col > c_halfx) else 0
    e_right = np.max(np.where(e_col > e_halfx)) if np.any(e_col > e_halfx) else 0
    e_left = np.min(np.where(e_col > e_halfx)) if np.any(e_col > e_halfx) else 0
    print(f"\nFOV X-extent (col indices at 50% illumination):")
    print(f"  Counts: left={c_left}, right={c_right}, width={c_right-c_left}")
    print(f"  Exposure: left={e_left}, right={e_right}, width={e_right-e_left}")
    print(f"  Right edge shift: {c_right - e_right} pix = {(c_right-e_right)*pixel_arcsec:.1f}\"")
    print(f"  Left edge shift: {c_left - e_left} pix = {(c_left-e_left)*pixel_arcsec:.1f}\"")
    print(f"  Center shift: {((c_right+c_left)/2 - (e_right+e_left)/2):.1f} pix "
          f"= {((c_right+c_left)/2 - (e_right+e_left)/2)*pixel_arcsec:.1f}\"")

    # --- 6. Diagnostic plot ---
    fig, axes = plt.subplots(2, 3, figsize=(18, 12), constrained_layout=True)

    # Panel 1: Counts FOV boundary
    ax = axes[0, 0]
    ax.imshow(counts, origin="lower", norm=matplotlib.colors.PowerNorm(0.3, vmin=0, vmax=3))
    ax.contour(c_mask, levels=[0.5], colors="cyan", linewidths=0.8)
    ax.contour(e_mask, levels=[0.5], colors="red", linewidths=0.8, linestyles="--")
    ax.plot(cx, cy, "c+", ms=15, mew=2, label="counts centroid")
    ax.plot(ex, ey, "r+", ms=15, mew=2, label="expo centroid")
    ax.legend(fontsize=8)
    ax.set_title("Counts + FOV contours\n(cyan=counts, red=expo)")

    # Panel 2: Exposure map
    ax = axes[0, 1]
    ax.imshow(expo, origin="lower", cmap="inferno")
    ax.contour(c_mask, levels=[0.5], colors="cyan", linewidths=0.8)
    ax.contour(e_mask, levels=[0.5], colors="red", linewidths=0.8, linestyles="--")
    ax.set_title("Exposure map + FOV contours")

    # Panel 3: Rate = counts/exposure (the problematic image)
    rate = np.full_like(counts, np.nan)
    good = expo > 0
    rate[good] = counts[good] / expo[good]
    ax = axes[0, 2]
    finite_rate = rate[np.isfinite(rate)]
    if len(finite_rate):
        vmin, vmax = np.percentile(finite_rate, [1, 99])
    else:
        vmin, vmax = 0, 1
    ax.imshow(rate, origin="lower", vmin=vmin, vmax=vmax)
    ax.set_title(f"Rate (counts/expo)\nvmin={vmin:.2e}, vmax={vmax:.2e}")

    # Panel 4: Edge profiles (Y axis)
    ax = axes[1, 0]
    y_idx = np.arange(len(c_row))
    ax.plot(c_row / c_row.max(), y_idx, "c-", label="counts", alpha=0.7)
    ax.plot(e_row / e_row.max(), y_idx, "r-", label="exposure", alpha=0.7)
    ax.axhline(c_top, color="cyan", ls=":", lw=0.7)
    ax.axhline(c_bot, color="cyan", ls=":", lw=0.7)
    ax.axhline(e_top, color="red", ls=":", lw=0.7)
    ax.axhline(e_bot, color="red", ls=":", lw=0.7)
    ax.set_ylabel("Row index (Y)")
    ax.set_xlabel("Normalized illumination")
    ax.legend(fontsize=8)
    ax.set_title(f"Y edge profiles\nY-shift: {((c_top+c_bot)/2 - (e_top+e_bot)/2)*pixel_arcsec:.1f}\"")

    # Panel 5: Edge profiles (X axis)
    ax = axes[1, 1]
    x_idx = np.arange(len(c_col))
    ax.plot(x_idx, c_col / c_col.max(), "c-", label="counts", alpha=0.7)
    ax.plot(x_idx, e_col / e_col.max(), "r-", label="exposure", alpha=0.7)
    ax.axvline(c_right, color="cyan", ls=":", lw=0.7)
    ax.axvline(c_left, color="cyan", ls=":", lw=0.7)
    ax.axvline(e_right, color="red", ls=":", lw=0.7)
    ax.axvline(e_left, color="red", ls=":", lw=0.7)
    ax.set_xlabel("Col index (X)")
    ax.set_ylabel("Normalized illumination")
    ax.legend(fontsize=8)
    ax.set_title(f"X edge profiles\nX-shift: {((c_right+c_left)/2 - (e_right+e_left)/2)*pixel_arcsec:.1f}\"")

    # Panel 6: Rate with shifted exposure (test correction)
    ax = axes[1, 2]
    # Shift the exposure map by the measured offset to see if it fixes the streak
    shift_y = (c_top + c_bot) / 2 - (e_top + e_bot) / 2
    shift_x = (c_right + c_left) / 2 - (e_right + e_left) / 2
    expo_shifted = ndshift(expo, [shift_y, shift_x], order=1, mode="constant", cval=0.0)
    rate_fixed = np.full_like(counts, np.nan)
    good2 = expo_shifted > 0
    rate_fixed[good2] = counts[good2] / expo_shifted[good2]
    if len(finite_rate):
        ax.imshow(rate_fixed, origin="lower", vmin=vmin, vmax=vmax)
    ax.set_title(f"Rate with shifted expo\n(shift: dx={shift_x:.1f}px={shift_x*pixel_arcsec:.0f}\", "
                 f"dy={shift_y:.1f}px={shift_y*pixel_arcsec:.0f}\")")

    fig.suptitle("PN Counts vs Exposure Alignment Diagnostic", fontsize=14, fontweight="bold")
    png = outdir / "diag_pn_alignment.png"
    fig.savefig(png, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"\nDiagnostic plot saved: {png}")
    print("=" * 60)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
