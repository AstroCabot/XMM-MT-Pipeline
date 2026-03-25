#!/usr/bin/env python3
"""Build masked, background, and per-instrument rate products for QC mosaics.

Called by xmm_qc_build_mosaics.sh after evselect + eexpmap have produced
raw per-instrument count images and exposure maps.

Products created underneath --stillsky-root and --comet-root:

Still-sky products (per band, per instrument + EPIC):
  {inst}_source_mask.fits    — circular mask around detected field sources
  {inst}_counts_masked.fits  — counts with masked pixels zeroed
  {inst}_exp_masked.fits     — exposure with masked pixels zeroed
  {inst}_bkg.fits            — smoothed background estimate
  {inst}_rate.fits           — counts / exposure
  {inst}_rate_masked.fits    — masked counts / masked exposure

Comet-frame products (per band, per instrument):
  {inst}_rate.fits           — per-instrument rate (counts / exposure)
  {inst}_bkg.fits            — smoothed background in comet frame
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
from astropy.io import fits

from trail_mask_utils import read_srclist, read_track, paint_disk


# ---------------------------------------------------------------------------
# I/O helpers
# ---------------------------------------------------------------------------


def first_image_hdu(hdul: fits.HDUList) -> tuple[int, fits.ImageHDU | fits.PrimaryHDU]:
    for idx, hdu in enumerate(hdul):
        if hdu.data is not None and getattr(hdu.data, "ndim", 0) == 2:
            return idx, hdu
    raise RuntimeError("No 2D image HDU found")


def read_image(path: str) -> tuple[np.ndarray, fits.Header]:
    with fits.open(path) as hdul:
        _, hdu = first_image_hdu(hdul)
        data = np.asarray(hdu.data, dtype=float)
        header = hdu.header.copy()
    return data, header


def write_like(header: fits.Header, data: np.ndarray, out_path: str) -> None:
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    fits.PrimaryHDU(data=data, header=header).writeto(out_path, overwrite=True)


# ---------------------------------------------------------------------------
# Source masking (for still-sky images where sources are point-like)
# ---------------------------------------------------------------------------


def build_source_mask(
    shape: tuple[int, int],
    header: fits.Header,
    srclist_fits: str,
    mask_radius_arcsec: float,
) -> np.ndarray:
    """Build a boolean mask with discs around detected field sources.

    Uses the WCS in the image header (CRPIX/CRVAL/CDELT or CD matrix) to
    convert source RA/Dec to pixel positions.  Falls back to the custom
    still-sky mosaic header keywords (CTRRA/CTRDEC/BIN_AS/XMIN_AS/YMIN_AS)
    if standard WCS keywords are absent.
    """
    src = read_srclist(srclist_fits)
    mask = np.zeros(shape, dtype=bool)
    if len(src.ra) == 0:
        return mask

    # Try standard WCS from SAS evselect images first
    try:
        from astropy.wcs import WCS

        wcs = WCS(header, naxis=2)
        if wcs.has_celestial:
            from astropy.coordinates import SkyCoord
            import astropy.units as u

            coords = SkyCoord(ra=src.ra * u.deg, dec=src.dec * u.deg, frame="icrs")
            px, py = wcs.world_to_pixel(coords)
            cdelt = header.get("CDELT1", header.get("CD1_1", None))
            if cdelt is not None:
                r_px = mask_radius_arcsec / (abs(float(cdelt)) * 3600.0)
            else:
                r_px = mask_radius_arcsec / 4.0  # assume ~4"/px default
            for x, y in zip(px, py):
                if np.isfinite(x) and np.isfinite(y):
                    paint_disk(mask, float(x), float(y), r_px)
            return mask
    except Exception:
        pass

    # Fallback: custom mosaic header (CTRRA, BIN_AS, XMIN_AS, …)
    import math

    ctr_ra = header.get("CTRRA")
    ctr_dec = header.get("CTRDEC")
    bin_as = header.get("BIN_AS")
    xmin_as = header.get("XMIN_AS")
    ymin_as = header.get("YMIN_AS")
    if ctr_ra is None or ctr_dec is None or bin_as is None:
        return mask  # cannot locate sources without WCS info

    cos_dec = math.cos(math.radians(float(ctr_dec)))
    r_px = mask_radius_arcsec / float(bin_as)
    for ra_s, dec_s in zip(src.ra, src.dec):
        dx_as = (
            ((float(ra_s) - float(ctr_ra) + 180.0) % 360.0 - 180.0) * cos_dec * 3600.0
        )
        dy_as = (float(dec_s) - float(ctr_dec)) * 3600.0
        col = (dx_as - float(xmin_as)) / float(bin_as)
        row = (dy_as - float(ymin_as)) / float(bin_as)
        paint_disk(mask, col, row, r_px)
    return mask


# ---------------------------------------------------------------------------
# Smoothed background estimation
# ---------------------------------------------------------------------------


def gaussian_smooth(data: np.ndarray, sigma_px: float) -> np.ndarray:
    """Smooth a 2D image with a Gaussian kernel, ignoring NaN/zero pixels."""
    from scipy.ndimage import gaussian_filter

    valid = np.isfinite(data) & (data != 0)
    filled = np.where(valid, data, 0.0)
    weight = valid.astype(float)
    sm_data = gaussian_filter(filled, sigma=sigma_px, mode="constant", cval=0.0)
    sm_weight = gaussian_filter(weight, sigma=sigma_px, mode="constant", cval=0.0)
    with np.errstate(divide="ignore", invalid="ignore"):
        result = np.where(sm_weight > 0, sm_data / sm_weight, np.nan)
    return result


def estimate_background(
    counts: np.ndarray,
    exposure: np.ndarray,
    source_mask: np.ndarray,
    sigma_px: float = 20.0,
) -> np.ndarray:
    """Estimate a background map from source-masked rate, rescaled to counts.

    1. Compute rate = counts / exposure (masking zero-exposure pixels)
    2. Zero out source-masked pixels
    3. Gaussian-smooth the masked rate to fill in the gaps
    4. Multiply by exposure to return a background count estimate
    """
    with np.errstate(divide="ignore", invalid="ignore"):
        rate = np.where(exposure > 0, counts / exposure, 0.0)
    rate[source_mask] = 0.0
    smooth_rate = gaussian_smooth(rate, sigma_px)
    bkg_counts = smooth_rate * np.where(exposure > 0, exposure, 0.0)
    return bkg_counts


# ---------------------------------------------------------------------------
# Product builders
# ---------------------------------------------------------------------------


def process_band_dir(
    banddir: Path,
    instruments: list[str],
    srclist_fits: str | None,
    mask_radius_arcsec: float,
    bkg_sigma_px: float,
) -> dict[str, object]:
    """Process one band directory: build mask, background, and rate for each instrument."""
    stats: dict[str, object] = {"instruments": {}}

    for inst in instruments + ["EPIC"]:
        cpath = banddir / f"{inst}_counts.fits"
        epath = banddir / f"{inst}_exp.fits"
        if not cpath.exists() or not epath.exists():
            continue

        counts, chdr = read_image(str(cpath))
        expo, _ = read_image(str(epath))

        # Rate image
        with np.errstate(divide="ignore", invalid="ignore"):
            rate = np.where(expo > 0, counts / expo, np.nan)
        rate_hdr = chdr.copy()
        rate_hdr["BUNIT"] = "count / s"
        write_like(rate_hdr, rate, str(banddir / f"{inst}_rate.fits"))

        # Source mask
        if srclist_fits and Path(srclist_fits).exists():
            smask = build_source_mask(
                counts.shape, chdr, srclist_fits, mask_radius_arcsec
            )
        else:
            smask = np.zeros(counts.shape, dtype=bool)
        mask_hdr = chdr.copy()
        mask_hdr["BUNIT"] = "1"
        mask_hdr["HISTORY"] = "Source mask (1=masked)"
        write_like(
            mask_hdr, smask.astype(np.int16), str(banddir / f"{inst}_source_mask.fits")
        )

        # Masked counts & exposure
        masked_counts = counts.copy()
        masked_counts[smask] = 0.0
        write_like(chdr, masked_counts, str(banddir / f"{inst}_counts_masked.fits"))

        masked_expo = expo.copy()
        masked_expo[smask] = 0.0
        exp_hdr = chdr.copy()
        exp_hdr.setdefault("BUNIT", "s")
        write_like(exp_hdr, masked_expo, str(banddir / f"{inst}_exp_masked.fits"))

        # Background map
        bkg = estimate_background(counts, expo, smask, sigma_px=bkg_sigma_px)
        bkg_hdr = chdr.copy()
        bkg_hdr["BUNIT"] = "counts (background estimate)"
        bkg_hdr["HISTORY"] = "Smoothed background from source-masked rate * exposure"
        write_like(bkg_hdr, bkg, str(banddir / f"{inst}_bkg.fits"))

        # Masked rate
        with np.errstate(divide="ignore", invalid="ignore"):
            masked_rate = np.where(masked_expo > 0, masked_counts / masked_expo, np.nan)
        write_like(rate_hdr, masked_rate, str(banddir / f"{inst}_rate_masked.fits"))

        stats["instruments"][inst] = {
            "mask_fraction": float(np.mean(smask)),
            "finite_rate_pixels": int(np.count_nonzero(np.isfinite(rate))),
        }

    return stats


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--config", help="Pipeline config .env (unused, reserved)")
    ap.add_argument(
        "--stillsky-root", required=True, help="Root dir for still-sky products"
    )
    ap.add_argument(
        "--comet-root", required=True, help="Root dir for comet-frame products"
    )
    ap.add_argument("--detect-srclist", help="Field source list FITS for masking")
    ap.add_argument("--grid-json", help="Image grid JSON from xmm_event_geometry.py")
    ap.add_argument("--track", help="Comet track FITS")
    ap.add_argument("--mask-radius-arcsec", type=float, default=20.0)
    ap.add_argument(
        "--bkg-sigma-px",
        type=float,
        default=20.0,
        help="Gaussian smoothing sigma in pixels for background estimation",
    )
    return ap.parse_args()


def main() -> int:
    args = parse_args()
    stillsky_root = Path(args.stillsky_root)
    comet_root = Path(args.comet_root)
    instruments = ["PN", "M1", "M2"]
    srclist = (
        args.detect_srclist
        if args.detect_srclist and Path(args.detect_srclist).exists()
        else None
    )

    # Process still-sky bands
    for banddir in sorted(stillsky_root.iterdir()):
        if not banddir.is_dir():
            continue
        print(f"  stillsky/{banddir.name}: ", end="", flush=True)
        stats = process_band_dir(
            banddir, instruments, srclist, args.mask_radius_arcsec, args.bkg_sigma_px
        )
        n = len(stats.get("instruments", {}))
        print(f"{n} instruments processed")

    # Process comet-frame bands
    for banddir in sorted(comet_root.iterdir()):
        if not banddir.is_dir():
            continue
        print(f"  comet/{banddir.name}: ", end="", flush=True)
        stats = process_band_dir(
            banddir, instruments, srclist, args.mask_radius_arcsec, args.bkg_sigma_px
        )
        n = len(stats.get("instruments", {}))
        print(f"{n} instruments processed")

    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(1)
