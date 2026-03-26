#!/usr/bin/env python3
"""Post-processing of comet-frame and QC mosaic images.

Two CLI subcommands:
    image_postproc.py science  -- produce net-rate and background products
                                  for the science image stage
    image_postproc.py qc       -- build rate, masked, and background images
                                  for QC still-sky and comet mosaics

Replaces the old build_qc_mosaics.py with expanded science-product support.
"""
from __future__ import annotations
import argparse, json, math, sys
from pathlib import Path
import numpy as np
from astropy.io import fits
from trail_mask_tools import read_srclist, paint_disk

INSTS = ["PN", "M1", "M2", "EPIC"]


def first_image_hdu(hdul: fits.HDUList) -> tuple[int, fits.ImageHDU | fits.PrimaryHDU]:
    for idx, hdu in enumerate(hdul):
        if hdu.data is not None and getattr(hdu.data, "ndim", 0) == 2:
            return idx, hdu
    raise RuntimeError("No 2D image HDU found")


def read_image(path: str) -> tuple[np.ndarray, fits.Header]:
    with fits.open(path) as hdul:
        _, hdu = first_image_hdu(hdul)
        return np.asarray(hdu.data, dtype=float), hdu.header.copy()


def write_like(header: fits.Header, data: np.ndarray, out_path: str | Path) -> None:
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fits.PrimaryHDU(data=data, header=header).writeto(out_path, overwrite=True)


def gaussian_smooth(
    data: np.ndarray, sigma_px: float, zero_is_invalid: bool = False
) -> np.ndarray:
    """Weighted Gaussian smooth that ignores NaN / zero pixels."""
    from scipy.ndimage import gaussian_filter

    valid = np.isfinite(data) & (~zero_is_invalid | (data != 0))
    filled = np.where(valid, data, 0.0)
    weight = valid.astype(float)
    sm_data = gaussian_filter(filled, sigma=sigma_px, mode="constant", cval=0.0)
    sm_weight = gaussian_filter(weight, sigma=sigma_px, mode="constant", cval=0.0)
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where(sm_weight > 0, sm_data / sm_weight, np.nan)


def estimate_background_rate(
    counts: np.ndarray, exposure: np.ndarray, trail_mask: np.ndarray, sigma_px: float
) -> np.ndarray:
    """Estimate the smooth background rate, masking trail-contaminated pixels."""
    with np.errstate(divide="ignore", invalid="ignore"):
        rate = np.where(exposure > 0.0, counts / exposure, np.nan)
    masked_rate = np.array(rate, copy=True)
    masked_rate[trail_mask] = np.nan
    return gaussian_smooth(masked_rate, sigma_px)


def build_source_mask(
    shape: tuple[int, int],
    header: fits.Header,
    srclist_fits: str,
    mask_radius_arcsec: float,
) -> np.ndarray:
    src = read_srclist(srclist_fits)
    mask = np.zeros(shape, dtype=bool)
    if len(src.ra) == 0:
        return mask
    try:
        from astropy.wcs import WCS
        from astropy.coordinates import SkyCoord
        import astropy.units as u

        wcs = WCS(header, naxis=2)
        if wcs.has_celestial:
            coords = SkyCoord(ra=src.ra * u.deg, dec=src.dec * u.deg, frame="icrs")
            px, py = wcs.world_to_pixel(coords)
            cdelt = header.get("CDELT1", header.get("CD1_1", None))
            r_px = (
                mask_radius_arcsec / (abs(float(cdelt)) * 3600.0)
                if cdelt is not None
                else mask_radius_arcsec / 4.0
            )
            for x, y in zip(px, py):
                if np.isfinite(x) and np.isfinite(y):
                    paint_disk(mask, float(x), float(y), r_px)
            return mask
    except Exception:
        pass
    ctr_ra = header.get("CTRRA")
    ctr_dec = header.get("CTRDEC")
    bin_as = header.get("BIN_AS")
    xmin_as = header.get("XMIN_AS")
    ymin_as = header.get("YMIN_AS")
    if ctr_ra is None or ctr_dec is None or bin_as is None:
        return mask
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


def estimate_background(
    counts: np.ndarray,
    exposure: np.ndarray,
    source_mask: np.ndarray,
    sigma_px: float = 20.0,
) -> np.ndarray:
    with np.errstate(divide="ignore", invalid="ignore"):
        rate = np.where(exposure > 0, counts / exposure, 0.0)
    rate[source_mask] = 0.0
    smooth_rate = gaussian_smooth(rate, sigma_px, zero_is_invalid=True)
    return smooth_rate * np.where(exposure > 0, exposure, 0.0)


def process_qc_band_dir(
    banddir: Path,
    instruments: list[str],
    srclist_fits: str | None,
    mask_radius_arcsec: float,
    bkg_sigma_px: float,
) -> dict[str, object]:
    stats: dict[str, object] = {"instruments": {}}
    for inst in instruments + ["EPIC"]:
        cpath = banddir / f"{inst}_counts.fits"
        epath = banddir / f"{inst}_exp.fits"
        if not cpath.exists() or not epath.exists():
            continue
        counts, chdr = read_image(str(cpath))
        expo, _ = read_image(str(epath))
        with np.errstate(divide="ignore", invalid="ignore"):
            rate = np.where(expo > 0, counts / expo, np.nan)
        rate_hdr = chdr.copy()
        rate_hdr["BUNIT"] = "count / s"
        write_like(rate_hdr, rate, banddir / f"{inst}_rate.fits")
        trail_mask_path = banddir / "EPIC_trail_mask.fits"
        if trail_mask_path.exists():
            tmask, _ = read_image(str(trail_mask_path))
            smask = np.asarray(tmask > 0, dtype=bool)
            mask_note = "Stationary-source trail mask (1=masked)"
        elif srclist_fits and Path(srclist_fits).exists():
            smask = build_source_mask(
                counts.shape, chdr, srclist_fits, mask_radius_arcsec
            )
            mask_note = "Source mask (1=masked)"
        else:
            smask = np.zeros(counts.shape, dtype=bool)
            mask_note = "Empty mask (1=masked)"
        mask_hdr = chdr.copy()
        mask_hdr["BUNIT"] = "1"
        mask_hdr["HISTORY"] = mask_note
        write_like(
            mask_hdr, smask.astype(np.int16), banddir / f"{inst}_source_mask.fits"
        )
        masked_counts = counts.copy()
        masked_counts[smask] = 0.0
        write_like(chdr, masked_counts, banddir / f"{inst}_counts_masked.fits")
        masked_expo = expo.copy()
        masked_expo[smask] = 0.0
        exp_hdr = chdr.copy()
        exp_hdr.setdefault("BUNIT", "s")
        write_like(exp_hdr, masked_expo, banddir / f"{inst}_exp_masked.fits")
        bkg = estimate_background(counts, expo, smask, sigma_px=bkg_sigma_px)
        bkg_hdr = chdr.copy()
        bkg_hdr["BUNIT"] = "counts (background estimate)"
        bkg_hdr["HISTORY"] = "Smoothed background from source-masked rate * exposure"
        write_like(bkg_hdr, bkg, banddir / f"{inst}_bkg.fits")
        with np.errstate(divide="ignore", invalid="ignore"):
            masked_rate = np.where(masked_expo > 0, masked_counts / masked_expo, np.nan)
        write_like(rate_hdr, masked_rate, banddir / f"{inst}_rate_masked.fits")
        stats["instruments"][inst] = {
            "mask_fraction": float(np.mean(smask)),
            "finite_rate_pixels": int(np.count_nonzero(np.isfinite(rate))),
        }
    return stats


def make_science_products(banddir: Path, sigma_px: float) -> dict[str, object]:
    """Create background-subtracted (net) rate and count images per band.

    Produces for each instrument:  bkg_rate, bkg_counts, rate_clean,
    net_rate, net_counts FITS images.  Trail-mask pixels are replaced
    with the smooth background estimate before subtraction.
    """
    mask_path = banddir / "EPIC_trail_mask.fits"
    trail_mask = None
    if mask_path.exists():
        mask_data, _ = read_image(str(mask_path))
        trail_mask = np.asarray(mask_data > 0, dtype=bool)
    summary: dict[str, object] = {"banddir": str(banddir), "products": {}}
    for inst in INSTS:
        cpath = banddir / f"{inst}_counts.fits"
        epath = banddir / f"{inst}_exp.fits"
        if not (cpath.exists() and epath.exists()):
            continue
        counts, chdr = read_image(str(cpath))
        expo, _ = read_image(str(epath))
        if trail_mask is None:
            mask = np.zeros_like(counts, dtype=bool)
        else:
            if trail_mask.shape != counts.shape:
                raise RuntimeError(
                    f"Trail-mask shape mismatch for {inst}: {trail_mask.shape} vs {counts.shape}"
                )
            mask = np.asarray(trail_mask, dtype=bool)
        with np.errstate(divide="ignore", invalid="ignore"):
            rate = np.where(expo > 0.0, counts / expo, np.nan)
        bkg_rate = estimate_background_rate(counts, expo, mask, sigma_px=sigma_px)
        bkg_rate = np.where(np.isfinite(bkg_rate), bkg_rate, 0.0)
        bkg_rate = np.where(expo > 0.0, bkg_rate, np.nan)
        bkg_counts = bkg_rate * np.where(expo > 0.0, expo, 0.0)
        rate_clean = np.array(rate, copy=True)
        rate_clean[mask] = bkg_rate[mask]
        net_rate = np.where(expo > 0.0, rate_clean - bkg_rate, np.nan)
        net_counts = counts - bkg_counts
        net_counts[mask] = 0.0
        rate_hdr = chdr.copy()
        rate_hdr["BUNIT"] = "count / s"
        bkg_rate_hdr = chdr.copy()
        bkg_rate_hdr["BUNIT"] = "count / s"
        bkg_rate_hdr["HISTORY"] = "Smoothed comet-frame background-rate estimate"
        counts_hdr = chdr.copy()
        counts_hdr["BUNIT"] = "counts"
        write_like(bkg_rate_hdr, bkg_rate, banddir / f"{inst}_bkg_rate.fits")
        write_like(counts_hdr, bkg_counts, banddir / f"{inst}_bkg_counts.fits")
        write_like(rate_hdr, rate_clean, banddir / f"{inst}_rate_clean.fits")
        write_like(rate_hdr, net_rate, banddir / f"{inst}_net_rate.fits")
        write_like(counts_hdr, net_counts, banddir / f"{inst}_net_counts.fits")
        summary["products"][inst] = {
            "trail_mask_fraction": float(np.mean(mask)),
            "finite_rate_fraction": float(np.mean(np.isfinite(rate))),
            "median_background_rate": (
                float(np.nanmedian(bkg_rate)) if np.any(np.isfinite(bkg_rate)) else None
            ),
            "median_net_rate": (
                float(np.nanmedian(net_rate)) if np.any(np.isfinite(net_rate)) else None
            ),
        }
    return summary


def cmd_qc(args: argparse.Namespace) -> int:
    stillsky_root = Path(args.stillsky_root)
    comet_root = Path(args.comet_root)
    instruments = ["PN", "M1", "M2"]
    srclist = (
        args.detect_srclist
        if args.detect_srclist and Path(args.detect_srclist).exists()
        else None
    )
    for banddir in sorted(stillsky_root.iterdir()):
        if banddir.is_dir():
            print(
                f"  stillsky/{banddir.name}: {len(process_qc_band_dir(banddir, instruments, srclist, args.mask_radius_arcsec, args.bkg_sigma_px).get('instruments', {}))} instruments processed"
            )
    for banddir in sorted(comet_root.iterdir()):
        if banddir.is_dir():
            print(
                f"  comet/{banddir.name}: {len(process_qc_band_dir(banddir, instruments, srclist, args.mask_radius_arcsec, args.bkg_sigma_px).get('instruments', {}))} instruments processed"
            )
    return 0


def cmd_science(args: argparse.Namespace) -> int:
    banddir = Path(args.banddir)
    summary = make_science_products(banddir, sigma_px=float(args.sigma_px))
    if args.out_json:
        out = Path(args.out_json)
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(
            json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
    print(banddir)
    return 0


def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser()
    sp = ap.add_subparsers(dest="cmd", required=True)
    p = sp.add_parser("science")
    p.add_argument("--banddir", required=True)
    p.add_argument("--sigma-px", type=float, default=20.0)
    p.add_argument("--out-json")
    p.set_defaults(func=cmd_science)
    p = sp.add_parser("qc")
    p.add_argument("--config")
    p.add_argument("--stillsky-root", required=True)
    p.add_argument("--comet-root", required=True)
    p.add_argument("--detect-srclist")
    p.add_argument("--grid-json")
    p.add_argument("--track")
    p.add_argument("--mask-radius-arcsec", type=float, default=20.0)
    p.add_argument("--bkg-sigma-px", type=float, default=20.0)
    p.set_defaults(func=cmd_qc)
    return ap


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    return int(args.func(args))


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(1)
