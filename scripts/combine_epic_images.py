#!/usr/bin/env python3
from __future__ import annotations
import argparse
import json
import sys
from pathlib import Path
import numpy as np
from astropy.io import fits
from trail_mask_tools import build_trail_mask, read_srclist

def first_image_hdu(hdul: fits.HDUList) -> tuple[int, fits.ImageHDU | fits.PrimaryHDU]:
    for idx, hdu in enumerate(hdul):
        if hdu.data is not None and getattr(hdu.data, 'ndim', 0) == 2:
            return (idx, hdu)
    raise RuntimeError('No 2D image HDU found')

def read_image(path: str) -> tuple[np.ndarray, fits.Header]:
    with fits.open(path) as hdul:
        _, hdu = first_image_hdu(hdul)
        data = np.asarray(hdu.data, dtype=float)
        header = hdu.header.copy()
    return (data, header)

def write_like(template_header: fits.Header, data: np.ndarray, out_path: str) -> None:
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    fits.PrimaryHDU(data=data, header=template_header).writeto(out_path, overwrite=True)

def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser()
    ap.add_argument('--counts', nargs='+', required=True, help='Per-instrument count images')
    ap.add_argument('--exposure', nargs='+', required=True, help='Per-instrument exposure maps')
    ap.add_argument('--out-counts', required=True)
    ap.add_argument('--out-exposure', required=True)
    ap.add_argument('--out-rate', required=True)
    ap.add_argument('--grid-json')
    ap.add_argument('--track')
    ap.add_argument('--srclist')
    ap.add_argument('--mask-radius-arcsec', type=float, default=0.0)
    ap.add_argument('--out-mask')
    ap.add_argument('--out-rate-masked')
    ap.add_argument('--track-tmin-sec', type=float)
    ap.add_argument('--track-tmax-sec', type=float)
    ap.add_argument('--out-json')
    return ap.parse_args()

def main() -> int:
    args = parse_args()
    if len(args.counts) != len(args.exposure):
        raise RuntimeError('--counts and --exposure must have the same length')
    counts_sum = None
    exp_sum = None
    template_header = None
    for cpath, epath in zip(args.counts, args.exposure):
        cdat, chdr = read_image(cpath)
        edat, _ehdr = read_image(epath)
        if cdat.shape != edat.shape:
            raise RuntimeError(f'Shape mismatch between {cpath} and {epath}')
        if counts_sum is None:
            counts_sum = np.zeros_like(cdat, dtype=float)
            exp_sum = np.zeros_like(edat, dtype=float)
            template_header = chdr.copy()
        if cdat.shape != counts_sum.shape:
            raise RuntimeError('All images must share the same shape')
        counts_sum += cdat
        exp_sum += edat
    if counts_sum is None or exp_sum is None or template_header is None:
        raise RuntimeError('No images supplied')
    with np.errstate(divide='ignore', invalid='ignore'):
        rate = np.where(exp_sum > 0.0, counts_sum / exp_sum, np.nan)
    write_like(template_header, counts_sum, args.out_counts)
    exp_hdr = template_header.copy()
    exp_hdr.setdefault('BUNIT', 's')
    write_like(exp_hdr, exp_sum, args.out_exposure)
    rate_hdr = template_header.copy()
    rate_hdr['BUNIT'] = 'count / s'
    write_like(rate_hdr, rate, args.out_rate)
    out_summary = {'counts_inputs': [str(Path(p).resolve()) for p in args.counts], 'exposure_inputs': [str(Path(p).resolve()) for p in args.exposure], 'shape': list(counts_sum.shape), 'finite_exposure_fraction': float(np.mean(np.isfinite(exp_sum) & (exp_sum > 0.0)))}
    need_mask = bool(args.grid_json and args.track and args.srclist and args.out_mask and args.out_rate_masked and (args.mask_radius_arcsec > 0.0))
    if need_mask:
        src_rows = 0
        try:
            src_rows = len(read_srclist(args.srclist).ra)
        except Exception:
            src_rows = 0
        if src_rows > 0:
            grid = json.loads(Path(args.grid_json).read_text(encoding='utf-8'))
            mask = build_trail_mask(grid=grid, track_fits=args.track, srclist_fits=args.srclist, mask_radius_arcsec=float(args.mask_radius_arcsec), shape=rate.shape, t_min_sec=args.track_tmin_sec, t_max_sec=args.track_tmax_sec)
        else:
            mask = np.zeros_like(rate, dtype=bool)
        masked_rate = np.array(rate, copy=True)
        masked_rate[mask] = np.nan
        mask_hdr = template_header.copy()
        mask_hdr['BUNIT'] = '1'
        mask_hdr['HISTORY'] = 'Stationary-source trail mask in the comet frame'
        if args.track_tmin_sec is not None:
            mask_hdr['TMINSEC'] = float(args.track_tmin_sec)
        if args.track_tmax_sec is not None:
            mask_hdr['TMAXSEC'] = float(args.track_tmax_sec)
        write_like(mask_hdr, mask.astype(np.int16), args.out_mask)
        masked_hdr = rate_hdr.copy()
        masked_hdr['HISTORY'] = 'Pixels crossed by stationary-source trails set to NaN'
        write_like(masked_hdr, masked_rate, args.out_rate_masked)
        out_summary['trail_mask_fraction'] = float(np.mean(mask))
    else:
        out_summary['trail_mask_fraction'] = 0.0
    if args.out_json:
        out_path = Path(args.out_json)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        out_path.write_text(json.dumps(out_summary, indent=2, sort_keys=True) + '\n', encoding='utf-8')
    return 0
if __name__ == '__main__':
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f'ERROR: {exc}', file=sys.stderr)
        raise SystemExit(1)
