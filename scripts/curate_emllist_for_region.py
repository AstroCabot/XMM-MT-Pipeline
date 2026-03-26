#!/usr/bin/env python3
from __future__ import annotations
import argparse
import math
import sys
from pathlib import Path
import numpy as np
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.io import fits

def first_table(hdul: fits.HDUList, *, allow_empty: bool=False) -> fits.BinTableHDU:
    empty_candidate: fits.BinTableHDU | None = None
    for hdu in hdul[1:]:
        if not isinstance(hdu, fits.BinTableHDU):
            continue
        if hdu.data is not None and len(hdu.data) > 0:
            return hdu
        if allow_empty and empty_candidate is None:
            empty_candidate = hdu
    if allow_empty and empty_candidate is not None:
        return empty_candidate
    raise RuntimeError('No binary table found in input FITS file.')

def first_nonempty_table(hdul: fits.HDUList) -> fits.BinTableHDU:
    return first_table(hdul, allow_empty=False)

def find_col(names: list[str], candidates: list[str]) -> str | None:
    lut = {name.upper(): name for name in names}
    for cand in candidates:
        if cand.upper() in lut:
            return lut[cand.upper()]
    return None

def load_track(track_fits: str) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    with fits.open(track_fits) as hdul:
        hdu = first_nonempty_table(hdul)
        names = list(hdu.columns.names)
        mjd_col = find_col(names, ['MJD'])
        ra_col = find_col(names, ['RA'])
        dec_col = find_col(names, ['DEC'])
        if not (mjd_col and ra_col and dec_col):
            raise RuntimeError('Track FITS must contain MJD, RA, DEC columns.')
        mjd = np.asarray(hdu.data[mjd_col], dtype=float)
        ra = np.asarray(hdu.data[ra_col], dtype=float)
        dec = np.asarray(hdu.data[dec_col], dtype=float)
    return (mjd, ra, dec)

def deduplicate_positions(ra_deg: np.ndarray, dec_deg: np.ndarray, tol_arcsec: float) -> np.ndarray:
    keep = np.ones(len(ra_deg), dtype=bool)
    coords = SkyCoord(ra=ra_deg * u.deg, dec=dec_deg * u.deg, frame='icrs')
    for i in range(len(coords)):
        if not keep[i]:
            continue
        if i + 1 >= len(coords):
            break
        sep = coords[i].separation(coords[i + 1:]).to_value(u.arcsec)
        dup = np.where(sep <= tol_arcsec)[0]
        if dup.size:
            keep[i + 1 + dup] = False
    return keep

def write_ds9(ds9_path: str, ra_deg: np.ndarray, dec_deg: np.ndarray, radius_arcsec: float) -> None:
    with open(ds9_path, 'w', encoding='utf-8') as fh:
        fh.write('# Region file format: DS9 version 4.1\n')
        fh.write("global color=green dashlist=8 3 width=1 font='helvetica 10 normal' select=1 highlite=1\n")
        fh.write('fk5\n')
        for ra, dec in zip(ra_deg, dec_deg):
            fh.write(f'circle({ra:.8f},{dec:.8f},{radius_arcsec:.3f}")\n')

def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument('--input', required=True, help='Input emldetect/edetect_chain source-list FITS')
    ap.add_argument('--output', required=True, help='Output curated FITS file')
    ap.add_argument('--ds9', help='Optional DS9 region file for manual inspection')
    ap.add_argument('--min-det-ml', type=float, default=10.0, help='Minimum detection likelihood')
    ap.add_argument('--exclude-track-radius-arcsec', type=float, default=0.0, help='Exclude detections within this distance of the comet track')
    ap.add_argument('--track', help='Track FITS for comet position as a function of time')
    ap.add_argument('--dedup-tol-arcsec', type=float, default=1.0, help='Deduplicate sources closer than this angular distance')
    ap.add_argument('--ds9-radius-arcsec', type=float, default=20.0, help='Circle radius for optional DS9 output')
    args = ap.parse_args()
    with fits.open(args.input) as hdul:
        src_hdu = first_table(hdul, allow_empty=True)
        tab = src_hdu.data
        names = list(src_hdu.columns.names)
        hdr = src_hdu.header.copy()
    if tab is None or len(tab) == 0:
        cols = fits.ColDefs([fits.Column(name='SRCID', format='J', array=np.array([], dtype=np.int32)), fits.Column(name='RA', format='D', unit='deg', array=np.array([], dtype=np.float64)), fits.Column(name='DEC', format='D', unit='deg', array=np.array([], dtype=np.float64)), fits.Column(name='DET_ML', format='E', array=np.array([], dtype=np.float32)), fits.Column(name='EXTENT', format='E', unit='arcsec', array=np.array([], dtype=np.float32))])
        out_hdu = fits.BinTableHDU.from_columns(cols, name='SRCLIST')
        out_hdu.header['HISTORY'] = 'Input source list was empty'
        fits.HDUList([fits.PrimaryHDU(), out_hdu]).writeto(args.output, overwrite=True)
        if args.ds9:
            Path(args.ds9).parent.mkdir(parents=True, exist_ok=True)
            Path(args.ds9).write_text('# Region file format: DS9 version 4.1\nfk5\n', encoding='utf-8')
        print('Input source list is empty; wrote empty curated source list')
        return 0
    ra_col = find_col(names, ['RA'])
    dec_col = find_col(names, ['DEC'])
    detml_col = find_col(names, ['DET_ML', 'ML', 'ML_ID_SRC', 'ML_ID_BND'])
    inst_col = find_col(names, ['ID_INST'])
    band_col = find_col(names, ['ID_BAND'])
    extent_col = find_col(names, ['EXTENT'])
    srcid_col = find_col(names, ['SRCID', 'SRC_NUM', 'ML_ID_SRC'])
    if not (ra_col and dec_col):
        raise RuntimeError('Input source list does not contain RA/DEC columns.')
    mask = np.ones(len(tab), dtype=bool)
    if inst_col:
        mask &= np.asarray(tab[inst_col], dtype=int) == 0
    if band_col:
        mask &= np.asarray(tab[band_col], dtype=int) == 0
    if detml_col:
        mask &= np.asarray(tab[detml_col], dtype=float) >= args.min_det_ml
    rows = tab[mask]
    if len(rows) == 0:
        cols = fits.ColDefs([fits.Column(name='SRCID', format='J', array=np.array([], dtype=np.int32)), fits.Column(name='RA', format='D', unit='deg', array=np.array([], dtype=np.float64)), fits.Column(name='DEC', format='D', unit='deg', array=np.array([], dtype=np.float64)), fits.Column(name='DET_ML', format='E', array=np.array([], dtype=np.float32)), fits.Column(name='EXTENT', format='E', unit='arcsec', array=np.array([], dtype=np.float32))])
        out_hdu = fits.BinTableHDU.from_columns(cols, name='SRCLIST')
        out_hdu.header['HISTORY'] = 'No curated field sources remained after filtering'
        fits.HDUList([fits.PrimaryHDU(), out_hdu]).writeto(args.output, overwrite=True)
        if args.ds9:
            Path(args.ds9).write_text('# Region file format: DS9 version 4.1\nfk5\n', encoding='utf-8')
        print('Kept 0 sources')
        return 0
    ra = np.asarray(rows[ra_col], dtype=float)
    dec = np.asarray(rows[dec_col], dtype=float)
    detml = np.asarray(rows[detml_col], dtype=float) if detml_col else np.full(len(rows), np.nan)
    extent = np.asarray(rows[extent_col], dtype=float) if extent_col else np.zeros(len(rows), dtype=float)
    srcid = np.asarray(rows[srcid_col], dtype=int) if srcid_col else np.arange(1, len(rows) + 1, dtype=int)
    order = np.argsort(np.nan_to_num(detml, nan=-np.inf))[::-1]
    ra = ra[order]
    dec = dec[order]
    detml = detml[order]
    extent = extent[order]
    srcid = srcid[order]
    if args.track and args.exclude_track_radius_arcsec > 0:
        _mjd_t, ra_t, dec_t = load_track(args.track)
        track = SkyCoord(ra=ra_t * u.deg, dec=dec_t * u.deg, frame='icrs')
        src = SkyCoord(ra=ra * u.deg, dec=dec * u.deg, frame='icrs')
        min_sep = np.empty(len(src), dtype=float)
        for i, one_src in enumerate(src):
            min_sep[i] = float(np.min(one_src.separation(track).to_value(u.arcsec)))
        keep = min_sep >= args.exclude_track_radius_arcsec
        ra = ra[keep]
        dec = dec[keep]
        detml = detml[keep]
        extent = extent[keep]
        srcid = srcid[keep]
    if len(ra) == 0:
        keep_mask = np.zeros(0, dtype=bool)
    else:
        keep_mask = deduplicate_positions(ra, dec, tol_arcsec=args.dedup_tol_arcsec)
        ra = ra[keep_mask]
        dec = dec[keep_mask]
        detml = detml[keep_mask]
        extent = extent[keep_mask]
        srcid = srcid[keep_mask]
    cols = fits.ColDefs([fits.Column(name='SRCID', format='J', array=np.asarray(srcid, dtype=np.int32)), fits.Column(name='RA', format='D', unit='deg', array=np.asarray(ra, dtype=np.float64)), fits.Column(name='DEC', format='D', unit='deg', array=np.asarray(dec, dtype=np.float64)), fits.Column(name='DET_ML', format='E', array=np.asarray(detml, dtype=np.float32)), fits.Column(name='EXTENT', format='E', unit='arcsec', array=np.asarray(extent, dtype=np.float32))])
    out_hdu = fits.BinTableHDU.from_columns(cols, name='SRCLIST')
    out_hdu.header['HISTORY'] = 'Curated from emldetect/edetect_chain source list'
    out_hdu.header['HISTORY'] = f'Minimum DET_ML = {args.min_det_ml}'
    if inst_col:
        out_hdu.header['HISTORY'] = 'Kept summary rows with ID_INST==0'
    if band_col:
        out_hdu.header['HISTORY'] = 'Kept summary rows with ID_BAND==0'
    if args.track and args.exclude_track_radius_arcsec > 0:
        out_hdu.header['HISTORY'] = f'Excluded sources within {args.exclude_track_radius_arcsec:.3f} arcsec of comet track'
    out_hdu.header['HISTORY'] = f'Deduplicated positions within {args.dedup_tol_arcsec:.3f} arcsec'
    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fits.HDUList([fits.PrimaryHDU(), out_hdu]).writeto(out_path, overwrite=True)
    if args.ds9:
        Path(args.ds9).parent.mkdir(parents=True, exist_ok=True)
        write_ds9(args.ds9, ra, dec, args.ds9_radius_arcsec)
    print(f'Input rows: {len(tab)}')
    print(f'Curated sources kept: {len(ra)}')
    if len(ra):
        print(f'DET_ML range: {np.nanmin(detml):.3f} .. {np.nanmax(detml):.3f}')
    return 0
if __name__ == '__main__':
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f'ERROR: {exc}', file=sys.stderr)
        raise SystemExit(1)
