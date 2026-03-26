#!/usr/bin/env python3
from __future__ import annotations
import argparse
import json
import math
import sys
from pathlib import Path
import numpy as np
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.io import fits
from astropy.time import Time

def first_table(hdul: fits.HDUList) -> fits.BinTableHDU:
    for hdu in hdul[1:]:
        if isinstance(hdu, fits.BinTableHDU) and hdu.data is not None and len(hdu.data):
            return hdu
    raise RuntimeError('No non-empty binary table found')

def _find_col_case_insensitive(names: list[str], target: str) -> str:
    lut = {name.upper(): name for name in names}
    if target.upper() not in lut:
        raise RuntimeError(f'Missing required column {target!r}')
    return lut[target.upper()]

def _circular_mean_deg(values: np.ndarray) -> float:
    ang = np.radians(np.asarray(values, dtype=float))
    s = np.nanmean(np.sin(ang))
    c = np.nanmean(np.cos(ang))
    return float(np.degrees(np.arctan2(s, c)) % 360.0)

def _circular_median_deg(values: np.ndarray) -> float:
    vals = np.asarray(values, dtype=float) % 360.0
    if vals.size == 0:
        return math.nan
    diffs = (vals[:, None] - vals[None, :] + 180.0) % 360.0 - 180.0
    score = np.sum(np.abs(diffs), axis=1)
    return float(vals[int(np.argmin(score))])

def _linear_median(values: np.ndarray) -> float:
    vals = np.asarray(values, dtype=float)
    return float(np.nanmedian(vals)) if vals.size else math.nan

class TrackInterpolator:

    def __init__(self, t_sec: np.ndarray, ra_deg: np.ndarray, dec_deg: np.ndarray) -> None:
        sc = SkyCoord(ra=np.asarray(ra_deg, dtype=float) * u.deg, dec=np.asarray(dec_deg, dtype=float) * u.deg, frame='icrs')
        self.t_sec = np.asarray(t_sec, dtype=float)
        self.xyz = sc.cartesian.xyz.value.T.astype(float)

    def at(self, t_sec: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        ts = np.asarray(t_sec, dtype=float)
        if ts.size == 0:
            return (np.asarray([], dtype=float), np.asarray([], dtype=float))
        if ts.min() < self.t_sec.min() - 1e-06 or ts.max() > self.t_sec.max() + 1e-06:
            raise RuntimeError(f'ATTHK time grid extends outside comet track coverage: [{ts.min():.3f}, {ts.max():.3f}] vs [{self.t_sec.min():.3f}, {self.t_sec.max():.3f}]')
        x = np.interp(ts, self.t_sec, self.xyz[:, 0])
        y = np.interp(ts, self.t_sec, self.xyz[:, 1])
        z = np.interp(ts, self.t_sec, self.xyz[:, 2])
        n = np.sqrt(x * x + y * y + z * z)
        n[n == 0.0] = 1.0
        x /= n
        y /= n
        z /= n
        ra = np.degrees(np.arctan2(y, x)) % 360.0
        dec = np.degrees(np.arctan2(z, np.sqrt(x * x + y * y)))
        return (np.asarray(ra, dtype=float), np.asarray(dec, dtype=float))

def read_track(track_fits: str) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    with fits.open(track_fits) as hdul:
        hdu = first_table(hdul)
        tab = hdu.data
        hdr = hdu.header
        names = list(hdu.columns.names)
        mjd_col = _find_col_case_insensitive(names, 'MJD')
        ra_col = _find_col_case_insensitive(names, 'RA')
        dec_col = _find_col_case_insensitive(names, 'DEC')
        mjd = np.asarray(tab[mjd_col], dtype=float)
        ra = np.asarray(tab[ra_col], dtype=float)
        dec = np.asarray(tab[dec_col], dtype=float)
        mjdref = float(hdr.get('MJDREF', 0.0))
        timesys = str(hdr.get('TIMESYS', 'TT')).strip().upper()
    track_time = Time(mjd, format='mjd', scale='utc')
    ref_time = Time(mjdref, format='mjd', scale=timesys.lower())
    t_sec = np.asarray((track_time - ref_time).to_value(u.s), dtype=float)
    order = np.argsort(t_sec)
    return (t_sec[order], ra[order], dec[order])

def read_atthk(path: str) -> tuple[fits.HDUList, fits.BinTableHDU, np.ndarray, np.ndarray, np.ndarray]:
    hdul = fits.open(path)
    if 'ATTHK' not in hdul:
        hdul.close()
        raise RuntimeError(f'No ATTHK extension in {path}')
    hdu = hdul['ATTHK']
    names = list(hdu.columns.names)
    t_col = _find_col_case_insensitive(names, 'TIME')
    ra_col = _find_col_case_insensitive(names, 'AHFRA')
    dec_col = _find_col_case_insensitive(names, 'AHFDEC')
    times = np.asarray(hdu.data[t_col], dtype=float)
    ra = np.asarray(hdu.data[ra_col], dtype=float)
    dec = np.asarray(hdu.data[dec_col], dtype=float)
    return (hdul, hdu, times, ra, dec)

def update_summary_keywords(header: fits.Header, ra: np.ndarray, dec: np.ndarray) -> None:
    keymap = {'MAHFRA': _circular_median_deg(ra), 'MAHFDEC': _linear_median(dec), 'AAHFRA': _circular_mean_deg(ra), 'AAHFDEC': float(np.nanmean(dec)) if len(dec) else math.nan, 'MOMRA': _circular_median_deg(ra), 'MOMDEC': _linear_median(dec), 'AOMRA': _circular_mean_deg(ra), 'AOMDEC': float(np.nanmean(dec)) if len(dec) else math.nan}
    for key, value in keymap.items():
        if key in header and np.isfinite(value):
            header[key] = value

def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument('--input-atthk', required=True)
    ap.add_argument('--track', required=True)
    ap.add_argument('--output', required=True)
    ap.add_argument('--ref-ra', required=True, type=float)
    ap.add_argument('--ref-dec', required=True, type=float)
    ap.add_argument('--report-json')
    args = ap.parse_args()
    trk_t, trk_ra, trk_dec = read_track(args.track)
    interp = TrackInterpolator(trk_t, trk_ra, trk_dec)
    hdul, atthk, times, pnt_ra, pnt_dec = read_atthk(args.input_atthk)
    obj_ra, obj_dec = interp.at(times)
    obj = SkyCoord(ra=obj_ra * u.deg, dec=obj_dec * u.deg, frame='icrs')
    pnt = SkyCoord(ra=pnt_ra * u.deg, dec=pnt_dec * u.deg, frame='icrs')
    ref = SkyCoord(ra=float(args.ref_ra) * u.deg, dec=float(args.ref_dec) * u.deg, frame='icrs')
    dlon, dlat = obj.spherical_offsets_to(pnt)
    moved = ref.spherical_offsets_by(dlon, dlat)
    new_ra = np.asarray(moved.ra.deg, dtype=float)
    new_dec = np.asarray(moved.dec.deg, dtype=float)
    atthk.data[_find_col_case_insensitive(list(atthk.columns.names), 'AHFRA')] = new_ra
    atthk.data[_find_col_case_insensitive(list(atthk.columns.names), 'AHFDEC')] = new_dec
    update_summary_keywords(atthk.header, new_ra, new_dec)
    atthk.header['HISTORY'] = 'Comet-tracking ATTHK generated from original ATTHK + track'
    atthk.header['HISTORY'] = f'Reference attitude fixed at RA={float(args.ref_ra):.8f} DEC={float(args.ref_dec):.8f}'
    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    hdul.writeto(out_path, overwrite=True)
    hdul.close()
    if args.report_json:
        sep = SkyCoord(ra=new_ra * u.deg, dec=new_dec * u.deg).separation(SkyCoord(ra=pnt_ra * u.deg, dec=pnt_dec * u.deg)).to_value(u.arcsec)
        payload = {'input_atthk': str(Path(args.input_atthk).resolve()), 'track': str(Path(args.track).resolve()), 'output': str(out_path.resolve()), 'ref_ra_deg': float(args.ref_ra), 'ref_dec_deg': float(args.ref_dec), 'n_rows': int(len(times)), 'median_shift_arcsec': float(np.nanmedian(sep)) if len(sep) else 0.0, 'max_shift_arcsec': float(np.nanmax(sep)) if len(sep) else 0.0, 'time_min_s': float(np.nanmin(times)) if len(times) else math.nan, 'time_max_s': float(np.nanmax(times)) if len(times) else math.nan}
        report = Path(args.report_json)
        report.parent.mkdir(parents=True, exist_ok=True)
        report.write_text(json.dumps(payload, indent=2, sort_keys=True) + '\n', encoding='utf-8')
    print(out_path)
    return 0
if __name__ == '__main__':
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f'ERROR: {exc}', file=sys.stderr)
        raise SystemExit(1)
