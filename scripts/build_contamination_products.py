#!/usr/bin/env python3
from __future__ import annotations
import argparse
import csv
import json
import math
import shutil
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable
import numpy as np
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.io import fits
from astropy.time import Time

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
    raise RuntimeError('No binary table found')

def first_nonempty_table(hdul: fits.HDUList) -> fits.BinTableHDU:
    return first_table(hdul, allow_empty=False)

@dataclass(frozen=True)
class TrackData:
    t_sec: np.ndarray
    mjd_utc: np.ndarray
    ra: np.ndarray
    dec: np.ndarray
    obs_t0: float
    obs_t1: float
    mjdref: float
    timesys: str

class TrackInterpolator:

    def __init__(self, t_sec: np.ndarray, ra_deg: np.ndarray, dec_deg: np.ndarray) -> None:
        sc = SkyCoord(ra=np.asarray(ra_deg, dtype=float) * u.deg, dec=np.asarray(dec_deg, dtype=float) * u.deg, frame='icrs')
        self.t_sec = np.asarray(t_sec, dtype=float)
        self.xyz = sc.cartesian.xyz.value.T.astype(float)

    def at(self, t_sec: np.ndarray | float) -> tuple[np.ndarray, np.ndarray]:
        ts = np.atleast_1d(np.asarray(t_sec, dtype=float))
        if ts.size == 0:
            return (np.asarray([], dtype=float), np.asarray([], dtype=float))
        if ts.min() < self.t_sec.min() - 1e-09 or ts.max() > self.t_sec.max() + 1e-09:
            raise RuntimeError('Requested interpolation time is outside the track range')
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
        if np.isscalar(t_sec):
            return (np.asarray([ra[0]], dtype=float), np.asarray([dec[0]], dtype=float))
        return (np.asarray(ra, dtype=float), np.asarray(dec, dtype=float))

def read_track(track_fits: str) -> TrackData:
    with fits.open(track_fits) as hdul:
        hdu = first_nonempty_table(hdul)
        tab = hdu.data
        hdr = hdu.header
        mjd_utc = np.asarray(tab['MJD'], dtype=float)
        ra = np.asarray(tab['RA'], dtype=float)
        dec = np.asarray(tab['DEC'], dtype=float)
        mjdref = float(hdr.get('MJDREF', 0.0))
        timesys = str(hdr.get('TIMESYS', 'TT')).strip().upper()
        obs_t0 = float(hdr.get('OBS_T0', np.nan)) if 'OBS_T0' in hdr else math.nan
        obs_t1 = float(hdr.get('OBS_T1', np.nan)) if 'OBS_T1' in hdr else math.nan
    track_time = Time(mjd_utc, format='mjd', scale='utc')
    ref_time = Time(mjdref, format='mjd', scale=timesys.lower())
    t_sec = np.asarray((track_time - ref_time).to_value(u.s), dtype=float)
    if not np.isfinite(obs_t0):
        obs_t0 = float(t_sec[0])
    if not np.isfinite(obs_t1):
        obs_t1 = float(t_sec[-1])
    return TrackData(t_sec=t_sec, mjd_utc=mjd_utc, ra=ra, dec=dec, obs_t0=obs_t0, obs_t1=obs_t1, mjdref=mjdref, timesys=timesys)

def read_srclist(src_fits: str) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    with fits.open(src_fits) as hdul:
        hdu = first_table(hdul, allow_empty=True)
        names = [n.upper() for n in hdu.columns.names]
        if 'RA' not in names or 'DEC' not in names:
            raise RuntimeError('Curated source list must contain RA and DEC columns')
        data = hdu.data
        if data is None or len(data) == 0:
            return (np.asarray([], dtype=int), np.asarray([], dtype=float), np.asarray([], dtype=float))
        ra = np.asarray(data[hdu.columns.names[names.index('RA')]], dtype=float)
        dec = np.asarray(data[hdu.columns.names[names.index('DEC')]], dtype=float)
        if 'SRCID' in names:
            srcid = np.asarray(data[hdu.columns.names[names.index('SRCID')]], dtype=int)
        else:
            srcid = np.arange(1, len(ra) + 1, dtype=int)
    return (srcid, ra, dec)

def arcsec_offsets(comet_ra: np.ndarray, comet_dec: np.ndarray, src_ra: float, src_dec: float) -> tuple[np.ndarray, np.ndarray]:
    comet = SkyCoord(ra=np.asarray(comet_ra, dtype=float) * u.deg, dec=np.asarray(comet_dec, dtype=float) * u.deg, frame='icrs')
    src = SkyCoord(ra=float(src_ra) * u.deg, dec=float(src_dec) * u.deg, frame='icrs')
    dlon, dlat = comet.spherical_offsets_to(src)
    return (dlon.to_value(u.arcsec), dlat.to_value(u.arcsec))

@dataclass(frozen=True)
class CircleAperture:
    dx: float
    dy: float
    r: float

@dataclass(frozen=True)
class AnnulusAperture:
    dx: float
    dy: float
    rin: float
    rout: float
Aperture = CircleAperture | AnnulusAperture

@dataclass
class SourceReport:
    srcid: int
    min_sep_src_arcsec: float
    min_sep_bkg_arcsec: float
    samples_src_overlap: int
    samples_bkg_overlap: int
    samples_strict_overlap: int

def overlaps_circle(src_dx: np.ndarray, src_dy: np.ndarray, ap: CircleAperture, mask_r: float) -> np.ndarray:
    dist = np.hypot(src_dx - ap.dx, src_dy - ap.dy)
    return dist <= ap.r + mask_r

def overlaps_annulus(src_dx: np.ndarray, src_dy: np.ndarray, ap: AnnulusAperture, mask_r: float) -> np.ndarray:
    dist = np.hypot(src_dx - ap.dx, src_dy - ap.dy)
    return (dist + mask_r >= ap.rin) & (dist - mask_r <= ap.rout)

def overlaps(src_dx: np.ndarray, src_dy: np.ndarray, ap: Aperture, mask_r: float) -> np.ndarray:
    if isinstance(ap, CircleAperture):
        return overlaps_circle(src_dx, src_dy, ap, mask_r)
    return overlaps_annulus(src_dx, src_dy, ap, mask_r)

def good_runs_to_intervals(times: np.ndarray, good: np.ndarray, obs_t0: float, obs_t1: float, min_good_span_s: float) -> list[tuple[float, float]]:
    if len(times) == 0:
        return []
    if len(times) == 1:
        return [(float(obs_t0), float(obs_t1))] if bool(good[0]) else []
    starts: list[float] = []
    stops: list[float] = []
    inside = False
    cur_start = 0.0
    for i, is_good in enumerate(good):
        left = float(obs_t0) if i == 0 else float(0.5 * (times[i - 1] + times[i]))
        right = float(obs_t1) if i == len(times) - 1 else float(0.5 * (times[i] + times[i + 1]))
        if is_good and (not inside):
            inside = True
            cur_start = left
        if inside and (not is_good):
            inside = False
            starts.append(cur_start)
            stops.append(left)
        if inside and i == len(times) - 1:
            starts.append(cur_start)
            stops.append(right)
    intervals = [(float(a), float(b)) for a, b in zip(starts, stops) if b > a]
    if min_good_span_s > 0.0:
        intervals = [(a, b) for a, b in intervals if b - a >= float(min_good_span_s)]
    return intervals
POLICY_LABELS = {'full': 'No contamination GTI filtering; report-only timeline', 'src': 'Reject source-aperture contamination only', 'bkg': 'Reject background-aperture contamination only', 'strict': 'Reject either source- or background-aperture contamination'}

def write_gti(out_fits: str, intervals: Iterable[tuple[float, float]], mjdref: float, timesys: str, notes: list[str]) -> None:
    start = np.asarray([a for a, _ in intervals], dtype=np.float64)
    stop = np.asarray([b for _, b in intervals], dtype=np.float64)
    cols = fits.ColDefs([fits.Column(name='START', format='D', unit='s', array=start), fits.Column(name='STOP', format='D', unit='s', array=stop)])
    hdu = fits.BinTableHDU.from_columns(cols, name='GTI')
    hdu.header['MJDREF'] = (float(mjdref), 'Reference MJD for XMM TIME')
    hdu.header['TIMESYS'] = (timesys, 'Time scale')
    for note in notes:
        hdu.header['HISTORY'] = note[:72]
    fits.HDUList([fits.PrimaryHDU(), hdu]).writeto(out_fits, overwrite=True)

def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser()
    ap.add_argument('--track', required=True)
    ap.add_argument('--srclist', required=True)
    ap.add_argument('--output-dir')
    ap.add_argument('--output-gti')
    ap.add_argument('--report-csv')
    ap.add_argument('--summary-json')
    ap.add_argument('--science-policy', default='full', choices=['full', 'src', 'bkg', 'strict'])
    ap.add_argument('--sample-step-s', type=float, default=5.0)
    ap.add_argument('--mask-radius-arcsec', type=float, required=True)
    ap.add_argument('--min-good-span-s', type=float, default=0.0)
    ap.add_argument('--src-dx-arcsec', type=float, required=True)
    ap.add_argument('--src-dy-arcsec', type=float, required=True)
    ap.add_argument('--src-r-arcsec', type=float, required=True)
    ap.add_argument('--bkg-shape', choices=['circle', 'annulus'], default='annulus')
    ap.add_argument('--bkg-dx-arcsec', type=float, required=True)
    ap.add_argument('--bkg-dy-arcsec', type=float, required=True)
    ap.add_argument('--bkg-r-arcsec', type=float)
    ap.add_argument('--bkg-rin-arcsec', type=float)
    ap.add_argument('--bkg-rout-arcsec', type=float)
    return ap.parse_args()

def main() -> int:
    args = parse_args()
    if args.sample_step_s <= 0:
        raise RuntimeError('--sample-step-s must be > 0')
    if args.mask_radius_arcsec < 0:
        raise RuntimeError('--mask-radius-arcsec must be >= 0')
    outdir = Path(args.output_dir) if args.output_dir else None
    if outdir is None and args.output_gti:
        outdir = Path(args.output_gti).parent
    if outdir is None and args.report_csv:
        outdir = Path(args.report_csv).parent
    if outdir is None:
        raise RuntimeError('Supply --output-dir or legacy --output-gti/--report-csv')
    outdir.mkdir(parents=True, exist_ok=True)
    track = read_track(args.track)
    interp = TrackInterpolator(track.t_sec, track.ra, track.dec)
    srcid, src_ra, src_dec = read_srclist(args.srclist)
    if args.bkg_shape == 'circle':
        if args.bkg_r_arcsec is None:
            raise RuntimeError('bkg-shape=circle requires --bkg-r-arcsec')
        bkg_ap: Aperture = CircleAperture(args.bkg_dx_arcsec, args.bkg_dy_arcsec, args.bkg_r_arcsec)
    else:
        if args.bkg_rin_arcsec is None or args.bkg_rout_arcsec is None:
            raise RuntimeError('bkg-shape=annulus requires --bkg-rin-arcsec and --bkg-rout-arcsec')
        bkg_ap = AnnulusAperture(args.bkg_dx_arcsec, args.bkg_dy_arcsec, args.bkg_rin_arcsec, args.bkg_rout_arcsec)
    src_ap = CircleAperture(args.src_dx_arcsec, args.src_dy_arcsec, args.src_r_arcsec)
    times = np.arange(track.obs_t0, track.obs_t1 + 0.5 * args.sample_step_s, args.sample_step_s, dtype=float)
    times = times[times <= track.obs_t1 + 1e-06]
    comet_ra, comet_dec = interp.at(times)
    mjd_utc = (Time(track.mjdref, format='mjd', scale=track.timesys.lower()) + times * u.s).utc.mjd
    n_src_overlap = np.zeros(len(times), dtype=int)
    n_bkg_overlap = np.zeros(len(times), dtype=int)
    src_ids_by_sample = [list() for _ in range(len(times))]
    bkg_ids_by_sample = [list() for _ in range(len(times))]
    reports: list[SourceReport] = []
    for sid, ra, dec in zip(srcid, src_ra, src_dec):
        dx, dy = arcsec_offsets(comet_ra, comet_dec, float(ra), float(dec))
        bad_src = overlaps(dx, dy, src_ap, float(args.mask_radius_arcsec))
        bad_bkg = overlaps(dx, dy, bkg_ap, float(args.mask_radius_arcsec))
        n_src_overlap += bad_src.astype(int)
        n_bkg_overlap += bad_bkg.astype(int)
        for idx in np.flatnonzero(bad_src):
            src_ids_by_sample[int(idx)].append(int(sid))
        for idx in np.flatnonzero(bad_bkg):
            bkg_ids_by_sample[int(idx)].append(int(sid))
        dist_src = np.hypot(dx - src_ap.dx, dy - src_ap.dy)
        dist_bkg = np.hypot(dx - getattr(bkg_ap, 'dx'), dy - getattr(bkg_ap, 'dy'))
        reports.append(SourceReport(srcid=int(sid), min_sep_src_arcsec=float(np.min(dist_src)) if len(dist_src) else math.nan, min_sep_bkg_arcsec=float(np.min(dist_bkg)) if len(dist_bkg) else math.nan, samples_src_overlap=int(np.count_nonzero(bad_src)), samples_bkg_overlap=int(np.count_nonzero(bad_bkg)), samples_strict_overlap=int(np.count_nonzero(bad_src | bad_bkg))))
    good_masks = {'full': np.ones(len(times), dtype=bool), 'src': n_src_overlap == 0, 'bkg': n_bkg_overlap == 0, 'strict': (n_src_overlap == 0) & (n_bkg_overlap == 0)}
    gti_paths: dict[str, str] = {}
    good_time: dict[str, float] = {}
    n_intervals: dict[str, int] = {}
    for policy, good in good_masks.items():
        intervals = good_runs_to_intervals(times, good, track.obs_t0, track.obs_t1, float(args.min_good_span_s))
        path = outdir / f'gti_{policy}.fits'
        write_gti(out_fits=str(path), intervals=intervals, mjdref=track.mjdref, timesys=track.timesys, notes=[POLICY_LABELS[policy], f'Sampling step = {args.sample_step_s:.6f} s', f'Mask radius = {args.mask_radius_arcsec:.3f} arcsec'])
        gti_paths[policy] = str(path)
        good_time[policy] = float(sum((b - a for a, b in intervals)))
        n_intervals[policy] = int(len(intervals))
    science_gti = outdir / 'science_gti.fits'
    shutil.copy2(gti_paths[args.science_policy], science_gti)
    timeline_csv = outdir / 'contamination_timeline.csv'
    with timeline_csv.open('w', newline='', encoding='utf-8') as fh:
        writer = csv.writer(fh)
        writer.writerow(['time_s', 'mjd_utc', 'good_full', 'good_src', 'good_bkg', 'good_strict', 'n_src_overlap', 'n_bkg_overlap', 'src_ids', 'bkg_ids'])
        for t, mjd, g_full, g_src, g_bkg, g_strict, ns, nb, sids, bids in zip(times, mjd_utc, good_masks['full'], good_masks['src'], good_masks['bkg'], good_masks['strict'], n_src_overlap, n_bkg_overlap, src_ids_by_sample, bkg_ids_by_sample):
            writer.writerow([f'{float(t):.6f}', f'{float(mjd):.10f}', int(g_full), int(g_src), int(g_bkg), int(g_strict), int(ns), int(nb), ';'.join((str(v) for v in sids)), ';'.join((str(v) for v in bids))])
    report_csv = Path(args.report_csv) if args.report_csv else outdir / 'contamination_report.csv'
    with report_csv.open('w', newline='', encoding='utf-8') as fh:
        writer = csv.writer(fh)
        writer.writerow(['srcid', 'min_sep_src_arcsec', 'min_sep_bkg_arcsec', 'samples_src_overlap', 'samples_bkg_overlap', 'samples_strict_overlap'])
        for rep in sorted(reports, key=lambda r: (r.samples_strict_overlap, -r.min_sep_src_arcsec), reverse=True):
            writer.writerow([rep.srcid, f'{rep.min_sep_src_arcsec:.6f}', f'{rep.min_sep_bkg_arcsec:.6f}', rep.samples_src_overlap, rep.samples_bkg_overlap, rep.samples_strict_overlap])
    if args.output_gti:
        shutil.copy2(science_gti, args.output_gti)
    summary = {'track': str(Path(args.track).resolve()), 'srclist': str(Path(args.srclist).resolve()), 'science_policy': args.science_policy, 'sample_step_s': float(args.sample_step_s), 'mask_radius_arcsec': float(args.mask_radius_arcsec), 'observation_span_s': float(track.obs_t1 - track.obs_t0), 'n_sources': int(len(srcid)), 'good_time_s': good_time, 'n_intervals': n_intervals, 'gti_paths': gti_paths, 'science_gti': str(science_gti), 'timeline_csv': str(timeline_csv), 'report_csv': str(report_csv)}
    summary_path = Path(args.summary_json) if args.summary_json else outdir / 'contamination_summary.json'
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + '\n', encoding='utf-8')
    print(f'Sources considered: {len(srcid)}')
    print(f'Observation span (s): {track.obs_t1 - track.obs_t0:.3f}')
    for policy in ['full', 'src', 'bkg', 'strict']:
        print(f'Policy {policy:6s}: good exposure {good_time[policy]:.3f} s in {n_intervals[policy]} GTI intervals')
    print(f'Selected science policy: {args.science_policy}')
    print(f'science_gti.fits -> {science_gti}')
    return 0
if __name__ == '__main__':
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f'ERROR: {exc}', file=sys.stderr)
        raise SystemExit(1)
