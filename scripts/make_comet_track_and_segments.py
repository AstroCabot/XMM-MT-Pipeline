#!/usr/bin/env python3
"""Build the comet ephemeris track and motion-segmentation table.

Queries JPL Horizons (via astroquery) or reads a user-supplied
ephemeris, then splits the observation into segments whose midpoint
blur stays within a configurable threshold.
"""
from __future__ import annotations
import argparse
import csv
import math
import re
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable
import numpy as np
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.io import fits
from astropy.time import Time, TimeDelta

try:
    from astroquery.jplhorizons import Horizons
except Exception:
    Horizons = None


def read_event_time_info(event_file: str) -> tuple[float, float, float, str]:
    with fits.open(event_file) as hdul:
        evt = hdul["EVENTS"].data
        hdr = hdul["EVENTS"].header
        if len(evt) == 0:
            raise RuntimeError(f"EVENTS table is empty in {event_file}")
        tmin = float(np.min(evt["TIME"]))
        tmax = float(np.max(evt["TIME"]))
        mjdref = hdr.get("MJDREF")
        if mjdref is None:
            mjdref = float(hdr.get("MJDREFI", 0)) + float(hdr.get("MJDREFF", 0.0))
        timesys = str(hdr.get("TIMESYS", "TT")).strip().upper()
    return (tmin, tmax, float(mjdref), timesys)


def xmm_seconds_to_time(
    seconds: np.ndarray | float, mjdref: float, timesys: str
) -> Time:
    return Time(mjdref, format="mjd", scale=timesys.lower()) + TimeDelta(
        seconds, format="sec"
    )


def xmm_seconds_to_utc_iso(seconds: float, mjdref: float, timesys: str) -> str:
    return xmm_seconds_to_time(seconds, mjdref, timesys).utc.isot


def _find_col_case_insensitive(
    names: Iterable[str], candidates: Iterable[str]
) -> str | None:
    lut = {n.lower(): n for n in names}
    for cand in candidates:
        if cand.lower() in lut:
            return lut[cand.lower()]
    return None


def load_track_from_csv(
    csv_file: str,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    data = np.genfromtxt(csv_file, delimiter=",", names=True, dtype=None, encoding=None)
    if data.size == 0:
        raise RuntimeError(f"Track CSV is empty: {csv_file}")
    names = data.dtype.names or ()
    tcol = _find_col_case_insensitive(names, ["time_iso", "utc_iso", "isot", "time"])
    racol = _find_col_case_insensitive(names, ["ra_deg", "ra"])
    deccol = _find_col_case_insensitive(names, ["dec_deg", "dec"])
    dcol = _find_col_case_insensitive(names, ["delta_au", "delta"])
    if not (tcol and racol and deccol and dcol):
        raise RuntimeError(
            "Track CSV must contain columns time_iso, ra_deg, dec_deg, delta_au (case-insensitive)."
        )
    tt = Time(np.array(data[tcol]), format="isot", scale="utc")
    return (
        tt.mjd.astype(float),
        np.array(data[racol], float),
        np.array(data[deccol], float),
        np.array(data[dcol], float),
    )


def load_track_from_fits(
    fits_file: str,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    with fits.open(fits_file) as hdul:
        tab = None
        for hdu in hdul[1:]:
            if isinstance(hdu, fits.BinTableHDU) and len(hdu.data) > 0:
                tab = hdu.data
                break
        if tab is None:
            raise RuntimeError(f"No non-empty table found in {fits_file}")
        names = tab.columns.names
        mjdcol = _find_col_case_insensitive(names, ["MJD"])
        racol = _find_col_case_insensitive(names, ["RA"])
        deccol = _find_col_case_insensitive(names, ["DEC"])
        dcol = _find_col_case_insensitive(names, ["DELTA"])
        if not (mjdcol and racol and deccol and dcol):
            raise RuntimeError("Track FITS must contain columns MJD, RA, DEC, DELTA.")
        return (
            np.array(tab[mjdcol], float),
            np.array(tab[racol], float),
            np.array(tab[deccol], float),
            np.array(tab[dcol], float),
        )


_STEP_RE = re.compile("^\\s*([0-9]*\\.?[0-9]+)\\s*([A-Za-z]+)?\\s*$")


def _parse_fixed_step_seconds(step: str) -> tuple[float | None, str | None]:
    s = step.strip()
    m = _STEP_RE.match(s)
    if not m:
        return (None, None)
    value = float(m.group(1))
    unit = (m.group(2) or "").strip().lower()
    if unit == "":
        return (None, None)
    if unit in {"s", "sec", "secs", "second", "seconds"}:
        return (value, "s")
    if unit in {"m", "min", "mins", "minute", "minutes"}:
        return (value * 60.0, "m")
    if unit in {"h", "hr", "hrs", "hour", "hours"}:
        return (value * 3600.0, "h")
    if unit in {"d", "day", "days"}:
        return (value * 86400.0, "d")
    if unit in {"mo", "mon", "mons", "month", "months"}:
        return (None, "mo")
    if unit in {"y", "yr", "yrs", "year", "years"}:
        return (None, "y")
    return (None, None)


def _build_epoch_list(start_utc: str, stop_utc: str, step_seconds: float) -> np.ndarray:
    if step_seconds <= 0.0:
        raise RuntimeError("Requested Horizons cadence must be > 0 seconds.")
    t0 = Time(start_utc, format="isot", scale="utc")
    t1 = Time(stop_utc, format="isot", scale="utc")
    total = max(0.0, (t1 - t0).to_value(u.s))
    offsets = np.arange(0.0, total + 0.5 * step_seconds, step_seconds, dtype=float)
    times = t0 + TimeDelta(offsets, format="sec")
    if times[-1] < t1:
        times = Time(
            np.concatenate([times.jd, np.array([t1.jd])]), format="jd", scale="utc"
        )
    jd = np.unique(
        np.concatenate(
            [np.array([t0.jd]), np.array(times.jd, float), np.array([t1.jd])]
        )
    )
    return jd.astype(float)


def _parse_eph_table(eph) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    cols = eph.colnames
    racol = _find_col_case_insensitive(cols, ["RA", "RA_ICRF", "RA_app", "RA_ICRF_app"])
    deccol = _find_col_case_insensitive(cols, ["DEC", "DEC_ICRF", "DEC_app", "DEC_ICRF_app"])
    dcol = _find_col_case_insensitive(cols, ["delta"])
    tcol = _find_col_case_insensitive(cols, ["datetime_jd", "JD", "jd"])
    if not (racol and deccol and dcol and tcol):
        raise RuntimeError(
            f"Horizons output missing one of RA/DEC/delta/datetime_jd columns; got {cols}"
        )
    jd = np.array(eph[tcol], float)
    mjd = jd - 2400000.5
    ra = np.array(eph[racol], float)
    dec = np.array(eph[deccol], float)
    delta = np.array(eph[dcol], float)
    return (mjd, ra, dec, delta)


def _query_horizons_epoch_chunk(
    target_id: str, id_type: str, epoch_chunk_jd: np.ndarray
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    obj = Horizons(
        id=target_id,
        id_type=id_type,
        location="500@399",
        epochs=[float(x) for x in epoch_chunk_jd],
    )
    eph = obj.ephemerides(quantities="1,20", extra_precision=True)
    return _parse_eph_table(eph)


def query_horizons_track(
    target_id: str,
    id_type: str,
    start_utc: str,
    stop_utc: str,
    step: str,
    epoch_chunk_size: int = 40,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    if Horizons is None:
        raise RuntimeError(
            "astroquery is unavailable, but no external track file was supplied."
        )
    step_seconds, canonical_unit = _parse_fixed_step_seconds(step)
    use_explicit_epochs = bool(step_seconds is not None and step_seconds < 60.0)
    if step_seconds is not None and canonical_unit in {"m", "h", "d"}:
        base = {"m": 60.0, "h": 3600.0, "d": 86400.0}[canonical_unit]
        use_explicit_epochs = (
            use_explicit_epochs
            or abs(step_seconds / base - round(step_seconds / base)) > 1e-12
        )
    if use_explicit_epochs:
        epochs_jd = _build_epoch_list(
            start_utc=start_utc, stop_utc=stop_utc, step_seconds=float(step_seconds)
        )
        mjd_parts: list[np.ndarray] = []
        ra_parts: list[np.ndarray] = []
        dec_parts: list[np.ndarray] = []
        delta_parts: list[np.ndarray] = []
        for i in range(0, len(epochs_jd), int(epoch_chunk_size)):
            chunk = epochs_jd[i : i + int(epoch_chunk_size)]
            mjd_i, ra_i, dec_i, delta_i = _query_horizons_epoch_chunk(
                target_id=target_id, id_type=id_type, epoch_chunk_jd=chunk
            )
            mjd_parts.append(mjd_i)
            ra_parts.append(ra_i)
            dec_parts.append(dec_i)
            delta_parts.append(delta_i)
        mjd = np.concatenate(mjd_parts)
        ra = np.concatenate(ra_parts)
        dec = np.concatenate(dec_parts)
        delta = np.concatenate(delta_parts)
    else:
        stop_for_query = stop_utc
        query_step = step.strip()
        if step_seconds is not None and canonical_unit in {"m", "h", "d"}:
            t1 = Time(stop_utc, format="isot", scale="utc") + TimeDelta(
                float(step_seconds), format="sec"
            )
            stop_for_query = t1.isot
            qval = int(
                round(
                    step_seconds
                    / {"m": 60.0, "h": 3600.0, "d": 86400.0}[canonical_unit]
                )
            )
            query_step = f"{qval}{canonical_unit}"
        obj = Horizons(
            id=target_id,
            id_type=id_type,
            location="500@399",
            epochs={
                "start": start_utc.replace("T", " "),
                "stop": stop_for_query.replace("T", " "),
                "step": query_step,
            },
        )
        eph = obj.ephemerides(quantities="1,20")
        mjd, ra, dec, delta = _parse_eph_table(eph)
    order = np.argsort(mjd)
    mjd = np.asarray(mjd, dtype=float)[order]
    ra = np.asarray(ra, dtype=float)[order]
    dec = np.asarray(dec, dtype=float)[order]
    delta = np.asarray(delta, dtype=float)[order]
    keep = np.ones(len(mjd), dtype=bool)
    if len(mjd) > 1:
        keep[1:] = np.diff(mjd) > 0.0
    return (mjd[keep], ra[keep], dec[keep], delta[keep])


@dataclass(frozen=True)
class TrackInterpolator:
    t_sec: np.ndarray
    xyz: np.ndarray
    delta_au: np.ndarray

    @classmethod
    def from_track(
        cls,
        t_sec: np.ndarray,
        ra_deg: np.ndarray,
        dec_deg: np.ndarray,
        delta_au: np.ndarray,
    ) -> "TrackInterpolator":
        sc = SkyCoord(ra=ra_deg * u.deg, dec=dec_deg * u.deg, frame="icrs")
        xyz = sc.cartesian.xyz.value.T.astype(float)
        return cls(
            t_sec=np.array(t_sec, float), xyz=xyz, delta_au=np.array(delta_au, float)
        )

    def at(
        self, t_sec: float | np.ndarray
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        ts = np.atleast_1d(np.array(t_sec, dtype=float))
        x = np.interp(ts, self.t_sec, self.xyz[:, 0])
        y = np.interp(ts, self.t_sec, self.xyz[:, 1])
        z = np.interp(ts, self.t_sec, self.xyz[:, 2])
        norm = np.sqrt(x * x + y * y + z * z)
        norm[norm == 0.0] = 1.0
        x /= norm
        y /= norm
        z /= norm
        ra = np.degrees(np.arctan2(y, x)) % 360.0
        dec = np.degrees(np.arctan2(z, np.sqrt(x * x + y * y)))
        delta = np.interp(ts, self.t_sec, self.delta_au)
        if np.isscalar(t_sec):
            return (ra[0], dec[0], delta[0])
        return (ra, dec, delta)

    def max_offset_from_midpoint_arcsec(self, t0: float, t1: float) -> float:
        tm = 0.5 * (t0 + t1)
        ra_m, dec_m, _ = self.at(tm)
        mid = SkyCoord(ra=ra_m * u.deg, dec=dec_m * u.deg, frame="icrs")
        sel = (self.t_sec >= t0) & (self.t_sec <= t1)
        sample_t = self.t_sec[sel]
        sample_t = np.unique(np.concatenate([sample_t, np.array([t0, tm, t1], float)]))
        ra_s, dec_s, _ = self.at(sample_t)
        sc = SkyCoord(ra=ra_s * u.deg, dec=dec_s * u.deg, frame="icrs")
        sep = sc.separation(mid).to_value(u.arcsec)
        return float(np.max(sep)) if sep.size else 0.0


def build_segments(
    interp: TrackInterpolator,
    obs_tstart: float,
    obs_tstop: float,
    max_blur_arcsec: float,
    snap_s: float,
    min_seg_s: float,
) -> list[tuple[int, float, float, float, float, float, float]]:
    if not obs_tstop > obs_tstart:
        raise RuntimeError("Observation stop time must exceed start time.")
    if max_blur_arcsec <= 0:
        raise RuntimeError("max_blur_arcsec must be > 0.")
    if snap_s <= 0:
        dt = np.diff(interp.t_sec)
        snap_s = float(np.median(dt[dt > 0])) if np.any(dt > 0) else 1.0
    if min_seg_s <= 0:
        min_seg_s = snap_s
    cur = obs_tstart
    segid = 1
    segments: list[tuple[int, float, float, float, float, float, float]] = []
    tol = 1e-06
    while cur < obs_tstop - tol:
        max_end = obs_tstop
        nsteps = int(math.floor((max_end - cur) / snap_s + tol))
        candidates = [cur + k * snap_s for k in range(1, nsteps + 1)]
        if not candidates:
            candidates = [obs_tstop]
        chosen_end = None
        chosen_blur = None
        for end in candidates:
            if end - cur + tol < min_seg_s:
                continue
            blur = interp.max_offset_from_midpoint_arcsec(cur, end)
            if blur <= max_blur_arcsec + 1e-09:
                chosen_end = end
                chosen_blur = blur
            else:
                break
        if chosen_end is None:
            remaining = obs_tstop - cur
            if remaining < min_seg_s and segments:
                prev = segments[-1]
                new_end = obs_tstop
                new_mid = 0.5 * (prev[1] + new_end)
                ra_m, dec_m, delta_m = interp.at(new_mid)
                new_blur = interp.max_offset_from_midpoint_arcsec(prev[1], new_end)
                segments[-1] = (
                    prev[0],
                    prev[1],
                    new_end,
                    float(ra_m),
                    float(dec_m),
                    float(delta_m),
                    float(new_blur),
                )
                break
            if remaining < min_seg_s and (not segments):
                blur = interp.max_offset_from_midpoint_arcsec(cur, obs_tstop)
                tmid = 0.5 * (cur + obs_tstop)
                ra_mid, dec_mid, delta_mid = interp.at(tmid)
                segments.append(
                    (
                        segid,
                        cur,
                        obs_tstop,
                        float(ra_mid),
                        float(dec_mid),
                        float(delta_mid),
                        float(blur),
                    )
                )
                break
            test_end = min(cur + snap_s, obs_tstop)
            test_blur = interp.max_offset_from_midpoint_arcsec(cur, test_end)
            raise RuntimeError(
                f"Cannot satisfy segmentation constraints: one snapped interval of {test_end - cur:.3f} s already has midpoint blur {test_blur:.3f} arcsec, exceeding max_blur_arcsec={max_blur_arcsec:.3f}. Increase MAX_SEG_BLUR_ARCSEC or reduce SEGMENT_SNAP_S / LC_BIN_S."
            )
        tmid = 0.5 * (cur + chosen_end)
        ra_mid, dec_mid, delta_mid = interp.at(tmid)
        segments.append(
            (
                segid,
                cur,
                chosen_end,
                float(ra_mid),
                float(dec_mid),
                float(delta_mid),
                float(chosen_blur),
            )
        )
        segid += 1
        cur = chosen_end
    return segments


def write_track_fits(
    out_file: str,
    mjd: np.ndarray,
    ra: np.ndarray,
    dec: np.ndarray,
    delta: np.ndarray,
    obs_tstart: float,
    obs_tstop: float,
    mjdref: float,
    timesys: str,
    source_desc: str,
) -> None:
    cols = fits.ColDefs(
        [
            fits.Column(
                name="MJD", format="D", unit="d", array=np.array(mjd, dtype=np.float64)
            ),
            fits.Column(
                name="RA", format="E", unit="deg", array=np.array(ra, dtype=np.float32)
            ),
            fits.Column(
                name="DEC",
                format="E",
                unit="deg",
                array=np.array(dec, dtype=np.float32),
            ),
            fits.Column(
                name="DELTA",
                format="E",
                unit="AU",
                array=np.array(delta, dtype=np.float32),
            ),
        ]
    )
    hdu = fits.BinTableHDU.from_columns(cols, name="OBJTRACK")
    hdu.header["OBS_T0"] = (float(obs_tstart), "Obs start in XMM seconds")
    hdu.header["OBS_T1"] = (float(obs_tstop), "Obs stop in XMM seconds")
    hdu.header["MJDREF"] = (float(mjdref), "Reference MJD of XMM TIME")
    hdu.header["TIMESYS"] = (timesys, "Time scale of XMM TIME")
    hdu.header["ORIGINTRK"] = (source_desc[:68], "Track source description")
    hdul = fits.HDUList([fits.PrimaryHDU(), hdu])
    hdul.writeto(out_file, overwrite=True)


def write_segments_csv(
    out_file: str,
    segments: list[tuple[int, float, float, float, float, float, float]],
    mjdref: float,
    timesys: str,
) -> None:
    fields = [
        "segid",
        "tstart",
        "tstop",
        "ra_mid_deg",
        "dec_mid_deg",
        "delta_mid_au",
        "max_mid_offset_arcsec",
        "utc_start",
        "utc_stop",
    ]
    with open(out_file, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(fields)
        for segid, t0, t1, ra, dec, delta, blur in segments:
            w.writerow(
                [
                    segid,
                    f"{t0:.6f}",
                    f"{t1:.6f}",
                    f"{ra:.10f}",
                    f"{dec:.10f}",
                    f"{delta:.8f}",
                    f"{blur:.6f}",
                    xmm_seconds_to_utc_iso(t0, mjdref, timesys),
                    xmm_seconds_to_utc_iso(t1, mjdref, timesys),
                ]
            )


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--reference-event",
        required=True,
        help="EPIC event file used to define TIME range and MJDREF",
    )
    ap.add_argument(
        "--track-out", required=True, help="Output FITS file for movecalc track"
    )
    ap.add_argument(
        "--segments-out",
        required=True,
        help="Output CSV file for motion-limited segments",
    )
    ap.add_argument(
        "--track-input",
        help="Optional input comet track (CSV or FITS) with time/RA/DEC/DELTA",
    )
    ap.add_argument(
        "--target-id", help="Horizons target ID if no track-input is supplied"
    )
    ap.add_argument(
        "--id-type", default="smallbody", help="Horizons id_type (default: smallbody)"
    )
    ap.add_argument(
        "--step",
        default="1m",
        help="Track sampling, e.g. 30s, 1m, 5m. Sub-minute uses explicit epochs.",
    )
    ap.add_argument(
        "--max-blur-arcsec",
        type=float,
        required=True,
        help="Maximum allowed comet offset from segment midpoint",
    )
    ap.add_argument(
        "--snap-s",
        type=float,
        default=0.0,
        help="Snap segment boundaries to this cadence in seconds",
    )
    ap.add_argument(
        "--min-seg-s",
        type=float,
        default=0.0,
        help="Minimum segment duration in seconds",
    )
    return ap.parse_args()


def main() -> int:
    args = parse_args()
    obs_tstart, obs_tstop, mjdref, timesys = read_event_time_info(args.reference_event)
    utc_start = xmm_seconds_to_utc_iso(obs_tstart, mjdref, timesys)
    utc_stop = xmm_seconds_to_utc_iso(obs_tstop, mjdref, timesys)
    if args.track_input:
        track_path = Path(args.track_input)
        if not track_path.exists():
            raise RuntimeError(f"Track input not found: {track_path}")
        if track_path.suffix.lower() in {".fits", ".fit", ".fts"}:
            mjd, ra, dec, delta = load_track_from_fits(str(track_path))
        else:
            mjd, ra, dec, delta = load_track_from_csv(str(track_path))
        source_desc = f"external:{track_path.name}"
    else:
        if not args.target_id:
            raise RuntimeError("Either --track-input or --target-id must be supplied.")
        pad = TimeDelta(120.0, format="sec")
        padded_start = (Time(utc_start, format="isot", scale="utc") - pad).utc.isot
        padded_stop = (Time(utc_stop, format="isot", scale="utc") + pad).utc.isot
        mjd, ra, dec, delta = query_horizons_track(
            target_id=args.target_id,
            id_type=args.id_type,
            start_utc=padded_start,
            stop_utc=padded_stop,
            step=args.step,
        )
        source_desc = f"Horizons:{args.target_id}"
    track_time = Time(np.array(mjd, float), format="mjd", scale="utc")
    ref_time = Time(mjdref, format="mjd", scale=timesys.lower())
    t_sec = (track_time - ref_time).to_value(u.s)
    order = np.argsort(t_sec)
    t_sec = t_sec[order]
    mjd = np.array(mjd, float)[order]
    ra = np.array(ra, float)[order]
    dec = np.array(dec, float)[order]
    delta = np.array(delta, float)[order]
    if t_sec[0] > obs_tstart or t_sec[-1] < obs_tstop:
        raise RuntimeError(
            f"Track does not fully cover the observation. Track spans [{t_sec[0]:.3f}, {t_sec[-1]:.3f}] s but observation spans [{obs_tstart:.3f}, {obs_tstop:.3f}] s."
        )
    interp = TrackInterpolator.from_track(
        t_sec=t_sec, ra_deg=ra, dec_deg=dec, delta_au=delta
    )
    segments = build_segments(
        interp=interp,
        obs_tstart=obs_tstart,
        obs_tstop=obs_tstop,
        max_blur_arcsec=float(args.max_blur_arcsec),
        snap_s=float(args.snap_s),
        min_seg_s=float(args.min_seg_s),
    )
    Path(args.track_out).parent.mkdir(parents=True, exist_ok=True)
    Path(args.segments_out).parent.mkdir(parents=True, exist_ok=True)
    write_track_fits(
        out_file=args.track_out,
        mjd=mjd,
        ra=ra,
        dec=dec,
        delta=delta,
        obs_tstart=obs_tstart,
        obs_tstop=obs_tstop,
        mjdref=mjdref,
        timesys=timesys,
        source_desc=source_desc,
    )
    write_segments_csv(
        out_file=args.segments_out, segments=segments, mjdref=mjdref, timesys=timesys
    )
    print(f"Wrote track FITS:     {args.track_out}")
    print(f"Wrote segments CSV:   {args.segments_out}")
    print(f"Observation UTC span: {utc_start} -> {utc_stop}")
    print(f"Number of segments:   {len(segments)}")
    if segments:
        blurs = [s[6] for s in segments]
        durs = [s[2] - s[1] for s in segments]
        print(
            f"Segment duration s:   min={min(durs):.3f} med={np.median(durs):.3f} max={max(durs):.3f}"
        )
        print(
            f"Midpoint blur arcsec: min={min(blurs):.3f} med={np.median(blurs):.3f} max={max(blurs):.3f}"
        )
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(1)
