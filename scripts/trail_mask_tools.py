#!/usr/bin/env python3
"""Trail-mask construction and application for the comet-frame pipeline.

Merges the old make_trail_mask.py, apply_spatial_trail_mask.py, and
trail_mask_utils.py into one module with two CLI subcommands:

    trail_mask_tools.py make   -- build a 2-D stationary-source trail mask
    trail_mask_tools.py apply  -- filter events that fall inside the mask

Also exports helpers used by combine_epic_images.py, make_still_sky_mosaic.py,
and image_postproc.py:
    TrackData, SrcList, TrackInterpolator, read_track, read_srclist,
    build_trail_mask, paint_disk
"""
from __future__ import annotations
import argparse, json, math, sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any
import numpy as np
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.io import fits
from astropy.time import Time

# ── Data containers ──────────────────────────────────────────────────────────


@dataclass(frozen=True)
class TrackData:
    mjd: np.ndarray
    ra: np.ndarray
    dec: np.ndarray
    mjdref: float
    timesys: str
    obs_t0: float
    obs_t1: float
    t_sec: np.ndarray


@dataclass(frozen=True)
class SrcList:
    srcid: np.ndarray
    ra: np.ndarray
    dec: np.ndarray


# ── Track interpolation ──────────────────────────────────────────────────────


class TrackInterpolator:
    """Spherical linear interpolator for a comet ephemeris track.

    Converts (RA, Dec) to Cartesian, interpolates each component linearly
    in time, then re-normalises back to the unit sphere.  This avoids RA
    wrap-around artefacts near 0/360 deg.
    """

    def __init__(
        self, t_sec: np.ndarray, ra_deg: np.ndarray, dec_deg: np.ndarray
    ) -> None:
        sc = SkyCoord(
            ra=np.asarray(ra_deg, dtype=float) * u.deg,
            dec=np.asarray(dec_deg, dtype=float) * u.deg,
            frame="icrs",
        )
        self._t = np.asarray(t_sec, dtype=float)
        self._xyz = sc.cartesian.xyz.value.T.astype(float)

    def at(self, t_sec: np.ndarray | float) -> tuple[np.ndarray, np.ndarray]:
        ts = np.atleast_1d(np.asarray(t_sec, dtype=float))
        if ts.size == 0:
            return (np.asarray([], dtype=float), np.asarray([], dtype=float))
        if ts.min() < self._t.min() - 1e-09 or ts.max() > self._t.max() + 1e-09:
            raise RuntimeError("Requested time is outside the track range")
        x = np.interp(ts, self._t, self._xyz[:, 0])
        y = np.interp(ts, self._t, self._xyz[:, 1])
        z = np.interp(ts, self._t, self._xyz[:, 2])
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


# ── FITS I/O helpers ─────────────────────────────────────────────────────────


def first_table(hdul: fits.HDUList, *, allow_empty: bool = False) -> fits.BinTableHDU:
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
    raise RuntimeError("No binary table found")


def first_nonempty_table(hdul: fits.HDUList) -> fits.BinTableHDU:
    return first_table(hdul, allow_empty=False)


def read_track(track_fits: str) -> TrackData:
    with fits.open(track_fits) as hdul:
        hdu = first_nonempty_table(hdul)
        tab = hdu.data
        hdr = hdu.header
        mjd = np.asarray(tab["MJD"], dtype=float)
        ra = np.asarray(tab["RA"], dtype=float)
        dec = np.asarray(tab["DEC"], dtype=float)
        mjdref = float(hdr.get("MJDREF", 0.0))
        timesys = str(hdr.get("TIMESYS", "TT")).strip().upper()
        obs_t0 = float(hdr.get("OBS_T0", np.nan)) if "OBS_T0" in hdr else math.nan
        obs_t1 = float(hdr.get("OBS_T1", np.nan)) if "OBS_T1" in hdr else math.nan
    track_time = Time(mjd, format="mjd", scale="utc")
    ref_time = Time(mjdref, format="mjd", scale=timesys.lower())
    t_sec = np.asarray((track_time - ref_time).to_value(u.s), dtype=float)
    if not np.isfinite(obs_t0):
        obs_t0 = float(t_sec[0])
    if not np.isfinite(obs_t1):
        obs_t1 = float(t_sec[-1])
    return TrackData(
        mjd=mjd,
        ra=ra,
        dec=dec,
        mjdref=mjdref,
        timesys=timesys,
        obs_t0=obs_t0,
        obs_t1=obs_t1,
        t_sec=t_sec,
    )


def read_srclist(src_fits: str) -> SrcList:
    with fits.open(src_fits) as hdul:
        hdu = first_table(hdul, allow_empty=True)
        names_upper = [name.upper() for name in hdu.columns.names]
        if "RA" not in names_upper or "DEC" not in names_upper:
            raise RuntimeError("Source list must contain RA and DEC columns")
        data = hdu.data
        if data is None or len(data) == 0:
            return SrcList(
                srcid=np.asarray([], dtype=int),
                ra=np.asarray([], dtype=float),
                dec=np.asarray([], dtype=float),
            )
        ra = np.asarray(data[hdu.columns.names[names_upper.index("RA")]], dtype=float)
        dec = np.asarray(data[hdu.columns.names[names_upper.index("DEC")]], dtype=float)
        if "SRCID" in names_upper:
            srcid = np.asarray(
                data[hdu.columns.names[names_upper.index("SRCID")]], dtype=int
            )
        else:
            srcid = np.arange(1, len(ra) + 1, dtype=int)
    return SrcList(srcid=srcid, ra=ra, dec=dec)


# ── Mask painting primitives ──────────────────────────────────────────────────


def paint_disk(mask: np.ndarray, col: float, row: float, r_px: float) -> None:
    """Set mask pixels within *r_px* of (col, row) to True."""
    if not np.isfinite(col) or not np.isfinite(row):
        return
    ny, nx = mask.shape
    xmin = max(0, int(math.floor(col - r_px - 1.0)))
    xmax = min(nx - 1, int(math.ceil(col + r_px + 1.0)))
    ymin = max(0, int(math.floor(row - r_px - 1.0)))
    ymax = min(ny - 1, int(math.ceil(row + r_px + 1.0)))
    if xmin > xmax or ymin > ymax:
        return
    yy, xx = np.mgrid[ymin : ymax + 1, xmin : xmax + 1]
    dist2 = (xx + 0.5 - col) ** 2 + (yy + 0.5 - row) ** 2
    mask[ymin : ymax + 1, xmin : xmax + 1] |= dist2 <= r_px * r_px


def paint_segment(
    mask: np.ndarray, c0: float, r0: float, c1: float, r1: float, r_px: float
) -> None:
    """Paint a capsule (thick line) between two pixel positions."""
    if not (
        np.isfinite(c0) and np.isfinite(r0) and np.isfinite(c1) and np.isfinite(r1)
    ):
        return
    seg_len = float(math.hypot(c1 - c0, r1 - r0))
    nstep = max(1, int(math.ceil(seg_len / 0.5)))
    for frac in np.linspace(0.0, 1.0, nstep + 1):
        col = (1.0 - frac) * c0 + frac * c1
        row = (1.0 - frac) * r0 + frac * r1
        paint_disk(mask, col, row, r_px)


def clip_track_to_time_range(
    track: TrackData, t_min_sec: float | None = None, t_max_sec: float | None = None
) -> TrackData:
    if t_min_sec is None and t_max_sec is None:
        return track
    t_lo = track.obs_t0 if t_min_sec is None else float(t_min_sec)
    t_hi = track.obs_t1 if t_max_sec is None else float(t_max_sec)
    if t_hi < t_lo:
        raise RuntimeError("track time window has t_max_sec < t_min_sec")
    t_lo = max(float(track.t_sec.min()), t_lo)
    t_hi = min(float(track.t_sec.max()), t_hi)
    if (
        t_lo <= float(track.t_sec.min()) + 1e-09
        and t_hi >= float(track.t_sec.max()) - 1e-09
    ):
        return track
    interp = TrackInterpolator(track.t_sec, track.ra, track.dec)
    keep = (track.t_sec > t_lo) & (track.t_sec < t_hi)
    times = np.concatenate(
        [
            np.asarray([t_lo], dtype=float),
            np.asarray(track.t_sec[keep], dtype=float),
            np.asarray([t_hi], dtype=float),
        ]
    )
    ra, dec = interp.at(times)
    mjd_times = (
        Time(track.mjdref, format="mjd", scale=track.timesys.lower()) + times * u.s
    )
    mjd_utc = mjd_times.utc.mjd
    return TrackData(
        mjd=np.asarray(mjd_utc, dtype=float),
        ra=np.asarray(ra, dtype=float),
        dec=np.asarray(dec, dtype=float),
        mjdref=track.mjdref,
        timesys=track.timesys,
        obs_t0=t_lo,
        obs_t1=t_hi,
        t_sec=np.asarray(times, dtype=float),
    )


# ── Core mask construction ────────────────────────────────────────────────────


def build_trail_mask(
    grid: dict[str, Any],
    track_fits: str,
    srclist_fits: str,
    mask_radius_arcsec: float,
    shape: tuple[int, int],
    t_min_sec: float | None = None,
    t_max_sec: float | None = None,
) -> np.ndarray:
    """Build a boolean trail mask in the comet frame.

    For each detected field source, traces its apparent path across the
    comet-frame image grid over the observation and paints a disk/capsule
    of radius *mask_radius_arcsec* at every track sample.
    """
    track = clip_track_to_time_range(
        read_track(track_fits), t_min_sec=t_min_sec, t_max_sec=t_max_sec
    )
    src = read_srclist(srclist_fits)
    center_x = float(grid["center_x_phys"])
    center_y = float(grid["center_y_phys"])
    x_min = float(grid["x_min_phys"])
    y_min = float(grid["y_min_phys"])
    bin_phys = float(grid["bin_phys"])
    x_arcsec_per_phys = float(grid["x_arcsec_per_phys"])
    y_arcsec_per_phys = float(grid["y_arcsec_per_phys"])
    r_phys = float(mask_radius_arcsec) / float(grid["scale_arcsec_per_phys"])
    r_px = r_phys / bin_phys
    comet = SkyCoord(ra=track.ra * u.deg, dec=track.dec * u.deg, frame="icrs")
    mask = np.zeros(shape, dtype=bool)
    for ra_s, dec_s in zip(src.ra, src.dec):
        src_coord = SkyCoord(
            ra=float(ra_s) * u.deg, dec=float(dec_s) * u.deg, frame="icrs"
        )
        dlon, dlat = comet.spherical_offsets_to(src_coord)
        dx = dlon.to_value(u.arcsec)
        dy = dlat.to_value(u.arcsec)
        x_phys = center_x + dx / x_arcsec_per_phys
        y_phys = center_y + dy / y_arcsec_per_phys
        cols = (x_phys - x_min) / bin_phys
        rows = (y_phys - y_min) / bin_phys
        if len(cols) == 1:
            paint_disk(mask, cols[0], rows[0], r_px)
            continue
        for i in range(len(cols) - 1):
            paint_segment(mask, cols[i], rows[i], cols[i + 1], rows[i + 1], r_px)
    return mask


def _read_mask(path: str) -> np.ndarray:
    with fits.open(path) as hdul:
        for hdu in hdul:
            if hdu.data is not None and getattr(hdu.data, "ndim", 0) == 2:
                return np.asarray(hdu.data, dtype=bool)
    raise RuntimeError(f"No 2D mask image found in {path}")


def _filter_events(
    event_path: str, mask: np.ndarray, grid: dict[str, float]
) -> tuple[fits.HDUList, int, int, int]:
    with fits.open(event_path) as hdul:
        if "EVENTS" not in hdul:
            raise RuntimeError(f"No EVENTS extension in {event_path}")
        evt = hdul["EVENTS"]
        data = evt.data
        if data is None:
            raise RuntimeError(f"EVENTS table is empty in {event_path}")
        x = np.asarray(data["X"], dtype=float)
        y = np.asarray(data["Y"], dtype=float)
        col = np.floor((x - float(grid["x_min_phys"])) / float(grid["bin_phys"]))
        row = np.floor((y - float(grid["y_min_phys"])) / float(grid["bin_phys"]))
        col_i = col.astype(int, copy=False)
        row_i = row.astype(int, copy=False)
        in_bounds = (
            np.isfinite(col)
            & np.isfinite(row)
            & (col_i >= 0)
            & (col_i < mask.shape[1])
            & (row_i >= 0)
            & (row_i < mask.shape[0])
        )
        masked = np.zeros(len(data), dtype=bool)
        masked[in_bounds] = mask[row_i[in_bounds], col_i[in_bounds]]
        keep = ~masked
        out_hdus = []
        for hdu in hdul:
            if hdu.name == "EVENTS":
                new_hdu = fits.BinTableHDU(
                    data=data[keep], header=hdu.header.copy(), name=hdu.name
                )
                new_hdu.header["HISTORY"] = (
                    "Stationary-source trail-mask events removed"
                )
                new_hdu.header["TRAILMSK"] = (
                    True,
                    "Stationary-source trail mask applied",
                )
                new_hdu.header["TRMSKREM"] = (
                    int(masked.sum()),
                    "Events removed by trail mask",
                )
                out_hdus.append(new_hdu)
            else:
                out_hdus.append(hdu.copy())
    return (
        fits.HDUList(out_hdus),
        int(len(data)),
        int(masked.sum()),
        int(np.count_nonzero(in_bounds)),
    )


# ── CLI subcommands ──────────────────────────────────────────────────────────


def cmd_make(args: argparse.Namespace) -> int:
    """Subcommand: build the 2-D trail mask FITS image."""
    grid = json.loads(Path(args.grid_json).read_text(encoding="utf-8"))
    nx = int(grid["nx"])
    ny = int(grid["ny"])
    if nx <= 0 or ny <= 0:
        raise RuntimeError("grid.json has non-positive nx/ny")
    try:
        src_rows = len(read_srclist(args.srclist).ra)
    except Exception:
        src_rows = 0
    if args.mask_radius_arcsec <= 0.0 or src_rows == 0:
        mask = np.zeros((ny, nx), dtype=np.int16)
    else:
        mask = build_trail_mask(
            grid=grid,
            track_fits=args.track,
            srclist_fits=args.srclist,
            mask_radius_arcsec=float(args.mask_radius_arcsec),
            shape=(ny, nx),
            t_min_sec=args.track_tmin_sec,
            t_max_sec=args.track_tmax_sec,
        ).astype(np.int16)
    hdr = fits.Header()
    hdr["BUNIT"] = ("1", "1=masked trail pixel")
    hdr["XMINPHY"] = float(grid["x_min_phys"])
    hdr["XMAXPHY"] = float(grid["x_max_phys"])
    hdr["YMINPHY"] = float(grid["y_min_phys"])
    hdr["YMAXPHY"] = float(grid["y_max_phys"])
    hdr["BINPHYS"] = float(grid["bin_phys"])
    hdr["SCALEAS"] = float(grid["scale_arcsec_per_phys"])
    hdr["MASKRAD"] = float(args.mask_radius_arcsec)
    if args.track_tmin_sec is not None:
        hdr["TMINSEC"] = float(args.track_tmin_sec)
    if args.track_tmax_sec is not None:
        hdr["TMAXSEC"] = float(args.track_tmax_sec)
    hdr["HISTORY"] = "Stationary-source trail mask in the comet frame"
    out_path = Path(args.out_mask)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fits.PrimaryHDU(data=mask, header=hdr).writeto(out_path, overwrite=True)
    if args.out_json:
        payload = {
            "grid_json": str(Path(args.grid_json).resolve()),
            "track": str(Path(args.track).resolve()),
            "srclist": str(Path(args.srclist).resolve()),
            "mask_radius_arcsec": float(args.mask_radius_arcsec),
            "shape": [int(ny), int(nx)],
            "n_masked_pixels": int(mask.sum()),
            "masked_fraction": float(np.mean(mask > 0)),
            "n_sources": int(src_rows),
            "track_tmin_sec": (
                None if args.track_tmin_sec is None else float(args.track_tmin_sec)
            ),
            "track_tmax_sec": (
                None if args.track_tmax_sec is None else float(args.track_tmax_sec)
            ),
        }
        out_json = Path(args.out_json)
        out_json.parent.mkdir(parents=True, exist_ok=True)
        out_json.write_text(
            json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
    print(out_path)
    return 0


def cmd_apply(args: argparse.Namespace) -> int:
    """Subcommand: remove events falling inside a pre-built trail mask."""
    grid = json.loads(Path(args.grid_json).read_text(encoding="utf-8"))
    mask = _read_mask(args.mask)
    hdul, n_in, n_removed, n_in_bounds = _filter_events(args.event, mask, grid)
    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    hdul.writeto(out_path, overwrite=True)
    if args.report_json:
        payload = {
            "event": str(Path(args.event).resolve()),
            "mask": str(Path(args.mask).resolve()),
            "output": str(out_path.resolve()),
            "input_rows": int(n_in),
            "removed_rows": int(n_removed),
            "kept_rows": int(n_in - n_removed),
            "removed_fraction": float(n_removed / n_in) if n_in else 0.0,
            "in_bounds_rows": int(n_in_bounds),
        }
        report_path = Path(args.report_json)
        report_path.parent.mkdir(parents=True, exist_ok=True)
        report_path.write_text(
            json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
    print(out_path)
    return 0


def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser()
    sp = ap.add_subparsers(dest="cmd", required=True)
    p = sp.add_parser("make")
    p.add_argument("--grid-json", required=True)
    p.add_argument("--track", required=True)
    p.add_argument("--srclist", required=True)
    p.add_argument("--mask-radius-arcsec", type=float, required=True)
    p.add_argument("--out-mask", required=True)
    p.add_argument("--out-json")
    p.add_argument("--track-tmin-sec", type=float)
    p.add_argument("--track-tmax-sec", type=float)
    p.set_defaults(func=cmd_make)
    p = sp.add_parser("apply")
    p.add_argument("--event", required=True)
    p.add_argument("--mask", required=True)
    p.add_argument("--grid-json", required=True)
    p.add_argument("--output", required=True)
    p.add_argument("--report-json")
    p.set_defaults(func=cmd_apply)
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
