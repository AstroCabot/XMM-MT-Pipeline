#!/usr/bin/env python3
"""Shared helpers for stationary-source trail masks in the comet frame.

These helpers are intentionally lightweight and pure-Python so they can be used
both by the science pipeline and by QC scripts without requiring SAS.
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Any

import numpy as np
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.io import fits
from astropy.time import Time


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


class TrackInterpolator:
    """Interpolate RA/Dec on the unit sphere using Cartesian linear blends."""

    def __init__(self, t_sec: np.ndarray, ra_deg: np.ndarray, dec_deg: np.ndarray) -> None:
        sc = SkyCoord(ra=np.asarray(ra_deg, dtype=float) * u.deg, dec=np.asarray(dec_deg, dtype=float) * u.deg, frame="icrs")
        self._t = np.asarray(t_sec, dtype=float)
        self._xyz = sc.cartesian.xyz.value.T.astype(float)

    def at(self, t_sec: np.ndarray | float) -> tuple[np.ndarray, np.ndarray]:
        ts = np.atleast_1d(np.asarray(t_sec, dtype=float))
        if ts.size == 0:
            return np.asarray([], dtype=float), np.asarray([], dtype=float)
        if ts.min() < self._t.min() - 1e-9 or ts.max() > self._t.max() + 1e-9:
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
            return np.asarray([ra[0]], dtype=float), np.asarray([dec[0]], dtype=float)
        return np.asarray(ra, dtype=float), np.asarray(dec, dtype=float)


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
            srcid = np.asarray(data[hdu.columns.names[names_upper.index("SRCID")]], dtype=int)
        else:
            srcid = np.arange(1, len(ra) + 1, dtype=int)
    return SrcList(srcid=srcid, ra=ra, dec=dec)


def paint_disk(mask: np.ndarray, col: float, row: float, r_px: float) -> None:
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
    mask[ymin : ymax + 1, xmin : xmax + 1] |= dist2 <= (r_px * r_px)


def paint_segment(mask: np.ndarray, c0: float, r0: float, c1: float, r1: float, r_px: float) -> None:
    if not (np.isfinite(c0) and np.isfinite(r0) and np.isfinite(c1) and np.isfinite(r1)):
        return
    seg_len = float(math.hypot(c1 - c0, r1 - r0))
    nstep = max(1, int(math.ceil(seg_len / 0.5)))
    for frac in np.linspace(0.0, 1.0, nstep + 1):
        col = (1.0 - frac) * c0 + frac * c1
        row = (1.0 - frac) * r0 + frac * r1
        paint_disk(mask, col, row, r_px)


def clip_track_to_time_range(track: TrackData, t_min_sec: float | None = None, t_max_sec: float | None = None) -> TrackData:
    """Clip a track to a time window, preserving interpolated endpoints.

    The returned track still follows the original great-circle-interpolated path
    but is restricted to the requested interval. If the interval fully contains
    the observation span, the original track object is returned.
    """
    if t_min_sec is None and t_max_sec is None:
        return track
    t_lo = track.obs_t0 if t_min_sec is None else float(t_min_sec)
    t_hi = track.obs_t1 if t_max_sec is None else float(t_max_sec)
    if t_hi < t_lo:
        raise RuntimeError("track time window has t_max_sec < t_min_sec")
    t_lo = max(float(track.t_sec.min()), t_lo)
    t_hi = min(float(track.t_sec.max()), t_hi)
    if t_lo <= float(track.t_sec.min()) + 1e-9 and t_hi >= float(track.t_sec.max()) - 1e-9:
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
    mjd_times = Time(track.mjdref, format="mjd", scale=track.timesys.lower()) + times * u.s
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


def build_trail_mask(
    grid: dict[str, Any],
    track_fits: str,
    srclist_fits: str,
    mask_radius_arcsec: float,
    shape: tuple[int, int],
    t_min_sec: float | None = None,
    t_max_sec: float | None = None,
) -> np.ndarray:
    track = clip_track_to_time_range(read_track(track_fits), t_min_sec=t_min_sec, t_max_sec=t_max_sec)
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
        src_coord = SkyCoord(ra=float(ra_s) * u.deg, dec=float(dec_s) * u.deg, frame="icrs")
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
