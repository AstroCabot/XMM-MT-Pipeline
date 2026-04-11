#!/usr/bin/env python3
"""Comprehensive multi-panel QC check suite for the comet pipeline.

Produces diagnostic PNG plots and a Markdown report for:
  clean, exposures, track, detect, image, contam, lcurve, spectrum,
  mosaics, spotcheck, and region-support.

Usage:  xmm_comet_quick_checks.py --config <env> <check> [<check> ...]
"""
from __future__ import annotations
import argparse
import csv
import json
import math
import os
import re
import shutil
import shlex
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.io import fits
from astropy.time import Time, TimeDelta
from astropy.visualization import simple_norm
from astropy.wcs import WCS
from matplotlib.colors import PowerNorm

try:
    from astroquery.jplhorizons import Horizons
except Exception:
    Horizons = None


def load_shell_env(path: str) -> dict[str, str]:
    quoted = shlex.quote(str(path))
    cmd = ["bash", "-lc", f"set -a; source {quoted} >/dev/null 2>&1; env -0"]
    proc = subprocess.run(cmd, capture_output=True, check=True)
    env: dict[str, str] = {}
    for chunk in proc.stdout.split(b"\x00"):
        if not chunk or b"=" not in chunk:
            continue
        key_b, val_b = chunk.split(b"=", 1)
        key = key_b.decode("utf-8", errors="replace")
        val = val_b.decode("utf-8", errors="replace")
        env[key] = val
    return env


def env_float(env: dict[str, str], key: str, default: float) -> float:
    val = env.get(key, "")
    if val == "":
        return float(default)
    return float(val)


def env_int(env: dict[str, str], key: str, default: int) -> int:
    val = env.get(key, "")
    if val == "":
        return int(default)
    return int(val)


def parse_image_bands(bands: str) -> list[tuple[str, int, int]]:
    out: list[tuple[str, int, int]] = []
    for entry in bands.split(";"):
        entry = entry.strip()
        if not entry:
            continue
        parts = entry.split(":")
        if len(parts) != 3:
            raise RuntimeError(f"Bad IMAGE_BANDS entry: {entry}")
        label, pimin, pimax = parts
        out.append((label.strip(), int(pimin), int(pimax)))
    return out


def read_list(path: Path) -> list[str]:
    if not path.exists():
        return []
    return [
        ln.strip() for ln in path.read_text(encoding="utf-8").splitlines() if ln.strip()
    ]


def fits_rows(path: Path) -> int | None:
    if not path.exists():
        return None
    try:
        with fits.open(path) as hdul:
            if "EVENTS" in hdul:
                data = hdul["EVENTS"].data
                return 0 if data is None else int(len(data))
            for hdu in hdul[1:]:
                if isinstance(hdu, fits.BinTableHDU) and hdu.data is not None:
                    return int(len(hdu.data))
    except Exception:
        return None
    return None


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


@dataclass
class TrackData:
    mjd: np.ndarray
    ra: np.ndarray
    dec: np.ndarray
    delta: np.ndarray
    mjdref: float
    timesys: str
    obs_t0: float
    obs_t1: float

    @property
    def t_sec(self) -> np.ndarray:
        track_time = Time(self.mjd, format="mjd", scale="utc")
        ref_time = Time(self.mjdref, format="mjd", scale=self.timesys.lower())
        return np.asarray((track_time - ref_time).to_value(u.s), dtype=float)


@dataclass
class SrcList:
    srcid: np.ndarray
    ra: np.ndarray
    dec: np.ndarray
    detml: np.ndarray | None = None


@dataclass
class LightCurve:
    name: str
    time: np.ndarray
    rate: np.ndarray
    error: np.ndarray
    fracexp: np.ndarray | None
    header: fits.Header


@dataclass
class Spectrum:
    channel: np.ndarray
    counts: np.ndarray
    exposure: float | None
    backscal: float | None
    header: fits.Header


@dataclass
class GTI:
    start: np.ndarray
    stop: np.ndarray
    mjdref: float | None = None
    timesys: str | None = None


def read_track(path: str) -> TrackData:
    with fits.open(path) as hdul:
        hdu = first_nonempty_table(hdul)
        tab = hdu.data
        hdr = hdu.header
        delta = (
            np.asarray(tab["DELTA"], dtype=float)
            if "DELTA" in hdu.columns.names
            else np.full(len(tab), np.nan)
        )
        return TrackData(
            mjd=np.asarray(tab["MJD"], dtype=float),
            ra=np.asarray(tab["RA"], dtype=float),
            dec=np.asarray(tab["DEC"], dtype=float),
            delta=delta,
            mjdref=float(hdr.get("MJDREF", 0.0)),
            timesys=str(hdr.get("TIMESYS", "TT")),
            obs_t0=float(hdr.get("OBS_T0", np.nan)),
            obs_t1=float(hdr.get("OBS_T1", np.nan)),
        )


def read_motion_segments(csv_path: str) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    with open(csv_path, "r", newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            row2: dict[str, Any] = dict(row)
            for key in [
                "segid",
                "tstart",
                "tstop",
                "ra_mid_deg",
                "dec_mid_deg",
                "delta_mid_au",
                "max_mid_offset_arcsec",
            ]:
                if key in row2 and row2[key] not in {"", None}:
                    row2[key] = float(row2[key])
            rows.append(row2)
    return rows


def read_srclist(path: str) -> SrcList:
    with fits.open(path) as hdul:
        hdu = first_table(hdul, allow_empty=True)
        tab = hdu.data
        names = {n.upper(): n for n in hdu.columns.names}
        if tab is None or len(tab) == 0:
            return SrcList(
                srcid=np.asarray([], dtype=int),
                ra=np.asarray([], dtype=float),
                dec=np.asarray([], dtype=float),
                detml=np.asarray([], dtype=float),
            )
        srcid = np.asarray(tab[names.get("SRCID", hdu.columns.names[0])], dtype=int)
        ra = np.asarray(tab[names["RA"]], dtype=float)
        dec = np.asarray(tab[names["DEC"]], dtype=float)
        detml = (
            np.asarray(tab[names["DET_ML"]], dtype=float) if "DET_ML" in names else None
        )
        return SrcList(srcid=srcid, ra=ra, dec=dec, detml=detml)


def read_image(path: str) -> tuple[np.ndarray, fits.Header]:
    with fits.open(path) as hdul:
        for hdu in hdul:
            if hdu.data is not None and getattr(hdu.data, "ndim", 0) == 2:
                return (np.asarray(hdu.data, dtype=float), hdu.header.copy())
    raise RuntimeError(f"No 2D image found in {path}")


def read_gti(path: str) -> GTI:
    with fits.open(path) as hdul:
        hdu = first_nonempty_table(hdul)
        tab = hdu.data
        hdr = hdu.header
        return GTI(
            start=np.asarray(tab["START"], dtype=float),
            stop=np.asarray(tab["STOP"], dtype=float),
            mjdref=float(hdr.get("MJDREF", np.nan)) if "MJDREF" in hdr else None,
            timesys=str(hdr.get("TIMESYS", "TT")) if "TIMESYS" in hdr else None,
        )


def read_lightcurve(path: str, name: str) -> LightCurve:
    with fits.open(path) as hdul:
        target = None
        for hdu in hdul[1:]:
            if isinstance(hdu, fits.BinTableHDU) and hdu.data is not None:
                cols = {c.upper() for c in hdu.columns.names}
                if {"TIME", "RATE", "ERROR"}.issubset(cols):
                    target = hdu
                    break
        if target is None:
            raise RuntimeError(f"No RATE table found in {path}")
        tab = target.data
        fracexp = (
            np.asarray(tab["FRACEXP"], dtype=float)
            if "FRACEXP" in target.columns.names
            else None
        )
        return LightCurve(
            name=name,
            time=np.asarray(tab["TIME"], dtype=float),
            rate=np.asarray(tab["RATE"], dtype=float),
            error=np.asarray(tab["ERROR"], dtype=float),
            fracexp=fracexp,
            header=target.header.copy(),
        )


def read_spectrum(path: str) -> Spectrum:
    with fits.open(path) as hdul:
        target = None
        for hdu in hdul[1:]:
            if isinstance(hdu, fits.BinTableHDU) and hdu.data is not None:
                cols = {c.upper() for c in hdu.columns.names}
                if "COUNTS" in cols:
                    target = hdu
                    break
        if target is None:
            raise RuntimeError(f"No COUNTS table found in {path}")
        tab = target.data
        names = {n.upper(): n for n in target.columns.names}
        chan_col = names.get("CHANNEL", list(names.values())[0])
        exposure = target.header.get("EXPOSURE")
        backscal = target.header.get("BACKSCAL")
        return Spectrum(
            channel=np.asarray(tab[chan_col], dtype=float),
            counts=np.asarray(tab[names["COUNTS"]], dtype=float),
            exposure=float(exposure) if exposure is not None else None,
            backscal=float(backscal) if backscal is not None else None,
            header=target.header.copy(),
        )


def savefig(fig: plt.Figure, path: str, tight_kw: dict | None = None) -> None:
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    # Skip tight_layout when the figure uses constrained_layout (incompatible).
    engine = fig.get_layout_engine()
    uses_constrained = (
        engine is not None and type(engine).__name__ == "ConstrainedLayoutEngine"
    )
    if not uses_constrained:
        fig.tight_layout(**tight_kw or {})
    fig.savefig(path, dpi=180, bbox_inches="tight", pad_inches=0.25)
    plt.close(fig)


def circle_mask_from_extent(
    shape: tuple[int, int],
    extent: tuple[float, float, float, float],
    dx: float,
    dy: float,
    r: float,
) -> np.ndarray:
    ny, nx = shape
    xmin, xmax, ymin, ymax = extent
    xs = np.linspace(xmin, xmax, nx, endpoint=False) + 0.5 * (xmax - xmin) / nx
    ys = np.linspace(ymin, ymax, ny, endpoint=False) + 0.5 * (ymax - ymin) / ny
    xx, yy = np.meshgrid(xs, ys)
    return (xx - dx) ** 2 + (yy - dy) ** 2 <= r * r


def annulus_mask_from_extent(
    shape: tuple[int, int],
    extent: tuple[float, float, float, float],
    dx: float,
    dy: float,
    rin: float,
    rout: float,
) -> np.ndarray:
    ny, nx = shape
    xmin, xmax, ymin, ymax = extent
    xs = np.linspace(xmin, xmax, nx, endpoint=False) + 0.5 * (xmax - xmin) / nx
    ys = np.linspace(ymin, ymax, ny, endpoint=False) + 0.5 * (ymax - ymin) / ny
    xx, yy = np.meshgrid(xs, ys)
    rr2 = (xx - dx) ** 2 + (yy - dy) ** 2
    return (rr2 >= rin * rin) & (rr2 <= rout * rout)


def add_apertures(ax: plt.Axes, env: dict[str, str]) -> None:
    from matplotlib.patches import Circle

    src_dx = env_float(env, "SRC_DX_ARCSEC", 0.0)
    src_dy = env_float(env, "SRC_DY_ARCSEC", 0.0)
    src_r = env_float(env, "SRC_R_ARCSEC", 60.0)
    ax.add_patch(Circle((src_dx, src_dy), src_r, fill=False, lw=1.5, ls="-"))
    bkg_mode = env.get("BKG_MODE", "annulus").strip().lower()
    bkg_dx = env_float(env, "BKG_DX_ARCSEC", 0.0)
    bkg_dy = env_float(env, "BKG_DY_ARCSEC", 0.0)
    if bkg_mode == "circle":
        bkg_r = env_float(env, "BKG_R_ARCSEC", 150.0)
        ax.add_patch(Circle((bkg_dx, bkg_dy), bkg_r, fill=False, lw=1.2, ls="--"))
    else:
        rin = env_float(env, "BKG_RIN_ARCSEC", 90.0)
        rout = env_float(env, "BKG_ROUT_ARCSEC", 150.0)
        ax.add_patch(Circle((bkg_dx, bkg_dy), rin, fill=False, lw=1.2, ls="--"))
        ax.add_patch(Circle((bkg_dx, bkg_dy), rout, fill=False, lw=1.2, ls="--"))


def nice_norm(data: np.ndarray, stretch: str = "sqrt"):
    finite = np.isfinite(data)
    if not np.any(finite):
        return None
    vals = data[finite]
    try:
        return simple_norm(vals, stretch=stretch, min_percent=1.0, max_percent=99.5)
    except Exception:
        return None


def split_visible_segments(mask: np.ndarray) -> list[tuple[int, int]]:
    segs: list[tuple[int, int]] = []
    i = 0
    n = len(mask)
    while i < n:
        while i < n and (not mask[i]):
            i += 1
        if i >= n:
            break
        j = i + 1
        while j < n and mask[j]:
            j += 1
        segs.append((i, j))
        i = j
    return segs


def gti_mask_from_samples(times: np.ndarray, gti: GTI) -> np.ndarray:
    mask = np.zeros(len(times), dtype=bool)
    for start, stop in zip(gti.start, gti.stop):
        mask |= (times >= start) & (times <= stop)
    return mask


def sky_offsets_arcsec(
    track_ra: np.ndarray, track_dec: np.ndarray, ref_ra: float, ref_dec: float
) -> tuple[np.ndarray, np.ndarray]:
    track = SkyCoord(
        ra=np.asarray(track_ra) * u.deg, dec=np.asarray(track_dec) * u.deg, frame="icrs"
    )
    ref = SkyCoord(ra=float(ref_ra) * u.deg, dec=float(ref_dec) * u.deg, frame="icrs")
    dlon, dlat = ref.spherical_offsets_to(track)
    return (dlon.to_value(u.arcsec), dlat.to_value(u.arcsec))


def source_offsets_from_track(
    track_ra: np.ndarray, track_dec: np.ndarray, src_ra: float, src_dec: float
) -> tuple[np.ndarray, np.ndarray]:
    comet = SkyCoord(
        ra=np.asarray(track_ra) * u.deg, dec=np.asarray(track_dec) * u.deg, frame="icrs"
    )
    src = SkyCoord(ra=float(src_ra) * u.deg, dec=float(src_dec) * u.deg, frame="icrs")
    dlon, dlat = comet.spherical_offsets_to(src)
    return (dlon.to_value(u.arcsec), dlat.to_value(u.arcsec))


def overlaps_circle(dist: np.ndarray, radius: float, mask_r: float) -> np.ndarray:
    return dist <= radius + mask_r


def overlaps_annulus(
    dist: np.ndarray, rin: float, rout: float, mask_r: float
) -> np.ndarray:
    return (dist + mask_r >= rin) & (dist - mask_r <= rout)


_HSTEP_RE = re.compile("^\\s*([0-9]*\\.?[0-9]+)\\s*([A-Za-z]+)?\\s*$")


def _qc_horizons_explicit_epochs(
    start_utc: str, stop_utc: str, step_seconds: float
) -> np.ndarray:
    t0 = Time(start_utc, format="isot", scale="utc")
    t1 = Time(stop_utc, format="isot", scale="utc")
    total = max(0.0, (t1 - t0).to_value(u.s))
    offsets = np.arange(0.0, total + 0.5 * step_seconds, step_seconds, dtype=float)
    times = t0 + TimeDelta(offsets, format="sec")
    if times[-1] < t1:
        times = Time(
            np.concatenate([times.jd, np.array([t1.jd])]), format="jd", scale="utc"
        )
    return np.unique(
        np.concatenate(
            [np.array([t0.jd]), np.array(times.jd, float), np.array([t1.jd])]
        )
    ).astype(float)


def compare_with_horizons(
    env: dict[str, str], track: TrackData
) -> tuple[float, float] | None:
    target_id = env.get("TARGET_ID", "").strip()
    if not target_id or Horizons is None:
        return None
    id_type = env.get("TARGET_ID_TYPE", "smallbody").strip() or "smallbody"
    utc = Time(track.mjd, format="mjd", scale="utc").utc
    step_sec = np.median(np.diff(track.t_sec)) if len(track.t_sec) > 1 else 60.0
    step_sec = max(1.0, float(step_sec))
    try:
        epochs_jd = _qc_horizons_explicit_epochs(utc[0].isot, utc[-1].isot, step_sec)
        mjd_parts = []
        ra_parts = []
        dec_parts = []
        for i in range(0, len(epochs_jd), 80):
            obj = Horizons(
                id=target_id,
                id_type=id_type,
                location="500@399",
                epochs=[float(x) for x in epochs_jd[i : i + 80]],
            )
            eph = obj.ephemerides()
            jd = np.asarray(eph["datetime_jd"], dtype=float)
            ra_parts.append(np.asarray(eph["RA"], dtype=float))
            dec_parts.append(np.asarray(eph["DEC"], dtype=float))
            mjd_parts.append(jd - 2400000.5)
        mjd_h = np.concatenate(mjd_parts)
        ra_h = np.concatenate(ra_parts)
        dec_h = np.concatenate(dec_parts)
    except Exception:
        return None
    sc_h = SkyCoord(ra=ra_h * u.deg, dec=dec_h * u.deg, frame="icrs")
    xyz = sc_h.cartesian.xyz.value.T.astype(float)
    t_h = np.asarray(
        (
            Time(mjd_h, format="mjd", scale="utc")
            - Time(track.mjdref, format="mjd", scale=track.timesys.lower())
        ).to_value(u.s),
        dtype=float,
    )
    t = track.t_sec
    if t.min() < t_h.min() or t.max() > t_h.max():
        return None
    x = np.interp(t, t_h, xyz[:, 0])
    y = np.interp(t, t_h, xyz[:, 1])
    z = np.interp(t, t_h, xyz[:, 2])
    n = np.sqrt(x * x + y * y + z * z)
    n[n == 0] = 1.0
    x /= n
    y /= n
    z /= n
    ra_i = np.degrees(np.arctan2(y, x)) % 360.0
    dec_i = np.degrees(np.arctan2(z, np.sqrt(x * x + y * y)))
    sep = (
        SkyCoord(ra=track.ra * u.deg, dec=track.dec * u.deg)
        .separation(SkyCoord(ra=ra_i * u.deg, dec=dec_i * u.deg))
        .to_value(u.arcsec)
    )
    return (float(np.nanmedian(sep)), float(np.nanmax(sep)))


def _copy_png(src: Path, dst: Path) -> None:
    if src.resolve() == dst.resolve():
        return
    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(src, dst)


def _detect_mosaic_entries(
    env: dict[str, str], workdir: Path
) -> list[tuple[str, Path]]:
    out: list[tuple[str, Path]] = []
    spec = env.get("DETECT_QC_BANDS", "") or "soft:200:1000;hard:1000:12000"
    for label, _pimin, _pimax in parse_image_bands(spec):
        path = workdir / "detect" / f"EPIC_{label}_mosaic.fits"
        if path.exists():
            out.append((label, path))
    fallback = workdir / "detect" / "EPIC_hard_mosaic.fits"
    if not out and fallback.exists():
        out.append(("hard", fallback))
    return out


def _extent_from_detect_header(
    header: fits.Header, shape: tuple[int, int]
) -> tuple[tuple[float, float, float, float] | None, float | None, float | None]:
    if all(
        (
            k in header
            for k in ["XMIN_AS", "XMAX_AS", "YMIN_AS", "YMAX_AS", "CTRRA", "CTRDEC"]
        )
    ):
        extent = (
            float(header["XMIN_AS"]),
            float(header["XMAX_AS"]),
            float(header["YMIN_AS"]),
            float(header["YMAX_AS"]),
        )
        return (extent, float(header["CTRRA"]), float(header["CTRDEC"]))
    try:
        wcs = WCS(header)
        if wcs.has_celestial:
            ny, nx = shape
            pts = np.array(
                [[0.0, 0.0], [nx - 1.0, ny - 1.0], [0.5 * (nx - 1.0), 0.5 * (ny - 1.0)]]
            )
            sky = wcs.celestial.pixel_to_world(pts[:, 0], pts[:, 1])
            ctr_ra = float(np.asarray(sky.ra.deg)[-1])
            ctr_dec = float(np.asarray(sky.dec.deg)[-1])
            corners = wcs.celestial.pixel_to_world(
                [0.0, nx - 1.0, nx - 1.0, 0.0], [0.0, 0.0, ny - 1.0, ny - 1.0]
            )
            dx, dy = sky_offsets_arcsec(
                np.asarray(corners.ra.deg), np.asarray(corners.dec.deg), ctr_ra, ctr_dec
            )
            extent = (
                float(np.nanmin(dx)),
                float(np.nanmax(dx)),
                float(np.nanmin(dy)),
                float(np.nanmax(dy)),
            )
            return (extent, ctr_ra, ctr_dec)
    except Exception:
        pass
    return (None, None, None)


def _read_contam_timeline(path: Path) -> dict[str, np.ndarray]:
    with path.open("r", newline="", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh))
    if not rows:
        raise RuntimeError(f"Empty contamination timeline: {path}")
    out: dict[str, np.ndarray] = {}
    cols_float = ["time_s", "mjd_utc", "n_src_overlap", "n_bkg_overlap"]
    cols_bool = ["good_full", "good_src", "good_bkg", "good_strict"]
    for name in cols_float:
        out[name] = np.asarray([float(r[name]) for r in rows], dtype=float)
    for name in cols_bool:
        out[name] = np.asarray([int(r[name]) > 0 for r in rows], dtype=bool)
    out["src_ids"] = np.asarray([r.get("src_ids", "") for r in rows], dtype=object)
    out["bkg_ids"] = np.asarray([r.get("bkg_ids", "") for r in rows], dtype=object)
    return out


def _interp_track_to_times(
    track: TrackData, times: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    sc = SkyCoord(ra=track.ra * u.deg, dec=track.dec * u.deg)
    xyz = sc.cartesian.xyz.value.T.astype(float)
    t_sec = track.t_sec
    x = np.interp(times, t_sec, xyz[:, 0])
    y = np.interp(times, t_sec, xyz[:, 1])
    z = np.interp(times, t_sec, xyz[:, 2])
    n = np.sqrt(x * x + y * y + z * z)
    n[n == 0] = 1.0
    x /= n
    y /= n
    z /= n
    trk_ra = np.degrees(np.arctan2(y, x)) % 360.0
    trk_dec = np.degrees(np.arctan2(z, np.sqrt(x * x + y * y)))
    return (trk_ra.astype(float), trk_dec.astype(float))


def _read_report_rows(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open("r", newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def _instrument_rate_image(banddir: Path, inst: str) -> np.ndarray | None:
    cpath = banddir / f"{inst}_counts.fits"
    epath = banddir / f"{inst}_exp.fits"
    if not (cpath.exists() and epath.exists()):
        return None
    counts, _ = read_image(str(cpath))
    expo, _ = read_image(str(epath))
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where(expo > 0, counts / expo, np.nan)


def _read_events_header(path: Path) -> dict[str, Any]:
    info: dict[str, Any] = {}
    if not path.exists():
        return info
    try:
        with fits.open(path) as hdu:
            for ext in hdu:
                if ext.name == "EVENTS":
                    info["ontime"] = float(ext.header.get("ONTIME", 0))
                    info["livetime"] = float(ext.header.get("LIVETIME", 0))
                    info["filter"] = str(ext.header.get("FILTER", "")).strip()
                    info["tstart"] = float(ext.header.get("TSTART", 0))
                    info["tstop"] = float(ext.header.get("TSTOP", 0))
                    info["rows"] = len(ext.data) if ext.data is not None else 0
                    break
    except Exception:
        pass
    return info


def _read_espfilt_lc(espdir: Path) -> tuple[np.ndarray, np.ndarray] | None:
    candidates = sorted(espdir.glob("*fovlc.fits"))
    if not candidates:
        return None
    try:
        with fits.open(candidates[0]) as hdu:
            for ext in hdu:
                if ext.data is not None and hasattr(ext, "columns"):
                    cols = ext.columns.names
                    if "TIME" in cols and "RATE" in cols:
                        return (ext.data["TIME"].copy(), ext.data["RATE"].copy())
    except Exception:
        pass
    return None


def _read_espfilt_gti(espdir: Path) -> list[tuple[float, float]]:
    candidates = sorted(espdir.glob("*gti.fits"))
    intervals: list[tuple[float, float]] = []
    if not candidates:
        return intervals
    try:
        with fits.open(candidates[0]) as hdu:
            for ext in hdu:
                if ext.data is not None and hasattr(ext, "columns"):
                    cols = ext.columns.names
                    if "START" in cols and "STOP" in cols:
                        for s, e in zip(ext.data["START"], ext.data["STOP"]):
                            intervals.append((float(s), float(e)))
                        break
    except Exception:
        pass
    return intervals


def check_clean(env: dict[str, str], outdir: Path) -> dict[str, Any]:
    workdir = Path(env["WORKDIR"])
    inst_map = {"PN": "EPN", "M1": "EMOS1", "M2": "EMOS2"}
    skip_espfilt = env.get("CLEAN_SKIP_ESPFILT", "no").strip().lower() in (
        "1",
        "y",
        "yes",
        "true",
    )
    stats: list[dict[str, Any]] = []
    all_exposures: list[dict[str, Any]] = []
    any_clean = False
    for inst in ("PN", "M1", "M2"):
        raw_files = read_list(workdir / "repro" / "manifest" / f"{inst}_raw.txt")
        clean_files = read_list(workdir / "clean" / f"{inst}_clean_files.txt")
        clean_set = {Path(f).stem.replace(".clean", "") for f in clean_files}
        rows = []
        inst_raw_ontime = 0.0
        inst_clean_ontime = 0.0
        inst_closed_ontime = 0.0
        inst_flare_ontime = 0.0
        for rf in raw_files:
            rfp = Path(rf)
            raw_info = _read_events_header(rfp)
            raw_ontime = raw_info.get("ontime", 0.0)
            raw_filter = raw_info.get("filter", "")
            raw_rows = raw_info.get("rows", 0)
            stem = rfp.stem
            clean_path = workdir / "clean" / inst / f"{stem}.clean.fits"
            clean_info = _read_events_header(clean_path)
            clean_ontime = clean_info.get("ontime", 0.0)
            clean_rows = clean_info.get("rows", 0)
            has_clean = clean_path.exists()
            espdir = workdir / "clean" / inst / f"{stem}.clean.espfilt"
            per_exp_gti = workdir / "clean" / inst / f"{stem}.clean.flare_gti.fits"
            per_exp_gti_time = 0.0
            if per_exp_gti.exists():
                try:
                    with fits.open(per_exp_gti) as hdu:
                        for ext in hdu:
                            if ext.data is not None and hasattr(ext, "columns"):
                                cols = ext.columns.names
                                if "START" in cols and "STOP" in cols:
                                    for s, e in zip(
                                        ext.data["START"], ext.data["STOP"]
                                    ):
                                        dt = float(e) - float(s)
                                        if 0 < dt < 1000000000.0:
                                            per_exp_gti_time += dt
                                    break
                except Exception:
                    per_exp_gti_time = 0.0
            is_closed = raw_filter.lower() == "closed"
            if is_closed:
                status = "closed"
                inst_closed_ontime += raw_ontime
            elif has_clean:
                status = "cleaned"
                if per_exp_gti_time > 0:
                    inst_clean_ontime += per_exp_gti_time
                    inst_flare_ontime += max(0, raw_ontime - per_exp_gti_time)
                elif per_exp_gti.exists():
                    inst_flare_ontime += raw_ontime
                elif skip_espfilt:
                    inst_clean_ontime += clean_ontime
                    inst_flare_ontime += max(0, raw_ontime - clean_ontime)
                else:
                    inst_flare_ontime += raw_ontime
                rows.append(clean_rows)
            else:
                status = "skipped"
                inst_flare_ontime += raw_ontime
            parts = stem.split("_")
            short = parts[3] if len(parts) > 3 else stem
            inst_raw_ontime += raw_ontime
            if per_exp_gti_time > 0:
                display_clean_ontime = per_exp_gti_time
            elif per_exp_gti.exists() or not skip_espfilt:
                display_clean_ontime = 0.0
            else:
                display_clean_ontime = clean_ontime
            all_exposures.append(
                {
                    "inst": inst,
                    "short": short,
                    "status": status,
                    "filter": raw_filter,
                    "raw_ontime": raw_ontime,
                    "clean_ontime": display_clean_ontime if has_clean else 0.0,
                    "raw_rows": raw_rows,
                    "clean_rows": clean_rows if has_clean else 0,
                    "espdir": espdir,
                }
            )
        merged_rows = fits_rows(workdir / "clean" / f"{inst}_clean_merged.fits")
        stats.append(
            {
                "inst": inst,
                "n_raw": len(raw_files),
                "n_clean": len(clean_files),
                "skipped": max(0, len(raw_files) - len(clean_files)),
                "rows": rows,
                "merged_rows": int(merged_rows) if merged_rows is not None else 0,
                "raw_ontime": inst_raw_ontime,
                "clean_ontime": inst_clean_ontime,
                "closed_ontime": inst_closed_ontime,
                "flare_ontime": inst_flare_ontime,
            }
        )
        any_clean = any_clean or bool(clean_files)
    if not any_clean:
        return {"status": "missing", "reason": "No clean event lists found"}
    fig, axes = plt.subplots(2, 2, figsize=(10.5, 8.2), constrained_layout=True)
    insts = [s["inst"] for s in stats]
    x = np.arange(len(insts), dtype=float)
    raw_counts = np.asarray([s["n_raw"] for s in stats], dtype=float)
    clean_counts = np.asarray([s["n_clean"] for s in stats], dtype=float)
    merged = np.asarray([s["merged_rows"] for s in stats], dtype=float)
    ax = axes[0, 0]
    ax.bar(x - 0.18, raw_counts, width=0.35, label="raw event lists")
    ax.bar(x + 0.18, clean_counts, width=0.35, label="clean event lists")
    ax.set_xticks(x)
    ax.set_xticklabels(insts)
    ax.set_ylabel("Count")
    ax.set_title("Raw vs cleaned event-list multiplicity")
    ax.legend(loc="best", fontsize=8)
    ax = axes[0, 1]
    for s in stats:
        rr = np.asarray(sorted(s["rows"], reverse=True), dtype=float)
        if rr.size == 0:
            continue
        ax.plot(
            np.arange(1, len(rr) + 1), rr, marker="o", linestyle="-", label=s["inst"]
        )
    ax.set_xlabel("Clean event list rank")
    ax.set_ylabel("Rows after cleaning")
    ax.set_title("Per-exposure cleaned event rows")
    ax.legend(loc="best", fontsize=8)
    ax.set_yscale("log")
    ax = axes[1, 0]
    ax.bar(x, merged)
    ax.set_xticks(x)
    ax.set_xticklabels(insts)
    ax.set_ylabel("Rows")
    ax.set_title("Merged clean event-list rows")
    ax.set_yscale("log")
    ax = axes[1, 1]
    lines = []
    for s in stats:
        med = float(np.median(s["rows"])) if s["rows"] else math.nan
        lines.append(
            f"{s['inst']}: raw={s['n_raw']} clean={s['n_clean']} skipped={s['skipped']} merged_rows={s['merged_rows']}"
        )
        if np.isfinite(med):
            lines.append(f"    median clean rows={med:.1f}")
    ax.text(0.02, 0.98, "\n".join(lines), va="top", family="monospace", fontsize=8)
    ax.set_axis_off()
    ax.set_title("Cleaning summary")
    png = outdir / "00_clean.png"
    savefig(fig, str(png))
    s_exposures = [
        e for e in all_exposures if e["short"].startswith("S") and e["espdir"].exists()
    ]
    n_lc = min(len(s_exposures), 3)
    n_rows_fig2 = 2 + (1 if n_lc > 0 else 0)
    fig2 = plt.figure(figsize=(13, 4.5 * n_rows_fig2), constrained_layout=True)
    if n_lc > 0:
        from matplotlib.gridspec import GridSpec

        gs = GridSpec(n_rows_fig2, max(n_lc, 2), figure=fig2)
        axes2_00 = fig2.add_subplot(gs[0, : max(n_lc, 2) // 2])
        axes2_01 = fig2.add_subplot(gs[0, max(n_lc, 2) // 2 :])
        axes2_10 = fig2.add_subplot(gs[1, :])
        rate_axes = [fig2.add_subplot(gs[2, i]) for i in range(n_lc)]
    else:
        gs = GridSpec(n_rows_fig2, 2, figure=fig2)
        axes2_00 = fig2.add_subplot(gs[0, 0])
        axes2_01 = fig2.add_subplot(gs[0, 1])
        axes2_10 = fig2.add_subplot(gs[1, :])
        rate_axes = []
    ax = axes2_00
    clean_t = np.array([s["clean_ontime"] for s in stats]) / 1000.0
    flare_t = np.array([s["flare_ontime"] for s in stats]) / 1000.0
    closed_t = np.array([s["closed_ontime"] for s in stats]) / 1000.0
    ax.bar(x, clean_t, width=0.55, label="Good time", color="steelblue")
    ax.bar(
        x, flare_t, width=0.55, bottom=clean_t, label="Flare-filtered", color="tomato"
    )
    ax.bar(
        x,
        closed_t,
        width=0.55,
        bottom=clean_t + flare_t,
        label="Closed filter",
        color="gray",
        alpha=0.6,
    )
    ax.set_xticks(x)
    ax.set_xticklabels(insts)
    ax.set_ylabel("Time (ks)")
    ax.set_title("Exposure time budget")
    ax.legend(loc="best", fontsize=8)
    for i, s in enumerate(stats):
        if s["raw_ontime"] > 0:
            pct = 100 * s["clean_ontime"] / s["raw_ontime"]
            ax.text(
                i,
                (clean_t[i] + flare_t[i] + closed_t[i]) * 1.02,
                f"{pct:.0f}%",
                ha="center",
                va="bottom",
                fontsize=9,
                fontweight="bold",
            )
    ax = axes2_01
    sorted_exp = sorted(all_exposures, key=lambda e: (e["inst"], -e["raw_ontime"]))
    labels = []
    raw_vals = []
    clean_vals = []
    colors = []
    inst_colors = {"PN": "C0", "M1": "C1", "M2": "C2"}
    for e in sorted_exp:
        labels.append(f"{e['inst']} {e['short']}")
        raw_vals.append(e["raw_ontime"] / 1000.0)
        clean_vals.append(e["clean_ontime"] / 1000.0)
        colors.append(inst_colors.get(e["inst"], "C3"))
    y_pos = np.arange(len(labels))
    ax.barh(y_pos, raw_vals, height=0.7, color="lightgray", label="Raw ONTIME")
    for i, e in enumerate(sorted_exp):
        c = inst_colors.get(e["inst"], "C3")
        if e["status"] == "closed":
            ax.barh(i, raw_vals[i], height=0.7, color="gray", alpha=0.5)
        elif e["clean_ontime"] > 0:
            ax.barh(i, clean_vals[i], height=0.7, color=c, alpha=0.8)
    ax.set_yticks(y_pos)
    ax.set_yticklabels(labels, fontsize=8)
    ax.set_xlabel("ONTIME (ks)")
    ax.set_title("Per-exposure time: raw (gray) vs clean (colour)")
    ax.invert_yaxis()
    ax = axes2_10
    tlines = []
    tlines.append(
        f"{'Inst':<4} {'Exp':<6} {'Status':<8} {'Filter':<8} {'Raw(ks)':>8} {'Clean(ks)':>9} {'Retained':>8} {'Raw rows':>10} {'Clean rows':>10}"
    )
    tlines.append("-" * 95)
    for e in sorted_exp:
        raw_ks = f"{e['raw_ontime'] / 1000.0:.1f}"
        clean_ks = (
            f"{e['clean_ontime'] / 1000.0:.1f}" if e["status"] == "cleaned" else "--"
        )
        if e["status"] == "cleaned" and e["raw_ontime"] > 0:
            pct = f"{100 * e['clean_ontime'] / e['raw_ontime']:.0f}%"
        else:
            pct = "--"
        cr = str(e["clean_rows"]) if e["status"] == "cleaned" else "--"
        tlines.append(
            f"{e['inst']:<4} {e['short']:<6} {e['status']:<8} {e['filter']:<8} {raw_ks:>8} {clean_ks:>9} {pct:>8} {e['raw_rows']:>10} {cr:>10}"
        )
    tlines.append("")
    for s in stats:
        if s["raw_ontime"] > 0:
            pct = 100 * s["clean_ontime"] / s["raw_ontime"]
            tlines.append(
                f"{s['inst']}: {s['raw_ontime'] / 1000.0:.1f} ks raw -> {s['clean_ontime'] / 1000.0:.1f} ks clean = {pct:.0f}% retained (flare: {s['flare_ontime'] / 1000.0:.1f} ks, closed: {s['closed_ontime'] / 1000.0:.1f} ks)"
            )
    ax.text(
        0.01,
        0.99,
        "\n".join(tlines),
        va="top",
        family="monospace",
        fontsize=7.5,
        transform=ax.transAxes,
    )
    ax.set_axis_off()
    ax.set_title("Per-exposure cleaning detail", fontsize=10)
    if n_lc > 0:
        for i_lc, exp in enumerate(s_exposures[:n_lc]):
            ax = rate_axes[i_lc]
            lc = _read_espfilt_lc(exp["espdir"])
            gtis = _read_espfilt_gti(exp["espdir"])
            if lc is not None:
                times, rates = lc
                t0 = times.min()
                t_ks = (times - t0) / 1000.0
                ax.plot(t_ks, rates, linewidth=0.3, color="gray", alpha=0.7)
                bin_s = 100
                n_bins = max(1, int((times.max() - times.min()) / bin_s))
                if n_bins > 5:
                    edges = np.linspace(times.min(), times.max(), n_bins + 1)
                    idx = np.digitize(times, edges) - 1
                    idx = np.clip(idx, 0, n_bins - 1)
                    binned = np.zeros(n_bins)
                    counts = np.zeros(n_bins)
                    for j in range(len(times)):
                        binned[idx[j]] += rates[j]
                        counts[idx[j]] += 1
                    mask = counts > 0
                    binned[mask] /= counts[mask]
                    t_mid = 0.5 * (edges[:-1] + edges[1:])
                    ax.plot((t_mid - t0) / 1000.0, binned, linewidth=0.8, color="k")
                for gs, ge in gtis:
                    ax.axvspan(
                        (gs - t0) / 1000.0,
                        (ge - t0) / 1000.0,
                        color="steelblue",
                        alpha=0.15,
                    )
                ax.set_xlabel("Time since start (ks)")
                ax.set_ylabel("FOV rate (ct/s)")
                title_pct = ""
                if exp["raw_ontime"] > 0 and exp["clean_ontime"] > 0:
                    title_pct = f" — {100 * exp['clean_ontime'] / exp['raw_ontime']:.0f}% retained"
                ax.set_title(
                    f"{exp['inst']} {exp['short']} espfilt rate{title_pct}", fontsize=9
                )
            else:
                ax.text(
                    0.5,
                    0.5,
                    "No rate curve",
                    ha="center",
                    va="center",
                    transform=ax.transAxes,
                )
                ax.set_title(f"{exp['inst']} {exp['short']}", fontsize=9)
    png2 = outdir / "00_clean_time.png"
    savefig(fig2, str(png2))
    retained_pcts = []
    for s in stats:
        if s["raw_ontime"] > 0:
            retained_pcts.append(
                f"{s['inst']}:{100 * s['clean_ontime'] / s['raw_ontime']:.0f}%"
            )
    return {
        "status": "ok",
        "png": png.name,
        "extra_pngs": [png2.name],
        "per_instrument": ", ".join(
            (f"{s['inst']}:{s['n_clean']}/{s['n_raw']}" for s in stats)
        ),
        "time_retained": ", ".join(retained_pcts),
    }


def _read_exposure_tsv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    import csv as _csv

    rows: list[dict[str, str]] = []
    with path.open("r", encoding="utf-8") as fh:
        reader = _csv.DictReader(fh, delimiter="\t")
        for row in reader:
            rows.append(dict(row))
    return rows


def _quick_look_image(
    evt_path: str, pi_min: int = 200, pi_max: int = 12000
) -> np.ndarray | None:
    try:
        with fits.open(evt_path, memmap=True) as hdul:
            ext = hdul["EVENTS"] if "EVENTS" in hdul else hdul[1]
            data = ext.data
            if data is None or len(data) == 0:
                return None
            pi = np.asarray(data["PI"], dtype=float)
            mask = np.isfinite(pi) & (pi >= pi_min) & (pi <= pi_max)
            if not np.any(mask):
                return None
            if "X" in data.columns.names and "Y" in data.columns.names:
                x = np.asarray(data["X"][mask], dtype=float)
                y = np.asarray(data["Y"][mask], dtype=float)
            elif "DETX" in data.columns.names and "DETY" in data.columns.names:
                x = np.asarray(data["DETX"][mask], dtype=float)
                y = np.asarray(data["DETY"][mask], dtype=float)
            else:
                return None
            valid = np.isfinite(x) & np.isfinite(y)
            x, y = (x[valid], y[valid])
            if x.size < 10:
                return None
            nbins = min(256, max(32, int(np.sqrt(x.size) / 2)))
            img, _, _ = np.histogram2d(y, x, bins=nbins)
            return img.astype(float)
    except Exception:
        return None


def _event_rate_histogram(
    evt_path: str, nbins: int = 100
) -> tuple[np.ndarray, np.ndarray] | None:
    try:
        with fits.open(evt_path, memmap=True) as hdul:
            ext = hdul["EVENTS"] if "EVENTS" in hdul else hdul[1]
            data = ext.data
            if data is None or len(data) == 0:
                return None
            t = np.asarray(data["TIME"], dtype=float)
            t = t[np.isfinite(t)]
            if t.size < 2:
                return None
            t0, t1 = (float(t.min()), float(t.max()))
            if t1 - t0 < 1.0:
                return None
            edges = np.linspace(t0, t1, nbins + 1)
            counts, _ = np.histogram(t, bins=edges)
            dt = float(edges[1] - edges[0])
            rate = counts.astype(float) / dt
            midpoints = 0.5 * (edges[:-1] + edges[1:]) - t0
            return (midpoints, rate)
    except Exception:
        return None


def check_exposures(env: dict[str, str], outdir: Path) -> dict[str, Any]:
    workdir = Path(env["WORKDIR"])
    tsv_path = workdir / "detect" / "all_pseudoexposures.tsv"
    exposures = _read_exposure_tsv(tsv_path)
    if not exposures:
        return {
            "status": "missing",
            "reason": "No all_pseudoexposures.tsv found; re-run detect stage",
        }
    for exp in exposures:
        for k in ("events", "gti_rows"):
            try:
                exp[k] = int(exp.get(k, 0))
            except (ValueError, TypeError):
                exp[k] = 0
        for k in ("tstart", "tstop", "ontime", "livetime"):
            try:
                exp[k] = float(exp.get(k, "nan"))
            except (ValueError, TypeError):
                exp[k] = float("nan")
    accepted = [e for e in exposures if e.get("status") == "accepted"]
    skipped = [e for e in exposures if e.get("status") == "skipped"]
    skipped_with_events = [e for e in skipped if e["events"] > 0]
    n_ql = len(skipped_with_events)
    n_ql_rows = max(0, (n_ql + 3) // 4)
    n_rows = 2 + n_ql_rows
    fig_height = 4.5 * n_rows
    fig, axes = plt.subplots(
        n_rows, 4, figsize=(18.0, fig_height), constrained_layout=True
    )
    if n_rows == 1:
        axes = axes[np.newaxis, :]
    ax_tl = fig.add_subplot(n_rows, 2, 1)
    inst_colors = {
        "EMOS1": "C0",
        "EMOS2": "C1",
        "EPN": "C2",
        "M1": "C0",
        "M2": "C1",
        "PN": "C2",
    }
    # Per-CCD ontime bar chart grouped by instrument.  Pseudoexposures are
    # per-CCD spatial subsets sharing the same tstart/tstop, so a timeline
    # bar chart is misleading — instead show the GTI ontime for each CCD.
    ccd_labels: list[str] = []
    ccd_ontime: list[float] = []
    ccd_events: list[int] = []
    ccd_colors: list[str] = []
    ccd_status: list[str] = []
    for e in exposures:
        ccd_labels.append(e.get("filename", "?").replace(".ds", ""))
        ccd_ontime.append(e["ontime"] if np.isfinite(e["ontime"]) else 0.0)
        ccd_events.append(int(e["events"]))
        ccd_colors.append(inst_colors.get(e.get("instrument", "?"), "gray"))
        ccd_status.append(e.get("status", "?"))
    y_pos = np.arange(len(ccd_labels))
    bars = ax_tl.barh(
        y_pos,
        ccd_ontime,
        color=ccd_colors,
        alpha=0.8,
        edgecolor="black",
        linewidth=0.5,
    )
    # Hatch skipped entries
    for i, (bar, st) in enumerate(zip(bars, ccd_status)):
        if st == "skipped":
            bar.set_alpha(0.3)
            bar.set_hatch("///")
    ax_tl.set_yticks(y_pos)
    ax_tl.set_yticklabels(ccd_labels, fontsize=6)
    ax_tl.set_xlabel("GTI ontime (s)", fontsize=9)
    ax_tl.set_title("Per-CCD ontime (pseudoexposures)", fontsize=10)
    ax_tl.invert_yaxis()
    from matplotlib.patches import Patch

    legend_els = [
        Patch(facecolor="C0", alpha=0.8, label="EMOS1 accepted"),
        Patch(facecolor="C1", alpha=0.8, label="EMOS2 accepted"),
        Patch(facecolor="C2", alpha=0.8, label="EPN accepted"),
        Patch(facecolor="C0", alpha=0.3, hatch="///", label="skipped"),
    ]
    ax_tl.legend(handles=legend_els, loc="lower right", fontsize=7)
    ax_meta = fig.add_subplot(n_rows, 2, 2)
    table_lines = [
        f"{'filename':<26s} {'status':<8s} {'inst':<5s} {'submode':<10s} {'datamode':<10s} {'filter':<7s} {'events':>8s} {'gti':>4s} {'ontime':>9s}",
        "-" * 90,
    ]
    for e in exposures:
        ontime_str = f"{e['ontime']:.1f}" if np.isfinite(e["ontime"]) else "?"
        table_lines.append(
            f"{e.get('filename', '?'):<26s} {e.get('status', '?'):<8s} {e.get('instrument', '?'):<5s} {e.get('submode', ''):<10s} {e.get('datamode', ''):<10s} {e.get('filter', ''):<7s} {e['events']:>8d} {e['gti_rows']:>4d} {ontime_str:>9s}"
        )
    ax_meta.text(
        0.01,
        0.99,
        "\n".join(table_lines),
        va="top",
        ha="left",
        family="monospace",
        fontsize=6,
        transform=ax_meta.transAxes,
    )
    ax_meta.set_axis_off()
    ax_meta.set_title("Pseudo-exposure metadata", fontsize=10)
    for ax in axes[0]:
        ax.set_visible(False)
    if skipped_with_events:
        n_rate = min(4, len(skipped_with_events))
        for idx in range(4):
            ax = axes[1][idx]
            if idx < n_rate:
                e = skipped_with_events[idx]
                evt_path = e.get("path", "")
                rh = _event_rate_histogram(evt_path, nbins=80)
                if rh is not None:
                    midpts, rate = rh
                    ax.step(midpts, rate, where="mid", linewidth=0.8)
                    ax.set_xlabel("Time (s from start)")
                    ax.set_ylabel("Rate (ct/s)")
                    ax.set_title(
                        f"{e.get('filename', '?')}\nevents={e['events']}", fontsize=8
                    )
                else:
                    ax.text(
                        0.5,
                        0.5,
                        f"{e.get('filename', '?')}\nno TIME data",
                        ha="center",
                        va="center",
                        fontsize=8,
                    )
                    ax.set_axis_off()
            else:
                ax.set_axis_off()
        if len(skipped_with_events) > 4:
            axes[1][3].text(
                0.5,
                0.02,
                f"+{len(skipped_with_events) - 3} more skipped exposures with events",
                ha="center",
                fontsize=7,
                transform=axes[1][3].transAxes,
            )
    else:
        for ax in axes[1]:
            ax.set_axis_off()
        axes[1][0].text(
            0.5,
            0.5,
            "No skipped exposures with events > 0",
            ha="center",
            va="center",
            fontsize=10,
        )
    for row_idx in range(n_ql_rows):
        for col_idx in range(4):
            ax = axes[2 + row_idx][col_idx]
            ql_idx = row_idx * 4 + col_idx
            if ql_idx < n_ql:
                e = skipped_with_events[ql_idx]
                img = _quick_look_image(e.get("path", ""))
                if img is not None:
                    im = ax.imshow(
                        img, origin="lower", norm=nice_norm(img, stretch="sqrt")
                    )
                    plt.colorbar(im, ax=ax, pad=0.01, shrink=0.8)
                    ax.set_title(
                        f"{e.get('filename', '?')}\n{e.get('instrument', '?')} {e.get('submode', '')} {e.get('datamode', '')}",
                        fontsize=7,
                    )
                else:
                    ax.text(
                        0.5,
                        0.5,
                        f"{e.get('filename', '?')}\nno image data",
                        ha="center",
                        va="center",
                        fontsize=8,
                    )
                    ax.set_axis_off()
            else:
                ax.set_axis_off()
    fig.suptitle(
        f"Pseudo-exposure diagnostics: {len(accepted)} accepted, {len(skipped)} skipped ({len(skipped_with_events)} with events)",
        y=1.0,
        fontsize=12,
    )
    png = outdir / "00b_exposures.png"
    savefig(fig, str(png))
    return {
        "status": "ok",
        "png": png.name,
        "n_total": len(exposures),
        "n_accepted": len(accepted),
        "n_skipped": len(skipped),
        "n_skipped_with_events": len(skipped_with_events),
    }


def check_track(
    env: dict[str, str], outdir: Path, use_horizons: bool
) -> dict[str, Any]:
    workdir = Path(env["WORKDIR"])
    track_path = workdir / "track" / "comet_track.fits"
    if not track_path.exists():
        return {"status": "missing", "reason": str(track_path)}
    track = read_track(str(track_path))
    ref_ra = float(track.ra[len(track.ra) // 2])
    ref_dec = float(track.dec[len(track.dec) // 2])
    dx, dy = sky_offsets_arcsec(track.ra, track.dec, ref_ra, ref_dec)
    total_motion = 0.0
    if len(track.ra) > 1:
        sep = (
            SkyCoord(ra=track.ra[:-1] * u.deg, dec=track.dec[:-1] * u.deg)
            .separation(SkyCoord(ra=track.ra[1:] * u.deg, dec=track.dec[1:] * u.deg))
            .to_value(u.arcsec)
        )
        total_motion = float(np.sum(sep))
    duration = float(track.obs_t1 - track.obs_t0)
    step = float(np.median(np.diff(track.t_sec))) if len(track.t_sec) > 1 else math.nan
    seg_csv = workdir / "track" / "motion_segments.csv"
    segs = read_motion_segments(str(seg_csv)) if seg_csv.exists() else []
    horiz_stats = compare_with_horizons(env, track) if use_horizons else None
    src_path = workdir / "detect" / "field_sources_all.fits"
    if not src_path.exists():
        src_path = workdir / "detect" / "field_sources_curated.fits"
    nearest_source = math.nan
    if src_path.exists():
        src = read_srclist(str(src_path))
        if len(src.ra):
            comet = SkyCoord(ra=track.ra * u.deg, dec=track.dec * u.deg)
            field = SkyCoord(ra=src.ra * u.deg, dec=src.dec * u.deg)
            nearest_source = float(
                np.min(comet[:, None].separation(field[None, :]).to_value(u.arcsec))
            )
    fig, axes = plt.subplots(2, 2, figsize=(10.5, 8.0))
    t_ks = (track.t_sec - track.t_sec[0]) / 1000.0
    ax = axes[0, 0]
    ax.plot(dx, dy)
    ax.plot(dx[0], dy[0], marker="o", linestyle="None")
    ax.plot(dx[-1], dy[-1], marker="s", linestyle="None")
    ax.set_xlabel("East offset from reference (arcsec)")
    ax.set_ylabel("North offset from reference (arcsec)")
    ax.set_title("Stored comet track on the sky")
    ax.set_aspect("equal", adjustable="box")
    ax = axes[0, 1]
    ax.plot(t_ks, dx, label="dRA cosDec")
    ax.plot(t_ks, dy, label="dDec")
    ax.set_xlabel("Time from start (ks)")
    ax.set_ylabel("Offset (arcsec)")
    ax.set_title("Track offsets vs time")
    ax.legend(loc="best", fontsize=8)
    ax = axes[1, 0]
    if segs:
        tmid = 0.5 * (
            np.array([s["tstart"] for s in segs]) + np.array([s["tstop"] for s in segs])
        )
        tmid = (tmid - track.obs_t0) / 1000.0
        blur = np.array([s["max_mid_offset_arcsec"] for s in segs], dtype=float)
        dur = np.array([s["tstop"] - s["tstart"] for s in segs], dtype=float)
        ax.plot(tmid, blur, label="segment blur")
        ax2 = ax.twinx()
        ax2.plot(tmid, dur, linestyle="--", label="segment duration")
        ax.set_xlabel("Segment midpoint from start (ks)")
        ax.set_ylabel("Max midpoint blur (arcsec)")
        ax2.set_ylabel("Duration (s)")
        ax.set_title("Motion-segmentation diagnostics")
    else:
        ax.text(0.5, 0.5, "motion_segments.csv not found", ha="center", va="center")
        ax.set_axis_off()
    ax = axes[1, 1]
    text_lines = [
        f"samples            : {len(track.mjd)}",
        f"obs span           : {duration:.1f} s",
        (
            f"median track step  : {step:.3f} s"
            if np.isfinite(step)
            else "median track step  : n/a"
        ),
        f"total sky motion   : {total_motion:.2f} arcsec",
        (
            f"nearest detected src: {nearest_source:.2f} arcsec"
            if np.isfinite(nearest_source)
            else "nearest detected src: n/a"
        ),
    ]
    if horiz_stats is not None:
        med_sep, max_sep = horiz_stats
        text_lines.extend(
            [
                f"Horizons median Δ  : {med_sep:.4f} arcsec",
                f"Horizons max Δ     : {max_sep:.4f} arcsec",
            ]
        )
    elif use_horizons:
        text_lines.append("Horizons compare    : unavailable")
    if segs:
        blur = [float(s["max_mid_offset_arcsec"]) for s in segs]
        text_lines.extend(
            [
                f"segments           : {len(segs)}",
                f"max seg blur       : {max(blur):.3f} arcsec",
            ]
        )
    ax.text(0.02, 0.98, "\n".join(text_lines), va="top", family="monospace", fontsize=9)
    ax.set_title("Track summary")
    ax.set_axis_off()
    png = outdir / "01_track.png"
    savefig(fig, str(png))
    return {
        "status": "ok",
        "png": png.name,
        "samples": len(track.mjd),
        "obs_span_s": duration,
        "median_track_step_s": step,
        "total_motion_arcsec": total_motion,
        "nearest_detected_source_arcsec": nearest_source,
        "horizons_median_arcsec": horiz_stats[0] if horiz_stats is not None else None,
        "horizons_max_arcsec": horiz_stats[1] if horiz_stats is not None else None,
    }


def check_detect(env: dict[str, str], outdir: Path) -> dict[str, Any]:
    workdir = Path(env["WORKDIR"])
    entries = _detect_mosaic_entries(env, workdir)
    if not entries:
        return {"status": "missing", "reason": "No detect diagnostic mosaics found"}
    src_all = None
    src_cur = None
    all_path = workdir / "detect" / "field_sources_all.fits"
    cur_path = workdir / "detect" / "field_sources_curated.fits"
    if all_path.exists():
        src_all = read_srclist(str(all_path))
    if cur_path.exists():
        src_cur = read_srclist(str(cur_path))
    track = (
        read_track(str(workdir / "track" / "comet_track.fits"))
        if (workdir / "track" / "comet_track.fits").exists()
        else None
    )
    fig, axes = plt.subplots(
        1,
        len(entries),
        figsize=(6.0 * len(entries), 5.2),
        squeeze=False,
        constrained_layout=True,
    )
    axes = axes[0]
    visible_stats: list[str] = []
    nearest_all = math.nan
    nearest_cur = math.nan
    for ax, (label, path) in zip(axes, entries):
        data, header = read_image(str(path))
        extent, ctr_ra, ctr_dec = _extent_from_detect_header(header, data.shape)
        norm = nice_norm(data, stretch="sqrt")
        if extent is not None:
            im = ax.imshow(data, origin="lower", extent=extent, norm=norm)
            ax.set_xlabel("East offset (arcsec)")
            ax.set_ylabel("North offset (arcsec)")
            if ctr_ra is not None and ctr_dec is not None:
                if src_all is not None and len(src_all.ra):
                    dx_a, dy_a = sky_offsets_arcsec(
                        src_all.ra, src_all.dec, ctr_ra, ctr_dec
                    )
                    ax.scatter(
                        dx_a,
                        dy_a,
                        s=18,
                        facecolors="none",
                        edgecolors="white",
                        linewidths=0.7,
                        label="all detected sources",
                    )
                    if track is not None:
                        comet = SkyCoord(ra=track.ra * u.deg, dec=track.dec * u.deg)
                        field = SkyCoord(ra=src_all.ra * u.deg, dec=src_all.dec * u.deg)
                        nearest_all = float(
                            np.min(
                                comet[:, None]
                                .separation(field[None, :])
                                .to_value(u.arcsec)
                            )
                        )
                if src_cur is not None and len(src_cur.ra):
                    dx_c, dy_c = sky_offsets_arcsec(
                        src_cur.ra, src_cur.dec, ctr_ra, ctr_dec
                    )
                    ax.scatter(
                        dx_c,
                        dy_c,
                        s=16,
                        marker="s",
                        facecolors="none",
                        edgecolors="orange",
                        linewidths=0.8,
                        label="track-excluded list",
                    )
                    if track is not None:
                        comet = SkyCoord(ra=track.ra * u.deg, dec=track.dec * u.deg)
                        field = SkyCoord(ra=src_cur.ra * u.deg, dec=src_cur.dec * u.deg)
                        nearest_cur = float(
                            np.min(
                                comet[:, None]
                                .separation(field[None, :])
                                .to_value(u.arcsec)
                            )
                        )
                if track is not None:
                    dx_t, dy_t = sky_offsets_arcsec(
                        track.ra, track.dec, ctr_ra, ctr_dec
                    )
                    visible = (
                        np.isfinite(dx_t)
                        & np.isfinite(dy_t)
                        & (dx_t >= extent[0])
                        & (dx_t <= extent[1])
                        & (dy_t >= extent[2])
                        & (dy_t <= extent[3])
                    )
                    visible_frac = float(np.mean(visible)) if len(visible) else math.nan
                    n_visible = int(np.count_nonzero(visible))
                    visible_stats.append(
                        f"{label}:{n_visible}/{len(visible)} ({visible_frac:.1%})"
                    )
                    first = True
                    for i0, i1 in split_visible_segments(visible):
                        ax.plot(
                            dx_t[i0:i1],
                            dy_t[i0:i1],
                            linewidth=1.2,
                            label="visible comet track" if first else None,
                        )
                        first = False
                    if np.any(visible):
                        idxs = np.flatnonzero(visible)
                        ax.scatter(
                            dx_t[idxs[0]],
                            dy_t[idxs[0]],
                            marker="o",
                            s=22,
                            label="visible start",
                        )
                        ax.scatter(
                            dx_t[idxs[-1]],
                            dy_t[idxs[-1]],
                            marker="s",
                            s=22,
                            label="visible end",
                        )
        else:
            im = ax.imshow(data, origin="lower", norm=norm)
            ax.set_xlabel("X")
            ax.set_ylabel("Y")
        plt.colorbar(im, ax=ax, pad=0.01, shrink=0.9, label="counts")
        ax.set_title(f"{label} still-sky mosaic")
        if ax is axes[0]:
            ax.legend(loc="upper right", fontsize=7)
    summary_lines = [
        f"bands rendered      : {', '.join((label for label, _ in entries))}",
        f"all detected srcs   : {(len(src_all.ra) if src_all is not None else 0)}",
        f"track-excluded srcs : {(len(src_cur.ra) if src_cur is not None else 0)}",
    ]
    if visible_stats:
        summary_lines.append("visible track frac   : " + "; ".join(visible_stats))
    if np.isfinite(nearest_all):
        summary_lines.append(f"nearest all-source sep: {nearest_all:.2f} arcsec")
    if np.isfinite(nearest_cur):
        summary_lines.append(f"nearest curated sep  : {nearest_cur:.2f} arcsec")
    fig.text(
        0.01, -0.04, "\n".join(summary_lines), family="monospace", fontsize=8, va="top"
    )
    fig.suptitle("Still-sky detect diagnostics", y=1.02)
    png = outdir / "02_detect.png"
    savefig(fig, str(png))
    return {
        "status": "ok",
        "png": png.name,
        "bands": ", ".join((label for label, _ in entries)),
        "n_detected_sources": len(src_all.ra) if src_all is not None else 0,
        "n_track_excluded_sources": len(src_cur.ra) if src_cur is not None else 0,
        "nearest_all_source_arcsec": nearest_all,
        "nearest_curated_source_arcsec": nearest_cur,
    }


def check_image(
    env: dict[str, str], outdir: Path, band_label: str | None = None
) -> dict[str, Any]:
    workdir = Path(env["WORKDIR"])
    bands = parse_image_bands(
        env.get("IMAGE_BANDS", "soft:200:1000;broad:300:2000;hard:1000:12000")
    )
    if band_label:
        bands = [b for b in bands if b[0] == band_label]
    if not bands:
        return {"status": "missing", "reason": "No matching IMAGE_BANDS entries"}
    grid_path = workdir / "images" / "grid.json"
    if not grid_path.exists():
        return {"status": "missing", "reason": str(grid_path)}
    grid = json.loads(grid_path.read_text(encoding="utf-8"))
    center_x = float(grid["center_x_phys"])
    center_y = float(grid["center_y_phys"])
    x_min = float(grid["x_min_phys"])
    x_max = float(grid["x_max_phys"])
    y_min = float(grid["y_min_phys"])
    y_max = float(grid["y_max_phys"])
    scale = float(grid["scale_arcsec_per_phys"])
    extent = (
        (x_min - center_x) * scale,
        (x_max - center_x) * scale,
        (y_min - center_y) * scale,
        (y_max - center_y) * scale,
    )
    rendered: list[Path] = []
    first_stats: dict[str, Any] | None = None
    for label, _pimin, _pimax in bands:
        banddir = workdir / "images" / label
        counts_path = banddir / "EPIC_counts.fits"
        exp_path = banddir / "EPIC_exp.fits"
        rate_path = banddir / "EPIC_rate.fits"
        clean_path = banddir / "EPIC_rate_clean.fits"
        net_path = banddir / "EPIC_net_rate.fits"
        masked_path = banddir / "EPIC_rate_masked.fits"
        mask_path = banddir / "EPIC_trail_mask.fits"
        if not (counts_path.exists() and exp_path.exists() and rate_path.exists()):
            continue
        counts, _ = read_image(str(counts_path))
        expo, _ = read_image(str(exp_path))
        rate, _ = read_image(str(rate_path))
        clean = None
        net = None
        masked = None
        trail_mask = None
        if clean_path.exists():
            clean, _ = read_image(str(clean_path))
        if net_path.exists():
            net, _ = read_image(str(net_path))
        if masked_path.exists():
            masked, _ = read_image(str(masked_path))
        if mask_path.exists():
            trail_mask, _ = read_image(str(mask_path))
        src_mask = circle_mask_from_extent(
            rate.shape,
            extent,
            env_float(env, "SRC_DX_ARCSEC", 0.0),
            env_float(env, "SRC_DY_ARCSEC", 0.0),
            env_float(env, "SRC_R_ARCSEC", 60.0),
        )
        if env.get("BKG_MODE", "annulus").strip().lower() == "circle":
            bkg_mask = circle_mask_from_extent(
                rate.shape,
                extent,
                env_float(env, "BKG_DX_ARCSEC", 0.0),
                env_float(env, "BKG_DY_ARCSEC", 0.0),
                env_float(env, "BKG_R_ARCSEC", 150.0),
            )
        else:
            bkg_mask = annulus_mask_from_extent(
                rate.shape,
                extent,
                env_float(env, "BKG_DX_ARCSEC", 0.0),
                env_float(env, "BKG_DY_ARCSEC", 0.0),
                env_float(env, "BKG_RIN_ARCSEC", 90.0),
                env_float(env, "BKG_ROUT_ARCSEC", 150.0),
            )
        src_exp_median = (
            float(np.nanmedian(expo[src_mask])) if np.any(src_mask) else math.nan
        )
        bkg_exp_median = (
            float(np.nanmedian(expo[bkg_mask])) if np.any(bkg_mask) else math.nan
        )
        src_rate_median = (
            float(np.nanmedian(rate[src_mask])) if np.any(src_mask) else math.nan
        )
        bkg_rate_median = (
            float(np.nanmedian(rate[bkg_mask])) if np.any(bkg_mask) else math.nan
        )
        zero_exp_frac = float(np.mean(~np.isfinite(expo) | (expo <= 0)))
        masked_frac = float(
            np.mean(
                trail_mask > 0
                if trail_mask is not None
                else np.zeros_like(rate, dtype=bool)
            )
        )
        fig, axes = plt.subplots(2, 4, figsize=(14.0, 7.8), constrained_layout=True)
        top_panels = [
            (counts, "EPIC counts", "sqrt"),
            (expo, "EPIC exposure", "sqrt"),
            (
                clean if clean is not None else rate,
                "EPIC clean rate" if clean is not None else "EPIC rate",
                "rate",
            ),
            (
                net if net is not None else masked if masked is not None else rate,
                (
                    "EPIC net rate"
                    if net is not None
                    else (
                        "EPIC masked rate"
                        if masked is not None
                        else "EPIC rate (repeat)"
                    )
                ),
                "symlog",
            ),
        ]
        for ax, (data, title, stretch) in zip(axes[0], top_panels):
            if stretch == "rate":
                # Rate image: use sqrt stretch, but limit scale to well-exposed
                # pixels to avoid FOV-edge artifacts dominating.
                exp_thresh = np.nanmax(expo) * 0.1 if np.any(np.isfinite(expo)) else 0
                good = np.isfinite(data) & (expo > exp_thresh)
                if np.any(good):
                    vmin = 0.0
                    vmax = float(np.nanpercentile(data[good], 99.5))
                    vmax = max(vmax, 1e-30)
                else:
                    vmin, vmax = (0.0, 1.0)
                from matplotlib.colors import PowerNorm

                im = ax.imshow(
                    data,
                    origin="lower",
                    extent=extent,
                    norm=PowerNorm(gamma=0.5, vmin=vmin, vmax=vmax),
                )
            elif stretch == "symlog":
                exp_thresh = np.nanmax(expo) * 0.1 if np.any(np.isfinite(expo)) else 0
                good = np.isfinite(data) & (expo > exp_thresh)
                if np.any(good):
                    absmax = np.nanpercentile(np.abs(data[good]), 99.0)
                    absmax = max(absmax, 1e-30)
                else:
                    absmax = 1.0
                im = ax.imshow(
                    data,
                    origin="lower",
                    extent=extent,
                    vmin=-absmax,
                    vmax=absmax,
                    cmap="RdBu_r",
                )
            elif stretch == "linear":
                finite = np.isfinite(data)
                if np.any(finite):
                    vmin, vmax = np.nanpercentile(data[finite], [1.0, 99.0])
                else:
                    vmin, vmax = (None, None)
                im = ax.imshow(
                    data, origin="lower", extent=extent, vmin=vmin, vmax=vmax
                )
            else:
                im = ax.imshow(
                    data,
                    origin="lower",
                    extent=extent,
                    norm=nice_norm(data, stretch=stretch),
                )
            plt.colorbar(im, ax=ax, pad=0.01, shrink=0.8)
            add_apertures(ax, env)
            ax.set_xlabel("East offset (arcsec)")
            ax.set_ylabel("North offset (arcsec)")
            ax.set_title(title)
            if trail_mask is not None and title.endswith("rate"):
                try:
                    ax.contour(
                        (trail_mask > 0).astype(float),
                        levels=[0.5],
                        origin="lower",
                        extent=extent,
                        linewidths=0.6,
                    )
                except Exception:
                    pass
        for ax, inst in zip(axes[1][:3], ["PN", "M1", "M2"]):
            inst_net_path = banddir / f"{inst}_net_rate.fits"
            if inst_net_path.exists():
                inst_rate, _ = read_image(str(inst_net_path))
                inst_title = f"{inst} net rate"
            else:
                inst_rate = _instrument_rate_image(banddir, inst)
                inst_title = f"{inst} rate"
            if inst_rate is None:
                ax.text(0.5, 0.5, f"{inst} missing", ha="center", va="center")
                ax.set_axis_off()
                continue
            finite = np.isfinite(inst_rate)
            # Per-instrument exposure for masking FOV edges
            inst_exp_path = banddir / f"{inst}_exp.fits"
            inst_expo = None
            if inst_exp_path.exists():
                inst_expo, _ = read_image(str(inst_exp_path))
            if "net" in inst_title.lower():
                if inst_expo is not None and np.any(np.isfinite(inst_expo)):
                    exp_thr = np.nanmax(inst_expo) * 0.1
                    good = finite & (inst_expo > exp_thr)
                else:
                    good = finite
                if np.any(good):
                    absmax = np.nanpercentile(np.abs(inst_rate[good]), 99.0)
                    absmax = max(absmax, 1e-30)
                else:
                    absmax = 1.0
                im = ax.imshow(
                    inst_rate,
                    origin="lower",
                    extent=extent,
                    vmin=-absmax,
                    vmax=absmax,
                    cmap="RdBu_r",
                )
            else:
                if np.any(finite):
                    vmin, vmax = np.nanpercentile(inst_rate[finite], [1.0, 99.0])
                else:
                    vmin, vmax = (None, None)
                im = ax.imshow(
                    inst_rate, origin="lower", extent=extent, vmin=vmin, vmax=vmax
                )
            plt.colorbar(im, ax=ax, pad=0.01, shrink=0.8)
            add_apertures(ax, env)
            ax.set_xlabel("East offset (arcsec)")
            ax.set_ylabel("North offset (arcsec)")
            ax.set_title(inst_title)
        ax = axes[1][3]
        if trail_mask is not None:
            im = ax.imshow(
                (trail_mask > 0).astype(float),
                origin="lower",
                extent=extent,
                vmin=0.0,
                vmax=1.0,
            )
            plt.colorbar(im, ax=ax, pad=0.01, shrink=0.8)
            add_apertures(ax, env)
            ax.set_xlabel("East offset (arcsec)")
            ax.set_ylabel("North offset (arcsec)")
            ax.set_title("Trail mask")
        else:
            text = f"band={label}\nsrc median exp  : {src_exp_median:.4g}\nbkg median exp  : {bkg_exp_median:.4g}\nsrc median rate : {src_rate_median:.4g}\nbkg median rate : {bkg_rate_median:.4g}\nzero-exp frac   : {zero_exp_frac:.3%}"
            ax.text(0.02, 0.98, text, va="top", family="monospace", fontsize=8)
            ax.set_axis_off()
            ax.set_title("Image summary")
        fig.suptitle(f"Comet-frame image diagnostics ({label})", y=1.02)
        out_png = outdir / f"04_image_{label}.png"
        savefig(fig, str(out_png))
        rendered.append(out_png)
        if first_stats is None:
            first_stats = {
                "status": "ok",
                "png": out_png.name,
                "band": label,
                "src_median_exposure": src_exp_median,
                "bkg_median_exposure": bkg_exp_median,
                "src_median_rate": src_rate_median,
                "bkg_median_rate": bkg_rate_median,
                "zero_exposure_fraction": zero_exp_frac,
                "trail_mask_fraction": masked_frac,
            }
    if not rendered or first_stats is None:
        return {"status": "missing", "reason": "No image diagnostics could be rendered"}
    if len(rendered) > 1:
        first_stats["extra_pngs"] = [p.name for p in rendered[1:]]
    return first_stats


def check_contam(env: dict[str, str], outdir: Path) -> dict[str, Any]:
    workdir = Path(env["WORKDIR"])
    timeline_path = workdir / "contam" / "contamination_timeline.csv"
    report_csv = workdir / "contam" / "contamination_report.csv"
    summary_path = workdir / "contam" / "contamination_summary.json"
    track_path = workdir / "track" / "comet_track.fits"
    if not timeline_path.exists() or not track_path.exists():
        return {
            "status": "missing",
            "reason": "Need contamination_timeline.csv and track/comet_track.fits",
        }
    timeline = _read_contam_timeline(timeline_path)
    track = read_track(str(track_path))
    summary = (
        json.loads(summary_path.read_text(encoding="utf-8"))
        if summary_path.exists()
        else {}
    )
    report_rows = _read_report_rows(report_csv)
    science_policy = str(summary.get("science_policy", "unknown"))
    good_time_map = (
        summary.get("good_time_s", {})
        if isinstance(summary.get("good_time_s"), dict)
        else {}
    )
    selected_good = (
        float(good_time_map.get(science_policy, float("nan")))
        if good_time_map
        else float("nan")
    )
    src_list_path = Path(
        summary.get("srclist", workdir / "detect" / "field_sources_all.fits")
    )
    if not src_list_path.exists():
        src_list_path = workdir / "detect" / "field_sources_curated.fits"
    src = read_srclist(str(src_list_path)) if src_list_path.exists() else None
    src_lookup = (
        {
            int(sid): (float(ra), float(dec))
            for sid, ra, dec in zip(src.srcid, src.ra, src.dec)
        }
        if src is not None
        else {}
    )
    top_rows = sorted(
        report_rows,
        key=lambda r: int(r.get("samples_strict_overlap", "0") or 0),
        reverse=True,
    )[:5]
    times = timeline["time_s"]
    t_ks = (times - times[0]) / 1000.0
    trk_ra, trk_dec = _interp_track_to_times(track, times)
    fig, axes = plt.subplots(
        3, 1, figsize=(10.5, 8.5), sharex=True, constrained_layout=True
    )
    ax = axes[0]
    src_r = env_float(env, "SRC_R_ARCSEC", 60.0)
    mask_r = env_float(env, "FIELD_SOURCE_MASK_R_ARCSEC", 20.0)
    for row in top_rows:
        sid = int(row.get("srcid", "0") or 0)
        if sid not in src_lookup:
            continue
        ra_s, dec_s = src_lookup[sid]
        dx, dy = source_offsets_from_track(trk_ra, trk_dec, ra_s, dec_s)
        dist = np.hypot(
            dx - env_float(env, "SRC_DX_ARCSEC", 0.0),
            dy - env_float(env, "SRC_DY_ARCSEC", 0.0),
        )
        ax.plot(t_ks, dist, label=f"src {sid}")
    ax.axhline(
        src_r + mask_r, linestyle="--", linewidth=1.0, label="source overlap limit"
    )
    ax.set_ylabel("Dist. to source center (arcsec)")
    ax.set_title("Top contamination-driving sources")
    if top_rows:
        ax.legend(loc="best", fontsize=7, ncol=2)
    ax = axes[1]
    ax.step(t_ks, timeline["n_src_overlap"], where="mid", label="source overlaps")
    ax.step(t_ks, timeline["n_bkg_overlap"], where="mid", label="background overlaps")
    ax.set_ylabel("Overlapping sources")
    ax.set_title("Overlap counts vs time")
    ax.legend(loc="best", fontsize=8)
    ax = axes[2]
    policies = ["full", "src", "bkg", "strict"]
    for idx, pol in enumerate(policies):
        y = idx + 0.8 * timeline[f"good_{pol}"]
        ax.step(t_ks, y, where="mid", label=pol)
    ax.set_yticks(np.arange(len(policies)) + 0.4)
    ax.set_yticklabels(policies)
    ax.set_ylim(-0.2, len(policies))
    ax.set_xlabel("Time from start (ks)")
    ax.set_ylabel("GTI policy")
    ax.set_title("Contamination GTI policies")
    text = f"selected policy={science_policy}"
    if np.isfinite(selected_good):
        frac = selected_good / max(float(track.obs_t1 - track.obs_t0), 1.0)
        text += f" | good exposure={selected_good:.1f}s ({frac:.1%})"
    ax.text(
        0.99,
        0.04,
        text,
        ha="right",
        va="bottom",
        transform=ax.transAxes,
        family="monospace",
        fontsize=8,
    )
    png = outdir / "03_contam.png"
    savefig(fig, str(png))
    return {
        "status": "ok",
        "png": png.name,
        "science_policy": science_policy,
        "good_exposure_s": selected_good,
        "top_report_rows": top_rows,
    }


def check_lcurve(env: dict[str, str], outdir: Path) -> dict[str, Any]:
    workdir = Path(env["WORKDIR"]) / "lcurve"

    def pick_lc(name: str) -> tuple[Path | None, str | None]:
        candidates = [
            (workdir / f"{name}_corr_abs.fits", "abs"),
            (workdir / f"{name}_corr_relonly.fits", "relonly"),
            (workdir / f"{name}_corr.fits", "legacy"),
        ]
        for path, mode in candidates:
            if path.exists():
                return (path, mode)
        return (None, None)

    paths = {}
    modes = {}
    lcs: list[LightCurve] = []
    for name in ("PN", "M1", "M2"):
        path, mode = pick_lc(name)
        if path is not None:
            paths[name] = path
            modes[name] = mode
            lcs.append(read_lightcurve(str(path), name))
    if not lcs:
        return {"status": "missing", "reason": "no corrected light curves found"}
    fig, axes = plt.subplots(
        2, 1, figsize=(10.0, 7.0), sharex=True, constrained_layout=True
    )
    base_t0 = None
    for lc in lcs:
        if base_t0 is None:
            base_t0 = lc.time[0]
        t_ks = (lc.time - base_t0) / 1000.0
        axes[0].errorbar(
            t_ks, lc.rate, yerr=lc.error, fmt="none", elinewidth=0.6, capsize=0
        )
        label = lc.name if lc.name not in modes else f"{lc.name} ({modes[lc.name]})"
        axes[0].step(t_ks, lc.rate, where="mid", label=label)
        if lc.fracexp is not None:
            axes[1].step(t_ks, lc.fracexp, where="mid", label=label)
    axes[0].set_ylabel("Corrected rate (counts/s)")
    axes[0].set_title("Per-instrument corrected light curves")
    if lcs:
        axes[0].legend(loc="best", fontsize=8)
    axes[1].set_xlabel("Time from first bin (ks)")
    axes[1].set_ylabel("FRACEXP")
    axes[1].set_ylim(-0.05, 1.05)
    axes[1].set_title("Fractional exposure per time bin")
    axes[1].legend(loc="best", fontsize=8)
    png = outdir / "05_lcurve.png"
    savefig(fig, str(png))
    n_abs = sum((1 for m in modes.values() if m == "abs"))
    n_relonly = sum((1 for m in modes.values() if m == "relonly"))
    out: dict[str, Any] = {
        "status": "ok",
        "png": png.name,
        "n_instruments": len(lcs),
        "instrument_modes": ", ".join((f"{k}:{v}" for k, v in sorted(modes.items()))),
        "n_abs_cameras": n_abs,
        "n_relonly_cameras": n_relonly,
    }
    return out


def check_spectrum(env: dict[str, str], outdir: Path) -> dict[str, Any]:
    specdir = Path(env["WORKDIR"]) / "spectra"
    inst_data: list[tuple[str, Any, Any]] = []
    for inst in ("PN", "M1", "M2"):
        grp_path = specdir / f"{inst}_src_spec_grp.fits"
        src_path = grp_path if grp_path.exists() else specdir / f"{inst}_src_spec.fits"
        bkg_path = specdir / f"{inst}_bkg_spec.fits"
        if src_path.exists() and bkg_path.exists():
            inst_data.append(
                (inst, read_spectrum(str(src_path)), read_spectrum(str(bkg_path)))
            )
    if not inst_data:
        return {
            "status": "missing",
            "reason": "no per-instrument source/background spectra found",
        }
    n_inst = len(inst_data)
    fig, axes = plt.subplots(
        2,
        n_inst,
        figsize=(5.0 * n_inst, 7.5),
        sharex=True,
        squeeze=False,
        constrained_layout=True,
    )
    results: dict[str, Any] = {"status": "ok", "instruments": {}}
    for col, (inst, src, bkg) in enumerate(inst_data):
        if len(src.channel) != len(bkg.channel):
            results["instruments"][inst] = {"error": "channel length mismatch"}
            continue
        scale = 1.0
        if (
            src.exposure
            and bkg.exposure
            and src.backscal
            and bkg.backscal
            and (src.backscal != 0)
            and (bkg.backscal != 0)
        ):
            scale = src.exposure / bkg.exposure * (src.backscal / bkg.backscal)
        scaled_bkg = bkg.counts * scale
        net = src.counts - scaled_bkg
        axes[0, col].step(src.channel, src.counts, where="mid", label="source")
        axes[0, col].step(src.channel, scaled_bkg, where="mid", label="scaled bkg")
        axes[0, col].set_yscale("log")
        axes[0, col].set_ylabel("Counts / bin")
        axes[0, col].set_title(f"{inst} spectrum")
        axes[0, col].legend(loc="best", fontsize=7)
        axes[1, col].step(src.channel, net, where="mid")
        axes[1, col].axhline(0.0, linestyle="--", linewidth=1.0)
        axes[1, col].set_xlabel("Channel")
        axes[1, col].set_ylabel("Approx. net counts / bin")
        axes[1, col].set_title(f"{inst} source - scaled bkg")
        results["instruments"][inst] = {
            "src_total_counts": float(np.nansum(src.counts)),
            "scaled_bkg_total_counts": float(np.nansum(scaled_bkg)),
            "approx_net_total_counts": float(np.nansum(net)),
            "background_scale_used": scale,
        }
    png = outdir / "06_spectrum.png"
    savefig(fig, str(png))
    results["png"] = png.name
    return results


def _first_band(env: dict[str, str]) -> str:
    for entry in env.get("IMAGE_BANDS", "broad:300:2000").split(";"):
        entry = entry.strip()
        if entry:
            return entry.split(":")[0].strip()
    return "broad"


def _parse_seginfo(path: Path) -> dict[str, str]:
    out: dict[str, str] = {}
    if not path.exists():
        return out
    for line in path.read_text(encoding="utf-8").splitlines():
        if "=" in line:
            k, v = line.split("=", 1)
            out[k.strip()] = v.strip()
    return out


def check_pointings(
    env: dict[str, str],
    outdir: Path,
    band_label: str | None = None,
    inst_filter: str | None = None,
) -> dict[str, Any]:
    """Per-time-slice comet-frame mosaic — shows each temporal chunk before stacking."""
    import plot_utils as qpu

    workdir = Path(env["WORKDIR"])
    bands = parse_image_bands(
        env.get("IMAGE_BANDS", "soft:200:1000;broad:300:2000;hard:1000:12000")
    )
    label = band_label or _first_band(env)
    bands = [b for b in bands if b[0] == label]
    if not bands:
        return {"status": "missing", "reason": "No matching IMAGE_BANDS"}
    _, pimin, pimax = bands[0]

    grid_path = workdir / "images" / "grid.json"
    if not grid_path.exists():
        return {"status": "missing", "reason": "grid.json not found"}
    grid = qpu.load_grid(grid_path)
    extent = qpu.grid_extent(grid)
    bin_phys = int(grid["bin_phys"])
    x_min = int(grid["x_min_phys"])
    x_max = int(grid["x_max_phys"])
    y_min = int(grid["y_min_phys"])
    y_max = int(grid["y_max_phys"])
    nx = int(grid["nx"])
    ny = int(grid["ny"])

    # Collect comet-frame events (optionally filtered to one instrument)
    loop_insts = [inst_filter] if inst_filter else ["PN", "M1", "M2"]
    all_x, all_y, all_t, all_pi = [], [], [], []
    insts_used = []
    for inst in loop_insts:
        evt_path = workdir / "comet" / f"{inst}_comet.fits"
        if not evt_path.exists():
            continue
        try:
            with fits.open(str(evt_path), memmap=True) as hdul:
                data = hdul["EVENTS"].data
                if data is None or len(data) == 0:
                    continue
                pi = np.asarray(data["PI"], dtype=float)
                sel = (pi >= pimin) & (pi <= pimax)
                all_x.append(np.asarray(data["X"][sel], dtype=float))
                all_y.append(np.asarray(data["Y"][sel], dtype=float))
                all_t.append(np.asarray(data["TIME"][sel], dtype=float))
                all_pi.append(pi[sel])
                insts_used.append(inst)
        except Exception:
            continue
    if not all_x:
        return {"status": "missing", "reason": "No comet-frame events found"}

    x = np.concatenate(all_x)
    y = np.concatenate(all_y)
    t = np.concatenate(all_t)

    # Decide number of time slices — aim for ~12-20 panels
    obs_dur = t.max() - t.min()
    n_slices = min(20, max(6, int(obs_dur / 5000)))  # ~5 ks per slice
    t_edges = np.linspace(t.min(), t.max(), n_slices + 1)
    slice_dur_ks = (t_edges[1] - t_edges[0]) / 1000.0

    # Build per-slice images
    images = []
    slice_labels = []
    for i in range(n_slices):
        sel = (t >= t_edges[i]) & (t < t_edges[i + 1])
        xi = x[sel]
        yi = y[sel]
        img, _, _ = np.histogram2d(
            yi,
            xi,
            bins=[
                np.arange(y_min, y_max + bin_phys, bin_phys),
                np.arange(x_min, x_max + bin_phys, bin_phys),
            ],
        )
        images.append(img)
        dt_ks = (t_edges[i] - t.min()) / 1000.0
        slice_labels.append(f"t={dt_ks:.1f} ks")

    # Build cumulative stack
    cumulative = np.zeros_like(images[0])
    for img in images:
        cumulative += img

    # Determine layout: slices + 1 cumulative panel
    n_total = n_slices + 1
    ncols = 5
    nrows = int(np.ceil(n_total / ncols))

    # Shared colour scale from cumulative 99th percentile
    finite_cum = cumulative[cumulative > 0]
    vmax_cum = float(np.percentile(finite_cum, 99.5)) if len(finite_cum) > 0 else 1.0
    # Per-slice vmax: scale by fraction of total
    vmax_slice = max(vmax_cum / n_slices * 3, 1.0)

    inst_tag = f" [{inst_filter}]" if inst_filter else ""
    fig, axes = qpu.make_figure(
        nrows,
        ncols,
        figscale=3.2,
        suptitle=f"Per-time-slice comet-frame mosaic ({label}, {slice_dur_ks:.1f} ks bins){inst_tag}",
    )

    for i in range(n_slices):
        r, c = divmod(i, ncols)
        ax = axes[r, c]
        norm = PowerNorm(gamma=0.5, vmin=0, vmax=vmax_slice)
        qpu.imshow_arcsec(
            ax,
            images[i],
            extent,
            norm=norm,
            cmap=qpu.COUNTS_CMAP,
            title=slice_labels[i],
            cbar_label="counts",
            show_xlabel=(r == nrows - 1),
            show_ylabel=(c == 0),
        )
        qpu.add_apertures(ax, env, color="cyan", lw_src=0.8, lw_bkg=0.6)

    # Cumulative panel
    r_cum, c_cum = divmod(n_slices, ncols)
    ax_cum = axes[r_cum, c_cum]
    norm_cum = PowerNorm(gamma=0.5, vmin=0, vmax=vmax_cum)
    qpu.imshow_arcsec(
        ax_cum,
        cumulative,
        extent,
        norm=norm_cum,
        cmap=qpu.COUNTS_CMAP,
        title="CUMULATIVE",
        cbar_label="counts",
        show_xlabel=(r_cum == nrows - 1),
        show_ylabel=(c_cum == 0),
    )
    qpu.add_apertures(ax_cum, env, color="cyan", lw_src=1.2, lw_bkg=0.8)

    # Hide unused panels
    for j in range(n_slices + 1, nrows * ncols):
        r, c = divmod(j, ncols)
        axes[r, c].set_axis_off()

    inst_suffix = f"_{inst_filter}" if inst_filter else ""
    png = outdir / f"09_pointings_{label}{inst_suffix}.png"
    qpu.savefig(fig, png)

    return {
        "status": "ok",
        "png": png.name,
        "n_slices": n_slices,
        "slice_duration_ks": round(slice_dur_ks, 2),
        "total_events": int(len(x)),
        "instruments": insts_used,
    }


def check_pointings_per_inst(
    env: dict[str, str], outdir: Path, band_label: str | None = None
) -> dict[str, Any]:
    """Run check_pointings separately for each detector."""
    combined: dict[str, Any] = {"status": "ok", "extra_pngs": []}
    for inst in ["PN", "M1", "M2"]:
        res = check_pointings(env, outdir, band_label=band_label, inst_filter=inst)
        if res.get("status") == "ok" and res.get("png"):
            combined["extra_pngs"].append(res["png"])
            combined[f"{inst}_events"] = res.get("total_events", 0)
    if not combined["extra_pngs"]:
        return {"status": "missing", "reason": "No per-instrument pointings produced"}
    combined["png"] = combined["extra_pngs"][0]
    return combined


def check_source_filter(env: dict[str, str], outdir: Path) -> dict[str, Any]:
    """Diagnostics for sky-frame source filtering: events removed per source, spatial footprint."""
    import plot_utils as qpu

    workdir = Path(env["WORKDIR"])
    grid_path = workdir / "images" / "grid.json"

    # Read source filter reports
    reports = {}
    total_input = 0
    total_removed = 0
    for inst in ["PN", "M1", "M2"]:
        rpath = workdir / "clean" / f"{inst}_srcfilt.json"
        if rpath.exists():
            rpt = json.loads(rpath.read_text(encoding="utf-8"))
            reports[inst] = rpt
            total_input += rpt.get("n_input_events", 0)
            total_removed += rpt.get("n_removed_events", 0)
    if not reports:
        return {"status": "missing", "reason": "No srcfilt.json reports found"}

    # Read source list
    srclist_path = workdir / "detect" / "field_sources_all.fits"
    if not srclist_path.exists():
        srclist_path = workdir / "detect" / "field_sources_curated.fits"
    src = read_srclist(str(srclist_path)) if srclist_path.exists() else None

    # Aggregate per-source removals
    all_removals: dict[str, int] = {}
    for inst, rpt in reports.items():
        for srcid, count in rpt.get("per_source_removed", {}).items():
            all_removals[srcid] = all_removals.get(srcid, 0) + int(count)
    sorted_srcs = sorted(all_removals.items(), key=lambda kv: kv[1], reverse=True)

    # Layout: 2×2
    fig, axes = qpu.make_figure(
        2, 2, figscale=5.0, suptitle="Source filtering diagnostics"
    )

    # Panel 1: Sky-frame source positions (top-left)
    ax_sky = axes[0, 0]
    if src is not None and len(src.ra) > 0:
        track_path = workdir / "track" / "comet_track.fits"
        if track_path.exists():
            track = read_track(str(track_path))
            ref_ra = env_float(env, "COMET_REF_RA", float(np.mean(track.ra)))
            ref_dec = env_float(env, "COMET_REF_DEC", float(np.mean(track.dec)))
            # Source positions relative to reference
            src_dx, src_dy = sky_offsets_arcsec(src.ra, src.dec, ref_ra, ref_dec)
            # Track positions relative to reference
            trk_dx, trk_dy = sky_offsets_arcsec(track.ra, track.dec, ref_ra, ref_dec)
            ax_sky.plot(trk_dx, trk_dy, "c-", lw=1.5, alpha=0.6, label="Comet track")
            # Size by removal count
            sizes = []
            for s in src.srcid:
                sizes.append(all_removals.get(str(s), 0))
            sizes = np.array(sizes, dtype=float)
            sizes = 10 + 200 * sizes / max(sizes.max(), 1)
            sc = ax_sky.scatter(
                src_dx,
                src_dy,
                s=sizes,
                c=src.detml if src.detml is not None else "red",
                cmap="plasma",
                alpha=0.7,
                edgecolors="k",
                linewidths=0.3,
            )
            if src.detml is not None:
                qpu.matched_colorbar(ax_sky, sc, label="DET_ML")
            ax_sky.legend(fontsize=7, loc="upper right")
        else:
            ax_sky.scatter(src.ra, src.dec, s=20, c="red", alpha=0.7)
        ax_sky.set_title("Source positions (sky frame)", fontsize=9, fontweight="bold")
        ax_sky.set_xlabel("East offset (arcsec)", fontsize=8)
        ax_sky.set_ylabel("North offset (arcsec)", fontsize=8)
        ax_sky.tick_params(labelsize=7)
        ax_sky.set_aspect("equal")
    else:
        ax_sky.text(0.5, 0.5, "No sources", ha="center", va="center")
        ax_sky.set_axis_off()

    # Panel 2: Bar chart of removals per source (top-right)
    ax_bar = axes[0, 1]
    if sorted_srcs:
        top_n = min(15, len(sorted_srcs))
        top_ids = [s[0] for s in sorted_srcs[:top_n]]
        top_counts = [s[1] for s in sorted_srcs[:top_n]]
        bars = ax_bar.barh(
            range(top_n), top_counts, color=qpu.INST_COLORS["EPIC"], alpha=0.8
        )
        ax_bar.set_yticks(range(top_n))
        ax_bar.set_yticklabels([f"src {s}" for s in top_ids], fontsize=7)
        ax_bar.invert_yaxis()
        ax_bar.set_xlabel("Events removed (EPIC total)", fontsize=8)
        ax_bar.set_title("Events removed per source", fontsize=9, fontweight="bold")
        ax_bar.tick_params(labelsize=7)
        ax_bar.tick_params(axis="y", pad=1)
    else:
        ax_bar.text(0.5, 0.5, "No removals", ha="center", va="center")
        ax_bar.set_axis_off()

    # Panel 3: Per-instrument bar chart (bottom-left)
    ax_inst = axes[1, 0]
    inst_names = []
    inst_input = []
    inst_removed = []
    for inst in ["PN", "M1", "M2"]:
        if inst in reports:
            inst_names.append(inst)
            inst_input.append(reports[inst].get("n_input_events", 0))
            inst_removed.append(reports[inst].get("n_removed_events", 0))
    if inst_names:
        x_pos = np.arange(len(inst_names))
        pct = [100.0 * r / max(i, 1) for r, i in zip(inst_removed, inst_input)]
        colors = [qpu.INST_COLORS.get(n, "gray") for n in inst_names]
        ax_inst.bar(x_pos, pct, color=colors, alpha=0.8, edgecolor="k", linewidth=0.5)
        ax_inst.set_xticks(x_pos)
        ax_inst.set_xticklabels(inst_names, fontsize=9)
        ax_inst.set_ylabel("Events removed (%)", fontsize=8)
        ax_inst.set_title(
            "Removal fraction per instrument", fontsize=9, fontweight="bold"
        )
        ax_inst.tick_params(labelsize=7)
        for i, (p, r) in enumerate(zip(pct, inst_removed)):
            ax_inst.text(
                i, p + 0.1, f"{p:.2f}%\n({r})", ha="center", va="bottom", fontsize=7
            )
    else:
        ax_inst.text(0.5, 0.5, "No data", ha="center", va="center")
        ax_inst.set_axis_off()

    # Panel 4: Summary text (bottom-right)
    ax_txt = axes[1, 1]
    lines = [
        f"Total input events : {total_input:,}",
        f"Total removed      : {total_removed:,}",
        f"Removal fraction   : {100*total_removed/max(total_input,1):.3f}%",
        f"N sources          : {len(all_removals)}",
        f"Mask radius        : {reports[list(reports)[0]].get('mask_radius_arcsec', '?')}\"",
        "",
    ]
    for inst in ["PN", "M1", "M2"]:
        if inst in reports:
            r = reports[inst]
            lines.append(
                f"{inst}: {r.get('n_removed_events',0):,} / {r.get('n_input_events',0):,} "
                f"({100*r.get('n_removed_events',0)/max(r.get('n_input_events',1),1):.2f}%)"
            )
    qpu.summary_text_panel(ax_txt, lines, fontsize=8.5)
    ax_txt.set_title("Summary", fontsize=9, fontweight="bold")

    png = outdir / "03b_source_filter.png"
    qpu.savefig(fig, png)

    return {
        "status": "ok",
        "png": png.name,
        "total_input_events": total_input,
        "total_removed_events": total_removed,
        "removal_fraction": round(total_removed / max(total_input, 1), 6),
        "n_sources_masked": len(all_removals),
    }


def check_image_per_inst(
    env: dict[str, str], outdir: Path, band_label: str | None = None
) -> dict[str, Any]:
    """Per-detector 2×2 image breakdown: counts / exposure / clean rate / net rate."""
    import plot_utils as qpu

    workdir = Path(env["WORKDIR"])
    bands = parse_image_bands(
        env.get("IMAGE_BANDS", "soft:200:1000;broad:300:2000;hard:1000:12000")
    )
    if band_label:
        bands = [b for b in bands if b[0] == band_label]
    else:
        bands = bands[:1]  # default: first band only
    if not bands:
        return {"status": "missing", "reason": "No matching IMAGE_BANDS entries"}

    grid_path = workdir / "images" / "grid.json"
    if not grid_path.exists():
        return {"status": "missing", "reason": "grid.json not found"}
    grid = qpu.load_grid(grid_path)
    extent = qpu.grid_extent(grid)

    rendered: list[Path] = []
    first_stats: dict[str, Any] | None = None

    for label, _pimin, _pimax in bands:
        banddir = workdir / "images" / label
        for inst in ["PN", "M1", "M2"]:
            counts_path = banddir / f"{inst}_counts.fits"
            exp_path = banddir / f"{inst}_exp.fits"
            if not counts_path.exists():
                continue

            panels: list[tuple[str, np.ndarray | None, str]] = []

            # Counts
            if counts_path.exists():
                data, _ = read_image(str(counts_path))
                panels.append(("Counts", data, "counts"))
            # Exposure
            expo = None
            if exp_path.exists():
                expo, _ = read_image(str(exp_path))
                panels.append(("Exposure (s)", expo, "expo"))
            # Clean rate
            clean_path = banddir / f"{inst}_rate_clean.fits"
            rate_path = banddir / f"{inst}_rate.fits"
            if clean_path.exists():
                clean, _ = read_image(str(clean_path))
                panels.append(("Clean rate", clean, "rate"))
            elif rate_path.exists():
                rate, _ = read_image(str(rate_path))
                panels.append(("Rate", rate, "rate"))
            # Net rate
            net_path = banddir / f"{inst}_net_rate.fits"
            if net_path.exists():
                net, _ = read_image(str(net_path))
                panels.append(("Net rate", net, "net"))
            else:
                # Compute rate from counts/exp as fallback
                if expo is not None and counts_path.exists():
                    with np.errstate(divide="ignore", invalid="ignore"):
                        fallback_rate = np.where(expo > 0, data / expo, np.nan)
                    panels.append(("Rate (counts/exp)", fallback_rate, "rate"))

            if not panels:
                continue

            nrows_p = 2
            ncols_p = 2
            fig, axes = qpu.make_figure(
                nrows_p,
                ncols_p,
                figscale=4.5,
                suptitle=f"{inst} per-detector diagnostics ({label})",
            )

            for idx, (ptitle, pdata, ptype) in enumerate(panels[:4]):
                r, c = divmod(idx, ncols_p)
                ax = axes[r, c]
                if pdata is None:
                    ax.text(0.5, 0.5, f"{ptitle}: N/A", ha="center", va="center")
                    ax.set_axis_off()
                    continue

                if ptype == "counts":
                    norm = qpu.rate_norm(pdata)
                    qpu.imshow_arcsec(
                        ax,
                        pdata,
                        extent,
                        norm=norm,
                        cmap=qpu.COUNTS_CMAP,
                        title=ptitle,
                        cbar_label="counts",
                        show_xlabel=(r == nrows_p - 1),
                        show_ylabel=(c == 0),
                    )
                elif ptype == "expo":
                    norm = qpu.rate_norm(pdata, pct=99.5)
                    qpu.imshow_arcsec(
                        ax,
                        pdata,
                        extent,
                        norm=norm,
                        cmap=qpu.EXPO_CMAP,
                        title=ptitle,
                        cbar_label="s",
                        show_xlabel=(r == nrows_p - 1),
                        show_ylabel=(c == 0),
                    )
                elif ptype == "net":
                    # Symmetric scale for net rate
                    if expo is not None:
                        masked = qpu.mask_low_exposure(pdata, expo)
                    else:
                        masked = pdata
                    norm = qpu.symmetric_norm(masked)
                    qpu.imshow_arcsec(
                        ax,
                        masked,
                        extent,
                        norm=norm,
                        cmap=qpu.NET_CMAP,
                        title=ptitle,
                        cbar_label="ct/s/px",
                        show_xlabel=(r == nrows_p - 1),
                        show_ylabel=(c == 0),
                    )
                else:  # rate
                    if expo is not None:
                        masked = qpu.mask_low_exposure(pdata, expo)
                    else:
                        masked = pdata
                    norm = qpu.rate_norm(masked, exposure=expo)
                    qpu.imshow_arcsec(
                        ax,
                        masked,
                        extent,
                        norm=norm,
                        cmap=qpu.RATE_CMAP,
                        title=ptitle,
                        cbar_label="ct/s/px",
                        show_xlabel=(r == nrows_p - 1),
                        show_ylabel=(c == 0),
                    )

                qpu.add_apertures(ax, env, color="cyan")

            # Hide unused panels
            for idx in range(len(panels), nrows_p * ncols_p):
                r, c = divmod(idx, ncols_p)
                axes[r, c].set_axis_off()

            png = outdir / f"04d_image_{inst}_{label}.png"
            qpu.savefig(fig, png)
            rendered.append(png)

            if first_stats is None:
                # Source aperture stats for this instrument
                src_mask = circle_mask_from_extent(
                    panels[0][1].shape,
                    extent,
                    env_float(env, "SRC_DX_ARCSEC", 0.0),
                    env_float(env, "SRC_DY_ARCSEC", 0.0),
                    env_float(env, "SRC_R_ARCSEC", 60.0),
                )
                src_med = (
                    float(np.nanmedian(expo[src_mask]))
                    if expo is not None and np.any(src_mask)
                    else math.nan
                )
                first_stats = {
                    "status": "ok",
                    "png": png.name,
                    "band": label,
                    "instruments_rendered": [],
                    "first_inst_src_median_exposure": src_med,
                }
            first_stats["instruments_rendered"].append(inst)

    if not rendered or first_stats is None:
        return {"status": "missing", "reason": "No per-instrument images found"}
    if len(rendered) > 1:
        first_stats["extra_pngs"] = [p.name for p in rendered[1:]]
    return first_stats


def check_radial_profile(
    env: dict[str, str], outdir: Path, band_label: str | None = None
) -> dict[str, Any]:
    """Azimuthally-averaged radial brightness profile in the comet frame."""
    import plot_utils as qpu

    workdir = Path(env["WORKDIR"])
    label = band_label or _first_band(env)
    banddir = workdir / "images" / label

    grid_path = workdir / "images" / "grid.json"
    if not grid_path.exists():
        return {"status": "missing", "reason": "grid.json not found"}
    grid = qpu.load_grid(grid_path)
    extent = qpu.grid_extent(grid)

    counts_path = banddir / "EPIC_counts.fits"
    exp_path = banddir / "EPIC_exp.fits"
    bkg_path = banddir / "EPIC_bkg_rate.fits"
    clean_path = banddir / "EPIC_rate_clean.fits"
    if not (counts_path.exists() and exp_path.exists()):
        return {"status": "missing", "reason": "EPIC counts/exp not found"}

    counts, _ = read_image(str(counts_path))
    expo, _ = read_image(str(exp_path))
    bkg_rate = None
    if bkg_path.exists():
        bkg_rate, _ = read_image(str(bkg_path))
    clean_rate = None
    if clean_path.exists():
        clean_rate, _ = read_image(str(clean_path))

    # Build radial distance array (arcsec from center)
    ny, nx_img = counts.shape
    e_min, e_max, n_min, n_max = extent
    xs = (
        np.linspace(e_min, e_max, nx_img, endpoint=False)
        + 0.5 * (e_max - e_min) / nx_img
    )
    ys = np.linspace(n_min, n_max, ny, endpoint=False) + 0.5 * (n_max - n_min) / ny
    xx, yy = np.meshgrid(xs, ys)
    rr = np.sqrt(xx**2 + yy**2)

    # Radial bins
    max_r = min(abs(e_min), abs(e_max), abs(n_min), abs(n_max)) * 0.95
    r_edges = np.arange(0, max_r, 10.0)
    if len(r_edges) < 3:
        return {"status": "missing", "reason": "Image too small for radial profile"}
    r_centers = 0.5 * (r_edges[:-1] + r_edges[1:])

    # Compute azimuthal profiles
    with np.errstate(divide="ignore", invalid="ignore"):
        rate = np.where(expo > 0, counts / expo, np.nan)

    profile_rate = np.full(len(r_centers), np.nan)
    profile_err = np.full(len(r_centers), np.nan)
    profile_bkg = np.full(len(r_centers), np.nan)
    profile_counts = np.full(len(r_centers), np.nan)
    profile_expo = np.full(len(r_centers), np.nan)

    for i in range(len(r_centers)):
        ann = (rr >= r_edges[i]) & (rr < r_edges[i + 1])
        valid = ann & np.isfinite(rate) & (expo > 0)
        if np.sum(valid) < 2:
            continue
        total_cts = np.nansum(counts[valid])
        total_exp = np.nansum(expo[valid])
        profile_counts[i] = total_cts
        profile_expo[i] = total_exp
        profile_rate[i] = total_cts / total_exp if total_exp > 0 else np.nan
        profile_err[i] = np.sqrt(total_cts) / total_exp if total_exp > 0 else np.nan
        if bkg_rate is not None:
            bkg_vals = bkg_rate[valid]
            profile_bkg[i] = (
                float(np.nanmean(bkg_vals[np.isfinite(bkg_vals)]))
                if np.any(np.isfinite(bkg_vals))
                else np.nan
            )

    # Azimuthal sectors (8 sectors)
    theta = np.degrees(np.arctan2(yy, xx))  # -180 to 180
    sector_edges = np.linspace(-180, 180, 9)
    sector_labels = ["E", "NE", "N", "NW", "W", "SW", "S", "SE"]
    sector_profiles = {}
    for s in range(8):
        sec_mask = (theta >= sector_edges[s]) & (theta < sector_edges[s + 1])
        sprof = np.full(len(r_centers), np.nan)
        for i in range(len(r_centers)):
            ann = (rr >= r_edges[i]) & (rr < r_edges[i + 1]) & sec_mask
            valid = ann & np.isfinite(rate) & (expo > 0)
            if np.sum(valid) < 1:
                continue
            sprof[i] = np.nansum(counts[valid]) / np.nansum(expo[valid])
        sector_profiles[sector_labels[s]] = sprof

    # Plot: 3-panel layout
    fig, axes = qpu.make_figure(
        1, 3, figscale=4.0, suptitle=f"Radial profile ({label})"
    )

    # Panel 1: Rate vs radius with background
    ax1 = axes[0, 0]
    good = np.isfinite(profile_rate)
    ax1.errorbar(
        r_centers[good],
        profile_rate[good],
        yerr=profile_err[good],
        fmt="o-",
        ms=2,
        lw=1.0,
        color="k",
        label="Rate",
        capsize=1,
    )
    if bkg_rate is not None:
        good_b = np.isfinite(profile_bkg)
        ax1.plot(
            r_centers[good_b], profile_bkg[good_b], "r--", lw=1.2, label="Background"
        )
    ax1.set_xlabel("Radius (arcsec)", fontsize=9)
    ax1.set_ylabel("Rate (ct/s/pixel)", fontsize=9)
    ax1.set_title("Azimuthal average", fontsize=9, fontweight="bold")
    ax1.legend(fontsize=7)
    ax1.tick_params(labelsize=7)
    ax1.set_xlim(0, max_r)

    # Panel 2: Cumulative enclosed counts
    ax2 = axes[0, 1]
    cum_counts = np.nancumsum(np.where(np.isfinite(profile_counts), profile_counts, 0))
    ax2.plot(r_centers, cum_counts, "b-", lw=1.5)
    ax2.set_xlabel("Radius (arcsec)", fontsize=9)
    ax2.set_ylabel("Enclosed counts", fontsize=9)
    ax2.set_title("Growth curve", fontsize=9, fontweight="bold")
    ax2.tick_params(labelsize=7)
    ax2.set_xlim(0, max_r)
    # Mark 50% and 90%
    if cum_counts[-1] > 0:
        for frac, ls in [(0.5, ":"), (0.9, "--")]:
            target = frac * cum_counts[-1]
            idx = np.searchsorted(cum_counts, target)
            if idx < len(r_centers):
                ax2.axhline(target, color="gray", ls=ls, lw=0.8)
                ax2.axvline(r_centers[idx], color="gray", ls=ls, lw=0.8)
                ax2.text(
                    r_centers[idx] + 5,
                    target,
                    f"{int(frac*100)}%",
                    fontsize=7,
                    va="bottom",
                )

    # Panel 3: Sector profiles
    ax3 = axes[0, 2]
    sector_colors = plt.cm.hsv(np.linspace(0, 1, 8, endpoint=False))
    for (slabel, sprof), sc in zip(sector_profiles.items(), sector_colors):
        good_s = np.isfinite(sprof) & (r_centers < max_r * 0.7)
        if np.any(good_s):
            ax3.plot(
                r_centers[good_s],
                sprof[good_s],
                "-",
                color=sc,
                lw=1.0,
                alpha=0.8,
                label=slabel,
            )
    ax3.set_xlabel("Radius (arcsec)", fontsize=9)
    ax3.set_ylabel("Rate (ct/s/pixel)", fontsize=9)
    ax3.set_title("Sector profiles (8 azimuths)", fontsize=9, fontweight="bold")
    ax3.legend(fontsize=6, ncol=2, loc="upper right")
    ax3.tick_params(labelsize=7)
    ax3.set_xlim(0, max_r * 0.7)

    png = outdir / f"04b_radial_{label}.png"
    qpu.savefig(fig, png)

    # Compute summary stats
    src_r = env_float(env, "SRC_R_ARCSEC", 60.0)
    src_mask = r_centers <= src_r
    src_total = float(np.nansum(profile_counts[src_mask])) if np.any(src_mask) else 0

    return {
        "status": "ok",
        "png": png.name,
        "band": label,
        "max_radius_arcsec": round(float(max_r), 1),
        "n_radial_bins": len(r_centers),
        "source_counts_within_r": round(src_total, 1),
        "total_enclosed_counts": (
            round(float(cum_counts[-1]), 1) if len(cum_counts) else 0
        ),
    }


def check_frames(
    env: dict[str, str], outdir: Path, band_label: str | None = None
) -> dict[str, Any]:
    """Side-by-side sky-frame vs comet-frame comparison."""
    import plot_utils as qpu

    workdir = Path(env["WORKDIR"])
    label = band_label or _first_band(env)
    rendered: list[Path] = []

    grid_path = workdir / "images" / "grid.json"
    if not grid_path.exists():
        return {"status": "missing", "reason": "grid.json not found"}
    grid = qpu.load_grid(grid_path)
    comet_extent = qpu.grid_extent(grid)

    # Still-sky image
    ss_dir = workdir / "qc" / "mosaics" / "stillsky" / label
    ss_rate = ss_dir / "EPIC_rate.fits"

    # Comet-frame images
    cm_dir = workdir / "images" / label
    cm_rate = cm_dir / "EPIC_rate.fits"
    cm_clean = cm_dir / "EPIC_rate_clean.fits"
    cm_exp = cm_dir / "EPIC_exp.fits"

    if not cm_rate.exists():
        return {"status": "missing", "reason": f"No comet-frame rate image for {label}"}

    # Read all available images
    panels = []

    if ss_rate.exists():
        data_ss, hdr_ss = read_image(str(ss_rate))
        ss_extent_info = _extent_from_detect_header(hdr_ss, data_ss.shape)
        ss_extent = ss_extent_info[0]
        panels.append(("Still-sky rate", data_ss, ss_extent, "sky"))

    data_cm, _ = read_image(str(cm_rate))
    panels.append(("Comet-frame rate", data_cm, comet_extent, "comet"))

    if cm_clean.exists():
        data_cl, _ = read_image(str(cm_clean))
        panels.append(("Comet clean rate", data_cl, comet_extent, "comet"))

    ncols = len(panels)
    fig, axes = qpu.make_figure(
        1, ncols, figscale=4.5, suptitle=f"Frame comparison ({label})"
    )

    for idx, (title, data, ext, frame_type) in enumerate(panels):
        ax = axes[0, idx]
        if ext is None:
            ax.text(0.5, 0.5, "No WCS", ha="center", va="center")
            ax.set_axis_off()
            continue

        expo_data = None
        if frame_type == "comet" and cm_exp.exists():
            expo_data, _ = read_image(str(cm_exp))
        norm = qpu.rate_norm(data, exposure=expo_data)
        qpu.imshow_arcsec(
            ax,
            data,
            ext,
            norm=norm,
            cmap=qpu.RATE_CMAP,
            title=title,
            cbar_label="ct/s/px",
            show_xlabel=True,
            show_ylabel=(idx == 0),
        )
        if frame_type == "comet":
            qpu.add_apertures(ax, env, color="cyan")

        # Overlay track on sky frame
        if frame_type == "sky":
            track_path = workdir / "track" / "comet_track.fits"
            if track_path.exists():
                try:
                    track = read_track(str(track_path))
                    ref_ra = ss_extent_info[1]
                    ref_dec = ss_extent_info[2]
                    if ref_ra is not None and ref_dec is not None:
                        trk_dx, trk_dy = sky_offsets_arcsec(
                            track.ra, track.dec, ref_ra, ref_dec
                        )
                        ax.plot(trk_dx, trk_dy, "c-", lw=1.5, alpha=0.8)
                except Exception:
                    pass

    qpu.add_scalebar(axes[0, ncols - 1], 120.0, "2'", color="white")

    png = outdir / f"04c_frames_{label}.png"
    qpu.savefig(fig, png)
    rendered.append(png)

    return {
        "status": "ok",
        "png": png.name,
        "band": label,
        "n_panels": ncols,
    }


def check_spotcheck(
    env: dict[str, str], outdir: Path, band: str | None = None
) -> dict[str, Any]:
    workdir = Path(env["WORKDIR"])
    band = band or env.get("QC_SPOTCHECK_BAND", "") or _first_band(env)
    root = workdir / "qc" / "spotchecks" / band
    if not root.exists():
        return {
            "status": "missing",
            "reason": f"Spot-check directory not found: {root}",
        }
    segdirs = sorted(
        [p for p in root.iterdir() if p.is_dir() and p.name.startswith("seg")]
    )
    if not segdirs:
        return {
            "status": "missing",
            "reason": f"No seg* directories found under {root}",
        }
    panels = []
    for segdir in segdirs:
        img = segdir / "EPIC_rate.fits"
        masked = segdir / "EPIC_rate_masked.fits"
        if not img.exists():
            continue
        info = _parse_seginfo(segdir / "segment_info.txt")
        panels.append((segdir.name, img, masked if masked.exists() else None, info))
    if not panels:
        return {
            "status": "missing",
            "reason": f"No EPIC_rate.fits files found under {root}",
        }
    n = len(panels)
    nrows = (
        2 if any((masked is not None for _seg, _img, masked, _info in panels)) else 1
    )
    fig, axes = plt.subplots(
        nrows, n, figsize=(4.8 * n, 4.4 * nrows), squeeze=False, constrained_layout=True
    )
    summary = []
    for col, (segname, img_path, masked_path, info) in enumerate(panels):
        data, _ = read_image(str(img_path))
        ax = axes[0, col]
        im = ax.imshow(data, origin="lower", norm=nice_norm(data))
        plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
        segid = info.get("segid", segname.replace("seg", ""))
        utc0 = info.get("utc_start", "")
        utc1 = info.get("utc_stop", "")
        ax.set_title(f"seg {segid} unmasked\n{utc0} → {utc1}")
        ax.set_xlabel("X")
        ax.set_ylabel("Y")
        finite = np.isfinite(data)
        summary_row = {
            "segment": segname,
            "image": str(img_path),
            "finite_fraction": float(np.mean(finite)) if finite.size else None,
            "median_rate": (
                float(np.nanmedian(data[finite])) if np.any(finite) else None
            ),
        }
        if nrows == 2:
            ax2 = axes[1, col]
            data2, _ = (
                read_image(str(masked_path))
                if masked_path is not None
                else (data, None)
            )
            im2 = ax2.imshow(data2, origin="lower", norm=nice_norm(data2))
            plt.colorbar(im2, ax=ax2, fraction=0.046, pad=0.04)
            ax2.set_title(f"seg {segid} masked")
            ax2.set_xlabel("X")
            ax2.set_ylabel("Y")
            finite2 = np.isfinite(data2)
            summary_row["masked_image"] = (
                str(masked_path) if masked_path is not None else None
            )
            summary_row["masked_finite_fraction"] = (
                float(np.mean(finite2)) if finite2.size else None
            )
            summary_row["masked_median_rate"] = (
                float(np.nanmedian(data2[finite2])) if np.any(finite2) else None
            )
        summary.append(summary_row)
    fig.suptitle(f"Comet-frame spot-check rate images ({band})", y=1.02)
    png = outdir / f"07_spotchecks_{band}.png"
    savefig(fig, str(png))
    json_out = png.with_suffix(".json")
    json_out.write_text(
        json.dumps({"band": band, "panels": summary}, indent=2) + "\n", encoding="utf-8"
    )
    print(png)
    return {"status": "ok", "png": png.name, "band": band, "n_panels": n}


def _expmap_stats(data: np.ndarray, mask: np.ndarray) -> dict[str, float | int]:
    vals = data[mask]
    finite = np.isfinite(vals)
    vals = vals[finite]
    if vals.size == 0:
        return {
            "pixels": int(mask.sum()),
            "finite_pixels": 0,
            "median": float("nan"),
            "mean": float("nan"),
            "min": float("nan"),
            "max": float("nan"),
            "zero_fraction": float("nan"),
        }
    return {
        "pixels": int(mask.sum()),
        "finite_pixels": int(vals.size),
        "median": float(np.nanmedian(vals)),
        "mean": float(np.nanmean(vals)),
        "min": float(np.nanmin(vals)),
        "max": float(np.nanmax(vals)),
        "zero_fraction": float(np.mean(vals <= 0.0)),
    }


def check_region_support(
    env: dict[str, str], outdir: Path, band: str | None = None, inst: str = "EPIC"
) -> dict[str, Any]:
    workdir = Path(env["WORKDIR"])
    bands = parse_image_bands(env.get("IMAGE_BANDS", "broad:300:2000"))
    if not bands:
        return {"status": "error", "reason": "No IMAGE_BANDS configured"}
    label = band or bands[0][0]
    grid_path = workdir / "images" / "grid.json"
    if not grid_path.exists():
        return {"status": "missing", "reason": f"Grid file not found: {grid_path}"}
    grid = json.loads(grid_path.read_text(encoding="utf-8"))
    center_x = float(grid["center_x_phys"])
    center_y = float(grid["center_y_phys"])
    x_min = float(grid["x_min_phys"])
    x_max = float(grid["x_max_phys"])
    y_min = float(grid["y_min_phys"])
    y_max = float(grid["y_max_phys"])
    scale = float(grid["scale_arcsec_per_phys"])
    extent = (
        (x_min - center_x) * scale,
        (x_max - center_x) * scale,
        (y_min - center_y) * scale,
        (y_max - center_y) * scale,
    )
    if inst == "EPIC":
        path = workdir / "images" / label / "EPIC_exp.fits"
    else:
        path = workdir / "images" / label / f"{inst}_exp.fits"
    if not path.exists():
        return {"status": "missing", "reason": f"Exposure map not found: {path}"}
    data, _ = read_image(str(path))
    src_mask = circle_mask_from_extent(
        data.shape,
        extent,
        env_float(env, "SRC_DX_ARCSEC", 0.0),
        env_float(env, "SRC_DY_ARCSEC", 0.0),
        env_float(env, "SRC_R_ARCSEC", 60.0),
    )
    if env.get("BKG_MODE", "annulus").strip().lower() == "circle":
        bkg_mask = circle_mask_from_extent(
            data.shape,
            extent,
            env_float(env, "BKG_DX_ARCSEC", 0.0),
            env_float(env, "BKG_DY_ARCSEC", 0.0),
            env_float(env, "BKG_R_ARCSEC", 150.0),
        )
    else:
        bkg_mask = annulus_mask_from_extent(
            data.shape,
            extent,
            env_float(env, "BKG_DX_ARCSEC", 0.0),
            env_float(env, "BKG_DY_ARCSEC", 0.0),
            env_float(env, "BKG_RIN_ARCSEC", 90.0),
            env_float(env, "BKG_ROUT_ARCSEC", 150.0),
        )
    out = {
        "band": label,
        "instrument": inst,
        "expmap": str(path),
        "source": _expmap_stats(data, src_mask),
        "background": _expmap_stats(data, bkg_mask),
    }
    json_out = outdir / f"region_support_{inst}.json"
    json_out.write_text(
        json.dumps(out, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(out, indent=2, sort_keys=True))
    return {
        "status": "ok",
        "instrument": inst,
        "band": label,
        "src_median_exp": out["source"].get("median"),
        "bkg_median_exp": out["background"].get("median"),
        "src_zero_fraction": out["source"].get("zero_fraction"),
        "bkg_zero_fraction": out["background"].get("zero_fraction"),
    }


def _mosaic_extent(header: fits.Header, shape: tuple[int, int]):
    return _extent_from_detect_header(header, shape)


def _find_mosaic_bands(root: Path) -> list[str]:
    if not root.is_dir():
        return []
    return sorted((d.name for d in root.iterdir() if d.is_dir()))


def check_mosaics(
    env: dict[str, str], outdir: Path, band_label: str | None = None
) -> dict[str, Any]:
    workdir = Path(env["WORKDIR"])
    mosaics_root = workdir / "qc" / "mosaics"
    stillsky_root = mosaics_root / "stillsky"
    comet_root = mosaics_root / "comet"
    if not mosaics_root.exists():
        return {"status": "missing", "reason": "No qc/mosaics directory"}
    instruments = ["PN", "M1", "M2"]
    rendered: list[Path] = []
    first_stats: dict[str, Any] | None = None
    ss_bands = _find_mosaic_bands(stillsky_root)
    if band_label:
        ss_bands = [b for b in ss_bands if b == band_label]

    def _safe_imshow(ax, data, extent, **extra_kw):
        """imshow with aspect='equal' and matched colorbar."""
        kw = {"origin": "lower", "aspect": "equal"}
        if extent is not None:
            kw["extent"] = extent
        kw.update(extra_kw)
        im = ax.imshow(data, **kw)
        ax.figure.colorbar(im, ax=ax, pad=0.01, fraction=0.046)
        return im

    for band in ss_bands:
        banddir = stillsky_root / band
        epic_counts = banddir / "EPIC_counts.fits"
        epic_masked = banddir / "EPIC_counts_masked.fits"
        epic_bkg = banddir / "EPIC_bkg.fits"
        epic_exp = banddir / "EPIC_exp.fits"
        epic_exp_masked = banddir / "EPIC_exp_masked.fits"
        epic_rate = banddir / "EPIC_rate.fits"
        epic_rate_masked = banddir / "EPIC_rate_masked.fits"
        epic_mask = banddir / "EPIC_source_mask.fits"
        if not epic_counts.exists():
            continue
        counts_data, counts_hdr = read_image(str(epic_counts))
        extent_info = _mosaic_extent(counts_hdr, counts_data.shape)
        extent = extent_info[0]
        fig, axes = plt.subplots(2, 4, figsize=(20.0, 11.0), constrained_layout=True)
        top_items = [
            (epic_counts, "EPIC counts", "sqrt"),
            (epic_masked, "EPIC masked counts", "sqrt"),
            (epic_bkg, "EPIC background", "sqrt"),
            (epic_exp, "EPIC exposure", "sqrt"),
        ]
        for ax, (fpath, title, stretch) in zip(axes[0], top_items):
            if fpath.exists():
                data, hdr = read_image(str(fpath))
                extra = {}
                if stretch == "sqrt":
                    extra["norm"] = nice_norm(data, stretch="sqrt")
                else:
                    finite = np.isfinite(data)
                    if np.any(finite):
                        vmin, vmax = np.nanpercentile(data[finite], [1.0, 99.0])
                        extra["vmin"], extra["vmax"] = (vmin, vmax)
                _safe_imshow(ax, data, extent, **extra)
            else:
                ax.text(0.5, 0.5, "not available", ha="center", va="center")
                ax.set_axis_off()
            ax.set_title(title, fontsize=9)
            if extent is not None:
                ax.set_xlabel("East offset (arcsec)", fontsize=8)
                ax.set_ylabel("North offset (arcsec)", fontsize=8)
        for ax, inst in zip(axes[1][:3], instruments):
            rpath = banddir / f"{inst}_rate.fits"
            if rpath.exists():
                data, _ = read_image(str(rpath))
                extra = {}
                finite = np.isfinite(data)
                if np.any(finite):
                    vmin, vmax = np.nanpercentile(data[finite], [1.0, 99.0])
                    extra["vmin"], extra["vmax"] = (vmin, vmax)
                _safe_imshow(ax, data, extent, **extra)
            else:
                ax.text(0.5, 0.5, f"{inst} missing", ha="center", va="center")
                ax.set_axis_off()
            ax.set_title(f"{inst} rate (still-sky)", fontsize=9)
            if extent is not None:
                ax.set_xlabel("East offset (arcsec)", fontsize=8)
                ax.set_ylabel("North offset (arcsec)", fontsize=8)
        ax = axes[1][3]
        if epic_mask.exists():
            data, _ = read_image(str(epic_mask))
            _safe_imshow(ax, data, extent, vmin=0.0, vmax=1.0)
            ax.set_title("Source mask", fontsize=9)
            if extent is not None:
                ax.set_xlabel("East offset (arcsec)", fontsize=8)
                ax.set_ylabel("North offset (arcsec)", fontsize=8)
        else:
            ax.text(0.5, 0.5, "no mask", ha="center", va="center")
            ax.set_axis_off()
            ax.set_title("Source mask", fontsize=9)
        fig.suptitle(
            f"Still-sky mosaic diagnostics ({band})", fontsize=13, fontweight="bold"
        )
        png = outdir / f"08_stillsky_mosaic_{band}.png"
        savefig(fig, str(png))
        rendered.append(png)
        if first_stats is None:
            first_stats = {
                "status": "ok",
                "png": png.name,
                "band": band,
                "frame": "stillsky",
            }
    cm_bands = _find_mosaic_bands(comet_root)
    if band_label:
        cm_bands = [b for b in cm_bands if b == band_label]
    for band in cm_bands:
        banddir = comet_root / band
        grid_path = workdir / "images" / "grid.json"
        comet_extent = None
        if grid_path.exists():
            grid = json.loads(grid_path.read_text(encoding="utf-8"))
            cx = float(grid["center_x_phys"])
            cy = float(grid["center_y_phys"])
            xmin = float(grid["x_min_phys"])
            xmax = float(grid["x_max_phys"])
            ymin = float(grid["y_min_phys"])
            ymax = float(grid["y_max_phys"])
            scale = float(grid["scale_arcsec_per_phys"])
            comet_extent = (
                (xmin - cx) * scale,
                (xmax - cx) * scale,
                (ymin - cy) * scale,
                (ymax - cy) * scale,
            )
        fig, axes = plt.subplots(2, 4, figsize=(20.0, 11.0), constrained_layout=True)
        for ax, inst in zip(axes[0][:3], instruments):
            rpath = banddir / f"{inst}_rate.fits"
            if rpath.exists():
                data, _ = read_image(str(rpath))
                extra = {}
                finite = np.isfinite(data)
                if np.any(finite):
                    vmin, vmax = np.nanpercentile(data[finite], [1.0, 99.0])
                    extra["vmin"], extra["vmax"] = (vmin, vmax)
                _safe_imshow(ax, data, comet_extent, **extra)
                add_apertures(ax, env)
            else:
                ax.text(0.5, 0.5, f"{inst} missing", ha="center", va="center")
                ax.set_axis_off()
            ax.set_title(f"{inst} rate (comet)", fontsize=9)
            if comet_extent is not None:
                ax.set_xlabel("East offset (arcsec)", fontsize=8)
                ax.set_ylabel("North offset (arcsec)", fontsize=8)
        ax = axes[0][3]
        epic_rate_cm = banddir / "EPIC_rate.fits"
        if epic_rate_cm.exists():
            data, _ = read_image(str(epic_rate_cm))
            extra = {}
            finite = np.isfinite(data)
            if np.any(finite):
                vmin, vmax = np.nanpercentile(data[finite], [1.0, 99.0])
                extra["vmin"], extra["vmax"] = (vmin, vmax)
            _safe_imshow(ax, data, comet_extent, **extra)
            add_apertures(ax, env)
        else:
            ax.text(0.5, 0.5, "EPIC missing", ha="center", va="center")
            ax.set_axis_off()
        ax.set_title("EPIC rate (comet)", fontsize=9)
        if comet_extent is not None:
            ax.set_xlabel("East offset (arcsec)", fontsize=8)
            ax.set_ylabel("North offset (arcsec)", fontsize=8)
        for ax, inst in zip(axes[1][:3], instruments):
            bpath = banddir / f"{inst}_bkg.fits"
            if bpath.exists():
                data, _ = read_image(str(bpath))
                _safe_imshow(
                    ax, data, comet_extent, norm=nice_norm(data, stretch="sqrt")
                )
                add_apertures(ax, env)
            else:
                ax.text(0.5, 0.5, f"{inst} bkg missing", ha="center", va="center")
                ax.set_axis_off()
            ax.set_title(f"{inst} background (comet)", fontsize=9)
            if comet_extent is not None:
                ax.set_xlabel("East offset (arcsec)", fontsize=8)
                ax.set_ylabel("North offset (arcsec)", fontsize=8)
        ax = axes[1][3]
        masked_cm = banddir / "EPIC_rate_masked.fits"
        if masked_cm.exists():
            data, _ = read_image(str(masked_cm))
            extra = {}
            finite = np.isfinite(data)
            if np.any(finite):
                vmin, vmax = np.nanpercentile(data[finite], [1.0, 99.0])
                extra["vmin"], extra["vmax"] = (vmin, vmax)
            _safe_imshow(ax, data, comet_extent, **extra)
            add_apertures(ax, env)
        else:
            ax.text(0.5, 0.5, "EPIC masked rate missing", ha="center", va="center")
            ax.set_axis_off()
        ax.set_title("EPIC masked rate (comet)", fontsize=9)
        if comet_extent is not None:
            ax.set_xlabel("East offset (arcsec)", fontsize=8)
            ax.set_ylabel("North offset (arcsec)", fontsize=8)
        fig.suptitle(
            f"Comet-frame mosaic diagnostics ({band})", fontsize=13, fontweight="bold"
        )
        png = outdir / f"08_comet_mosaic_{band}.png"
        savefig(fig, str(png))
        rendered.append(png)
        if first_stats is None:
            first_stats = {
                "status": "ok",
                "png": png.name,
                "band": band,
                "frame": "comet",
            }
    if not rendered or first_stats is None:
        return {"status": "missing", "reason": "No mosaic products found"}
    _copy_png(rendered[0], outdir / "08_mosaics.png")
    if len(rendered) > 1:
        first_stats["extra_pngs"] = [p.name for p in rendered[1:]]
    first_stats["total_panels"] = len(rendered)
    return first_stats


def md_kv(lines: list[str], key: str, val: Any) -> None:
    if val is None:
        return
    if isinstance(val, float):
        if math.isnan(val):
            return
        lines.append(f"- {key}: {val:.6g}")
    else:
        lines.append(f"- {key}: {val}")


def write_report(
    path: Path, env: dict[str, str], results: dict[str, dict[str, Any]]
) -> None:
    lines: list[str] = []
    lines.append("# XMM comet pipeline quick-check report")
    lines.append("")
    lines.append(f"- WORKDIR: `{env.get('WORKDIR', '')}`")
    lines.append(f"- Generated by: `xmm_comet_quick_checks.py`")
    lines.append("")
    order = [
        ("clean", "Cleaning / flare filtering"),
        ("exposures", "Pseudo-exposure diagnostics"),
        ("track", "Track / ephemeris"),
        ("detect", "Still-sky detection"),
        ("source-filter", "Source filtering diagnostics"),
        ("contam", "Contamination GTI"),
        ("image", "Comet-frame imaging"),
        ("image-per-inst", "Per-detector imaging"),
        ("radial-profile", "Radial brightness profile"),
        ("frames", "Frame comparison (sky vs comet)"),
        ("lcurve", "Light curves"),
        ("spectrum", "Spectrum"),
        ("pointings", "Per-time-slice mosaic"),
        ("pointings-per-inst", "Per-time-slice mosaic (per detector)"),
        ("mosaics", "QC mosaics (still-sky + comet)"),
        ("spotcheck", "Comet-frame spot-checks"),
        ("region-support", "Region exposure support"),
    ]
    for key, title in order:
        res = results.get(key)
        if not res:
            continue
        lines.append(f"## {title}")
        lines.append("")
        lines.append(f"- status: {res.get('status', 'unknown')}")
        if res.get("status") != "ok":
            md_kv(lines, "reason", res.get("reason"))
            lines.append("")
            continue
        for k, v in res.items():
            if k in {"status", "png", "extra_pngs", "top_report_rows"}:
                continue
            md_kv(lines, k, v)
        if res.get("png"):
            png = res["png"]
            lines.append("")
            lines.append(f"![{title}]({png})")
        for extra_png in res.get("extra_pngs", []):
            lines.append("")
            lines.append(f"![{title}]({extra_png})")
        if res.get("top_report_rows"):
            lines.append("")
            lines.append("Top contamination-report rows:")
            lines.append("")
            for row in res["top_report_rows"]:
                lines.append(
                    f"- srcid={row.get('srcid')} min_sep_src={row.get('min_sep_src_arcsec')} min_sep_bkg={row.get('min_sep_bkg_arcsec')} src_overlap={row.get('samples_src_overlap')} bkg_overlap={row.get('samples_bkg_overlap')} strict_overlap={row.get('samples_strict_overlap')}"
                )
        lines.append("")
    path.write_text("\n".join(lines).rstrip() + "\n", encoding="utf-8")


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True, help="Pipeline config.env file")
    ap.add_argument(
        "--outdir", default="", help="Output QA directory (default: WORKDIR/qc)"
    )
    ap.add_argument(
        "--band",
        help="Image band label for image/spotcheck/region-support checks (default: first IMAGE_BANDS entry)",
    )
    ap.add_argument(
        "--inst",
        default="EPIC",
        choices=["EPIC", "PN", "M1", "M2"],
        help="Instrument for region-support check",
    )
    ap.add_argument(
        "--horizons",
        action="store_true",
        help="Compare stored track against a fresh Horizons query when possible",
    )
    ap.add_argument(
        "checks",
        nargs="+",
        choices=[
            "clean",
            "exposures",
            "track",
            "detect",
            "source-filter",
            "image",
            "image-per-inst",
            "radial-profile",
            "frames",
            "contam",
            "lcurve",
            "spectrum",
            "pointings",
            "pointings-per-inst",
            "mosaics",
            "spotcheck",
            "region-support",
            "all",
        ],
    )
    return ap.parse_args()


def main() -> int:
    args = parse_args()
    env = load_shell_env(args.config)
    if "WORKDIR" not in env or not env["WORKDIR"].strip():
        raise RuntimeError("WORKDIR is not defined in the config env file")
    outdir = Path(args.outdir) if args.outdir else Path(env["WORKDIR"]) / "qc"
    outdir.mkdir(parents=True, exist_ok=True)
    checks = args.checks
    if "all" in checks:
        checks = [
            "clean",
            "exposures",
            "track",
            "detect",
            "source-filter",
            "contam",
            "image",
            "image-per-inst",
            "radial-profile",
            "frames",
            "lcurve",
            "spectrum",
            "pointings",
            "pointings-per-inst",
            "mosaics",
            "spotcheck",
            "region-support",
        ]
    results: dict[str, dict[str, Any]] = {}
    for check in checks:
        if check == "clean":
            results[check] = check_clean(env, outdir)
        elif check == "exposures":
            results[check] = check_exposures(env, outdir)
        elif check == "track":
            results[check] = check_track(env, outdir, use_horizons=args.horizons)
        elif check == "detect":
            results[check] = check_detect(env, outdir)
        elif check == "source-filter":
            results[check] = check_source_filter(env, outdir)
        elif check == "image":
            results[check] = check_image(env, outdir, band_label=args.band)
        elif check == "image-per-inst":
            results[check] = check_image_per_inst(env, outdir, band_label=args.band)
        elif check == "radial-profile":
            results[check] = check_radial_profile(env, outdir, band_label=args.band)
        elif check == "frames":
            results[check] = check_frames(env, outdir, band_label=args.band)
        elif check == "contam":
            results[check] = check_contam(env, outdir)
        elif check == "lcurve":
            results[check] = check_lcurve(env, outdir)
        elif check == "spectrum":
            results[check] = check_spectrum(env, outdir)
        elif check == "pointings":
            results[check] = check_pointings(env, outdir, band_label=args.band)
        elif check == "pointings-per-inst":
            results[check] = check_pointings_per_inst(env, outdir, band_label=args.band)
        elif check == "mosaics":
            results[check] = check_mosaics(env, outdir, band_label=args.band)
        elif check == "spotcheck":
            results[check] = check_spotcheck(env, outdir, band=args.band)
        elif check == "region-support":
            results[check] = check_region_support(
                env, outdir, band=args.band, inst=args.inst
            )
        else:
            raise RuntimeError(f"Unknown check: {check}")
    report_path = outdir / "qc_report.md"
    write_report(report_path, env, results)
    print(f"Wrote QA outputs to: {outdir}")
    print(f"Report: {report_path}")
    for name, res in results.items():
        print(f"{name}: {res.get('status')}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(1)
