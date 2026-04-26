#!/usr/bin/env python3
"""Helpers for the reduction_v9 comet-frame stack stage."""

from __future__ import annotations

import argparse
import json
import math
import re
import sys
from pathlib import Path
from typing import Iterable


def die(message: str) -> None:
    raise SystemExit(message)


def event_hdu(hdul):
    if "EVENTS" in hdul:
        return hdul["EVENTS"]
    for hdu in hdul:
        if getattr(hdu, "data", None) is not None and hasattr(hdu, "columns"):
            return hdu
    raise RuntimeError("No EVENTS-like table HDU found")


def first_image_hdu(hdul):
    for hdu in hdul:
        data = getattr(hdu, "data", None)
        if data is not None and getattr(data, "ndim", 0) >= 2:
            return hdu
    raise RuntimeError("No image HDU found")


def read_manifest(manifest_arg: str) -> list[Path]:
    manifest = Path(manifest_arg)
    if not manifest.is_file():
        die(f"Manifest not found: {manifest}")
    paths: list[Path] = []
    seen: set[Path] = set()
    for line in manifest.read_text(encoding="utf-8").splitlines():
        text = line.strip()
        if not text:
            continue
        path = Path(text)
        if path not in seen:
            paths.append(path)
            seen.add(path)
    if not paths:
        die(f"No paths listed in manifest: {manifest}")
    return paths


def parse_grid_env(path_arg: str) -> dict[str, float | int]:
    path = Path(path_arg)
    if not path.is_file():
        die(f"Grid env not found: {path}")
    values: dict[str, str] = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if not line or line.startswith("#") or "=" not in line:
            continue
        key, value = line.split("=", 1)
        values[key.strip()] = value.strip().strip("'").strip('"')
    try:
        return {
            "bin": float(values["MAP_GRID_BIN_PHYS"]),
            "x_min": float(values["MAP_GRID_X_MIN"]),
            "x_max": float(values["MAP_GRID_X_MAX"]),
            "y_min": float(values["MAP_GRID_Y_MIN"]),
            "y_max": float(values["MAP_GRID_Y_MAX"]),
            "nx": int(values["MAP_GRID_NX"]),
            "ny": int(values["MAP_GRID_NY"]),
        }
    except KeyError as exc:
        die(f"Grid env missing key {exc.args[0]}: {path}")


def read_image(path_arg: str):
    import numpy as np
    from astropy.io import fits

    path = Path(path_arg)
    if not path.is_file():
        die(f"Image not found: {path}")
    with fits.open(path, memmap=False) as hdul:
        hdu = first_image_hdu(hdul)
        data = np.asarray(hdu.data, dtype=float)
        while data.ndim > 2:
            data = data[0]
        return data, hdu.header.copy()


def write_image(template_header, data, out_arg: str, bunit: str | None = None) -> None:
    from astropy.io import fits

    out = Path(out_arg)
    out.parent.mkdir(parents=True, exist_ok=True)
    header = template_header.copy()
    if bunit:
        header["BUNIT"] = bunit
    fits.PrimaryHDU(data=data, header=header).writeto(out, overwrite=True)


def _find_col_case_insensitive(
    names: Iterable[str], candidates: Iterable[str]
) -> str | None:
    lut = {name.lower(): name for name in names}
    for cand in candidates:
        if cand.lower() in lut:
            return lut[cand.lower()]
    return None


def read_event_time_info(path_arg: str) -> tuple[float, float, float, str]:
    import numpy as np
    from astropy.io import fits

    path = Path(path_arg)
    if not path.is_file():
        die(f"Event file not found: {path}")
    with fits.open(path, memmap=False) as hdul:
        hdu = event_hdu(hdul)
        hdr = hdu.header
        data = hdu.data
        have_times = False
        if data is not None and len(data) > 0 and "TIME" in hdu.columns.names:
            times = np.asarray(data["TIME"], dtype=float)
            good = np.isfinite(times)
            if np.any(good):
                tmin = float(times[good].min())
                tmax = float(times[good].max())
                have_times = True
        if not have_times:
            if "TSTART" not in hdr or "TSTOP" not in hdr:
                raise RuntimeError(f"Could not determine time range for {path}")
            tmin = float(hdr["TSTART"])
            tmax = float(hdr["TSTOP"])
        mjdref = hdr.get("MJDREF")
        if mjdref is None:
            mjdref = float(hdr.get("MJDREFI", 0.0)) + float(
                hdr.get("MJDREFF", 0.0)
            )
        timesys = str(hdr.get("TIMESYS", "TT")).strip().upper()
        return tmin, tmax, float(mjdref), timesys


def event_time_mid(path_arg: str) -> tuple[float, float, str]:
    tmin, tmax, mjdref, timesys = read_event_time_info(path_arg)
    return 0.5 * (tmin + tmax), mjdref, timesys


def observation_window(manifest_arg: str) -> tuple[list[Path], float, float, float, str]:
    paths = read_manifest(manifest_arg)
    tmin = math.inf
    tmax = -math.inf
    mjdref = None
    timesys = None
    for path in paths:
        one_min, one_max, one_mjdref, one_timesys = read_event_time_info(str(path))
        tmin = min(tmin, one_min)
        tmax = max(tmax, one_max)
        if mjdref is None:
            mjdref = one_mjdref
            timesys = one_timesys
        else:
            if abs(one_mjdref - mjdref) > 1.0e-9:
                raise RuntimeError(
                    f"MJDREF mismatch between event files: {path} ({one_mjdref})"
                )
            if one_timesys != timesys:
                raise RuntimeError(
                    f"TIMESYS mismatch between event files: {path} ({one_timesys})"
                )
    if not math.isfinite(tmin) or not math.isfinite(tmax):
        raise RuntimeError(f"Could not determine observation window from {manifest_arg}")
    return paths, float(tmin), float(tmax), float(mjdref), str(timesys)


def xmm_seconds_to_time(seconds, mjdref: float, timesys: str):
    from astropy.time import Time, TimeDelta

    return Time(mjdref, format="mjd", scale=timesys.lower()) + TimeDelta(
        seconds, format="sec"
    )


def _parse_fixed_step_seconds(step: str) -> tuple[float | None, str | None]:
    step_re = re.compile(r"^\s*([0-9]*\.?[0-9]+)\s*([A-Za-z]+)?\s*$")
    match = step_re.match(step.strip())
    if not match:
        return None, None
    value = float(match.group(1))
    unit = (match.group(2) or "").strip().lower()
    if unit == "":
        return None, None
    if unit in {"s", "sec", "secs", "second", "seconds"}:
        return value, "s"
    if unit in {"m", "min", "mins", "minute", "minutes"}:
        return value * 60.0, "m"
    if unit in {"h", "hr", "hrs", "hour", "hours"}:
        return value * 3600.0, "h"
    if unit in {"d", "day", "days"}:
        return value * 86400.0, "d"
    return None, None


def _build_epoch_list(start_utc: str, stop_utc: str, step_seconds: float):
    import numpy as np
    from astropy.time import Time, TimeDelta
    import astropy.units as u

    if step_seconds <= 0.0:
        raise RuntimeError("Track step must be > 0 seconds")
    t0 = Time(start_utc, format="isot", scale="utc")
    t1 = Time(stop_utc, format="isot", scale="utc")
    total = max(0.0, (t1 - t0).to_value(u.s))
    offsets = np.arange(0.0, total + 0.5 * step_seconds, step_seconds, dtype=float)
    times = t0 + TimeDelta(offsets, format="sec")
    if times[-1] < t1:
        times = Time(
            np.concatenate([times.jd, np.array([t1.jd], dtype=float)]),
            format="jd",
            scale="utc",
        )
    jd = np.unique(
        np.concatenate(
            [
                np.array([t0.jd], dtype=float),
                np.asarray(times.jd, float),
                np.array([t1.jd], dtype=float),
            ]
        )
    )
    return jd.astype(float)


def _parse_eph_table(eph):
    import numpy as np

    cols = eph.colnames
    racol = _find_col_case_insensitive(cols, ["RA_app", "RA_ICRF_app", "RA"])
    deccol = _find_col_case_insensitive(cols, ["DEC_app", "DEC_ICRF_app", "DEC"])
    dcol = _find_col_case_insensitive(cols, ["delta"])
    tcol = _find_col_case_insensitive(cols, ["datetime_jd", "JD", "jd"])
    if not (racol and deccol and dcol and tcol):
        raise RuntimeError(
            f"Horizons output missing RA/DEC/delta/time columns; got {cols}"
        )
    jd = np.asarray(eph[tcol], dtype=float)
    return (
        jd - 2400000.5,
        np.asarray(eph[racol], dtype=float),
        np.asarray(eph[deccol], dtype=float),
        np.asarray(eph[dcol], dtype=float),
    )


def query_horizons_track(
    target_id: str,
    id_type: str,
    start_utc: str,
    stop_utc: str,
    step: str,
):
    import numpy as np

    try:
        from astroquery.jplhorizons import Horizons
    except Exception as exc:  # pragma: no cover
        raise RuntimeError(
            "astroquery.jplhorizons is unavailable; supply track_input instead"
        ) from exc

    step_seconds, canonical_unit = _parse_fixed_step_seconds(step)
    use_explicit_epochs = bool(step_seconds is not None and step_seconds < 60.0)
    if step_seconds is not None and canonical_unit in {"m", "h", "d"}:
        base = {"m": 60.0, "h": 3600.0, "d": 86400.0}[canonical_unit]
        use_explicit_epochs = use_explicit_epochs or (
            abs(step_seconds / base - round(step_seconds / base)) > 1.0e-12
        )

    if use_explicit_epochs:
        epochs_jd = _build_epoch_list(start_utc, stop_utc, float(step_seconds))
        mjd_parts: list[np.ndarray] = []
        ra_parts: list[np.ndarray] = []
        dec_parts: list[np.ndarray] = []
        delta_parts: list[np.ndarray] = []
        chunk = 40
        for start in range(0, len(epochs_jd), chunk):
            one_epochs = [float(x) for x in epochs_jd[start : start + chunk]]
            obj = Horizons(
                id=target_id,
                id_type=id_type,
                location="500@399",
                epochs=one_epochs,
            )
            eph = obj.ephemerides(quantities="2,20", extra_precision=True)
            mjd, ra, dec, delta = _parse_eph_table(eph)
            mjd_parts.append(mjd)
            ra_parts.append(ra)
            dec_parts.append(dec)
            delta_parts.append(delta)
        return (
            np.concatenate(mjd_parts),
            np.concatenate(ra_parts),
            np.concatenate(dec_parts),
            np.concatenate(delta_parts),
        )

    epochs = {
        "start": start_utc,
        "stop": stop_utc,
        "step": step.strip(),
    }
    obj = Horizons(
        id=target_id,
        id_type=id_type,
        location="500@399",
        epochs=epochs,
    )
    eph = obj.ephemerides(quantities="2,20", extra_precision=True)
    return _parse_eph_table(eph)


def load_track_input(track_input: str):
    import numpy as np
    from astropy.io import fits
    from astropy.time import Time

    path = Path(track_input)
    if not path.is_file():
        die(f"Track input not found: {path}")
    if path.suffix.lower() in {".csv", ".txt"}:
        data = np.genfromtxt(path, delimiter=",", names=True, dtype=None, encoding=None)
        if data.size == 0:
            raise RuntimeError(f"Track CSV is empty: {path}")
        names = data.dtype.names or ()
        tcol = _find_col_case_insensitive(names, ["time_iso", "utc_iso", "isot", "time"])
        racol = _find_col_case_insensitive(names, ["ra_deg", "ra"])
        deccol = _find_col_case_insensitive(names, ["dec_deg", "dec"])
        dcol = _find_col_case_insensitive(names, ["delta_au", "delta"])
        if not (tcol and racol and deccol and dcol):
            raise RuntimeError(
                "Track CSV must contain time_iso, ra_deg, dec_deg, delta_au columns"
            )
        tt = Time(np.asarray(data[tcol]), format="isot", scale="utc")
        return (
            np.asarray(tt.mjd, dtype=float),
            np.asarray(data[racol], dtype=float),
            np.asarray(data[deccol], dtype=float),
            np.asarray(data[dcol], dtype=float),
        )

    with fits.open(path, memmap=False) as hdul:
        tab = None
        for hdu in hdul[1:]:
            if getattr(hdu, "data", None) is not None and len(hdu.data) > 0:
                tab = hdu.data
                break
        if tab is None:
            raise RuntimeError(f"No non-empty table found in {path}")
        names = tab.dtype.names or ()
        mjdcol = _find_col_case_insensitive(names, ["MJD"])
        racol = _find_col_case_insensitive(names, ["RA"])
        deccol = _find_col_case_insensitive(names, ["DEC"])
        dcol = _find_col_case_insensitive(names, ["DELTA"])
        if not (mjdcol and racol and deccol and dcol):
            raise RuntimeError("Track FITS must contain MJD, RA, DEC, DELTA columns")
        return (
            np.asarray(tab[mjdcol], dtype=float),
            np.asarray(tab[racol], dtype=float),
            np.asarray(tab[deccol], dtype=float),
            np.asarray(tab[dcol], dtype=float),
        )


def write_track_fits(
    out_arg: str,
    mjd,
    ra,
    dec,
    delta,
    source_desc: str,
    obs_tmin_sec: float,
    obs_tmax_sec: float,
    mjdref: float,
    timesys: str,
) -> None:
    import numpy as np
    from astropy.io import fits

    out = Path(out_arg)
    out.parent.mkdir(parents=True, exist_ok=True)
    cols = [
        fits.Column(name="MJD", format="D", array=np.asarray(mjd, dtype=float)),
        fits.Column(name="RA", format="D", array=np.asarray(ra, dtype=float)),
        fits.Column(name="DEC", format="D", array=np.asarray(dec, dtype=float)),
        fits.Column(name="DELTA", format="D", array=np.asarray(delta, dtype=float)),
    ]
    hdu = fits.BinTableHDU.from_columns(cols, name="OBJTRACK")
    hdu.header["ORIGINTRK"] = (source_desc[:68], "Track source")
    hdu.header["MJDREF"] = (float(mjdref), "XMM event-file MJD reference")
    hdu.header["TIMESYS"] = (str(timesys), "Time scale of source event list")
    hdu.header["TSTART"] = (float(obs_tmin_sec), "Observation start in XMM seconds")
    hdu.header["TSTOP"] = (float(obs_tmax_sec), "Observation stop in XMM seconds")
    fits.HDUList([fits.PrimaryHDU(), hdu]).writeto(out, overwrite=True)


def write_ref_env(out_arg: str, ref_ra: float, ref_dec: float) -> None:
    out = Path(out_arg)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(
        f'COMET_REF_RA="{ref_ra:.10f}"\nCOMET_REF_DEC="{ref_dec:.10f}"\n',
        encoding="utf-8",
    )


def interpolate_track(mjd, ra, dec, when_mjd: float) -> tuple[float, float]:
    import numpy as np

    mjd_arr = np.asarray(mjd, dtype=float)
    order = np.argsort(mjd_arr)
    mjd_arr = mjd_arr[order]
    ra_arr = np.unwrap(np.deg2rad(np.asarray(ra, dtype=float)[order]))
    dec_arr = np.asarray(dec, dtype=float)[order]
    if when_mjd < mjd_arr[0] or when_mjd > mjd_arr[-1]:
        raise RuntimeError(
            f"Requested time {when_mjd:.8f} outside track range "
            f"[{mjd_arr[0]:.8f}, {mjd_arr[-1]:.8f}]"
        )
    ra_deg = float(np.rad2deg(np.interp(when_mjd, mjd_arr, ra_arr)))
    dec_deg = float(np.interp(when_mjd, mjd_arr, dec_arr))
    return ra_deg % 360.0, dec_deg


def build_track_command(argv: list[str]) -> int:
    import numpy as np

    ap = argparse.ArgumentParser(prog="stack_tools.py build-track")
    ap.add_argument("--event-manifest", required=True)
    ap.add_argument("--track-out", required=True)
    ap.add_argument("--ref-env-out", required=True)
    ap.add_argument("--summary-out")
    ap.add_argument("--track-input", default="")
    ap.add_argument("--target-id", default="")
    ap.add_argument("--id-type", default="smallbody")
    ap.add_argument("--step", default="30s")
    ap.add_argument("--ref-ra", default="")
    ap.add_argument("--ref-dec", default="")
    args = ap.parse_args(argv)

    _paths, tmin, tmax, mjdref, timesys = observation_window(args.event_manifest)
    start_utc = xmm_seconds_to_time(tmin, mjdref, timesys).utc.isot
    stop_utc = xmm_seconds_to_time(tmax, mjdref, timesys).utc.isot

    if args.track_input.strip():
        mjd, ra, dec, delta = load_track_input(args.track_input.strip())
        source_desc = f"file:{Path(args.track_input).name}"
    else:
        if not args.target_id.strip():
            die("Set target_id in config.json or provide track_input before running stack")
        mjd, ra, dec, delta = query_horizons_track(
            target_id=args.target_id.strip(),
            id_type=(args.id_type or "smallbody").strip(),
            start_utc=start_utc,
            stop_utc=stop_utc,
            step=args.step,
        )
        source_desc = f"Horizons:{args.target_id.strip()}"

    mjd = np.asarray(mjd, dtype=float)
    ra = np.asarray(ra, dtype=float)
    dec = np.asarray(dec, dtype=float)
    delta = np.asarray(delta, dtype=float)
    order = np.argsort(mjd)
    mjd = mjd[order]
    ra = ra[order]
    dec = dec[order]
    delta = delta[order]
    keep = np.ones(len(mjd), dtype=bool)
    if len(mjd) > 1:
        keep[1:] = np.diff(mjd) > 0
    mjd = mjd[keep]
    ra = ra[keep]
    dec = dec[keep]
    delta = delta[keep]

    obs_mjd_min = mjdref + tmin / 86400.0
    obs_mjd_max = mjdref + tmax / 86400.0
    if len(mjd) == 0 or mjd[0] > obs_mjd_min or mjd[-1] < obs_mjd_max:
        raise RuntimeError(
            "Track does not cover the full observation window; "
            "provide a denser/longer track_input or a valid Horizons target_id"
        )
    if len(mjd) == 1:
        mjd = np.array([obs_mjd_min, obs_mjd_max], dtype=float)
        ra = np.array([ra[0], ra[0]], dtype=float)
        dec = np.array([dec[0], dec[0]], dtype=float)
        delta = np.array([delta[0], delta[0]], dtype=float)

    mid = len(mjd) // 2
    ref_ra = float(args.ref_ra) if str(args.ref_ra).strip() else float(ra[mid])
    ref_dec = float(args.ref_dec) if str(args.ref_dec).strip() else float(dec[mid])

    write_track_fits(
        args.track_out,
        mjd,
        ra,
        dec,
        delta,
        source_desc,
        tmin,
        tmax,
        mjdref,
        timesys,
    )
    write_ref_env(args.ref_env_out, ref_ra, ref_dec)

    if args.summary_out:
        summary = {
            "source": source_desc,
            "target_id": args.target_id.strip(),
            "track_input": args.track_input.strip(),
            "n_track_rows": int(len(mjd)),
            "observation_tstart": float(tmin),
            "observation_tstop": float(tmax),
            "observation_start_utc": start_utc,
            "observation_stop_utc": stop_utc,
            "track_mjd_min": float(mjd[0]),
            "track_mjd_max": float(mjd[-1]),
            "ref_ra_deg": ref_ra,
            "ref_dec_deg": ref_dec,
            "step": args.step,
        }
        out = Path(args.summary_out)
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")

    print(f"Wrote {Path(args.track_out).resolve()}")
    print(f"Wrote {Path(args.ref_env_out).resolve()}")
    if args.summary_out:
        print(f"Wrote {Path(args.summary_out).resolve()}")
    return 0


def mask_events_by_grid_command(argv: list[str]) -> int:
    import numpy as np
    from astropy.io import fits

    if len(argv) != 4:
        die("mask-events-by-grid EVENT MASK GRID_ENV OUT")
    event_arg, mask_arg, grid_env_arg, out_arg = argv
    grid = parse_grid_env(grid_env_arg)
    mask, _mask_header = read_image(mask_arg)
    if mask.shape != (int(grid["ny"]), int(grid["nx"])):
        raise RuntimeError(
            f"Mask shape {mask.shape} does not match grid {(grid['ny'], grid['nx'])}"
        )

    event_path = Path(event_arg)
    if not event_path.is_file():
        die(f"Event file not found: {event_path}")
    out_path = Path(out_arg)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    with fits.open(event_path, memmap=False) as hdul:
        evt_ref = event_hdu(hdul)
        evt_idx = next(idx for idx, hdu in enumerate(hdul) if hdu is evt_ref)
        evt_hdu = hdul[evt_idx]
        data = evt_hdu.data
        if data is None or len(data) == 0:
            keep = np.zeros(0, dtype=bool)
        else:
            x = np.asarray(data["X"], dtype=float)
            y = np.asarray(data["Y"], dtype=float)
            ix = np.floor((x - float(grid["x_min"])) / float(grid["bin"])).astype(int)
            iy = np.floor((y - float(grid["y_min"])) / float(grid["bin"])).astype(int)
            valid = (
                np.isfinite(x)
                & np.isfinite(y)
                & (ix >= 0)
                & (ix < int(grid["nx"]))
                & (iy >= 0)
                & (iy < int(grid["ny"]))
            )
            keep = np.zeros(len(data), dtype=bool)
            if np.any(valid):
                keep[valid] = mask[iy[valid], ix[valid]] > 0

        out_hdus = []
        for idx, hdu in enumerate(hdul):
            if idx == evt_idx:
                filtered = data[keep] if data is not None else data
                out_hdu = fits.BinTableHDU(
                    data=filtered,
                    header=hdu.header.copy(),
                    name=hdu.name,
                )
            elif isinstance(hdu, fits.PrimaryHDU):
                out_hdu = fits.PrimaryHDU(
                    data=None if hdu.data is None else np.array(hdu.data, copy=True),
                    header=hdu.header.copy(),
                )
            elif isinstance(hdu, fits.ImageHDU):
                out_hdu = fits.ImageHDU(
                    data=None if hdu.data is None else np.array(hdu.data, copy=True),
                    header=hdu.header.copy(),
                    name=hdu.name,
                )
            elif isinstance(hdu, fits.BinTableHDU):
                out_hdu = fits.BinTableHDU(
                    data=None if hdu.data is None else hdu.data.copy(),
                    header=hdu.header.copy(),
                    name=hdu.name,
                )
            else:
                out_hdu = hdu.copy()
            out_hdus.append(out_hdu)
        fits.HDUList(out_hdus).writeto(out_path, overwrite=True)

    print(0 if data is None else int(keep.sum()))
    return 0


def set_event_nominal_command(argv: list[str]) -> int:
    from astropy.io import fits

    if len(argv) != 3:
        die("set-event-nominal EVENT RA DEC")
    event_arg, ra_arg, dec_arg = argv
    ra = float(ra_arg)
    dec = float(dec_arg)
    with fits.open(event_arg, mode="update", memmap=False) as hdul:
        for hdu in hdul:
            header = hdu.header
            if "RA_NOM" in header:
                header["RA_NOM"] = ra
            if "DEC_NOM" in header:
                header["DEC_NOM"] = dec
            if "RA_PNT" in header:
                header["RA_PNT"] = ra
            if "DEC_PNT" in header:
                header["DEC_PNT"] = dec
        hdul.flush()
    return 0


def measure_shift_command(argv: list[str]) -> int:
    import numpy as np
    import astropy.units as u
    from astropy.coordinates import SkyCoord
    from astropy.io import fits

    if len(argv) != 5:
        die("measure-shift SKY_EVENT COMET_EVENT TRACK_FITS REF_RA REF_DEC")
    sky_arg, comet_arg, track_arg, ref_ra_arg, ref_dec_arg = argv
    ref_ra = float(ref_ra_arg)
    ref_dec = float(ref_dec_arg)

    dx = dy = 0.0
    method = "events"
    n_used = 0
    t_mid, mjdref, _timesys = event_time_mid(sky_arg)

    with fits.open(sky_arg, memmap=False) as sky_hdul, fits.open(
        comet_arg, memmap=False
    ) as comet_hdul:
        sky_evt = event_hdu(sky_hdul).data
        comet_evt = event_hdu(comet_hdul).data
        n_sky = 0 if sky_evt is None else len(sky_evt)
        n_comet = 0 if comet_evt is None else len(comet_evt)
        if n_sky > 0 and n_comet > 0:
            n = min(n_sky, n_comet)
            sx = np.asarray(sky_evt["X"][:n], dtype=float)
            sy = np.asarray(sky_evt["Y"][:n], dtype=float)
            cx = np.asarray(comet_evt["X"][:n], dtype=float)
            cy = np.asarray(comet_evt["Y"][:n], dtype=float)
            good = np.isfinite(sx) & np.isfinite(sy) & np.isfinite(cx) & np.isfinite(cy)
            if np.any(good):
                dx = float(np.median(cx[good] - sx[good]))
                dy = float(np.median(cy[good] - sy[good]))
                n_used = int(np.count_nonzero(good))
            else:
                method = "track"
        else:
            method = "track"

    if method == "track":
        with fits.open(track_arg, memmap=False) as hdul:
            tab = None
            for hdu in hdul[1:]:
                if getattr(hdu, "data", None) is not None and len(hdu.data) > 0:
                    tab = hdu.data
                    break
            if tab is None:
                raise RuntimeError(f"No non-empty track table found in {track_arg}")
            mjd = np.asarray(tab["MJD"], dtype=float)
            ra = np.asarray(tab["RA"], dtype=float)
            dec = np.asarray(tab["DEC"], dtype=float)
        when_mjd = mjdref + t_mid / 86400.0
        comet_ra, comet_dec = interpolate_track(mjd, ra, dec, when_mjd)
        ref = SkyCoord(ra=ref_ra * u.deg, dec=ref_dec * u.deg, frame="icrs")
        comet = SkyCoord(ra=comet_ra * u.deg, dec=comet_dec * u.deg, frame="icrs")
        dlon, dlat = ref.spherical_offsets_to(comet)
        arcsec_per_phys = 0.05
        dx = -float(dlon.to_value(u.arcsec)) / arcsec_per_phys
        dy = -float(dlat.to_value(u.arcsec)) / arcsec_per_phys
        n_used = 0

    print(f"{dx:.6f}\t{dy:.6f}\t{method}\t{n_used}\t{t_mid:.6f}")
    return 0


def shift_image_grid_command(argv: list[str]) -> int:
    import numpy as np

    if len(argv) != 7:
        die("shift-image-grid IMAGE SRC_GRID_ENV DST_GRID_ENV DX DY TEMPLATE OUT")
    image_arg, src_grid_arg, dst_grid_arg, dx_arg, dy_arg, template_arg, out_arg = (
        argv[0],
        argv[1],
        argv[2],
        argv[3],
        argv[4],
        argv[5],
        argv[6],
    )
    source, src_header = read_image(image_arg)
    src_grid = parse_grid_env(src_grid_arg)
    dst_grid = parse_grid_env(dst_grid_arg)
    dx = float(dx_arg)
    dy = float(dy_arg)
    _template_data, template_header = read_image(template_arg)

    if source.shape != (int(src_grid["ny"]), int(src_grid["nx"])):
        raise RuntimeError(
            f"Source image shape {source.shape} does not match source grid "
            f"{(src_grid['ny'], src_grid['nx'])}"
        )
    if abs(float(src_grid["bin"]) - float(dst_grid["bin"])) > 1.0e-9:
        raise RuntimeError(
            "shift-image-grid currently requires identical source/destination bin sizes"
        )

    src_bin = float(src_grid["bin"])
    x_src = (
        float(dst_grid["x_min"])
        + (np.arange(int(dst_grid["nx"]), dtype=float) + 0.5) * float(dst_grid["bin"])
        - dx
    )
    y_src = (
        float(dst_grid["y_min"])
        + (np.arange(int(dst_grid["ny"]), dtype=float) + 0.5) * float(dst_grid["bin"])
        - dy
    )
    ix = np.floor((x_src - float(src_grid["x_min"])) / src_bin).astype(int)
    iy = np.floor((y_src - float(src_grid["y_min"])) / src_bin).astype(int)
    valid_x = (ix >= 0) & (ix < int(src_grid["nx"]))
    valid_y = (iy >= 0) & (iy < int(src_grid["ny"]))

    shifted = np.zeros((int(dst_grid["ny"]), int(dst_grid["nx"])), dtype=float)
    if np.any(valid_x) and np.any(valid_y):
        y_idx = np.where(valid_y)[0]
        x_idx = np.where(valid_x)[0]
        sampled = source[np.ix_(iy[y_idx], ix[x_idx])]
        shifted[np.ix_(y_idx, x_idx)] = np.where(np.isfinite(sampled), sampled, 0.0)

    bunit = src_header.get("BUNIT")
    write_image(template_header, shifted, out_arg, bunit=bunit)
    return 0


def combine_stack_command(argv: list[str]) -> int:
    import numpy as np

    if len(argv) != 3:
        die("combine-stack STACK_MANIFEST FINAL_DIR PRODUCTS_TSV")
    manifest_arg, final_dir_arg, products_arg = argv
    manifest = Path(manifest_arg)
    if not manifest.is_file():
        die(f"Stack manifest not found: {manifest}")
    lines = manifest.read_text(encoding="utf-8").splitlines()
    if not lines:
        die(f"Empty stack manifest: {manifest}")
    header = lines[0].split("\t")
    index = {name: idx for idx, name in enumerate(header)}
    for required in (
        "band",
        "shifted_counts",
        "shifted_exposure_vig",
        "shifted_background",
    ):
        if required not in index:
            die(f"Stack manifest missing column '{required}': {manifest}")

    rows = [line.split("\t") for line in lines[1:] if line.strip()]
    if not rows:
        die(f"No rows in stack manifest: {manifest}")

    final_dir = Path(final_dir_arg)
    final_dir.mkdir(parents=True, exist_ok=True)
    summary = ["band\tcounts\texposure_vig\tbackground\tresidual_counts\tcorrected"]

    bands = []
    for row in rows:
        band = row[index["band"]]
        if band not in bands:
            bands.append(band)

    for band in bands:
        band_rows = [row for row in rows if row[index["band"]] == band]
        counts_sum = exp_sum = bg_sum = None
        template_header = None
        for row in band_rows:
            counts, header = read_image(row[index["shifted_counts"]])
            exp, _ = read_image(row[index["shifted_exposure_vig"]])
            bg, _ = read_image(row[index["shifted_background"]])
            if counts.shape != exp.shape or counts.shape != bg.shape:
                raise RuntimeError(
                    f"Shape mismatch while combining stack products for {band}"
                )
            if counts_sum is None:
                counts_sum = np.zeros_like(counts, dtype=float)
                exp_sum = np.zeros_like(exp, dtype=float)
                bg_sum = np.zeros_like(bg, dtype=float)
                template_header = header
            counts_sum += np.nan_to_num(counts, nan=0.0, posinf=0.0, neginf=0.0)
            exp_sum += np.nan_to_num(exp, nan=0.0, posinf=0.0, neginf=0.0)
            bg_sum += np.nan_to_num(bg, nan=0.0, posinf=0.0, neginf=0.0)

        if counts_sum is None or exp_sum is None or bg_sum is None or template_header is None:
            continue

        residual = counts_sum - bg_sum
        corrected = np.full_like(residual, np.nan, dtype=float)
        good = np.isfinite(exp_sum) & (exp_sum > 0)
        corrected[good] = residual[good] / exp_sum[good]

        band_dir = final_dir / band
        band_dir.mkdir(parents=True, exist_ok=True)
        counts_out = band_dir / f"EPIC_{band}_stack_counts.fits"
        exp_out = band_dir / f"EPIC_{band}_stack_exposure_vig.fits"
        bg_out = band_dir / f"EPIC_{band}_stack_background.fits"
        residual_out = band_dir / f"EPIC_{band}_stack_residual_counts.fits"
        corrected_out = band_dir / f"EPIC_{band}_stack_corrected.fits"

        write_image(template_header, counts_sum, str(counts_out), bunit="counts")
        write_image(template_header, exp_sum, str(exp_out), bunit="s")
        write_image(template_header, bg_sum, str(bg_out), bunit="counts")
        write_image(template_header, residual, str(residual_out), bunit="counts")
        write_image(template_header, corrected, str(corrected_out), bunit="counts / s")

        summary.append(
            f"{band}\t{counts_out.resolve()}\t{exp_out.resolve()}\t{bg_out.resolve()}\t"
            f"{residual_out.resolve()}\t{corrected_out.resolve()}"
        )

    Path(products_arg).write_text("\n".join(summary) + "\n", encoding="utf-8")
    print(f"Wrote {Path(products_arg).resolve()}")
    return 0


def stack_qc_command(argv: list[str]) -> int:
    if len(argv) != 2:
        die("stack-qc STACK_PRODUCTS_TSV OUTDIR")
    products_arg, outdir_arg = argv
    products = Path(products_arg)
    if not products.is_file():
        die(f"Stack products table not found: {products}")
    lines = products.read_text(encoding="utf-8").splitlines()
    if not lines:
        die(f"Empty stack products table: {products}")
    header = lines[0].split("\t")
    index = {name: idx for idx, name in enumerate(header)}
    rows = [line.split("\t") for line in lines[1:] if line.strip()]
    if not rows:
        die(f"No stack products listed in {products}")

    import tools as qc_tools

    outdir = Path(outdir_arg)
    outdir.mkdir(parents=True, exist_ok=True)
    summary = ["band\tproduct\tpng"]
    for row in rows:
        band = row[index["band"]]
        counts = Path(row[index["counts"]])
        exp = Path(row[index["exposure_vig"]])
        bg = Path(row[index["background"]])
        residual = Path(row[index["residual_counts"]])
        corrected = Path(row[index["corrected"]])

        images = [
            ("stack_counts", counts, "linear_percentile", f"{band} stacked counts (shifted cheese maps)"),
            ("stack_exposure_vig", exp, "linear_percentile", f"{band} stacked vignetted exposure (shifted cheese maps)"),
            ("stack_background", bg, "linear_percentile", f"{band} stacked background (shifted cheese maps)"),
            ("stack_residual_counts", residual, "linear_diverging_percentile", f"{band} stacked residual counts (C-B)"),
            ("stack_corrected", corrected, "linear_diverging_percentile", f"{band} stacked corrected ((C-B)/E)"),
        ]
        for suffix, path, scale, label in images:
            image = qc_tools.fits_image(path)
            png = outdir / f"{band}_{suffix}.png"
            qc_tools._save_mosaic_png(png, image, scale, label)
            summary.append(f"{band}\t{suffix}\t{png.resolve()}")

    (outdir / "stack_qc_summary.tsv").write_text(
        "\n".join(summary) + "\n", encoding="utf-8"
    )
    return 0


def usage() -> str:
    return (
        "Usage: stack_tools.py "
        "build-track --event-manifest MANIFEST --track-out TRACK --ref-env-out ENV [--summary-out JSON] "
        "[--track-input PATH | --target-id ID] [--id-type TYPE] [--step 30s] [--ref-ra RA] [--ref-dec DEC] | "
        "mask-events-by-grid EVENT MASK GRID_ENV OUT | "
        "set-event-nominal EVENT RA DEC | "
        "measure-shift SKY_EVENT COMET_EVENT TRACK_FITS REF_RA REF_DEC | "
        "shift-image-grid IMAGE SRC_GRID_ENV DST_GRID_ENV DX DY TEMPLATE OUT | "
        "combine-stack STACK_MANIFEST FINAL_DIR PRODUCTS_TSV | "
        "stack-qc STACK_PRODUCTS_TSV OUTDIR"
    )


def main(argv: list[str]) -> int:
    if len(argv) < 2:
        die(usage())
    cmd = argv[1]
    args = argv[2:]
    if cmd == "build-track":
        return build_track_command(args)
    if cmd == "mask-events-by-grid":
        return mask_events_by_grid_command(args)
    if cmd == "set-event-nominal":
        return set_event_nominal_command(args)
    if cmd == "measure-shift":
        return measure_shift_command(args)
    if cmd == "shift-image-grid":
        return shift_image_grid_command(args)
    if cmd == "combine-stack":
        return combine_stack_command(args)
    if cmd == "stack-qc":
        return stack_qc_command(args)
    die(usage())


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
