#!/usr/bin/env python3
"""Build a simple still-sky diagnostic mosaic from one or more event lists.

This is intended for QC and source-overlay plots. It does not replace the SAS
source-detection products.
"""
from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS

from trail_mask_utils import read_srclist, read_track


def _event_radec(hdu):
    """Return RA, DEC arrays from an EVENTS HDU.

    Uses explicit RA/DEC columns if present, otherwise converts X/Y
    sky-pixel columns via the TAN-projection WCS in the header.
    """
    cols_upper = [c.upper() for c in hdu.columns.names]
    if "RA" in cols_upper and "DEC" in cols_upper:
        ra_col = hdu.columns.names[cols_upper.index("RA")]
        dec_col = hdu.columns.names[cols_upper.index("DEC")]
        return (
            np.asarray(hdu.data[ra_col], dtype=float),
            np.asarray(hdu.data[dec_col], dtype=float),
        )
    # Fall back to X/Y + WCS
    h = hdu.header
    # Find 1-based column indices for X and Y
    xcol = ycol = None
    for i in range(1, h.get("TFIELDS", 0) + 1):
        name = h.get(f"TTYPE{i}", "").strip().upper()
        if name == "X":
            xcol = i
        elif name == "Y":
            ycol = i
    if xcol is None or ycol is None:
        raise RuntimeError("EVENTS HDU has neither RA/DEC nor X/Y columns")
    # Build a minimal astropy WCS from the column keywords
    w = WCS(naxis=2)
    w.wcs.crpix = [h[f"TCRPX{xcol}"], h[f"TCRPX{ycol}"]]
    w.wcs.cdelt = [h[f"TCDLT{xcol}"], h[f"TCDLT{ycol}"]]
    w.wcs.crval = [h[f"TCRVL{xcol}"], h[f"TCRVL{ycol}"]]
    w.wcs.ctype = [h.get(f"TCTYP{xcol}", "RA---TAN"), h.get(f"TCTYP{ycol}", "DEC--TAN")]
    x = np.asarray(hdu.data["X"], dtype=float)
    y = np.asarray(hdu.data["Y"], dtype=float)
    ra, dec = w.all_pix2world(x, y, 1)  # 1-based pixel origin
    return ra.astype(float), dec.astype(float)


def wrap_delta_ra_deg(ra_deg: np.ndarray, center_ra_deg: float) -> np.ndarray:
    return (
        (np.asarray(ra_deg, dtype=float) - float(center_ra_deg) + 180.0) % 360.0
    ) - 180.0


def offsets_arcsec(
    ra_deg: np.ndarray, dec_deg: np.ndarray, center_ra_deg: float, center_dec_deg: float
) -> tuple[np.ndarray, np.ndarray]:
    cos_dec = math.cos(math.radians(center_dec_deg))
    dx = wrap_delta_ra_deg(ra_deg, center_ra_deg) * cos_dec * 3600.0
    dy = (np.asarray(dec_deg, dtype=float) - float(center_dec_deg)) * 3600.0
    return dx.astype(float), dy.astype(float)


def choose_center(
    events: list[str],
    track_path: str | None,
    center_ra: float | None,
    center_dec: float | None,
    pimin: int,
    pimax: int,
) -> tuple[float, float]:
    if center_ra is not None and center_dec is not None:
        return float(center_ra), float(center_dec)
    if track_path:
        track = read_track(track_path)
        mid = len(track.ra) // 2
        return float(track.ra[mid]), float(track.dec[mid])
    for path in events:
        with fits.open(path, memmap=True) as hdul:
            if "EVENTS" not in hdul:
                continue
            evt = hdul["EVENTS"]
            data = evt.data
            if data is None or len(data) == 0:
                continue
            pi = np.asarray(data["PI"], dtype=float)
            m = np.isfinite(pi) & (pi >= pimin) & (pi <= pimax)
            if not np.any(m):
                continue
            ra, dec = _event_radec(evt)
            return float(np.nanmedian(ra[m])), float(np.nanmedian(dec[m]))
    raise RuntimeError("Could not determine mosaic center")


def selected_offsets(
    path: str, center_ra: float, center_dec: float, pimin: int, pimax: int
) -> tuple[np.ndarray, np.ndarray, int]:
    with fits.open(path, memmap=True) as hdul:
        if "EVENTS" not in hdul:
            raise RuntimeError(f"No EVENTS extension in {path}")
        evt = hdul["EVENTS"]
        data = evt.data
        if data is None or len(data) == 0:
            return np.asarray([], dtype=float), np.asarray([], dtype=float), 0
        pi = np.asarray(data["PI"], dtype=float)
        ra, dec = _event_radec(evt)
    m = (
        np.isfinite(pi)
        & np.isfinite(ra)
        & np.isfinite(dec)
        & (pi >= pimin)
        & (pi <= pimax)
    )
    if not np.any(m):
        return np.asarray([], dtype=float), np.asarray([], dtype=float), 0
    dx, dy = offsets_arcsec(ra[m], dec[m], center_ra, center_dec)
    return dx, dy, int(np.count_nonzero(m))


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser()
    ap.add_argument("--events", nargs="+", required=True)
    ap.add_argument("--pimin", type=int, required=True)
    ap.add_argument("--pimax", type=int, required=True)
    ap.add_argument("--out-fits", required=True)
    ap.add_argument("--out-json")
    ap.add_argument("--track")
    ap.add_argument("--srclist")
    ap.add_argument("--center-ra", type=float)
    ap.add_argument("--center-dec", type=float)
    ap.add_argument("--bin-arcsec", type=float, default=6.0)
    ap.add_argument("--pad-arcmin", type=float, default=2.0)
    ap.add_argument("--radius-arcmin", type=float)
    ap.add_argument("--label", default="diagnostic")
    return ap.parse_args()


def main() -> int:
    args = parse_args()
    if args.bin_arcsec <= 0:
        raise RuntimeError("--bin-arcsec must be > 0")
    center_ra, center_dec = choose_center(
        events=args.events,
        track_path=args.track,
        center_ra=args.center_ra,
        center_dec=args.center_dec,
        pimin=args.pimin,
        pimax=args.pimax,
    )

    pad_arcsec = float(args.pad_arcmin) * 60.0
    all_stats: list[dict[str, object]] = []
    xmins: list[float] = []
    xmaxs: list[float] = []
    ymins: list[float] = []
    ymaxs: list[float] = []
    nevents_total = 0

    for path in args.events:
        dx, dy, nsel = selected_offsets(
            path, center_ra, center_dec, args.pimin, args.pimax
        )
        nevents_total += nsel
        if dx.size:
            xmins.append(float(np.nanmin(dx)))
            xmaxs.append(float(np.nanmax(dx)))
            ymins.append(float(np.nanmin(dy)))
            ymaxs.append(float(np.nanmax(dy)))
        all_stats.append(
            {"path": str(Path(path).resolve()), "selected_events": int(nsel)}
        )

    if args.track:
        track = read_track(args.track)
        dx_t, dy_t = offsets_arcsec(track.ra, track.dec, center_ra, center_dec)
        if dx_t.size:
            xmins.append(float(np.nanmin(dx_t)))
            xmaxs.append(float(np.nanmax(dx_t)))
            ymins.append(float(np.nanmin(dy_t)))
            ymaxs.append(float(np.nanmax(dy_t)))
    if args.srclist:
        try:
            src = read_srclist(args.srclist)
            dx_s, dy_s = offsets_arcsec(src.ra, src.dec, center_ra, center_dec)
            if dx_s.size:
                xmins.append(float(np.nanmin(dx_s)))
                xmaxs.append(float(np.nanmax(dx_s)))
                ymins.append(float(np.nanmin(dy_s)))
                ymaxs.append(float(np.nanmax(dy_s)))
        except Exception:
            pass

    if args.radius_arcmin is not None:
        half = float(args.radius_arcmin) * 60.0
        xmin = -half
        xmax = half
        ymin = -half
        ymax = half
    elif xmins and xmaxs and ymins and ymaxs:
        xmin = float(min(xmins)) - pad_arcsec
        xmax = float(max(xmaxs)) + pad_arcsec
        ymin = float(min(ymins)) - pad_arcsec
        ymax = float(max(ymaxs)) + pad_arcsec
        # The lists interleave x/y entries. A square image is easier to compare between bands.
        half = max(abs(xmin), abs(xmax), abs(ymin), abs(ymax))
        xmin, xmax, ymin, ymax = -half, half, -half, half
    else:
        xmin = ymin = -300.0
        xmax = ymax = 300.0

    bin_arcsec = float(args.bin_arcsec)
    xmin = math.floor(xmin / bin_arcsec) * bin_arcsec
    xmax = math.ceil(xmax / bin_arcsec) * bin_arcsec
    ymin = math.floor(ymin / bin_arcsec) * bin_arcsec
    ymax = math.ceil(ymax / bin_arcsec) * bin_arcsec
    if xmax <= xmin:
        xmax = xmin + bin_arcsec
    if ymax <= ymin:
        ymax = ymin + bin_arcsec

    xedges = np.arange(xmin, xmax + 0.5 * bin_arcsec, bin_arcsec, dtype=float)
    yedges = np.arange(ymin, ymax + 0.5 * bin_arcsec, bin_arcsec, dtype=float)
    image = np.zeros((len(yedges) - 1, len(xedges) - 1), dtype=np.float32)

    for path in args.events:
        dx, dy, _nsel = selected_offsets(
            path, center_ra, center_dec, args.pimin, args.pimax
        )
        if dx.size == 0:
            continue
        hist, _ye, _xe = np.histogram2d(dy, dx, bins=[yedges, xedges])
        image += hist.astype(np.float32)

    hdr = fits.Header()
    hdr["BUNIT"] = "counts"
    hdr["LABEL"] = str(args.label)
    hdr["PIMIN"] = int(args.pimin)
    hdr["PIMAX"] = int(args.pimax)
    hdr["CTRRA"] = (float(center_ra), "Mosaic centre RA in deg")
    hdr["CTRDEC"] = (float(center_dec), "Mosaic centre Dec in deg")
    hdr["XMIN_AS"] = (float(xmin), "Min east offset from centre in arcsec")
    hdr["XMAX_AS"] = (float(xmax), "Max east offset from centre in arcsec")
    hdr["YMIN_AS"] = (float(ymin), "Min north offset from centre in arcsec")
    hdr["YMAX_AS"] = (float(ymax), "Max north offset from centre in arcsec")
    hdr["BIN_AS"] = (float(bin_arcsec), "Image bin size in arcsec")
    hdr["NINPUTS"] = int(len(args.events))
    hdr["NEVTSEL"] = int(nevents_total)
    hdr["HISTORY"] = "Diagnostic still-sky mosaic built from selected event lists"

    out_fits = Path(args.out_fits)
    out_fits.parent.mkdir(parents=True, exist_ok=True)
    fits.PrimaryHDU(data=image, header=hdr).writeto(out_fits, overwrite=True)

    if args.out_json:
        payload = {
            "label": str(args.label),
            "pimin": int(args.pimin),
            "pimax": int(args.pimax),
            "center_ra_deg": float(center_ra),
            "center_dec_deg": float(center_dec),
            "x_min_arcsec": float(xmin),
            "x_max_arcsec": float(xmax),
            "y_min_arcsec": float(ymin),
            "y_max_arcsec": float(ymax),
            "bin_arcsec": float(bin_arcsec),
            "selected_events_total": int(nevents_total),
            "inputs": all_stats,
        }
        out_json = Path(args.out_json)
        out_json.parent.mkdir(parents=True, exist_ok=True)
        out_json.write_text(
            json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )

    print(out_fits)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:  # pragma: no cover
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(1)
