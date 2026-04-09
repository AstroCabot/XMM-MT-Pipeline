#!/usr/bin/env python3
"""Filter background field sources from XMM event lists in the sky frame.

For each event, checks whether the event's sky-frame (X, Y) position falls
within a circular exclusion zone around any detected field source.  Events
inside any exclusion zone are removed; all others are kept.

This should be run on the **sky-frame** (pre-attcalc) clean merged event
files, where background sources are stationary points.  The surviving events
are then suitable for comet-frame reprocessing via attcalc, producing
source-clean comet-frame images with minimal data loss (~few percent removed
instead of the ~72% lost with full-trail masking).

Usage:
    filter_field_sources.py \\
        --event  clean/PN_clean_merged.fits \\
        --srclist detect/field_sources_all.fits \\
        --mask-radius-arcsec 20.0 \\
        --output clean/PN_srcfilt.fits \\
        [--report-json clean/PN_srcfilt.json]
"""
from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

import numpy as np
from astropy.io import fits


# ── WCS helpers ──────────────────────────────────────────────────────────────

def _find_col_index(header: fits.Header, name: str) -> int:
    """Find the 1-based column index for *name* in a FITS table header."""
    tfields = int(header.get("TFIELDS", 0))
    for idx in range(1, tfields + 1):
        if str(header.get(f"TTYPE{idx}", "")).strip().upper() == name.upper():
            return idx
    raise RuntimeError(f"Could not find column {name!r} in EVENTS header")


def _read_evt_wcs(header: fits.Header) -> dict[str, float]:
    """Extract XMM-style event-list WCS from the EVENTS header."""
    x_idx = _find_col_index(header, "X")
    y_idx = _find_col_index(header, "Y")

    def _g(prefix: str, idx: int, default: float | None = None) -> float:
        key = f"{prefix}{idx}"
        if key in header:
            return float(header[key])
        if default is not None:
            return default
        raise RuntimeError(f"Missing WCS keyword {key}")

    return {
        "x_ref_phys": _g("TCRPX", x_idx),
        "y_ref_phys": _g("TCRPX", y_idx),
        "x_ref_ra":   _g("TCRVL", x_idx, 0.0),
        "y_ref_dec":  _g("TCRVL", y_idx, 0.0),
        "x_cdelt":    _g("TCDLT", x_idx),   # deg / phys
        "y_cdelt":    _g("TCDLT", y_idx),   # deg / phys
    }


def _radec_to_xy(
    wcs: dict[str, float],
    ra_deg: np.ndarray,
    dec_deg: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Convert RA, Dec (degrees) to XMM physical X, Y using a tangent-plane
    projection centred on the event-list WCS reference pixel.

    This is the inverse of the standard SAS sky-coordinate mapping:
        RA  = x_ref_ra  + (X - x_ref_phys) * x_cdelt / cos(dec_ref)
        Dec = y_ref_dec  + (Y - y_ref_phys) * y_cdelt
    (with x_cdelt negative for RA increasing to the east).
    """
    ra0  = wcs["x_ref_ra"]
    dec0 = wcs["y_ref_dec"]
    cos_dec0 = math.cos(math.radians(dec0))

    # Arc-second offsets on the tangent plane
    dra  = ((np.asarray(ra_deg, dtype=float) - ra0 + 180.0) % 360.0 - 180.0)
    ddec = np.asarray(dec_deg, dtype=float) - dec0

    # Convert degrees to physical pixels
    x = wcs["x_ref_phys"] + (dra * cos_dec0) / wcs["x_cdelt"]
    y = wcs["y_ref_phys"] + ddec / wcs["y_cdelt"]
    return x, y


# ── Source list reader ───────────────────────────────────────────────────────

def _read_srclist(path: str) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Read RA, DEC, SRCID from a source-list FITS file."""
    with fits.open(path) as hdul:
        for hdu in hdul[1:]:
            if not isinstance(hdu, fits.BinTableHDU):
                continue
            names_upper = [n.upper() for n in hdu.columns.names]
            if "RA" in names_upper and "DEC" in names_upper:
                data = hdu.data
                if data is None or len(data) == 0:
                    return (
                        np.empty(0, dtype=float),
                        np.empty(0, dtype=float),
                        np.empty(0, dtype=int),
                    )
                ra_col  = hdu.columns.names[names_upper.index("RA")]
                dec_col = hdu.columns.names[names_upper.index("DEC")]
                ra  = np.asarray(data[ra_col], dtype=float)
                dec = np.asarray(data[dec_col], dtype=float)
                if "SRCID" in names_upper:
                    sid = np.asarray(
                        data[hdu.columns.names[names_upper.index("SRCID")]],
                        dtype=int,
                    )
                else:
                    sid = np.arange(1, len(ra) + 1, dtype=int)
                return ra, dec, sid
    raise RuntimeError(f"No table with RA/DEC columns found in {path}")


# ── Core filtering ───────────────────────────────────────────────────────────

def filter_events(
    event_path: str,
    srclist_path: str,
    mask_radius_arcsec: float,
    output_path: str,
    report_json: str | None = None,
) -> dict:
    """Remove events within *mask_radius_arcsec* of any field source.

    Returns a summary dict with filtering statistics.
    """
    src_ra, src_dec, src_id = _read_srclist(srclist_path)
    n_sources = len(src_ra)

    with fits.open(event_path) as hdul:
        if "EVENTS" not in hdul:
            raise RuntimeError(f"No EVENTS extension in {event_path}")
        evt_hdu = hdul["EVENTS"]
        data = evt_hdu.data
        header = evt_hdu.header
        if data is None or len(data) == 0:
            raise RuntimeError(f"EVENTS table is empty in {event_path}")

        n_input = len(data)
        x = np.asarray(data["X"], dtype=float)
        y = np.asarray(data["Y"], dtype=float)

        # --- Convert source RA/Dec to physical X, Y in the event frame ---
        wcs = _read_evt_wcs(header)
        src_x, src_y = _radec_to_xy(wcs, src_ra, src_dec)

        # --- Mask radius in physical units ---
        scale = 0.5 * (abs(wcs["x_cdelt"]) + abs(wcs["y_cdelt"])) * 3600.0
        r_phys = mask_radius_arcsec / scale

        # --- Build per-event exclusion mask ---
        # For ~77 sources × ~10M events, a vectorised loop over sources
        # is fast enough (~few seconds) with numpy broadcasting chunks.
        excluded = np.zeros(n_input, dtype=bool)
        r_phys_sq = r_phys * r_phys
        per_source_counts = {}

        # Process in source-loop (77 iterations, each O(N) with numpy)
        for i in range(n_sources):
            dx = x - src_x[i]
            dy = y - src_y[i]
            inside = (dx * dx + dy * dy) <= r_phys_sq
            n_hit = int(np.count_nonzero(inside & ~excluded))
            excluded |= inside
            per_source_counts[int(src_id[i])] = n_hit

        n_removed = int(np.count_nonzero(excluded))
        n_kept = n_input - n_removed
        keep = ~excluded

        # --- Build output FITS ---
        out_hdus = []
        for hdu in hdul:
            if hdu.name == "EVENTS":
                new_hdu = fits.BinTableHDU(
                    data=data[keep],
                    header=hdu.header.copy(),
                    name="EVENTS",
                )
                new_hdu.header["HISTORY"] = (
                    f"Filtered {n_removed} events within {mask_radius_arcsec} "
                    f"arcsec of {n_sources} field sources"
                )
                new_hdu.header["SRCFILT"] = (
                    True,
                    "Sky-frame field-source filter applied",
                )
                new_hdu.header["NSRCFLT"] = (
                    n_sources,
                    "Number of field sources filtered",
                )
                new_hdu.header["NSRCREM"] = (
                    n_removed,
                    "Events removed by source filter",
                )
                new_hdu.header["SRCRAD"] = (
                    mask_radius_arcsec,
                    "[arcsec] Source exclusion radius",
                )
                out_hdus.append(new_hdu)
            else:
                out_hdus.append(hdu.copy())

    # --- Write output ---
    out_path = Path(output_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fits.HDUList(out_hdus).writeto(out_path, overwrite=True)

    summary = {
        "event_input": str(Path(event_path).resolve()),
        "srclist": str(Path(srclist_path).resolve()),
        "output": str(out_path.resolve()),
        "n_sources": n_sources,
        "mask_radius_arcsec": mask_radius_arcsec,
        "r_phys": float(r_phys),
        "n_input_events": n_input,
        "n_removed_events": n_removed,
        "n_kept_events": n_kept,
        "removed_fraction": float(n_removed / n_input) if n_input > 0 else 0.0,
        "per_source_removed": per_source_counts,
    }

    if report_json:
        rpath = Path(report_json)
        rpath.parent.mkdir(parents=True, exist_ok=True)
        rpath.write_text(
            json.dumps(summary, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )

    return summary


# ── CLI ──────────────────────────────────────────────────────────────────────

def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(
        description="Remove sky-frame events contaminated by field sources."
    )
    ap.add_argument(
        "--event", required=True,
        help="Input event FITS file (sky-frame X/Y, e.g. clean merged).",
    )
    ap.add_argument(
        "--srclist", required=True,
        help="FITS source list with RA and DEC columns.",
    )
    ap.add_argument(
        "--mask-radius-arcsec", type=float, required=True,
        help="Circular exclusion radius around each source (arcsec).",
    )
    ap.add_argument(
        "--output", required=True,
        help="Output filtered event FITS file.",
    )
    ap.add_argument(
        "--report-json",
        help="Optional JSON report with filtering statistics.",
    )
    return ap.parse_args()


def main() -> int:
    args = parse_args()
    summary = filter_events(
        event_path=args.event,
        srclist_path=args.srclist,
        mask_radius_arcsec=args.mask_radius_arcsec,
        output_path=args.output,
        report_json=args.report_json,
    )
    n_in = summary["n_input_events"]
    n_rm = summary["n_removed_events"]
    frac = summary["removed_fraction"]
    n_src = summary["n_sources"]
    print(
        f"Filtered {n_rm}/{n_in} events ({frac:.1%}) "
        f"within {args.mask_radius_arcsec} arcsec of {n_src} sources"
    )
    print(args.output)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(1)
