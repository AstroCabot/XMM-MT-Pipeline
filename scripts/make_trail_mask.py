#!/usr/bin/env python3
"""Build a common stationary-source trail mask on the comet-frame image grid."""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
from astropy.io import fits

from trail_mask_utils import build_trail_mask, read_srclist


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser()
    ap.add_argument("--grid-json", required=True)
    ap.add_argument("--track", required=True)
    ap.add_argument("--srclist", required=True)
    ap.add_argument("--mask-radius-arcsec", type=float, required=True)
    ap.add_argument("--out-mask", required=True)
    ap.add_argument("--out-json")
    ap.add_argument("--track-tmin-sec", type=float)
    ap.add_argument("--track-tmax-sec", type=float)
    return ap.parse_args()


def main() -> int:
    args = parse_args()
    grid = json.loads(Path(args.grid_json).read_text(encoding="utf-8"))
    nx = int(grid["nx"])
    ny = int(grid["ny"])
    if nx <= 0 or ny <= 0:
        raise RuntimeError("grid.json has non-positive nx/ny")

    src_rows = 0
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
            "track_tmin_sec": None if args.track_tmin_sec is None else float(args.track_tmin_sec),
            "track_tmax_sec": None if args.track_tmax_sec is None else float(args.track_tmax_sec),
        }
        out_json = Path(args.out_json)
        out_json.parent.mkdir(parents=True, exist_ok=True)
        out_json.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    print(out_path)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:  # pragma: no cover
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(1)
