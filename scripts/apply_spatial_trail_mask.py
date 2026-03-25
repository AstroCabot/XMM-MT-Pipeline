#!/usr/bin/env python3
"""Remove events that fall inside a stationary-source trail mask image."""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
from astropy.io import fits


def read_mask(path: str) -> np.ndarray:
    with fits.open(path) as hdul:
        for hdu in hdul:
            if hdu.data is not None and getattr(hdu.data, "ndim", 0) == 2:
                return np.asarray(hdu.data, dtype=bool)
    raise RuntimeError(f"No 2D mask image found in {path}")


def filter_events(event_path: str, mask: np.ndarray, grid: dict[str, float]) -> tuple[fits.HDUList, int, int, int]:
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

        out_hdus: list[fits.hdu.base.ExtensionHDU | fits.PrimaryHDU] = []
        for hdu in hdul:
            if hdu.name == "EVENTS":
                new_hdu = fits.BinTableHDU(data=data[keep], header=hdu.header.copy(), name=hdu.name)
                new_hdu.header["HISTORY"] = "Stationary-source trail-mask events removed"
                new_hdu.header["TRAILMSK"] = (True, "Stationary-source trail mask applied")
                new_hdu.header["TRMSKREM"] = (int(masked.sum()), "Events removed by trail mask")
                out_hdus.append(new_hdu)
            else:
                out_hdus.append(hdu.copy())

    return fits.HDUList(out_hdus), int(len(data)), int(masked.sum()), int(np.count_nonzero(in_bounds))


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser()
    ap.add_argument("--event", required=True)
    ap.add_argument("--mask", required=True)
    ap.add_argument("--grid-json", required=True)
    ap.add_argument("--output", required=True)
    ap.add_argument("--report-json")
    return ap.parse_args()


def main() -> int:
    args = parse_args()
    grid = json.loads(Path(args.grid_json).read_text(encoding="utf-8"))
    mask = read_mask(args.mask)
    hdul, n_in, n_removed, n_in_bounds = filter_events(args.event, mask, grid)
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
        report_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    print(out_path)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:  # pragma: no cover
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(1)
