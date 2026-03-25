#!/usr/bin/env python3
"""Utilities for comet-frame XMM event geometry.

Reads the X/Y table WCS from an XMM EVENTS extension and converts comet-centric
sky offsets in arcsec into event-list X/Y coordinates (physical sky pixels).

Supported subcommands:
  xy          : print centre and pixel scale information
  region      : emit evselect/especget region expressions in X/Y
  image-grid  : emit a common image grid for evselect/eexpmap, optionally to JSON
"""
from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path
from typing import Any

from astropy.io import fits


def _find_col_index(header: fits.Header, name: str) -> int:
    tfields = int(header.get("TFIELDS", 0))
    for idx in range(1, tfields + 1):
        if str(header.get(f"TTYPE{idx}", "")).strip().upper() == name.upper():
            return idx
    raise RuntimeError(f"Could not find column {name!r} in EVENTS header")


def read_xy_wcs(event_file: str) -> dict[str, float]:
    with fits.open(event_file) as hdul:
        if "EVENTS" not in hdul:
            raise RuntimeError(f"No EVENTS extension in {event_file}")
        hdr = hdul["EVENTS"].header
        x_idx = _find_col_index(hdr, "X")
        y_idx = _find_col_index(hdr, "Y")

        def get_num(prefix: str, idx: int, default: float | None = None) -> float:
            key = f"{prefix}{idx}"
            if key in hdr:
                return float(hdr[key])
            if default is not None:
                return float(default)
            raise RuntimeError(f"Missing required WCS keyword {key} in EVENTS header")

        xref = get_num("TCRPX", x_idx)
        yref = get_num("TCRPX", y_idx)
        xcdelt = get_num("TCDLT", x_idx)
        ycdelt = get_num("TCDLT", y_idx)
        xval = get_num("TCRVL", x_idx, 0.0)
        yval = get_num("TCRVL", y_idx, 0.0)
        xmin = float(hdr.get(f"TLMIN{x_idx}", math.floor(xref - 30000.0)))
        xmax = float(hdr.get(f"TLMAX{x_idx}", math.ceil(xref + 30000.0)))
        ymin = float(hdr.get(f"TLMIN{y_idx}", math.floor(yref - 30000.0)))
        ymax = float(hdr.get(f"TLMAX{y_idx}", math.ceil(yref + 30000.0)))

    x_arcsec_per_phys = xcdelt * 3600.0
    y_arcsec_per_phys = ycdelt * 3600.0
    if abs(abs(x_arcsec_per_phys) - abs(y_arcsec_per_phys)) > 1e-6:
        # EPIC sky coordinates should be square; tolerate small floating noise only.
        raise RuntimeError(
            "X and Y event scales differ unexpectedly: "
            f"{x_arcsec_per_phys} vs {y_arcsec_per_phys} arcsec/phys"
        )

    return {
        "x_ref_phys": xref,
        "y_ref_phys": yref,
        "x_ref_world": xval,
        "y_ref_world": yval,
        "x_arcsec_per_phys": x_arcsec_per_phys,
        "y_arcsec_per_phys": y_arcsec_per_phys,
        "scale_arcsec_per_phys": 0.5
        * (abs(x_arcsec_per_phys) + abs(y_arcsec_per_phys)),
        "x_tlmin": xmin,
        "x_tlmax": xmax,
        "y_tlmin": ymin,
        "y_tlmax": ymax,
    }


def offset_to_xy(
    wcs: dict[str, float], dx_east_arcsec: float, dy_north_arcsec: float
) -> tuple[float, float]:
    x = wcs["x_ref_phys"] + dx_east_arcsec / wcs["x_arcsec_per_phys"]
    y = wcs["y_ref_phys"] + dy_north_arcsec / wcs["y_arcsec_per_phys"]
    return float(x), float(y)


def radius_to_phys(wcs: dict[str, float], radius_arcsec: float) -> float:
    return float(radius_arcsec / wcs["scale_arcsec_per_phys"])


def build_region_expr(
    wcs: dict[str, float],
    shape: str,
    dx: float,
    dy: float,
    r: float | None = None,
    rin: float | None = None,
    rout: float | None = None,
) -> tuple[str, float, float]:
    x, y = offset_to_xy(wcs, dx, dy)
    if shape == "circle":
        if r is None:
            raise RuntimeError("circle shape requires r")
        rp = radius_to_phys(wcs, r)
        expr = f"((X,Y) IN circle({x:.6f},{y:.6f},{rp:.6f}))"
    elif shape == "annulus":
        if rin is None or rout is None:
            raise RuntimeError("annulus shape requires rin and rout")
        rinp = radius_to_phys(wcs, rin)
        routp = radius_to_phys(wcs, rout)
        expr = f"((X,Y) IN annulus({x:.6f},{y:.6f},{rinp:.6f},{routp:.6f}))"
    else:
        raise RuntimeError(f"Unsupported shape: {shape}")
    return expr, x, y


def image_grid(
    wcs: dict[str, float],
    radius_arcsec: float,
    bin_arcsec: float,
    center_dx: float = 0.0,
    center_dy: float = 0.0,
) -> dict[str, Any]:
    if radius_arcsec <= 0:
        raise RuntimeError("radius_arcsec must be > 0")
    if bin_arcsec <= 0:
        raise RuntimeError("bin_arcsec must be > 0")

    cx, cy = offset_to_xy(wcs, center_dx, center_dy)
    scale = wcs["scale_arcsec_per_phys"]
    bin_phys = round(bin_arcsec / scale)
    if bin_phys < 1:
        raise RuntimeError(
            f"Requested image bin ({bin_arcsec} arcsec) is smaller than one X/Y physical pixel ({scale} arcsec)"
        )

    half_width_phys = radius_arcsec / scale
    xmin = math.floor(cx - half_width_phys)
    xmax = math.ceil(cx + half_width_phys)
    ymin = math.floor(cy - half_width_phys)
    ymax = math.ceil(cy + half_width_phys)

    # Make the range align with the requested image binning.
    nx = max(1, int(math.ceil((xmax - xmin) / bin_phys)))
    ny = max(1, int(math.ceil((ymax - ymin) / bin_phys)))
    xmax = xmin + nx * bin_phys
    ymax = ymin + ny * bin_phys

    return {
        "center_dx_arcsec": float(center_dx),
        "center_dy_arcsec": float(center_dy),
        "center_x_phys": float(cx),
        "center_y_phys": float(cy),
        "x_min_phys": float(xmin),
        "x_max_phys": float(xmax),
        "y_min_phys": float(ymin),
        "y_max_phys": float(ymax),
        "bin_arcsec": float(bin_arcsec),
        "bin_phys": float(bin_phys),
        "nx": int(nx),
        "ny": int(ny),
        "radius_arcsec": float(radius_arcsec),
        "scale_arcsec_per_phys": float(scale),
        "x_arcsec_per_phys": float(wcs["x_arcsec_per_phys"]),
        "y_arcsec_per_phys": float(wcs["y_arcsec_per_phys"]),
        "x_ref_phys": float(wcs["x_ref_phys"]),
        "y_ref_phys": float(wcs["y_ref_phys"]),
    }


def cmd_xy(args: argparse.Namespace) -> int:
    wcs = read_xy_wcs(args.event)
    for key in [
        "x_ref_phys",
        "y_ref_phys",
        "x_arcsec_per_phys",
        "y_arcsec_per_phys",
        "scale_arcsec_per_phys",
        "x_tlmin",
        "x_tlmax",
        "y_tlmin",
        "y_tlmax",
    ]:
        print(f"{key.upper()}={wcs[key]}")
    return 0


def cmd_region(args: argparse.Namespace) -> int:
    wcs = read_xy_wcs(args.event)
    expr, x, y = build_region_expr(
        wcs=wcs,
        shape=args.shape,
        dx=float(args.dx),
        dy=float(args.dy),
        r=args.r,
        rin=args.rin,
        rout=args.rout,
    )
    print(f"REGION_EXPR={json.dumps(expr)}")
    print(f"CENTER_X={x:.6f}")
    print(f"CENTER_Y={y:.6f}")
    return 0


def cmd_image_grid(args: argparse.Namespace) -> int:
    wcs = read_xy_wcs(args.event)
    grid = image_grid(
        wcs=wcs,
        radius_arcsec=float(args.radius_arcsec),
        bin_arcsec=float(args.bin_arcsec),
        center_dx=float(args.center_dx),
        center_dy=float(args.center_dy),
    )
    if args.out_json:
        out = Path(args.out_json)
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(
            json.dumps(grid, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
    for key, value in grid.items():
        if isinstance(value, str):
            print(f"{key.upper()}={json.dumps(value)}")
        else:
            print(f"{key.upper()}={value}")
    return 0


def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)

    ap_xy = sub.add_parser("xy")
    ap_xy.add_argument("--event", required=True)
    ap_xy.set_defaults(func=cmd_xy)

    ap_reg = sub.add_parser("region")
    ap_reg.add_argument("--event", required=True)
    ap_reg.add_argument("--shape", choices=["circle", "annulus"], required=True)
    ap_reg.add_argument(
        "--dx", type=float, default=0.0, help="Offset east of comet center in arcsec"
    )
    ap_reg.add_argument(
        "--dy", type=float, default=0.0, help="Offset north of comet center in arcsec"
    )
    ap_reg.add_argument("--r", type=float, help="Circle radius in arcsec")
    ap_reg.add_argument("--rin", type=float, help="Annulus inner radius in arcsec")
    ap_reg.add_argument("--rout", type=float, help="Annulus outer radius in arcsec")
    ap_reg.set_defaults(func=cmd_region)

    ap_img = sub.add_parser("image-grid")
    ap_img.add_argument("--event", required=True)
    ap_img.add_argument("--radius-arcsec", type=float, required=True)
    ap_img.add_argument("--bin-arcsec", type=float, required=True)
    ap_img.add_argument("--center-dx", type=float, default=0.0)
    ap_img.add_argument("--center-dy", type=float, default=0.0)
    ap_img.add_argument("--out-json")
    ap_img.set_defaults(func=cmd_image_grid)

    return ap


def main() -> int:
    args = build_parser().parse_args()
    return int(args.func(args))


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:  # pragma: no cover
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(1)
