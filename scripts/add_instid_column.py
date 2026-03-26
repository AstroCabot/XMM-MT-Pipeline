#!/usr/bin/env python3
"""Add an INSTID column to a FITS event list.

Used to tag events by instrument before merging PN + MOS event lists.
"""
from __future__ import annotations
import argparse
import sys
from pathlib import Path
import numpy as np
from astropy.io import fits

INST_MAP = {"PN": 1, "M1": 2, "MOS1": 2, "M2": 3, "MOS2": 3}


def _copy_header_cards(src: fits.Header, dst: fits.Header) -> None:
    skip = {
        "XTENSION",
        "BITPIX",
        "NAXIS",
        "NAXIS1",
        "NAXIS2",
        "PCOUNT",
        "GCOUNT",
        "TFIELDS",
    }
    for key, value in src.items():
        if key not in skip:
            dst[key] = value


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True)
    ap.add_argument("--output", required=True)
    ap.add_argument("--inst", required=True, help="PN, M1, or M2")
    args = ap.parse_args()
    inst_key = args.inst.strip().upper()
    if inst_key not in INST_MAP:
        raise RuntimeError(f"Unsupported instrument code: {args.inst}")
    instid = INST_MAP[inst_key]
    with fits.open(args.input) as hdul:
        out_hdus = []
        for hdu in hdul:
            if isinstance(hdu, fits.BinTableHDU) and hdu.name == "EVENTS":
                hdr = hdu.header.copy()
                names = list(hdu.columns.names)
                nrows = len(hdu.data) if hdu.data is not None else 0
                arr = np.full(nrows, instid, dtype=np.int16)
                if "INSTID" in names:
                    cols = []
                    for col in hdu.columns:
                        if col.name == "INSTID":
                            cols.append(
                                fits.Column(name="INSTID", format="I", array=arr)
                            )
                        else:
                            cols.append(col)
                    new_evt = fits.BinTableHDU.from_columns(cols, name="EVENTS")
                else:
                    new_evt = fits.BinTableHDU.from_columns(
                        hdu.columns
                        + fits.ColDefs(
                            [fits.Column(name="INSTID", format="I", array=arr)]
                        ),
                        name="EVENTS",
                    )
                _copy_header_cards(hdr, new_evt.header)
                new_evt.header["INSTID"] = (instid, "Constant per-row instrument code")
                new_evt.header["HISTORY"] = f"Added INSTID={instid} ({inst_key})"
                out_hdus.append(new_evt)
            else:
                out_hdus.append(hdu.copy())
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    fits.HDUList(out_hdus).writeto(args.output, overwrite=True)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(1)
