#!/usr/bin/env python3
"""Sum corrected EPIC light curves onto a common time grid.

This creates a convenience combined light curve by summing the per-instrument
corrected net count rates and combining their uncertainties in quadrature.  It
is intended for low-resolution timing products after each input light curve has
already been corrected by `epiclccorr`.
"""
from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

import numpy as np
from astropy.io import fits


RATE_COLS = ["TIME", "RATE", "ERROR"]
OPTIONAL_COLS = ["BACKV", "BACKE", "FRACEXP"]


def find_rate_hdu(hdul: fits.HDUList) -> fits.BinTableHDU:
    for hdu in hdul[1:]:
        if isinstance(hdu, fits.BinTableHDU) and hdu.data is not None:
            names = {n.upper() for n in hdu.columns.names}
            if all(col in names for col in RATE_COLS):
                return hdu
    raise RuntimeError("No RATE-like binary table found")


class LightCurve:
    def __init__(self, path: str):
        self.path = path
        with fits.open(path) as hdul:
            self.primary = hdul[0].header.copy()
            hdu = find_rate_hdu(hdul)
            self.header = hdu.header.copy()
            tab = hdu.data
            self.time = np.asarray(tab["TIME"], dtype=float)
            self.rate = np.asarray(tab["RATE"], dtype=float)
            self.error = np.asarray(tab["ERROR"], dtype=float)
            self.backv = np.asarray(tab["BACKV"], dtype=float) if "BACKV" in hdu.columns.names else None
            self.backe = np.asarray(tab["BACKE"], dtype=float) if "BACKE" in hdu.columns.names else None
            self.fracexp = np.asarray(tab["FRACEXP"], dtype=float) if "FRACEXP" in hdu.columns.names else None
            self.extname = hdu.name if hdu.name else "RATE"
            self.timedel = float(hdu.header.get("TIMEDEL", self.primary.get("TIMEDEL", np.nan)))



def assert_same_grid(lcs: list[LightCurve], tol: float = 1e-6) -> None:
    base = lcs[0]
    for other in lcs[1:]:
        if len(other.time) != len(base.time):
            raise RuntimeError(f"Light curves have different lengths: {base.path} vs {other.path}")
        if not np.allclose(other.time, base.time, atol=tol, rtol=0.0):
            raise RuntimeError(f"Light curves are not aligned in TIME: {base.path} vs {other.path}")
        if np.isfinite(base.timedel) and np.isfinite(other.timedel):
            if abs(base.timedel - other.timedel) > tol:
                raise RuntimeError(f"Light curves have different TIMEDEL: {base.path} vs {other.path}")



def write_output(lcs: list[LightCurve], out_fits: str, out_csv: str | None) -> None:
    assert lcs
    assert_same_grid(lcs)
    base = lcs[0]

    rates = np.vstack([lc.rate for lc in lcs])
    errors = np.vstack([lc.error for lc in lcs])
    good = np.isfinite(rates) & np.isfinite(errors)

    total_rate = np.where(np.any(good, axis=0), np.nansum(np.where(good, rates, 0.0), axis=0), np.nan)
    total_error = np.where(np.any(good, axis=0), np.sqrt(np.nansum(np.where(good, errors, 0.0) ** 2, axis=0)), np.nan)
    ninstr = np.sum(good, axis=0).astype(np.int16)

    if any(lc.backv is not None for lc in lcs):
        backv_stack = np.vstack([np.asarray(lc.backv if lc.backv is not None else np.full_like(base.time, np.nan), dtype=float) for lc in lcs])
        total_backv = np.where(np.any(np.isfinite(backv_stack), axis=0), np.nansum(backv_stack, axis=0), np.nan)
    else:
        total_backv = None

    if any(lc.backe is not None for lc in lcs):
        backe_stack = np.vstack([np.asarray(lc.backe if lc.backe is not None else np.full_like(base.time, np.nan), dtype=float) for lc in lcs])
        total_backe = np.where(np.any(np.isfinite(backe_stack), axis=0), np.sqrt(np.nansum(backe_stack ** 2, axis=0)), np.nan)
    else:
        total_backe = None

    if any(lc.fracexp is not None for lc in lcs):
        frac_stack = np.vstack([np.asarray(lc.fracexp if lc.fracexp is not None else np.full_like(base.time, np.nan), dtype=float) for lc in lcs])
        total_frac = np.where(np.any(np.isfinite(frac_stack), axis=0), np.nanmean(frac_stack, axis=0), np.nan)
    else:
        total_frac = None

    cols = [
        fits.Column(name="TIME", format="D", unit="s", array=base.time.astype(np.float64)),
        fits.Column(name="RATE", format="E", unit="count / s", array=total_rate.astype(np.float32)),
        fits.Column(name="ERROR", format="E", unit="count / s", array=total_error.astype(np.float32)),
    ]
    if total_backv is not None:
        cols.append(fits.Column(name="BACKV", format="E", unit="count / s", array=total_backv.astype(np.float32)))
    if total_backe is not None:
        cols.append(fits.Column(name="BACKE", format="E", unit="count / s", array=total_backe.astype(np.float32)))
    if total_frac is not None:
        cols.append(fits.Column(name="FRACEXP", format="E", array=total_frac.astype(np.float32)))
    cols.append(fits.Column(name="NINST", format="I", array=ninstr))

    rate_hdu = fits.BinTableHDU.from_columns(cols, name=base.extname)
    for key, val in base.header.items():
        if key not in {"XTENSION", "BITPIX", "NAXIS", "NAXIS1", "NAXIS2", "PCOUNT", "GCOUNT", "TFIELDS"}:
            rate_hdu.header[key] = val
    rate_hdu.header["HISTORY"] = "RATE is the sum of epiclccorr-corrected instrument rates"
    rate_hdu.header["HISTORY"] = f"Combined from {len(lcs)} input light curves"
    rate_hdu.header["NINSTR"] = (len(lcs), "Number of contributing light curves")

    primary = fits.PrimaryHDU(header=base.primary)
    hdul = fits.HDUList([primary, rate_hdu])
    Path(out_fits).parent.mkdir(parents=True, exist_ok=True)
    hdul.writeto(out_fits, overwrite=True)

    if out_csv:
        Path(out_csv).parent.mkdir(parents=True, exist_ok=True)
        with open(out_csv, "w", newline="", encoding="utf-8") as fh:
            writer = csv.writer(fh)
            header = ["TIME", "RATE", "ERROR"]
            if total_backv is not None:
                header.append("BACKV")
            if total_backe is not None:
                header.append("BACKE")
            if total_frac is not None:
                header.append("FRACEXP")
            header.append("NINST")
            writer.writerow(header)
            for i in range(len(base.time)):
                row = [f"{base.time[i]:.6f}", f"{total_rate[i]:.8g}", f"{total_error[i]:.8g}"]
                if total_backv is not None:
                    row.append(f"{total_backv[i]:.8g}")
                if total_backe is not None:
                    row.append(f"{total_backe[i]:.8g}")
                if total_frac is not None:
                    row.append(f"{total_frac[i]:.8g}")
                row.append(int(ninstr[i]))
                writer.writerow(row)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--inputs", nargs="+", required=True, help="Corrected light curves from epiclccorr")
    ap.add_argument("--out-fits", required=True)
    ap.add_argument("--out-csv")
    args = ap.parse_args()

    lcs = [LightCurve(path) for path in args.inputs]
    write_output(lcs, args.out_fits, args.out_csv)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:  # pragma: no cover
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(1)
