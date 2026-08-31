from pathlib import Path

import numpy as np
from astropy.io import fits

from common import directory

# CAL-SRN-0399 B4--B6 overlap weights for 2.5--8.5 keV.
HARD_WEIGHTS = (("B4", 1499 / 1999), ("B5", 1.0), ("B6", 0.125))


def ccf_path():
    cif = directory("01_init") / "ccf.cif"
    with fits.open(cif, memmap=True) as hdul:
        rows = hdul["CALINDEX"].data
        selected = [row for row in rows if str(row["TYPEID"]).strip() == "FLARE"]
    if len(selected) != 1:
        raise RuntimeError(f"expected one FLARE calibration in {cif}")
    path = Path(str(selected[0]["FNAME"]).strip())
    if not path.is_file():
        raise FileNotFoundError(path)
    return path


def hard_template(detector, calibration=None):
    instrument = {"mos1": "EMOS1", "mos2": "EMOS2", "pn": "EPN"}[detector]
    calibration = Path(calibration or ccf_path())
    output = np.zeros((780, 780))
    with fits.open(calibration, memmap=True) as hdul:
        rows = hdul["FLAREMAP"].data
        index = {
            (
                str(row["INSTRUME"]).strip(),
                str(row["BAND"]).strip(),
                str(row["FILTER"]).strip(),
            ): position
            for position, row in enumerate(rows)
        }
        for band, weight in HARD_WEIGHTS:
            for filter_name in ("Thin1", "Medium", "Thick"):
                key = (instrument, band, filter_name)
                if key not in index:
                    raise RuntimeError(
                        f"missing FLARE map {instrument}/{filter_name}/{band}"
                    )
                current = np.asarray(rows[index[key]]["FLAREMAP"], float)
                if current.shape != (780, 780) or not np.all(np.isfinite(current)):
                    raise RuntimeError(
                        f"invalid FLARE map {instrument}/{filter_name}/{band}"
                    )
                output += weight * current
    return output


def replace_hard_map(path, detector, chips, support=None, calibration=None):
    calibration = Path(calibration or ccf_path())
    template = hard_template(detector, calibration)
    if "F" in chips.split():
        if support is None:
            raise RuntimeError(f"missing chip-mask support for {path}")
        with fits.open(support, memmap=True) as hdul:
            support_data = next(
                item.data
                for item in hdul
                if item.data is not None and item.data.ndim == 2
            )
            template = template * (np.asarray(support_data) != 0)
    total = float(np.sum(template, dtype=np.float64))
    with fits.open(path, mode="update", memmap=False) as hdul:
        hdu = next(
            item for item in hdul if item.data is not None and item.data.ndim == 2
        )
        target = float(np.nansum(hdu.data, dtype=np.float64))
        if target < 0 or not np.isfinite(target) or total <= 0:
            raise RuntimeError(f"invalid hard soft-proton normalization for {path}")
        corrected = np.asarray(
            np.zeros_like(template) if target == 0 else template * (target / total),
            dtype=hdu.data.dtype,
        )
        if not np.isclose(
            float(np.sum(corrected, dtype=np.float64)), target, rtol=2e-7, atol=1e-5
        ):
            raise RuntimeError(f"hard soft-proton count conservation failed for {path}")
        hdu.data[...] = corrected
        hdu.header["SPCORR"] = "CCF FLARE B4-B6"
        hdu.header["SPFILT"] = "ALL"
        hdu.header["SPB4W"] = HARD_WEIGHTS[0][1]
        hdu.header["SPB5W"] = HARD_WEIGHTS[1][1]
        hdu.header["SPB6W"] = HARD_WEIGHTS[2][1]
        hdu.header["SPCCF"] = calibration.name
