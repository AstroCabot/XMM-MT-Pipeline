import argparse
import json
import os
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
RUNTIME = ROOT / "work" / "runtime" / "analysis_06_solar_wind"
for name in ("tmp", "cache", "config", "matplotlib"):
    (RUNTIME / name).mkdir(parents=True, exist_ok=True)
os.environ.update(
    TMPDIR=str(RUNTIME / "tmp"),
    XDG_CACHE_HOME=str(RUNTIME / "cache"),
    XDG_CONFIG_HOME=str(RUNTIME / "config"),
    MPLCONFIGDIR=str(RUNTIME / "matplotlib"),
)

import numpy as np

import solar_wind

SHARED = json.loads((ROOT / "parameters.json").read_text())
FINAL_PRODUCTS = (
    "mapped_flux.tsv",
    "result.json",
    "lightcurve.tsv",
    "figure4.png",
    "figure4.pdf",
)
ION_FRACTION = SHARED["ion_fraction_context"]


def summarize(wind, rows, band_ev, pattern):
    rates = np.array([float(row["rate_ct_s"]) for row in rows])
    if np.any(~np.isfinite(rates)) or rates.min() <= 0:
        raise RuntimeError("light-curve rates must be positive")
    return {
        "observation_utc": list(solar_wind.OBSERVATION_UTC),
        "proton_flux": {"units": "cm^-2 s^-1", **wind["stats"]},
        "ace_quality": wind["ace_quality"],
        "f107": wind["f107"],
        "ion_fraction": {"adopted_range": ION_FRACTION["adopted_range"]},
        "heavy_ion_prior": SHARED["heavy_ion_prior"],
        "lightcurve": {
            "band_eV": list(band_ev),
            "pattern": pattern,
            "frames": len(rows),
            "min_rate_ct_s": float(rates.min()),
            "max_rate_ct_s": float(rates.max()),
            "dynamic_range": float(rates.max() / rates.min()),
        },
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--wind-only",
        action="store_true",
        help="map the monitors and cache F10.7 without Task 01",
    )
    args = parser.parse_args()
    for name in FINAL_PRODUCTS:
        (solar_wind.OUTPUT / name).unlink(missing_ok=True)
    wind = solar_wind.products()
    if args.wind_only:
        return

    import lightcurve

    rows = lightcurve.build(lightcurve.pn_frames())
    result = summarize(wind, rows, lightcurve.BAND_EV, lightcurve.PATTERN)
    solar_wind.OUTPUT.mkdir(parents=True, exist_ok=True)
    (solar_wind.OUTPUT / "result.json").write_text(json.dumps(result, indent=2) + "\n")
    from figure import draw

    draw(wind, rows, solar_wind.OUTPUT)


if __name__ == "__main__":
    main()
