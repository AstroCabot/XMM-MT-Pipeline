import sys

from _common import DATA, REPRO, restore_env, stage_dir, write_json

sys.path.insert(0, str(REPRO / "analysis_06_solar_wind"))

import run as run6
import solar_wind

restore_env()

BAND_EV = (200, 1000)
PATTERN = 0


def main():
    output = stage_dir(6)
    solar_wind.CACHE = DATA / "wind" / "cache"
    solar_wind.OUTPUT = output
    wind = solar_wind.products(output=output)
    rows = solar_wind.read_tsv(DATA / "wind" / "lightcurve.tsv")
    result = run6.summarize(wind, rows, BAND_EV, PATTERN)
    write_json(output / "result.json", result)
    from figure import draw

    draw(wind, rows, output)
    stats = result["proton_flux"]
    print(
        "stage 06: mapped proton flux medians "
        + ", ".join(
            f"{name} {stats[name]['median']:.2e}"
            for name in ("stereo_a", "solar_orbiter", "wind")
            if name in stats
        )
    )


if __name__ == "__main__":
    main()
