import sys

import numpy as np

from _common import DATA, REPRO, read_json, stage_dir

sys.path.insert(0, str(REPRO))

from analysis_02_pps_luminosity import run


def main():
    run.OUTPUT = stage_dir(2)
    archive = read_json(DATA / "pps" / "profile_inputs.json")
    samples = {}
    with np.load(DATA / "pps" / "profiles.npz") as packed:
        for band in archive["bands"]:
            low, high = map(int, band.split("-"))
            prefix = f"{low}_{high}_"
            samples[low, high] = {
                key[len(prefix) :]: packed[key]
                for key in packed.files
                if key.startswith(prefix)
            }
    ecfs = {
        tuple(map(int, band.split("-"))): value
        for band, value in archive["ecfs"].items()
    }
    result, profiles = run.calculate(samples, ecfs)
    for name, band in result["bands"].items():
        band["ecf_source"] = archive["ecf_sources"][name]
    result.update(archive["metadata"])
    run.write_outputs(result, profiles)


if __name__ == "__main__":
    main()
