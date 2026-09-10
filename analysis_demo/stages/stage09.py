import sys

from _common import OUT, REPRO, read_json, stage_dir, write_json

sys.path.insert(0, str(REPRO / "analysis_09_recoil"))

import physics
from run import load_draws


def main():
    output = stage_dir(9)
    draws, names = load_draws(OUT / "08" / "chain.npz")
    task08 = read_json(OUT / "08" / "result.json")
    parameters = read_json(REPRO / "analysis_09_recoil" / "parameters.json")
    result = physics.calculate(draws, names, task08, parameters)
    write_json(output / "result.json", result)
    print(f"stage 09: wrote {', '.join(result)}")


if __name__ == "__main__":
    main()
