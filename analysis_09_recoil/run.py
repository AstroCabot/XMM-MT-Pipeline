import argparse
import json
from pathlib import Path

import numpy as np
from physics import calculate

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent


def load_draws(path):
    with np.load(path, allow_pickle=False) as data:
        draws = np.asarray(data["physical_draws"], float)
        names = [str(name) for name in data["physical_names"]]
    return draws, names


def run(task08_path, chain_path, parameters_path, output_path):
    output_resolved = output_path.resolve()
    try:
        output_resolved.relative_to(ROOT.resolve())
    except ValueError as error:
        raise ValueError("output must remain inside analysis_repro") from error
    inputs = tuple(
        path.resolve() for path in (task08_path, chain_path, parameters_path)
    )
    if any(
        output_resolved.exists() and path.exists() and output_resolved.samefile(path)
        for path in inputs
    ):
        raise ValueError("output must differ from every input")
    output_path.unlink(missing_ok=True)
    task08 = json.loads(task08_path.read_text())
    parameters = json.loads(parameters_path.read_text())
    draws, names = load_draws(chain_path)
    result = calculate(draws, names, task08, parameters)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(result, indent=2) + "\n")
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--task08",
        type=Path,
        default=ROOT / "analysis_08_chemistry/outputs/result.json",
    )
    parser.add_argument(
        "--chain",
        type=Path,
        default=ROOT / "analysis_08_chemistry/outputs/chain.npz",
    )
    parser.add_argument("--parameters", type=Path, default=HERE / "parameters.json")
    parser.add_argument("--output", type=Path, default=HERE / "outputs/result.json")
    args = parser.parse_args()
    run(args.task08, args.chain, args.parameters, args.output)


if __name__ == "__main__":
    main()
