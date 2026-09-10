import sys

from _common import OUT, REPRO, read_json, stage_dir

sys.path.insert(0, str(REPRO))

from analysis_07_scaling.run import (
    calculate,
    require_current_profile,
    write_outputs,
)


def main():
    output = stage_dir(7)
    require_current_profile(
        task5=OUT / "05" / "result.json",
        task3=OUT / "03" / "result.json",
        morphology=OUT / "04",
    )
    task5 = read_json(OUT / "05" / "result.json")
    task3 = read_json(OUT / "03" / "result.json")
    task6 = read_json(OUT / "06" / "result.json")
    result, rows = calculate(task5, task3, task6)
    write_outputs(result, rows, output=output)
    print(
        "stage 07: wrote production-rate inversions "
        f"({len(rows)} rows) to {output.name}/scaling.tsv"
    )


if __name__ == "__main__":
    main()
