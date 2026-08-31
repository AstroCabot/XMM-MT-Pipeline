import argparse
import subprocess
from pathlib import Path

from common import CFG, OUTPUT, REPRO, complete, sas_environment
from init import snapshot

STAGES = ("init", "events", "sources", "background", "soft_protons", "comove", "figure")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "stage", choices=("all", "prune", *STAGES), nargs="?", default="all"
    )
    args = parser.parse_args()
    if args.stage == "prune":
        from prune import prune

        print(f"prune: removed {prune()} scratch items", flush=True)
        return
    selected = STAGES if args.stage == "all" else (args.stage,)
    subprocess.run(
        ["sasversion"],
        env=sas_environment(),
        check=True,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )
    if "soft_protons" in selected:
        __import__("soft_protons").xspec_module()
    sources = {name: Path(CFG[name]) for name in ("odf", "pps", "mt_pps")}
    before = {name: snapshot(path) for name, path in sources.items()}
    try:
        for name in selected:
            module = __import__(name)
            if complete(module.STAGE):
                print(f"{name}: complete", flush=True)
                continue
            module.run()
    finally:
        changed = [
            name for name, path in sources.items() if snapshot(path) != before[name]
        ]
        if changed:
            raise RuntimeError(
                f"original input changed during the run: {', '.join(changed)}"
            )
        for path in OUTPUT.rglob("*"):
            if path.is_symlink() and REPRO not in path.resolve().parents:
                raise RuntimeError(f"output link escapes analysis_repro: {path}")


if __name__ == "__main__":
    main()
