import json
import os
import sys
from pathlib import Path

sys.dont_write_bytecode = True

DEMO = Path(__file__).resolve().parent.parent
REPRO = DEMO.parent
DATA = DEMO / "data"
OUT = DEMO / "outputs"

_SCRATCH = ("TMPDIR", "XDG_CACHE_HOME", "XDG_CONFIG_HOME", "MPLCONFIGDIR")
_SAVED = {name: os.environ.get(name) for name in _SCRATCH}


def restore_env():
    for name, value in _SAVED.items():
        if value is None:
            os.environ.pop(name, None)
        else:
            os.environ[name] = value


def stage_dir(number):
    path = OUT / f"{number:02d}"
    path.mkdir(parents=True, exist_ok=True)
    return path


def write_json(path, payload):
    Path(path).write_text(json.dumps(payload, indent=2) + "\n")


def read_json(path):
    return json.loads(Path(path).read_text())
