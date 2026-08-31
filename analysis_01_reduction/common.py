import csv
import json
import os
import re
import shutil
import subprocess
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPRO = HERE.parent.resolve()
OUTPUT_PATH = HERE / "outputs"
if OUTPUT_PATH.is_symlink():
    raise RuntimeError(f"outputs is a symlink: {OUTPUT_PATH}")
OUTPUT = OUTPUT_PATH.resolve()
if OUTPUT == REPRO or REPRO not in OUTPUT.parents:
    raise RuntimeError(f"output escapes analysis_repro: {OUTPUT}")
BOOTSTRAP = OUTPUT / "tmp"
for name in ("tmp", "cache", "config", "matplotlib"):
    (BOOTSTRAP / name).mkdir(parents=True, exist_ok=True)
os.environ.update(
    TMPDIR=str(BOOTSTRAP / "tmp"),
    XDG_CACHE_HOME=str(BOOTSTRAP / "cache"),
    XDG_CONFIG_HOME=str(BOOTSTRAP / "config"),
    MPLCONFIGDIR=str(BOOTSTRAP / "matplotlib"),
)

from astropy.io import fits
from astropy.utils import iers
from astropy.wcs import WCS

# Use bundled leap-second data without probing the user cache.
iers.conf.auto_download = False
iers.conf.iers_leap_second_auto_url = ""
iers.conf.ietf_leap_second_auto_url = ""

PARAMS = json.loads((REPRO / "parameters.json").read_text())
CFG = json.loads((HERE / "settings.json").read_text())
CFG["obsid"] = PARAMS["obsid"]
CFG["comove"]["target"] = PARAMS["target"]
CFG["sources"] = {"radius_arcsec": PARAMS["source_mask_radius_arcsec"]}
for band in ("soft", "hard"):
    CFG["bands"][band]["pi"] = PARAMS[f"{band}_band_ev"]
CFG["espfilt"]["elow"], CFG["espfilt"]["ehigh"] = PARAMS["hard_band_ev"]


def sky_wcs(header):
    header = header.copy()
    if "RADECSYS" in header:
        header["RADESYS"] = header["RADECSYS"]
        del header["RADECSYS"]
    return WCS(header, fix=False).celestial


def setup_output():
    if REPRO not in OUTPUT.parents:
        raise RuntimeError(f"output escapes analysis_repro: {OUTPUT}")
    OUTPUT.mkdir(parents=True, exist_ok=True)
    scratch = OUTPUT / "tmp"
    for name in ("cache", "config", "matplotlib", "pfiles"):
        (scratch / name).mkdir(parents=True, exist_ok=True)
    os.environ.update(
        TMPDIR=str(scratch),
        XDG_CACHE_HOME=str(scratch / "cache"),
        XDG_CONFIG_HOME=str(scratch / "config"),
        MPLCONFIGDIR=str(scratch / "matplotlib"),
    )
    system_pfiles = os.environ.get("PFILES", "").split(";")[-1]
    os.environ["PFILES"] = f"{scratch / 'pfiles'};{system_pfiles}"
    return OUTPUT


def directory(name):
    root = setup_output()
    path = (root / name).resolve()
    if path != root and root not in path.parents:
        raise RuntimeError(f"path escapes analysis_repro: {path}")
    path.mkdir(parents=True, exist_ok=True)
    return path


def complete(name):
    return (directory(name) / ".done").exists()


def finish(name):
    (directory(name) / ".done").write_text("complete\n")


def valid_fits(path):
    path = Path(path)
    if not path.exists() or path.stat().st_size == 0:
        return False
    try:
        with fits.open(path):
            pass
    except OSError:
        return False
    return True


def sas_environment():
    command = (
        f'export HEADAS="{CFG["headas"]}" && '
        'source "$HEADAS/headas-init.sh" && '
        f'source "{CFG["sas"]}" && env -0'
    )
    result = subprocess.run(["bash", "-lc", command], capture_output=True)
    if result.returncode:
        raise RuntimeError(result.stderr.decode(errors="replace"))
    env = dict(
        item.split("=", 1) for item in result.stdout.decode().split("\0") if "=" in item
    )
    env.update(
        SAS_CCFPATH=CFG["ccf"],
        SAS_VERBOSITY="4",
        ATOMDB=CFG["atomdb"],
    )
    system_pfiles = env.get("PFILES", "").split(";")[-1]
    pfiles = setup_output() / "tmp" / "pfiles"
    pfiles.mkdir(exist_ok=True)
    env["PFILES"] = f"{pfiles};{system_pfiles}"
    init = setup_output() / "01_init"
    if (init / "ccf.cif").exists():
        env["SAS_CCF"] = str(init / "ccf.cif")
    summaries = sorted((init / "odf").glob("*SUM.SAS"))
    if summaries:
        env["SAS_ODF"] = str(init / "odf" / summaries[0].name)
    env.update(
        TMPDIR=os.environ["TMPDIR"],
        XDG_CACHE_HOME=os.environ["XDG_CACHE_HOME"],
        XDG_CONFIG_HOME=os.environ["XDG_CONFIG_HOME"],
    )
    return env


def sas(tool, cwd, environment=None, **parameters):
    cwd = Path(cwd).resolve()
    if cwd != OUTPUT and OUTPUT not in cwd.parents:
        raise RuntimeError(f"SAS working directory escapes analysis_repro: {cwd}")
    argv = [tool]
    for key, value in parameters.items():
        if isinstance(value, bool):
            value = "yes" if value else "no"
        elif isinstance(value, (list, tuple)):
            value = " ".join(map(str, value))
        argv.append(f"{key}={value}")
    logs = directory("logs")
    log = logs / f"{tool}-{time.time_ns()}.log"
    env = sas_environment()
    if environment:
        env.update(environment)
    system_pfiles = env.get("PFILES", "").split(";")[-1]
    pfiles = cwd / ".pfiles"
    pfiles.mkdir(exist_ok=True)
    env["PFILES"] = f"{pfiles};{system_pfiles}"
    try:
        with log.open("wb") as stream:
            stream.write((" ".join(argv) + "\n\n").encode())
            stream.flush()
            result = subprocess.run(
                argv, cwd=cwd, env=env, stdout=stream, stderr=subprocess.STDOUT
            )
    finally:
        shutil.rmtree(pfiles, ignore_errors=True)
    text = log.read_text(errors="replace")
    errors = [
        line
        for line in text.splitlines()
        if line.startswith("** ") and re.search(r":\s+(?:fatal\s+)?error\b", line, re.I)
    ]
    if result.returncode or errors:
        tail = "\n".join(text.splitlines()[-30:])
        raise RuntimeError(f"{tool} failed; see {log}\n{tail}")
    log.unlink()
    print(f"{tool}: {time.strftime('%H:%M:%S')}", flush=True)


def write_tsv(path, columns, rows):
    with Path(path).open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=columns, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def read_tsv(path):
    with Path(path).open(newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))
