import os
import re
import shutil
from pathlib import Path

from common import CFG, directory, finish, sas

STAGE = "01_init"


def snapshot(source):
    return {
        path.name: (path.stat().st_size, path.stat().st_mtime_ns)
        for path in source.iterdir()
        if path.is_file()
    }


def copy_odf(source, target, before):
    target.mkdir(parents=True, exist_ok=True)
    for name, (size, _) in before.items():
        src, dst = source / name, target / name
        if dst.exists() and dst.stat().st_size == size:
            continue
        partial = dst.with_suffix(dst.suffix + ".partial")
        shutil.copyfile(src, partial)
        os.replace(partial, dst)


def install_summary(target):
    source = Path(CFG["summary"])
    destination = target / source.name
    shutil.copyfile(source, destination)
    text = destination.read_text()
    text = re.sub(r"^PATH\s+.*$", f"PATH {target}/", text, count=1, flags=re.MULTILINE)
    destination.write_text(text)
    if not any(line == f"PATH {target}/" for line in text.splitlines()):
        raise RuntimeError("SUM.SAS PATH was not rewritten")
    return destination


def validate_summary(path):
    lines = path.read_text(errors="replace").splitlines()
    obsid = next(
        (line.split()[0] for line in lines if "Observation Identifier" in line), ""
    )
    if obsid != CFG["obsid"]:
        raise RuntimeError(f"summary observation is {obsid!r}")
    timing = (
        "Scheduled Start Time",
        "Scheduled End Time",
        "Actual Start Time",
        "Actual End Time",
    )
    for index, line in enumerate(lines):
        if line == "EXPOSURE":
            block = lines[index : index + 14]
            if any(not any(key in item for item in block) for key in timing):
                raise RuntimeError(f"truncated summary exposure at line {index + 1}")


def run():
    work = directory(STAGE)
    source = Path(CFG["odf"])
    before = snapshot(source)
    copy_odf(source, work / "odf", before)
    summary = install_summary(work / "odf")
    validate_summary(summary)

    ccf = work / "ccf.cif"
    if not ccf.exists():
        sas(
            "cifbuild",
            work,
            analysisdate=CFG["analysis_date"],
            withccfpath=False,
            category="XMMCCF",
            calindexset=ccf,
            fullpath=True,
            environment={"SAS_ODF": str(work / "odf")},
        )
    if snapshot(source) != before:
        raise RuntimeError("the original ODF changed during initialization")
    finish(STAGE)
