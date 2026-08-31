import csv
import os
import re
from pathlib import Path

from common import OUTPUT, REPRO

GATES = ("06_soft_protons", "07_comove", "08_figure")
WORKSPACES = (
    "01_init",
    "02_reprocess",
    "03_events",
    "04_sources",
    "05_background",
    "06_soft_protons",
)
MANIFESTS = (
    "02_reprocess/exposures.tsv",
    "03_events/frames.tsv",
    "03_events/selection.tsv",
    "04_sources/events.tsv",
    "04_sources/sources.tsv",
    "05_background/background.tsv",
    "06_soft_protons/fits.tsv",
    "06_soft_protons/soft_protons.tsv",
)
FILES = ("01_init/ccf.cif", "02_reprocess/atthk.fits")
RESUME_GLOBS = (
    "02_reprocess/emproc/*ImagingEvts.ds",
    "02_reprocess/epproc/*ImagingEvts.ds",
    "02_reprocess/epproc_oot/*ImagingEvts.ds",
    "02_reprocess/filter/*/*-gti.fits",
    "02_reprocess/filter/*/*-allevc.fits",
    "02_reprocess/filter/*/*-allevcoot.fits",
    "02_reprocess/filter/*/*-corevc.fits",
    "04_sources/sources.reg",
)
ERROR = re.compile(r"^\*\* .*:\s+(?:fatal\s+)?error\b", re.I | re.M)
SUCCESS = re.compile(
    r"\bended:\s+\d{4}-\d\d-\d\dT|^\*\* .+\s+Finished\s*$|"
    r"successfully written RSP_MATRIX data\s*$",
    re.I | re.M,
)


def absolute(path):
    return Path(os.path.abspath(path))


def under(path, root):
    path, root = absolute(path), absolute(root)
    return path != root and root in path.parents


def keep_path(path, root, keep):
    path = absolute(path)
    if not under(path, root):
        raise RuntimeError(f"retained path escapes outputs: {path}")
    if not path.exists() and not path.is_symlink():
        raise FileNotFoundError(path)
    while path != root:
        keep.add(path)
        path = path.parent


def manifest_paths(path, root):
    with path.open(newline="") as stream:
        rows = csv.DictReader(stream, delimiter="\t")
        for row in rows:
            for value in row.values():
                if not value or not (
                    Path(value).is_absolute() or "/" in value or "\\" in value
                ):
                    continue
                candidate = Path(value)
                if not candidate.is_absolute():
                    candidate = path.parent / candidate
                if under(candidate, root):
                    yield candidate


def walk(root):
    for parent, directories, files in os.walk(root, followlinks=False):
        parent = Path(parent)
        for name in list(directories):
            path = parent / name
            if path.is_symlink():
                directories.remove(name)
                yield path
            else:
                yield path
        for name in files:
            yield parent / name


def successful_log(path):
    text = path.read_text(errors="replace")
    return bool(SUCCESS.search(text) and not ERROR.search(text))


def remove(path, root):
    path, root = absolute(path), absolute(root)
    if not under(path, root):
        raise RuntimeError(f"deletion escapes outputs: {path}")
    parent = path.parent
    while parent != root:
        if parent.is_symlink():
            raise RuntimeError(f"deletion crosses symlink: {path}")
        parent = parent.parent
    if path.is_symlink() or path.is_file():
        path.unlink()
    elif path.is_dir():
        path.rmdir()


def prerequisites(background):
    with background.open(newline="") as stream:
        for row in csv.DictReader(stream, delimiter="\t"):
            image = Path(row["events_image"])
            yield image.with_name(image.name.replace("-fovimsky-", "-fovimdet-"))
            if row["det"] == "pn":
                source = Path(row["source_pha"])
                yield source.with_name(
                    source.name.replace("-fovt.pi", "-fovtootsub.pi")
                )


def prune(root=OUTPUT):
    root = Path(root)
    if root.is_symlink():
        raise RuntimeError(f"outputs is a symlink: {root}")
    root = root.resolve()
    if root == REPRO or REPRO not in root.parents:
        raise RuntimeError(f"prune target escapes analysis_repro: {root}")
    for stage in GATES:
        if not (root / stage / ".done").is_file():
            raise RuntimeError(f"cannot prune before {stage} completes")

    keep = set()
    for relative in MANIFESTS:
        manifest = root / relative
        keep_path(manifest, root, keep)
        for path in manifest_paths(manifest, root):
            keep_path(path, root, keep)
    for path in walk(root):
        if path.name == ".done":
            keep_path(path, root, keep)
    odf = root / "01_init" / "odf"
    keep_path(odf, root, keep)
    for path in walk(odf):
        keep_path(path, root, keep)
    for relative in FILES:
        keep_path(root / relative, root, keep)
    for pattern in RESUME_GLOBS:
        matches = tuple(root.glob(pattern))
        if not matches:
            raise FileNotFoundError(f"resume pattern has no matches: {pattern}")
        for path in matches:
            keep_path(path, root, keep)
    for path in prerequisites(root / "05_background" / "background.tsv"):
        keep_path(path, root, keep)

    candidates = set()
    for name in WORKSPACES:
        candidates.update(walk(root / name))
    bootstrap = root / "tmp"
    if bootstrap.exists():
        candidates.update(walk(bootstrap))
        candidates.add(bootstrap)
    logs = root / "logs"
    if logs.exists():
        candidates.update(
            path for path in logs.iterdir() if path.is_file() and successful_log(path)
        )
    for path in walk(root):
        if path.name == ".pfiles" and path.is_dir():
            candidates.update(walk(path))
            candidates.add(path)

    removed = 0
    leaves = [
        path
        for path in candidates
        if path not in keep and (path.is_file() or path.is_symlink())
    ]
    for path in leaves:
        remove(path, root)
        removed += 1
    directories = sorted(
        (path for path in candidates if path not in keep and path.is_dir()),
        key=lambda path: len(path.parts),
        reverse=True,
    )
    for path in directories:
        try:
            remove(path, root)
        except OSError:
            pass
    return removed
