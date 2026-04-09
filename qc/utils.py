"""Shared I/O helpers used by multiple QC modules.

Consolidates load_shell_env, read_list, fits_rows, and other small
utilities that were previously duplicated across checks.py and manifest.py.
"""
from __future__ import annotations

import shlex
import subprocess
from pathlib import Path

from astropy.io import fits


def load_shell_env(path: str) -> dict[str, str]:
    """Source a bash env file and return all exported variables."""
    quoted = shlex.quote(str(path))
    cmd = ["bash", "-lc", f"set -a; source {quoted} >/dev/null 2>&1; env -0"]
    proc = subprocess.run(cmd, capture_output=True, check=True)
    env: dict[str, str] = {}
    for chunk in proc.stdout.split(b"\x00"):
        if not chunk or b"=" not in chunk:
            continue
        key_b, val_b = chunk.split(b"=", 1)
        env[key_b.decode("utf-8", errors="replace")] = val_b.decode(
            "utf-8", errors="replace"
        )
    return env


def read_list(path: Path) -> list[str]:
    """Read a text file into a list of stripped non-empty lines."""
    if not path.exists():
        return []
    return [
        ln.strip()
        for ln in path.read_text(encoding="utf-8").splitlines()
        if ln.strip()
    ]


def fits_rows(path: Path) -> int | None:
    """Return the number of rows in the first binary-table extension."""
    if not path.exists():
        return None
    try:
        with fits.open(path) as hdul:
            if "EVENTS" in hdul:
                data = hdul["EVENTS"].data
                return 0 if data is None else int(len(data))
            for hdu in hdul[1:]:
                if isinstance(hdu, fits.BinTableHDU) and hdu.data is not None:
                    return int(len(hdu.data))
    except Exception:
        return None
    return None


def nonempty(path: Path) -> bool:
    """Return True if *path* exists and has non-zero size."""
    return path.exists() and path.stat().st_size > 0


def env_float(env: dict[str, str], key: str, default: float) -> float:
    val = env.get(key, "")
    return float(val) if val != "" else float(default)


def env_int(env: dict[str, str], key: str, default: int) -> int:
    val = env.get(key, "")
    return int(val) if val != "" else int(default)


def parse_image_bands(bands: str) -> list[tuple[str, int, int]]:
    """Parse an IMAGE_BANDS spec like ``'soft:200:1000;broad:300:2000'``."""
    out: list[tuple[str, int, int]] = []
    for entry in bands.split(";"):
        entry = entry.strip()
        if not entry:
            continue
        parts = entry.split(":")
        if len(parts) != 3:
            raise RuntimeError(f"Bad IMAGE_BANDS entry: {entry}")
        out.append((parts[0].strip(), int(parts[1]), int(parts[2])))
    return out


def first_band(env: dict[str, str]) -> str:
    """Return the label of the first configured image band."""
    for entry in env.get("IMAGE_BANDS", "broad:300:2000").split(";"):
        entry = entry.strip()
        if entry:
            return entry.split(":")[0].strip()
    return "broad"
