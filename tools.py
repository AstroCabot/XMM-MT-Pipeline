#!/usr/bin/env python3
"""Small helpers for reduction_v6/pipeline.sh."""

from __future__ import annotations

import json
import shlex
import sys
from pathlib import Path
from typing import Any


def die(message: str) -> None:
    raise SystemExit(message)


def load_config(path: str) -> dict[str, Any]:
    try:
        cfg = json.loads(Path(path).read_text(encoding="utf-8"))
    except FileNotFoundError:
        die(f"Config file not found: {path}")
    except json.JSONDecodeError as exc:
        die(f"Invalid JSON in {path}: {exc}")
    if not isinstance(cfg, dict):
        die("Config must be a JSON object")
    for key in ("workdir", "odfdir"):
        if not cfg.get(key):
            die(f"Missing required config key: {key}")
    return cfg


def yesno(value: Any, default: bool = False) -> str:
    if value is None:
        return "yes" if default else "no"
    if isinstance(value, bool):
        return "yes" if value else "no"
    return "yes" if str(value).strip().lower() in {"1", "y", "yes", "true"} else "no"


def quote(value: Any) -> str:
    return shlex.quote("" if value is None else str(value))


def normalize_detectors(value: Any) -> list[str]:
    text = "PN" if value in (None, "") else value
    tokens = [str(item) for item in text] if isinstance(text, list) else str(text).replace(",", " ").replace(";", " ").split()
    aliases = {
        "PN": ["PN"],
        "EPN": ["PN"],
        "M1": ["M1"],
        "MOS1": ["M1"],
        "EMOS1": ["M1"],
        "M2": ["M2"],
        "MOS2": ["M2"],
        "EMOS2": ["M2"],
        "ALL": ["PN", "M1", "M2"],
        "EPIC": ["PN", "M1", "M2"],
    }
    detectors: list[str] = []
    for token in tokens:
        key = token.upper()
        if key not in aliases:
            die(f"Unknown detector in config: {token}")
        for detector in aliases[key]:
            if detector not in detectors:
                detectors.append(detector)
    return detectors or ["PN"]


def emit(name: str, value: Any) -> None:
    print(f"{name}={quote(value)}")


def shell(path: str) -> None:
    cfg = load_config(path)
    values = {
        "WORKDIR": str(cfg["workdir"]).rstrip("/"),
        "ODFDIR": str(cfg["odfdir"]).rstrip("/"),
        "SHORTLINK_NAME": cfg.get("shortlink_name", ""),
        "KEEP_SHORTLINK": yesno(cfg.get("keep_shortlink"), True),
        "SHORTLINK_MAX_PATH": cfg.get("shortlink_max_path", 60),
        "SAS_SETUP_SCRIPT": cfg.get("sas_setup_script", ""),
        "SAS_CCFPATH_CONFIG": cfg.get("sas_ccfpath", ""),
        "SAS_VERBOSITY_CONFIG": cfg.get("sas_verbosity", ""),
        "DETECTORS": " ".join(normalize_detectors(cfg.get("detectors"))),
        "AHF_INPUT": cfg.get("ahf_input", ""),
        "QC_SOFT_HARD_SPLIT_EV": cfg.get("qc_soft_hard_split_ev", 1000),
        "CLEAN_GTI_ENABLED": yesno(cfg.get("clean_gti_enabled"), True),
        "CLEAN_GTI_RATE_CUT": cfg.get("clean_gti_rate_cut", 4.80001211),
        "CLEAN_GTI_TIMEBIN": cfg.get("clean_gti_timebin", 10),
        "CLEAN_GTI_PI_MIN": cfg.get("clean_gti_pi_min", 7000),
        "CLEAN_GTI_PI_MAX": cfg.get("clean_gti_pi_max", 15000),
        "CLEAN_GTI_PATTERN_MAX": cfg.get("clean_gti_pattern_max", 0),
        "LINK_ODF_CONSTITUENTS": yesno(cfg.get("link_odf_constituents"), True),
        "SKIP_ODF_LINK_PATTERNS": "|".join(cfg.get("skip_odf_link_patterns", [])),
    }
    for name, value in values.items():
        emit(name, value)


def rewrite_path(summary: str, odf_path: str) -> None:
    path = odf_path.rstrip("/") + "/"
    summary_path = Path(summary)
    old = summary_path.read_text(encoding="utf-8", errors="surrogateescape")
    lines = old.splitlines()
    for idx, line in enumerate(lines):
        if line.startswith("PATH "):
            lines[idx] = f"PATH {path}"
            break
    else:
        lines.append(f"PATH {path}")
    new = "\n".join(lines) + "\n"
    if new != old:
        summary_path.write_text(new, encoding="utf-8", errors="surrogateescape")


def manifest_events(manifest_arg: str, detectors_arg: str = "PN") -> list[Path]:
    source = Path(manifest_arg)
    if source.is_file():
        manifests = [source]
    else:
        manifests = []
        roots = [source]
        if (source / "manifest").is_dir():
            roots.append(source / "manifest")
        for detector in normalize_detectors(detectors_arg):
            for root in roots:
                candidates = [
                    root / f"{detector}_raw.txt",
                    root / f"{detector}_clean_files.txt",
                ]
                candidates.extend(sorted(root.glob(f"{detector}_*_clean_files.txt")))
                for path in candidates:
                    if path.is_file() and path not in manifests:
                        manifests.append(path)
    if not manifests:
        die(f"Missing event manifest(s): {manifest_arg}")

    paths: list[Path] = []
    seen: set[Path] = set()
    for manifest in manifests:
        for line in manifest.read_text(encoding="utf-8").splitlines():
            text = line.strip()
            if not text:
                continue
            path = Path(text)
            if path not in seen:
                seen.add(path)
                paths.append(path)
    return paths


def clean_bands(path: str) -> list[dict[str, Any]]:
    cfg = load_config(path)
    bands = cfg.get("clean_bands", [])
    if not isinstance(bands, list) or not bands:
        die("Missing non-empty clean_bands list in config")
    seen: set[str] = set()
    for band in bands:
        if not isinstance(band, dict):
            die("Each clean_bands entry must be an object")
        label = str(band.get("label", "")).strip()
        if not label:
            die("Each clean_bands entry needs a label")
        if label in seen:
            die(f"Duplicate clean band label: {label}")
        seen.add(label)
        ranges = band.get("pi_ranges", [])
        if not isinstance(ranges, list) or not ranges:
            die(f"Band {label} needs non-empty pi_ranges")
        for item in ranges:
            if not isinstance(item, list) or len(item) != 2:
                die(f"Bad pi_ranges entry in band {label}: {item}")
            int(item[0])
            int(item[1])
        for family in ("pn", "mos"):
            spec = band.get(family)
            if not isinstance(spec, dict):
                die(f"Band {label} is missing {family} filter settings")
            extra = spec.get("extra", [])
            if not isinstance(extra, list):
                die(f"Band {label} {family}.extra must be a list")
    return bands


def detector_family(detector: str) -> str:
    return "pn" if detector.upper() == "PN" else "mos"


def pi_expr(ranges: list[list[Any]]) -> str:
    parts = [f"(PI in [{int(lo)}:{int(hi)}])" for lo, hi in ranges]
    return parts[0] if len(parts) == 1 else "(" + "||".join(parts) + ")"


def pattern_expr(spec: dict[str, Any]) -> str:
    pmin = int(spec.get("pattern_min", 0))
    pmax = int(spec.get("pattern_max", pmin))
    if pmin == pmax:
        return f"(PATTERN=={pmin})"
    return f"((PATTERN>={pmin})&&(PATTERN<={pmax}))"


def band_expression(band: dict[str, Any], detector: str) -> str:
    family = detector_family(detector)
    spec = band.get(family, {})
    if not isinstance(spec, dict):
        die(f"Band {band.get('label', '')} is missing {family} filter settings")
    terms = [pi_expr(band["pi_ranges"]), pattern_expr(spec)]
    flag = str(spec.get("flag", "")).strip()
    if flag:
        terms.append(flag if flag.startswith("#") else f"({flag})")
    for extra in spec.get("extra", []):
        text = str(extra).strip()
        if text:
            terms.append(f"({text})")
    return "&&".join(terms)


def clean_expression_text(path: str, detector: str, label: str) -> str:
    for band in clean_bands(path):
        if band["label"] != label:
            continue
        return band_expression(band, detector)
    die(f"Unknown clean band label: {label}")


def clean_expression(path: str, detector: str, label: str) -> None:
    print(clean_expression_text(path, detector, label))


def clean_band_labels(path: str) -> None:
    for band in clean_bands(path):
        print(band["label"])


def clean_band_table(path: str) -> None:
    print("label\tdescription\tpn_expression\tmos_expression")
    for band in clean_bands(path):
        label = band["label"]
        desc = str(band.get("description", ""))
        pn = clean_expression_text(path, "PN", label)
        mos = clean_expression_text(path, "M1", label)
        print(f"{label}\t{desc}\t{pn}\t{mos}")


def event_hdu(hdul):
    return hdul["EVENTS"] if "EVENTS" in hdul else hdul[1]


def event_columns(path: Path):
    from astropy.io import fits

    hdul = fits.open(path, memmap=True)
    hdu = event_hdu(hdul)
    names = set(hdu.columns.names)
    for needed in ("X", "Y", "PI"):
        if needed not in names:
            hdul.close()
            die(f"{path} does not contain EVENTS column {needed}")
    return hdul, hdu.data["X"], hdu.data["Y"], hdu.data["PI"]


def valid_sky(x, y):
    import numpy as np

    return np.isfinite(x) & np.isfinite(y) & (abs(x) < 1.0e7) & (abs(y) < 1.0e7)


def event_bounds(events: list[Path]):
    import numpy as np

    xmin = ymin = np.inf
    xmax = ymax = -np.inf
    for path in events:
        hdul, x, y, _pi = event_columns(path)
        good = valid_sky(x, y)
        if np.any(good):
            xmin = min(xmin, float(x[good].min()))
            xmax = max(xmax, float(x[good].max()))
            ymin = min(ymin, float(y[good].min()))
            ymax = max(ymax, float(y[good].max()))
        hdul.close()
    if not np.isfinite([xmin, xmax, ymin, ymax]).all() or xmin == xmax or ymin == ymax:
        die("Could not determine a valid mosaic extent from event lists")

    data_bounds = (xmin, xmax, ymin, ymax)
    pad_x = max(0.08 * (xmax - xmin), 500.0)
    pad_y = max(0.08 * (ymax - ymin), 500.0)
    extent = (xmin - pad_x, xmax + pad_x, ymin - pad_y, ymax + pad_y)
    nx = 900
    ny = max(200, int(round(nx * (extent[3] - extent[2]) / (extent[1] - extent[0]))))
    return data_bounds, extent, nx, min(ny, 1200)


def fits_rows(path_arg: str) -> None:
    from astropy.io import fits

    path = Path(path_arg)
    if not path.is_file():
        print(-1)
        return
    with fits.open(path, memmap=True) as hdul:
        hdu = event_hdu(hdul)
        print(0 if hdu.data is None else len(hdu.data))


def write_gray_png(path: Path, image) -> None:
    import struct
    import zlib

    import numpy as np

    gray = np.asarray(image, dtype=np.uint8)
    height, width = gray.shape
    rows = b"".join(b"\x00" + gray[row].tobytes() for row in range(height))

    def chunk(tag: bytes, data: bytes) -> bytes:
        crc = zlib.crc32(tag + data) & 0xFFFFFFFF
        return struct.pack(">I", len(data)) + tag + data + struct.pack(">I", crc)

    png = (
        b"\x89PNG\r\n\x1a\n"
        + chunk(b"IHDR", struct.pack(">IIBBBBB", width, height, 8, 0, 0, 0, 0))
        + chunk(b"IDAT", zlib.compress(rows, level=6))
        + chunk(b"IEND", b"")
    )
    path.write_bytes(png)


def scaled_log_image(image):
    import numpy as np

    positive = image[image > 0]
    if not positive.size:
        return np.zeros_like(image, dtype=np.uint8)
    logged = np.zeros_like(image, dtype=float)
    logged[image > 0] = np.log10(image[image > 0])
    vmin = max(float(np.percentile(np.log10(positive), 1.0)), 0.0)
    vmax = max(float(np.percentile(np.log10(positive), 99.7)), vmin + 0.01)
    scaled = np.clip((logged - vmin) / (vmax - vmin), 0.0, 1.0)
    return np.flipud((255.0 * scaled).astype(np.uint8))


def event_mosaic(manifest_dir: str, outdir: str, detectors_arg: str = "PN", split_ev_arg: str = "1000") -> None:
    import numpy as np

    split_ev = float(split_ev_arg)
    detectors = normalize_detectors(detectors_arg)
    events = manifest_events(manifest_dir, " ".join(detectors))
    data_bounds, extent, nx, ny = event_bounds(events)

    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)
    summary = [
        f"detectors {' '.join(detectors)}",
        f"event_lists {len(events)}",
        f"data_bounds {data_bounds}",
        f"plot_extent {extent}",
        f"plot_pixels {nx} {ny}",
        f"soft_hard_split_ev {split_ev:g}",
    ]

    for tag, (lo, hi) in {"soft": (None, split_ev), "hard": (split_ev, None)}.items():
        image = np.zeros((ny, nx), dtype=float)
        total = 0
        for path in events:
            hdul, x, y, pi = event_columns(path)
            good = valid_sky(x, y)
            if lo is not None:
                good &= pi > lo
            if hi is not None:
                good &= pi < hi
            if np.any(good):
                hist, _xe, _ye = np.histogram2d(
                    x[good],
                    y[good],
                    bins=(nx, ny),
                    range=((extent[0], extent[1]), (extent[2], extent[3])),
                )
                image += hist.T
                total += int(good.sum())
            hdul.close()

        write_gray_png(out / f"{tag}_mosaic.png", scaled_log_image(image))
        summary.append(f"{tag}_events {total}")

    (out / "mosaic_summary.txt").write_text("\n".join(summary) + "\n", encoding="utf-8")


def main(argv: list[str]) -> int:
    cmd, args = (argv[1], argv[2:]) if len(argv) > 1 else ("", [])
    if cmd == "shell" and len(args) == 1:
        shell(*args)
    elif cmd == "rewrite-path" and len(args) == 2:
        rewrite_path(*args)
    elif cmd == "clean-band-labels" and len(args) == 1:
        clean_band_labels(*args)
    elif cmd == "clean-band-table" and len(args) == 1:
        clean_band_table(*args)
    elif cmd == "clean-expr" and len(args) == 3:
        clean_expression(*args)
    elif cmd == "fits-rows" and len(args) == 1:
        fits_rows(*args)
    elif cmd == "event-mosaic" and len(args) in {2, 3, 4}:
        event_mosaic(args[0], args[1], args[2] if len(args) >= 3 else "PN", args[3] if len(args) == 4 else "1000")
    else:
        die(
            "Usage: tools.py shell CONFIG | rewrite-path SUM.SAS PATH | "
            "clean-band-labels CONFIG | clean-band-table CONFIG | clean-expr CONFIG DETECTOR LABEL | "
            "fits-rows FITS | event-mosaic MANIFEST_DIR OUTDIR [DETECTORS] [SPLIT_EV]"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
