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
    source_bands = cfg.get("map_source_detection_bands", ["3000", "4000", "5000"])
    if isinstance(source_bands, list):
        source_bands_text = " ".join(str(item) for item in source_bands)
    else:
        source_bands_text = str(source_bands)
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
        "MAP_BIN_PHYS": cfg.get("map_bin_phys", 80),
        "MAP_PAD_FRAC": cfg.get("map_pad_frac", 0.08),
        "MAP_SOFT_HARD_SPLIT_EV": cfg.get("map_soft_hard_split_ev", cfg.get("qc_soft_hard_split_ev", 1000)),
        "MAP_PI_MIN": cfg.get("map_pi_min", 201),
        "MAP_PI_MAX": cfg.get("map_pi_max", 12000),
        "MAP_EEXPMAP_ATTREBIN": cfg.get("map_eexpmap_attrebin", "0.020626481"),
        "MAP_SOURCE_DETECTION_BANDS": source_bands_text,
        "MAP_SOURCE_DETECTION_IDBAND": cfg.get("map_source_detection_idband", 0),
        "MAP_EBOX_LIKEMIN_LOCAL": cfg.get("map_ebox_likemin_local", 8),
        "MAP_EBOX_BOXSIZE": cfg.get("map_ebox_boxsize", 5),
        "MAP_EBOX_NRUNS": cfg.get("map_ebox_nruns", 3),
        "MAP_BKG_FITMETHOD": cfg.get("map_bkg_fitmethod", "model"),
        "MAP_BKG_MLMIN": cfg.get("map_bkg_mlmin", 1),
        "MAP_BKG_SCUT": cfg.get("map_bkg_scut", 0.01),
        "MAP_BKG_NSPLINENODES": cfg.get("map_bkg_nsplinenodes", 12),
        "MAP_BKG_NFITRUN": cfg.get("map_bkg_nfitrun", 3),
        "MAP_BKG_EXCESSSIGMA": cfg.get("map_bkg_excesssigma", 4),
        "MAP_BKG_SNRMIN": cfg.get("map_bkg_snrmin", 30),
        "MAP_BKG_SMOOTHSIGMA": cfg.get("map_bkg_smoothsigma", 6),
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


def band_expression(band: dict[str, Any], detector: str, ranges: list[list[Any]] | None = None) -> str:
    family = detector_family(detector)
    spec = band.get(family, {})
    if not isinstance(spec, dict):
        die(f"Band {band.get('label', '')} is missing {family} filter settings")
    terms = [pi_expr(ranges if ranges is not None else band["pi_ranges"]), pattern_expr(spec)]
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


def band_pi_bounds(band: dict[str, Any]) -> tuple[int, int]:
    override = band.get("exposure_pi_range")
    ranges = override if override is not None else band["pi_ranges"]
    if isinstance(ranges, list) and len(ranges) == 2 and not isinstance(ranges[0], list):
        ranges = [ranges]
    if not isinstance(ranges, list) or not ranges:
        die(f"Band {band.get('label', '')} has bad exposure_pi_range")
    bounds: list[tuple[int, int]] = []
    for item in ranges:
        if not isinstance(item, list) or len(item) != 2:
            die(f"Band {band.get('label', '')} has bad PI range: {item}")
        bounds.append((int(item[0]), int(item[1])))
    return min(lo for lo, _hi in bounds), max(hi for _lo, hi in bounds)


def maps_band_table(path: str) -> None:
    print("label\tdescription\tpimin\tpimax")
    for band in map_bands(path):
        print(f"{band['label']}\t{band['description']}\t{band['pimin']}\t{band['pimax']}")


def maps_clean_labels(path: str, label: str) -> None:
    selected = None
    for band in map_bands(path):
        if band["label"] == label:
            selected = band
            break
    if selected is None:
        die(f"Unknown maps band label: {label}")

    lo, hi = int(selected["pimin"]), int(selected["pimax"])
    for band in clean_bands(path):
        if intersect_pi_ranges(band["pi_ranges"], lo, hi):
            print(band["label"])


def map_pi_limits(cfg: dict[str, Any]) -> tuple[int, int]:
    bands = clean_bands_from_config(cfg)
    pimin = min(int(lo) for band in bands for lo, _hi in band["pi_ranges"])
    pimax = max(int(hi) for band in bands for _lo, hi in band["pi_ranges"])
    return int(cfg.get("map_pi_min", pimin)), int(cfg.get("map_pi_max", pimax))


def clean_bands_from_config(cfg: dict[str, Any]) -> list[dict[str, Any]]:
    pathless = dict(cfg)
    tmp = Path("__unused_config_path__")
    # Reuse clean_bands validation without exposing a second schema path.
    bands = pathless.get("clean_bands", [])
    if not isinstance(bands, list) or not bands:
        die("Missing non-empty clean_bands list in config")
    return bands


def map_bands(path: str) -> list[dict[str, Any]]:
    cfg = load_config(path)
    pimin, pimax = map_pi_limits(cfg)
    split = int(cfg.get("map_soft_hard_split_ev", cfg.get("qc_soft_hard_split_ev", 1000)))
    if not (pimin <= split < pimax):
        die("map_soft_hard_split_ev must fall between map_pi_min and map_pi_max")
    return [
        {"label": "soft", "description": f"{pimin}-{split} eV", "pimin": pimin, "pimax": split},
        {"label": "hard", "description": f"{split + 1}-{pimax} eV", "pimin": split + 1, "pimax": pimax},
    ]


def intersect_pi_ranges(ranges: list[list[Any]], lo: int, hi: int) -> list[list[int]]:
    out: list[list[int]] = []
    for item in ranges:
        rlo, rhi = int(item[0]), int(item[1])
        ilo, ihi = max(rlo, lo), min(rhi, hi)
        if ilo <= ihi:
            out.append([ilo, ihi])
    return out


def maps_expression_text(path: str, detector: str, label: str) -> str:
    selected = None
    for band in map_bands(path):
        if band["label"] == label:
            selected = band
            break
    if selected is None:
        die(f"Unknown maps band label: {label}")

    terms: list[str] = []
    lo, hi = int(selected["pimin"]), int(selected["pimax"])
    for band in clean_bands(path):
        ranges = intersect_pi_ranges(band["pi_ranges"], lo, hi)
        if ranges:
            terms.append(f"({band_expression(band, detector, ranges)})")
    if not terms:
        die(f"No clean-band filters overlap maps band: {label}")
    return terms[0] if len(terms) == 1 else "(" + "||".join(terms) + ")"


def maps_expression(path: str, detector: str, label: str) -> None:
    print(maps_expression_text(path, detector, label))


def maps_grid(outdir_arg: str, bin_phys_arg: str, pad_frac_arg: str, manifest_dir: str, detectors_arg: str = "PN") -> None:
    import math
    import numpy as np

    events = manifest_events(manifest_dir, detectors_arg)
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
        die("Could not determine a valid map grid from clean event lists")

    bin_phys = float(bin_phys_arg)
    pad_frac = float(pad_frac_arg)
    if bin_phys <= 0:
        die("map_bin_phys must be positive")

    pad_x = max(pad_frac * (xmax - xmin), bin_phys)
    pad_y = max(pad_frac * (ymax - ymin), bin_phys)
    xlo = math.floor((xmin - pad_x) / bin_phys) * bin_phys
    xhi = math.ceil((xmax + pad_x) / bin_phys) * bin_phys
    ylo = math.floor((ymin - pad_y) / bin_phys) * bin_phys
    yhi = math.ceil((ymax + pad_y) / bin_phys) * bin_phys
    nx = int(round((xhi - xlo) / bin_phys))
    ny = int(round((yhi - ylo) / bin_phys))

    outdir = Path(outdir_arg)
    outdir.mkdir(parents=True, exist_ok=True)
    env = [
        "# Generated by reduction_v6/tools.py maps-grid",
        f"MAP_GRID_BIN_PHYS={bin_phys:g}",
        f"MAP_GRID_X_MIN={xlo:g}",
        f"MAP_GRID_X_MAX={xhi:g}",
        f"MAP_GRID_Y_MIN={ylo:g}",
        f"MAP_GRID_Y_MAX={yhi:g}",
        f"MAP_GRID_NX={nx}",
        f"MAP_GRID_NY={ny}",
    ]
    summary = [
        f"event_lists\t{len(events)}",
        f"data_bounds\t{xmin:g}\t{xmax:g}\t{ymin:g}\t{ymax:g}",
        f"grid_bounds\t{xlo:g}\t{xhi:g}\t{ylo:g}\t{yhi:g}",
        f"bin_phys\t{bin_phys:g}",
        f"pixels\t{nx}\t{ny}",
    ]
    (outdir / "grid.env").write_text("\n".join(env) + "\n", encoding="utf-8")
    (outdir / "grid_summary.tsv").write_text("\n".join(summary) + "\n", encoding="utf-8")


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


def fits_table_hdu(path_arg: str, preferred_arg: str = "SRCLIST") -> None:
    from astropy.io import fits

    path = Path(path_arg)
    if not path.is_file():
        die(f"Missing FITS file: {path}")
    preferred = str(preferred_arg).strip()
    with fits.open(path, memmap=True) as hdul:
        for hdu in hdul:
            if isinstance(hdu, (fits.BinTableHDU, fits.TableHDU)) and hdu.name == preferred:
                print(hdu.name)
                return
        for hdu in hdul:
            if isinstance(hdu, (fits.BinTableHDU, fits.TableHDU)) and str(hdu.name).strip():
                print(hdu.name)
                return
    die(f"No named table HDU found in {path}")


def fits_table_has_column(path_arg: str, table_arg: str, column_arg: str) -> None:
    from astropy.io import fits

    path = Path(path_arg)
    if not path.is_file():
        print("no")
        return
    table = str(table_arg).strip()
    column = str(column_arg).strip().upper()
    with fits.open(path, memmap=True) as hdul:
        try:
            hdu = hdul[table]
        except Exception:
            print("no")
            return
        names = {str(name).upper() for name in getattr(hdu.columns, "names", [])}
        print("yes" if column in names else "no")


def apply_image_mask(image_arg: str, mask_arg: str, out_arg: str, mode_arg: str = "image") -> None:
    import numpy as np
    from astropy.io import fits

    image_path = Path(image_arg)
    mask_path = Path(mask_arg)
    out_path = Path(out_arg)
    mode = str(mode_arg).strip().lower()
    if mode not in {"image", "mask"}:
        die("apply-image-mask mode must be image or mask")
    if not image_path.is_file():
        die(f"Missing image file: {image_path}")
    if not mask_path.is_file():
        die(f"Missing mask file: {mask_path}")

    def first_image(path: Path):
        hdul = fits.open(path, memmap=True)
        for idx, hdu in enumerate(hdul):
            if getattr(hdu, "data", None) is not None:
                arr = np.asarray(hdu.data)
                if arr.ndim >= 2 and arr.dtype.kind in {"i", "u", "f"}:
                    data = np.array(np.squeeze(arr), copy=True)
                    while data.ndim > 2:
                        data = data[0]
                    return hdul, idx, data
        hdul.close()
        die(f"No image data found in {path}")

    image_hdul, image_idx, image_data = first_image(image_path)
    mask_hdul, _mask_idx, mask_data = first_image(mask_path)
    try:
        if image_data.shape != mask_data.shape:
            die(f"Shape mismatch between {image_path} and {mask_path}")
        mask_bool = np.asarray(mask_data) > 0
        header = image_hdul[image_idx].header.copy()
        for key in ("XTENSION", "PCOUNT", "GCOUNT", "EXTNAME"):
            if key in header:
                del header[key]
        if mode == "mask":
            out_array = (((np.asarray(image_data) > 0) & mask_bool)).astype(np.uint8)
        else:
            out_array = np.where(mask_bool, np.asarray(image_data), 0)
            out_array = out_array.astype(image_data.dtype, copy=False)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        fits.PrimaryHDU(data=out_array, header=header).writeto(out_path, overwrite=True)
    finally:
        image_hdul.close()
        mask_hdul.close()


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


def scaled_linear_image(image):
    import numpy as np

    finite = image[np.isfinite(image)]
    if not finite.size:
        return np.zeros_like(image, dtype=np.uint8)
    vmin = float(np.percentile(finite, 1.0))
    vmax = float(np.percentile(finite, 99.0))
    if vmax <= vmin:
        vmax = vmin + 1.0
    scaled = np.clip((image - vmin) / (vmax - vmin), 0.0, 1.0)
    return np.flipud((255.0 * scaled).astype(np.uint8))


def fits_image(path: Path):
    import numpy as np
    from astropy.io import fits

    with fits.open(path, memmap=True) as hdul:
        data = None
        for hdu in hdul:
            if getattr(hdu, "data", None) is not None:
                arr = np.asarray(hdu.data)
                if arr.ndim >= 2 and arr.dtype.kind in {"i", "u", "f"}:
                    data = np.squeeze(arr).astype(float)
                    break
        if data is None:
            die(f"No image data found in {path}")
        while data.ndim > 2:
            data = data[0]
        return np.nan_to_num(data, nan=0.0, posinf=0.0, neginf=0.0)


def maps_qc(manifest_arg: str, outdir_arg: str) -> None:
    import numpy as np

    manifest = Path(manifest_arg)
    if not manifest.is_file():
        die(f"Maps manifest not found: {manifest}")
    rows = []
    lines = manifest.read_text(encoding="utf-8").splitlines()
    if not lines:
        die(f"Empty maps manifest: {manifest}")
    header = lines[0].split("\t")
    index = {name: idx for idx, name in enumerate(header)}
    for line in lines[1:]:
        parts = line.split("\t")
        if len(parts) >= len(header):
            rows.append(parts)
    if not rows:
        die(f"No map products listed in {manifest}")

    products = [
        ("counts", "sum", "log"),
        ("exposure_vig", "sum", "log"),
        ("exposure_unvig", "sum", "log"),
    ]
    outdir = Path(outdir_arg)
    outdir.mkdir(parents=True, exist_ok=True)
    summary = ["band\tproduct\tfiles\toperation\tpng"]

    bands = []
    for row in rows:
        band = row[index["band"]]
        if band not in bands:
            bands.append(band)

    for band in bands:
        band_rows = [row for row in rows if row[index["band"]] == band]
        for product, op, scale in products:
            acc = None
            used = 0
            for row in band_rows:
                path = Path(row[index[product]])
                if not path.is_file():
                    continue
                image = fits_image(path)
                if acc is None:
                    acc = np.zeros_like(image, dtype=float)
                if image.shape != acc.shape:
                    die(f"Shape mismatch in {product} for band {band}: {path}")
                if op == "max":
                    acc = np.maximum(acc, image)
                else:
                    acc += image
                used += 1
            if acc is None:
                continue
            png = outdir / f"{band}_{product}_mosaic.png"
            write_gray_png(png, scaled_linear_image(acc) if scale == "linear" else scaled_log_image(acc))
            summary.append(f"{band}\t{product}\t{used}\t{op}\t{png}")

    (outdir / "maps_qc_summary.tsv").write_text("\n".join(summary) + "\n", encoding="utf-8")


def cheese_qc(manifest_arg: str, outdir_arg: str) -> None:
    import numpy as np

    manifest = Path(manifest_arg)
    if not manifest.is_file():
        die(f"Cheese manifest not found: {manifest}")
    rows = []
    lines = manifest.read_text(encoding="utf-8").splitlines()
    if not lines:
        die(f"Empty cheese manifest: {manifest}")
    header = lines[0].split("\t")
    index = {name: idx for idx, name in enumerate(header)}
    for line in lines[1:]:
        parts = line.split("\t")
        if len(parts) >= len(header):
            rows.append(parts)
    if not rows:
        die(f"No cheese products listed in {manifest}")

    products = [
        ("cheese_mask", "max", "linear"),
        ("cheesed_counts", "sum", "log"),
        ("cheesed_exposure_vig", "sum", "log"),
    ]
    outdir = Path(outdir_arg)
    outdir.mkdir(parents=True, exist_ok=True)
    summary = ["band\tproduct\tfiles\toperation\tpng"]

    bands = []
    for row in rows:
        band = row[index["band"]]
        if band not in bands:
            bands.append(band)

    for band in bands:
        band_rows = [row for row in rows if row[index["band"]] == band]
        for product, op, scale in products:
            acc = None
            used = 0
            for row in band_rows:
                path = Path(row[index[product]])
                if not path.is_file():
                    continue
                image = fits_image(path)
                if acc is None:
                    acc = np.zeros_like(image, dtype=float)
                if image.shape != acc.shape:
                    die(f"Shape mismatch in {product} for band {band}: {path}")
                if op == "max":
                    acc = np.maximum(acc, image)
                else:
                    acc += image
                used += 1
            if acc is None:
                continue
            png = outdir / f"{band}_{product}_mosaic.png"
            write_gray_png(png, scaled_linear_image(acc) if scale == "linear" else scaled_log_image(acc))
            summary.append(f"{band}\t{product}\t{used}\t{op}\t{png}")

    (outdir / "cheese_qc_summary.tsv").write_text("\n".join(summary) + "\n", encoding="utf-8")


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
    elif cmd == "maps-band-table" and len(args) == 1:
        maps_band_table(*args)
    elif cmd == "maps-clean-labels" and len(args) == 2:
        maps_clean_labels(*args)
    elif cmd == "maps-expr" and len(args) == 3:
        maps_expression(*args)
    elif cmd == "maps-grid" and len(args) == 4:
        maps_grid(*args)
    elif cmd == "maps-grid" and len(args) == 5:
        maps_grid(args[0], args[1], args[2], args[3], args[4])
    elif cmd == "maps-qc" and len(args) == 2:
        maps_qc(*args)
    elif cmd == "cheese-qc" and len(args) == 2:
        cheese_qc(*args)
    elif cmd == "fits-rows" and len(args) == 1:
        fits_rows(*args)
    elif cmd == "fits-table-hdu" and len(args) in {1, 2}:
        fits_table_hdu(args[0], args[1] if len(args) == 2 else "SRCLIST")
    elif cmd == "fits-table-has-column" and len(args) == 3:
        fits_table_has_column(*args)
    elif cmd == "apply-image-mask" and len(args) in {3, 4}:
        apply_image_mask(args[0], args[1], args[2], args[3] if len(args) == 4 else "image")
    elif cmd == "event-mosaic" and len(args) in {2, 3, 4}:
        event_mosaic(args[0], args[1], args[2] if len(args) >= 3 else "PN", args[3] if len(args) == 4 else "1000")
    else:
        die(
            "Usage: tools.py shell CONFIG | rewrite-path SUM.SAS PATH | "
            "clean-band-labels CONFIG | clean-band-table CONFIG | clean-expr CONFIG DETECTOR LABEL | "
            "maps-band-table CONFIG | maps-clean-labels CONFIG MAP_LABEL | maps-expr CONFIG DETECTOR LABEL | "
            "maps-grid OUTDIR BIN_PHYS PAD_FRAC MANIFEST_DIR [DETECTORS] | "
            "maps-qc MAPS_MANIFEST OUTDIR | cheese-qc CHEESE_MANIFEST OUTDIR | "
            "fits-rows FITS | fits-table-hdu FITS [PREFERRED] | fits-table-has-column FITS TABLE COLUMN | "
            "apply-image-mask IMAGE MASK OUT [image|mask] | "
            "event-mosaic MANIFEST_DIR OUTDIR [DETECTORS] [SPLIT_EV]"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
