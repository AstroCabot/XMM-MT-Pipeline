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
    tokens = (
        [str(item) for item in text]
        if isinstance(text, list)
        else str(text).replace(",", " ").replace(";", " ").split()
    )
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
        "MAP_SOFT_HARD_SPLIT_EV": cfg.get(
            "map_soft_hard_split_ev", cfg.get("qc_soft_hard_split_ev", 1000)
        ),
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
        "MAPS_BACKGROUND_FRACTION": cfg.get("maps_background_fraction", 0.10),
        "MAPS_BACKGROUND_OUTSIDE_RADIUS_FRACTION": cfg.get(
            "maps_background_outside_radius_fraction", 0.5
        ),
        "CHEESE_EMIN": cfg.get("cheese_emin", 1000),
        "CHEESE_EMAX": cfg.get("cheese_emax", 4000),
        "CHEESE_MLMIN": cfg.get("cheese_mlmin", 10),
        "CHEESE_FLUX": cfg.get("cheese_flux", 0.2),
        "CHEESE_ECF_PN": cfg.get("cheese_ecf_pn", 3.2),
        "CHEESE_ECF_MOS": cfg.get("cheese_ecf_mos", 1.65),
        "CHEESE_EMASK_THRESH1": cfg.get("cheese_emask_thresh1", 0.3),
        "CHEESE_EMASK_THRESH2": cfg.get("cheese_emask_thresh2", 0.5),
        "CHEESE_REGION_ENERGYFRACTION": cfg.get("cheese_region_energyfraction", 0.9),
        "CHEESE_REGION_FIXEDRADIUS": cfg.get("cheese_region_fixedradius", 12),
        "CHEESE_REGION_BKGFRACTION": cfg.get("cheese_region_bkgfraction", 0.6),
        "CHEESE_SOURCE_RADIUS_PIX": cfg.get("cheese_source_radius_pix", 4),
        "CHEESE_ESPLINEMAP_NSPLINENODES": cfg.get("cheese_esplinemap_nsplinenodes", 20),
        "CHEESE_ESPLINEMAP_EXCESSSIGMA": cfg.get("cheese_esplinemap_excesssigma", 4),
        "CHEESE_ESPLINEMAP_NFITRUN": cfg.get("cheese_esplinemap_nfitrun", 3),
        "CHEESE_ESPLINEMAP_SNRMIN": cfg.get("cheese_esplinemap_snrmin", 30),
        "CHEESE_ESPLINEMAP_SMOOTHSIGMA": cfg.get("cheese_esplinemap_smoothsigma", 15),
        "CHEESE_ESPLINEMAP_SCUT": cfg.get("cheese_esplinemap_scut", 0.01),
        "CHEESE_ESPLINEMAP_MLMIN": cfg.get("cheese_esplinemap_mlmin", 1),
        "CHEESE_EBOXDETECT_NRUNS": cfg.get("cheese_eboxdetect_nruns", 3),
        "CHEESE_EBOXDETECT_BOXSIZE": cfg.get("cheese_eboxdetect_boxsize", 5),
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


def band_expression(
    band: dict[str, Any], detector: str, ranges: list[list[Any]] | None = None
) -> str:
    family = detector_family(detector)
    spec = band.get(family, {})
    if not isinstance(spec, dict):
        die(f"Band {band.get('label', '')} is missing {family} filter settings")
    terms = [
        pi_expr(ranges if ranges is not None else band["pi_ranges"]),
        pattern_expr(spec),
    ]
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
    if (
        isinstance(ranges, list)
        and len(ranges) == 2
        and not isinstance(ranges[0], list)
    ):
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
        print(
            f"{band['label']}\t{band['description']}\t{band['pimin']}\t{band['pimax']}"
        )


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


def cheese_clean_labels(path: str, emin: int, emax: int) -> None:
    """Emit clean-stage labels whose PI ranges overlap [emin, emax].

    Used by the cheese stage to pick the event files to merge into the
    cheese detection band (v8-style, default 1000-4000 eV).
    """
    lo, hi = int(emin), int(emax)
    if lo > hi:
        die(f"cheese PI range invalid: {lo} > {hi}")
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
    split = int(
        cfg.get("map_soft_hard_split_ev", cfg.get("qc_soft_hard_split_ev", 1000))
    )
    if not (pimin <= split < pimax):
        die("map_soft_hard_split_ev must fall between map_pi_min and map_pi_max")
    return [
        {
            "label": "soft",
            "description": f"{pimin}-{split} eV",
            "pimin": pimin,
            "pimax": split,
        },
        {
            "label": "hard",
            "description": f"{split + 1}-{pimax} eV",
            "pimin": split + 1,
            "pimax": pimax,
        },
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


def maps_grid(
    outdir_arg: str,
    bin_phys_arg: str,
    pad_frac_arg: str,
    manifest_dir: str,
    detectors_arg: str = "PN",
) -> None:
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
    (outdir / "grid_summary.tsv").write_text(
        "\n".join(summary) + "\n", encoding="utf-8"
    )


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
            if (
                isinstance(hdu, (fits.BinTableHDU, fits.TableHDU))
                and hdu.name == preferred
            ):
                print(hdu.name)
                return
        for hdu in hdul:
            if (
                isinstance(hdu, (fits.BinTableHDU, fits.TableHDU))
                and str(hdu.name).strip()
            ):
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


def combine_masks(detmask_arg: str, regmask_arg: str, out_arg: str) -> None:
    """Combine a detector validity mask and a region (source-exclusion) mask.

    The output is a binary uint8 mask where a pixel is 1 iff both inputs are
    > 0, matching v8's cheese combine-masks semantics. Shapes must match.
    """
    import numpy as np
    from astropy.io import fits

    det_path = Path(detmask_arg)
    reg_path = Path(regmask_arg)
    out_path = Path(out_arg)
    for p in (det_path, reg_path):
        if not p.is_file():
            die(f"Missing mask file: {p}")

    def first_image(path: Path):
        with fits.open(path, memmap=False) as hdul:
            for hdu in hdul:
                data = getattr(hdu, "data", None)
                if data is None:
                    continue
                arr = np.asarray(data)
                if arr.ndim >= 2 and arr.dtype.kind in {"i", "u", "f"}:
                    out = np.array(np.squeeze(arr), copy=True)
                    while out.ndim > 2:
                        out = out[0]
                    return out, hdu.header.copy()
        die(f"No image HDU in {path}")

    det_data, det_header = first_image(det_path)
    reg_data, _ = first_image(reg_path)
    if det_data.shape != reg_data.shape:
        die(
            f"Shape mismatch {det_path} {det_data.shape} vs {reg_path} {reg_data.shape}"
        )
    combined = ((det_data > 0) & (reg_data > 0)).astype(np.uint8)
    for key in ("XTENSION", "PCOUNT", "GCOUNT", "EXTNAME"):
        if key in det_header:
            del det_header[key]
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fits.PrimaryHDU(data=combined, header=det_header).writeto(out_path, overwrite=True)
    print(int(combined.sum()))


def fov_mask(exp_arg: str, out_arg: str) -> None:
    """Build a binary FOV mask from an exposure (or any image) FITS file.

    Output is a uint8 image with 1 where the input is finite and > 0
    (the per-frame field of view on the maps grid), 0 elsewhere. The
    output WCS is taken from the input image HDU so SAS region/regionmask
    can paint sky-coordinate regions on it without resampling.
    """
    import numpy as np
    from astropy.io import fits

    exp_path = Path(exp_arg)
    out_path = Path(out_arg)
    if not exp_path.is_file():
        die(f"Missing exposure image: {exp_path}")

    with fits.open(exp_path, memmap=False) as hdul:
        src_hdu = None
        for hdu in hdul:
            data = getattr(hdu, "data", None)
            if data is None:
                continue
            arr = np.asarray(data)
            if arr.ndim >= 2 and arr.dtype.kind in {"i", "u", "f"}:
                src_hdu = hdu
                break
        if src_hdu is None:
            die(f"No image HDU in {exp_path}")
        data = np.array(np.squeeze(src_hdu.data), copy=True)
        while data.ndim > 2:
            data = data[0]
        header = src_hdu.header.copy()

    mask = (np.isfinite(data) & (data > 0)).astype(np.uint8)
    for key in ("XTENSION", "PCOUNT", "GCOUNT", "EXTNAME", "BUNIT"):
        if key in header:
            del header[key]
    header["BUNIT"] = ("", "1=inside FOV, 0=outside")
    header["HISTORY"] = f"fov_mask from {exp_path.name} (>0)"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fits.PrimaryHDU(data=mask, header=header).writeto(out_path, overwrite=True)
    print(f"fov-mask: in_pixels={int(mask.size)} fov_pixels={int(mask.sum())}")


def mask_subtract(fov_arg: str, src_arg: str, out_arg: str) -> None:
    """Build a cheese mask = (fov > 0) AND (src == 0).

    Both inputs are read as 2D image arrays. The output preserves the FOV
    image's WCS header. Useful for combining a SAS regionmask source mask
    with the per-frame FOV mask to produce the final cheese mask.
    """
    import numpy as np
    from astropy.io import fits

    fov_path = Path(fov_arg)
    src_path = Path(src_arg)
    out_path = Path(out_arg)
    for p in (fov_path, src_path):
        if not p.is_file():
            die(f"Missing mask file: {p}")

    def first_image(path: Path):
        with fits.open(path, memmap=False) as hdul:
            for hdu in hdul:
                data = getattr(hdu, "data", None)
                if data is None:
                    continue
                arr = np.asarray(data)
                if arr.ndim >= 2 and arr.dtype.kind in {"i", "u", "f"}:
                    out = np.array(np.squeeze(arr), copy=True)
                    while out.ndim > 2:
                        out = out[0]
                    return out, hdu.header.copy()
        die(f"No image HDU in {path}")

    fov_data, fov_header = first_image(fov_path)
    src_data, _ = first_image(src_path)
    if fov_data.shape != src_data.shape:
        die(
            f"Shape mismatch {fov_path} {fov_data.shape} vs {src_path} {src_data.shape}"
        )
    out = ((fov_data > 0) & (src_data == 0)).astype(np.uint8)
    for key in ("XTENSION", "PCOUNT", "GCOUNT", "EXTNAME", "BUNIT"):
        if key in fov_header:
            del fov_header[key]
    fov_header["BUNIT"] = ("", "1=keep, 0=excluded source or outside FOV")
    fov_header["HISTORY"] = f"mask_subtract fov={fov_path.name} src={src_path.name}"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fits.PrimaryHDU(data=out, header=fov_header).writeto(out_path, overwrite=True)
    print(
        f"mask-subtract: fov_pixels={int((fov_data>0).sum())} "
        f"src_pixels={int((src_data>0).sum())} kept_pixels={int(out.sum())}"
    )


def cheese_make_mask(
    emllist_arg: str,
    fov_arg: str,
    output_arg: str,
    id_inst_arg: str,
    id_band_arg: str,
    mlmin_arg: str,
    fluxmin_e14_arg: str,
    radius_pix_arg: str,
) -> None:
    """Build a cheese mask by painting fixed-radius circles on a FOV image.

    Reads an emldetect SRCLIST, filters rows by ID_INST/ID_BAND/DET_ML/FLUX,
    and zeros out a circular region of `radius_pix` image pixels around each
    surviving source's (X_IMA, Y_IMA). The starting mask is the FOV image
    thresholded at >0 (typically the maps-stage exposure_vig FITS, which is
    the reliable per-frame footprint on the maps grid). Output is uint8 with
    1 = keep (background), 0 = exclude (outside FOV or inside source).

    This replaces the SAS region + regionmask + combine-masks chain, which
    silently fails when emldetect outputs lack BG_MAP / EXP_MAP, leaving
    base_mask = detmask (no source holes punched).
    """
    import numpy as np
    from astropy.io import fits

    emllist = Path(emllist_arg)
    fov_path = Path(fov_arg)
    output = Path(output_arg)
    id_inst = int(id_inst_arg)
    id_band = int(id_band_arg)
    mlmin = float(mlmin_arg)
    fluxmin = float(fluxmin_e14_arg) * 1e-14
    radius_pix = float(radius_pix_arg)

    if not fov_path.is_file():
        die(f"Cheese FOV reference image not found: {fov_path}")

    with fits.open(fov_path, memmap=False) as hdul:
        fov_hdu = None
        for hdu in hdul:
            if getattr(hdu, "data", None) is not None:
                arr = np.asarray(hdu.data)
                if arr.ndim >= 2 and arr.dtype.kind in {"i", "u", "f"}:
                    fov_hdu = hdu
                    break
        if fov_hdu is None:
            die(f"No image HDU in FOV reference: {fov_path}")
        fov_data = np.array(np.squeeze(fov_hdu.data), copy=True)
        while fov_data.ndim > 2:
            fov_data = fov_data[0]
        fov_header = fov_hdu.header.copy()

    ny, nx = fov_data.shape
    mask = (np.isfinite(fov_data) & (fov_data > 0)).astype(np.uint8)

    n_kept = 0
    n_total = 0
    if emllist.is_file():
        try:
            with fits.open(emllist, memmap=False) as hdul:
                table_hdu = None
                for hdu in hdul:
                    if isinstance(hdu, (fits.BinTableHDU, fits.TableHDU)):
                        if str(hdu.name).upper() == "SRCLIST":
                            table_hdu = hdu
                            break
                if table_hdu is None:
                    for hdu in hdul:
                        if isinstance(hdu, (fits.BinTableHDU, fits.TableHDU)):
                            table_hdu = hdu
                            break
                if table_hdu is not None:
                    cols = {c.upper() for c in table_hdu.data.names}
                    need = {"ID_INST", "ID_BAND", "X_IMA", "Y_IMA", "DET_ML", "FLUX"}
                    if need.issubset(cols):
                        t = table_hdu.data
                        inst = np.asarray(t["ID_INST"], dtype=int)
                        band = np.asarray(t["ID_BAND"], dtype=int)
                        ml = np.asarray(t["DET_ML"], dtype=float)
                        flux = np.asarray(t["FLUX"], dtype=float)
                        xs = np.asarray(t["X_IMA"], dtype=float)
                        ys = np.asarray(t["Y_IMA"], dtype=float)
                        n_total = int(((inst == id_inst) & (band == id_band)).sum())
                        keep = (
                            (inst == id_inst)
                            & (band == id_band)
                            & np.isfinite(ml)
                            & (ml >= mlmin)
                            & np.isfinite(flux)
                            & (flux >= fluxmin)
                            & np.isfinite(xs)
                            & np.isfinite(ys)
                        )
                        xs, ys = xs[keep], ys[keep]
                        n_kept = int(keep.sum())
                        if n_kept > 0:
                            r2 = radius_pix * radius_pix
                            for cx, cy in zip(xs, ys):
                                # X_IMA / Y_IMA are 1-indexed FITS pixel coords.
                                cx0 = float(cx) - 1.0
                                cy0 = float(cy) - 1.0
                                x0 = max(int(np.floor(cx0 - radius_pix)), 0)
                                x1 = min(int(np.ceil(cx0 + radius_pix)) + 1, nx)
                                y0 = max(int(np.floor(cy0 - radius_pix)), 0)
                                y1 = min(int(np.ceil(cy0 + radius_pix)) + 1, ny)
                                if x1 <= x0 or y1 <= y0:
                                    continue
                                sub_y = np.arange(y0, y1)[:, None]
                                sub_x = np.arange(x0, x1)[None, :]
                                circ = (sub_x - cx0) ** 2 + (sub_y - cy0) ** 2 <= r2
                                mask[y0:y1, x0:x1][circ] = 0
        except Exception as exc:  # noqa: BLE001
            print(
                f"cheese-mask: failed to read SRCLIST {emllist}: {exc}", file=sys.stderr
            )

    out_hdr = fov_header
    for key in ("XTENSION", "PCOUNT", "GCOUNT", "EXTNAME", "BUNIT"):
        if key in out_hdr:
            del out_hdr[key]
    out_hdr["BUNIT"] = ("", "1=keep background, 0=exclude")
    out_hdr["HISTORY"] = (
        f"cheese_make_mask id_inst={id_inst} id_band={id_band} "
        f"mlmin={mlmin:g} fluxmin={fluxmin:.3e} radius_pix={radius_pix:g} "
        f"sources_in_band={n_total} sources_used={n_kept}"
    )

    output.parent.mkdir(parents=True, exist_ok=True)
    fits.PrimaryHDU(data=mask, header=out_hdr).writeto(output, overwrite=True)
    fov_pix = int(((np.isfinite(fov_data) & (fov_data > 0))).sum())
    kept = int(mask.sum())
    print(
        f"cheese-mask: fov_pixels={fov_pix} kept_pixels={kept} "
        f"sources_in_band={n_total} sources_used={n_kept}"
    )


def apply_image_mask(
    image_arg: str, mask_arg: str, out_arg: str, mode_arg: str = "image"
) -> None:
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


def write_rgb_png(path: Path, image) -> None:
    import struct
    import zlib

    import numpy as np

    rgb = np.asarray(image, dtype=np.uint8)
    if rgb.ndim != 3 or rgb.shape[2] != 3:
        die(f"write_rgb_png expects (H, W, 3) array, got {rgb.shape}")
    height, width, _ = rgb.shape
    rows = b"".join(b"\x00" + rgb[row].tobytes() for row in range(height))

    def chunk(tag: bytes, data: bytes) -> bytes:
        crc = zlib.crc32(tag + data) & 0xFFFFFFFF
        return struct.pack(">I", len(data)) + tag + data + struct.pack(">I", crc)

    png = (
        b"\x89PNG\r\n\x1a\n"
        # color_type=2 (truecolor, 8-bit RGB)
        + chunk(b"IHDR", struct.pack(">IIBBBBB", width, height, 8, 2, 0, 0, 0))
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


def scaled_mask_image(image):
    """Render a 0/1-valued mask image as a full-range binary PNG.

    The generic percentile stretch produces a black image when >99 % of
    pixels are 1 (vmin == vmax == 1) or when the mask is essentially
    constant, so cheese mask mosaics were "black on black". Here we simply
    map mask > 0 -> 255 so kept regions are white and excluded regions are
    black.
    """
    import numpy as np

    arr = np.asarray(image)
    out = np.where(arr > 0, 255, 0).astype(np.uint8)
    return np.flipud(out)


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


def _save_mosaic_png(
    path: Path,
    image,
    scale: str,
    label: str,
) -> None:
    """Render a 2-D image to PNG with a side colorbar (matplotlib).

    scale is one of "log", "linear", "log_positive", or "linear_full".
    "linear_full" uses 0..max with no percentile clipping (preserves
    discrete step structure in flat per-frame mosaics).
    """
    import numpy as np
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm, Normalize
    from matplotlib.ticker import LogLocator, MaxNLocator

    arr = np.asarray(image, dtype=float)
    arr = np.nan_to_num(arr, nan=0.0, posinf=0.0, neginf=0.0)

    nonzero = arr[(arr != 0) & np.isfinite(arr)]
    is_log = False
    if scale in {"log", "log_positive"}:
        is_log = True
        positive = arr[arr > 0]
        if positive.size:
            vmin = float(np.percentile(positive, 10.0))
            vmax = float(np.percentile(positive, 90.0))
            if vmax <= vmin:
                vmax = vmin * 10.0 if vmin > 0 else 1.0
            norm = LogNorm(vmin=max(vmin, 1e-12), vmax=vmax, clip=True)
            display = np.where(arr > 0, arr, np.nan)
        else:
            norm = Normalize(vmin=0.0, vmax=1.0)
            display = arr
    elif scale == "linear_full":
        if nonzero.size:
            vmin = float(np.percentile(nonzero, 10.0))
            vmax = float(np.percentile(nonzero, 90.0))
        else:
            vmin, vmax = 0.0, 1.0
        if vmax <= vmin:
            vmax = vmin + 1.0
        norm = Normalize(vmin=vmin, vmax=vmax, clip=True)
        display = arr
    else:
        if nonzero.size:
            vmin = float(np.percentile(nonzero, 10.0))
            vmax = float(np.percentile(nonzero, 90.0))
            if vmax <= vmin:
                vmax = vmin + 1.0
        else:
            vmin, vmax = 0.0, 1.0
        norm = Normalize(vmin=vmin, vmax=vmax, clip=True)
        display = arr

    height, width = arr.shape
    aspect = width / max(height, 1)
    fig_h = 6.0
    fig_w = max(fig_h * aspect + 1.6, 4.0)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h), facecolor="black")
    ax.set_facecolor("black")
    cmap = plt.get_cmap("gray").copy()
    cmap.set_bad("black")
    im = ax.imshow(
        display, origin="lower", cmap=cmap, norm=norm, interpolation="nearest"
    )
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)
    cbar = fig.colorbar(im, ax=ax, fraction=0.04, pad=0.02)
    if is_log:
        cbar.ax.yaxis.set_major_locator(LogLocator(base=10.0, numticks=8))
        cbar.ax.yaxis.set_minor_locator(
            LogLocator(base=10.0, subs=tuple(np.arange(2, 10) * 0.1), numticks=80)
        )
    else:
        cbar.ax.yaxis.set_major_locator(MaxNLocator(nbins=6, steps=[1, 2, 5, 10]))
        cbar.ax.yaxis.set_minor_locator(MaxNLocator(nbins=24))
    cbar.update_ticks()
    cbar.ax.tick_params(which="major", colors="white", labelsize=8, length=4, width=0.8)
    cbar.ax.tick_params(which="minor", colors="white", length=2, width=0.5)
    cbar.outline.set_edgecolor("white")
    cbar.set_label(label, color="white", fontsize=9)
    fig.savefig(
        path,
        dpi=150,
        facecolor="black",
        edgecolor="black",
        bbox_inches="tight",
        pad_inches=0.1,
    )
    plt.close(fig)


def maps_qc(
    manifest_arg: str,
    outdir_arg: str,
    background_manifest_arg: str | None = None,
) -> None:
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

    # Optional flat-background lookup keyed by (inst, band, base) -> path.
    bg_paths: dict[tuple[str, str, str], Path] = {}
    if background_manifest_arg:
        bg_path = Path(background_manifest_arg)
        if bg_path.is_file():
            bg_lines = bg_path.read_text(encoding="utf-8").splitlines()
            if bg_lines:
                bg_header = bg_lines[0].split("\t")
                bg_idx = {name: i for i, name in enumerate(bg_header)}
                for line in bg_lines[1:]:
                    parts = line.split("\t")
                    if len(parts) < len(bg_header):
                        continue
                    key = (
                        parts[bg_idx["inst"]],
                        parts[bg_idx["band"]],
                        parts[bg_idx["base"]],
                    )
                    bg_paths[key] = Path(parts[bg_idx["background"]])

    products = [
        ("counts", "sum", "log", "counts/pixel (summed)"),
        ("exposure_vig", "sum", "log", "exposure_vig (s, summed)"),
        ("exposure_unvig", "sum", "log", "exposure_unvig (s, summed)"),
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
        for product, op, scale, label in products:
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
            _save_mosaic_png(png, acc, scale, f"{band} {label}")
            summary.append(f"{band}\t{product}\t{used}\t{op}\t{png}")

        # Exposure-corrected rate mosaic (counts / exposure_vig), no
        # background subtraction. Always emitted when both products exist.
        rate_counts = None
        rate_exp = None
        rate_used = 0
        for row in band_rows:
            counts_file = Path(row[index["counts"]])
            exp_file = Path(row[index["exposure_vig"]])
            if not counts_file.is_file() or not exp_file.is_file():
                continue
            counts_img = fits_image(counts_file)
            exp_img = fits_image(exp_file)
            if counts_img.shape != exp_img.shape:
                die(f"Counts/exposure shape mismatch for band {band}: {counts_file}")
            if rate_counts is None:
                rate_counts = np.zeros_like(counts_img, dtype=float)
                rate_exp = np.zeros_like(counts_img, dtype=float)
            rate_counts += counts_img
            rate_exp += exp_img
            rate_used += 1
        if rate_counts is not None:
            with np.errstate(divide="ignore", invalid="ignore"):
                rate = np.where(rate_exp > 0, rate_counts / rate_exp, 0.0)
            rate = np.nan_to_num(rate, nan=0.0, posinf=0.0, neginf=0.0)
            png = outdir / f"{band}_rate_mosaic.png"
            _save_mosaic_png(png, rate, "log_positive", f"{band} C/E (counts/s/pixel)")
            summary.append(f"{band}\trate\t{rate_used}\tC/E\t{png}")

        # Background + corrected mosaics if a background manifest is provided.
        if not bg_paths:
            continue
        counts_acc = None
        bg_acc = None
        exp_acc = None
        used_bg = 0
        for row in band_rows:
            key = (row[index["inst"]], band, row[index["base"]])
            bg_file = bg_paths.get(key)
            counts_file = Path(row[index["counts"]])
            exp_file = Path(row[index["exposure_vig"]])
            if (
                bg_file is None
                or not bg_file.is_file()
                or not counts_file.is_file()
                or not exp_file.is_file()
            ):
                continue
            counts_img = fits_image(counts_file)
            bg_img = fits_image(bg_file)
            exp_img = fits_image(exp_file)
            if counts_img.shape != bg_img.shape or counts_img.shape != exp_img.shape:
                die(f"Shape mismatch in background/corrected for band {band}: {key}")
            if counts_acc is None:
                counts_acc = np.zeros_like(counts_img, dtype=float)
                bg_acc = np.zeros_like(counts_img, dtype=float)
                exp_acc = np.zeros_like(counts_img, dtype=float)
            counts_acc += counts_img
            bg_acc += bg_img
            exp_acc += exp_img
            used_bg += 1
        if counts_acc is None:
            continue
        png = outdir / f"{band}_background_mosaic.png"
        _save_mosaic_png(
            png, bg_acc, "linear", f"{band} background (counts/pixel, summed)"
        )
        summary.append(f"{band}\tbackground\t{used_bg}\tsum\t{png}")
        with np.errstate(divide="ignore", invalid="ignore"):
            corrected = np.where(exp_acc > 0, (counts_acc - bg_acc) / exp_acc, 0.0)
        corrected = np.nan_to_num(corrected, nan=0.0, posinf=0.0, neginf=0.0)
        png = outdir / f"{band}_corrected_mosaic.png"
        _save_mosaic_png(
            png, corrected, "log_positive", f"{band} (C-B)/E (counts/s/pixel)"
        )
        summary.append(f"{band}\tcorrected\t{used_bg}\t(C-B)/E\t{png}")

    (outdir / "maps_qc_summary.tsv").write_text(
        "\n".join(summary) + "\n", encoding="utf-8"
    )


def maps_background(
    manifest_arg: str,
    outdir_arg: str,
    fraction_arg: str,
    outside_radius_fraction_arg: str,
) -> None:
    """Write flat constant-valued background maps for each maps manifest row.

    For each per-pointing image on the common maps grid, the constant is

        fraction * mean(counts[footprint & outside_circle])

    where:

    - ``footprint`` is the per-pointing exposure footprint (``exposure_vig
      > 0``). All non-footprint pixels are excluded from the mean and are
      written as zero in the output map.
    - ``outside_circle`` is the complement of a per-pointing source
      exclusion disk centred on the **pointing's footprint centroid**
      (i.e. the chip centre on the maps grid for this exposure, which is
      where the source sits for a tracked observation). Its radius is
      ``outside_radius_fraction * R_foot`` where ``R_foot`` is half the
      larger of the footprint's bounding-box dimensions for this
      pointing — so radius is measured relative to the chip FOV, not the
      common grid.

    The earlier implementation centred the exclusion circle at the
    *grid* centre with radius set by the *grid* half-dimension. For
    multi-pointing mosaics the per-pointing footprints lie hundreds of
    pixels off-grid-centre, which (a) caused the "background" sample to
    include the source for off-centre pointings and (b) chose a circle
    larger than a single pointing's FOV so the centred case produced an
    empty sample.
    """
    import numpy as np
    from astropy.io import fits

    try:
        fraction = float(fraction_arg)
    except (TypeError, ValueError):
        die(f"Invalid background fraction: {fraction_arg}")
    if not (fraction >= 0):
        die(f"Background fraction must be >= 0, got {fraction}")
    try:
        outside_radius_fraction = float(outside_radius_fraction_arg)
    except (TypeError, ValueError):
        die(f"Invalid outside-radius fraction: {outside_radius_fraction_arg}")
    if not (0.0 <= outside_radius_fraction):
        die(f"Outside-radius fraction must be >= 0, got {outside_radius_fraction}")

    manifest = Path(manifest_arg)
    if not manifest.is_file():
        die(f"Maps manifest not found: {manifest}")
    lines = manifest.read_text(encoding="utf-8").splitlines()
    if not lines:
        die(f"Empty maps manifest: {manifest}")
    header = lines[0].split("\t")
    index = {name: i for i, name in enumerate(header)}
    for required in ("inst", "band", "base", "counts", "exposure_vig"):
        if required not in index:
            die(f"Maps manifest missing column '{required}': {manifest}")

    outdir = Path(outdir_arg)
    outdir.mkdir(parents=True, exist_ok=True)
    out_manifest = outdir / "maps_background_manifest.tsv"
    out_lines = [
        "inst\tband\tbase\tbackground\tfraction\toutside_radius_fraction"
        "\tmean_counts\tconstant\tcenter_x\tcenter_y\tradius_pix\tn_sample"
    ]

    for line in lines[1:]:
        parts = line.split("\t")
        if len(parts) < len(header):
            continue
        inst = parts[index["inst"]]
        band = parts[index["band"]]
        base = parts[index["base"]]
        counts_path = Path(parts[index["counts"]])
        exp_path = Path(parts[index["exposure_vig"]])
        if not counts_path.is_file():
            die(f"Counts image missing: {counts_path}")

        counts_img = fits_image(counts_path)
        if not exp_path.is_file():
            die(f"Exposure image missing (needed for footprint): {exp_path}")
        exp_img = fits_image(exp_path)
        if exp_img.shape != counts_img.shape:
            die(f"Counts/exposure shape mismatch for {base}: {counts_path}")
        footprint = exp_img > 0

        # Per-pointing exclusion circle. Centre and radius are measured
        # from this pointing's exposure footprint — NOT the common grid.
        ny, nx = counts_img.shape
        if footprint.any():
            ys, xs = np.where(footprint)
            cy = float(ys.mean())
            cx = float(xs.mean())
            half_extent = 0.5 * max(ys.max() - ys.min(), xs.max() - xs.min())
            radius = outside_radius_fraction * half_extent
            yy, xx = np.ogrid[:ny, :nx]
            outside_circle = (yy - cy) ** 2 + (xx - cx) ** 2 >= radius * radius
            sample = footprint & outside_circle
            n_sample = int(sample.sum())
            if n_sample > 0:
                mean_counts = float(counts_img[sample].mean())
            else:
                # Circle swallowed the entire footprint; fall back to
                # the full footprint (rare, only if outside_radius_fraction
                # is set to 0 or a value > sqrt(2)).
                mean_counts = float(counts_img[footprint].mean())
                n_sample = int(footprint.sum())
        else:
            cy = (ny - 1) / 2.0
            cx = (nx - 1) / 2.0
            radius = 0.0
            mean_counts = 0.0
            n_sample = 0

        constant = fraction * mean_counts
        # Build a per-frame map: constant inside the exposure footprint,
        # zero outside (matches the FOV of every other map product).
        bg_image = np.zeros_like(counts_img, dtype=np.float32)
        bg_image[footprint] = constant

        # Preserve the WCS/header from the counts image so the background
        # aligns with downstream products.
        with fits.open(counts_path, memmap=True) as hdul:
            src_header = None
            for hdu in hdul:
                if getattr(hdu, "data", None) is not None and hdu.data.ndim >= 2:
                    src_header = hdu.header.copy()
                    break
        new_header = fits.Header()
        if src_header is not None:
            for key in src_header:
                if key in {"SIMPLE", "BITPIX", "NAXIS", "NAXIS1", "NAXIS2", "EXTEND"}:
                    continue
                try:
                    new_header[key] = src_header[key]
                except (ValueError, KeyError):
                    continue

        bg_dir = outdir / "background" / inst / base
        bg_dir.mkdir(parents=True, exist_ok=True)
        bg_path = bg_dir / f"{base}_{band}_background_flat.fits"
        hdu = fits.PrimaryHDU(data=bg_image, header=new_header)
        hdu.header["BG_FRAC"] = (fraction, "maps_background fraction of mean counts")
        hdu.header["BG_RFRAC"] = (
            outside_radius_fraction,
            "exclusion radius / per-pointing half-extent",
        )
        hdu.header["BG_CX"] = (cx, "exclusion circle centre, pixel x")
        hdu.header["BG_CY"] = (cy, "exclusion circle centre, pixel y")
        hdu.header["BG_RAD"] = (radius, "exclusion circle radius, pixels")
        hdu.header["BG_NSAMP"] = (n_sample, "n pixels in background sample")
        hdu.header["BG_MEAN"] = (
            mean_counts,
            "mean counts/pixel in background sample",
        )
        hdu.header["BG_CONST"] = (constant, "flat background value (counts/pixel)")
        hdu.header["BG_SRC"] = (str(counts_path), "source counts image")
        hdu.writeto(bg_path, overwrite=True)

        out_lines.append(
            f"{inst}\t{band}\t{base}\t{bg_path.resolve()}"
            f"\t{fraction:.6g}\t{outside_radius_fraction:.6g}"
            f"\t{mean_counts:.6g}\t{constant:.6g}"
            f"\t{cx:.2f}\t{cy:.2f}\t{radius:.2f}\t{n_sample}"
        )

    out_manifest.write_text("\n".join(out_lines) + "\n", encoding="utf-8")
    print(f"Wrote {out_manifest}")


def cheese_qc(
    manifest_arg: str, outdir_arg: str, maps_manifest_arg: str | None = None
) -> None:
    import numpy as np
    from astropy.io import fits

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

    # Optional lookup: (inst, band, base) -> (unmasked maps counts path,
    # unmasked maps exposure_vig path). The counts QC mosaic prefers the
    # original counts so the red overlay is painted on top of the un-excised
    # image, and the exposure_vig is used to gate the red overlay so that
    # pixels OUTSIDE a frame's FOV are never marked as "excluded" by that
    # frame.
    maps_counts: dict[tuple[str, str, str], Path] = {}
    maps_exp_vig: dict[tuple[str, str, str], Path] = {}
    if maps_manifest_arg:
        mm_path = Path(maps_manifest_arg)
        if mm_path.is_file():
            mm_lines = mm_path.read_text(encoding="utf-8").splitlines()
            if mm_lines:
                mm_header = mm_lines[0].split("\t")
                mm_idx = {name: idx for idx, name in enumerate(mm_header)}
                for line in mm_lines[1:]:
                    parts = line.split("\t")
                    if len(parts) < len(mm_header):
                        continue
                    key = (
                        parts[mm_idx["inst"]],
                        parts[mm_idx["band"]],
                        parts[mm_idx["base"]],
                    )
                    maps_counts[key] = Path(parts[mm_idx["counts"]])
                    if "exposure_vig" in mm_idx:
                        maps_exp_vig[key] = Path(parts[mm_idx["exposure_vig"]])

    outdir = Path(outdir_arg)
    outdir.mkdir(parents=True, exist_ok=True)
    summary = ["band\tproduct\tfiles\toperation\tpng"]

    bands = []
    for row in rows:
        band = row[index["band"]]
        if band not in bands:
            bands.append(band)

    # Simple mosaic helper.
    def mosaic(paths: list[Path], op: str):
        acc = None
        used = 0
        for p in paths:
            if not p.is_file():
                continue
            image = fits_image(p)
            if acc is None:
                acc = np.zeros_like(image, dtype=float)
            if image.shape != acc.shape:
                die(f"Shape mismatch in mosaic input: {p}")
            if op == "max":
                acc = np.maximum(acc, image)
            else:
                acc += image
            used += 1
        return acc, used

    for band in bands:
        band_rows = [row for row in rows if row[index["band"]] == band]

        # cheese_mask mosaic (binary white/black on black background).
        mask_paths = [Path(r[index["cheese_mask"]]) for r in band_rows]
        mask_acc, mask_used = mosaic(mask_paths, "max")
        if mask_acc is not None:
            png = outdir / f"{band}_cheese_mask_mosaic.png"
            write_gray_png(png, scaled_mask_image(mask_acc))
            summary.append(f"{band}\tcheese_mask\t{mask_used}\tmax\t{png}")

        # counts mosaic with red overlay on excluded pixels. We mosaic the
        # original (pre-mask) maps counts when available, and paint semi-
        # transparent red wherever the cheese mask excludes pixels in any
        # contributing frame (the "holes" the cheese stage cut). The
        # exclusion is gated by each frame's original exposure_vig > 0 so
        # pixels outside a frame's FOV are NOT counted as excluded by that
        # frame (without this gating, a frame would "exclude" every pixel
        # it never observed, swamping the overlay with red everywhere).
        counts_paths: list[Path] = []
        excl_specs: list[tuple[Path, Path | None]] = []
        for row in band_rows:
            key = (row[index["inst"]], band, row[index["base"]])
            orig = maps_counts.get(key)
            if orig is not None and orig.is_file():
                counts_paths.append(orig)
            else:
                # Fallback: use cheesed_counts, which still sums correctly but
                # shows the excision as zeros under the red overlay.
                counts_paths.append(Path(row[index["cheesed_counts"]]))
            excl_specs.append((Path(row[index["cheese_mask"]]), maps_exp_vig.get(key)))

        counts_acc, counts_used = mosaic(counts_paths, "sum")
        if counts_acc is not None:
            # For the overlay, count how many contributing frames exclude
            # each pixel (mask==0) AND actually cover that pixel
            # (exposure_vig > 0). Using the unmasked maps exposure_vig
            # avoids painting red over pixels outside the frame's FOV.
            excl_acc = None
            for mask_path, fov_path in excl_specs:
                if not mask_path.is_file():
                    continue
                m = fits_image(mask_path)
                if m.shape != counts_acc.shape:
                    die(f"Mask/counts shape mismatch for band {band}: {mask_path}")
                if fov_path is not None and fov_path.is_file():
                    fov = fits_image(fov_path)
                    if fov.shape != counts_acc.shape:
                        die(
                            f"Exp_vig/counts shape mismatch for band {band}: {fov_path}"
                        )
                    fov_ok = np.isfinite(fov) & (fov > 0)
                else:
                    # No FOV reference -> fall back to "every pixel in image".
                    fov_ok = np.ones_like(m, dtype=bool)
                ex = (fov_ok & (m <= 0)).astype(float)
                excl_acc = ex if excl_acc is None else excl_acc + ex

            gray = scaled_log_image(counts_acc)  # already flipped vertically
            rgb = np.stack([gray, gray, gray], axis=-1).astype(float)
            if excl_acc is not None:
                excl_norm = np.clip(excl_acc / max(excl_acc.max(), 1.0), 0.0, 1.0)
                excl_norm = np.flipud(excl_norm)
                alpha = 0.5 * excl_norm  # 0 where no frame excludes, 0.5 peak
                red = np.zeros_like(rgb)
                red[..., 0] = 255.0
                alpha3 = alpha[..., None]
                rgb = (1.0 - alpha3) * rgb + alpha3 * red
            rgb_u8 = np.clip(rgb, 0.0, 255.0).astype(np.uint8)
            png = outdir / f"{band}_cheesed_counts_mosaic.png"
            write_rgb_png(png, rgb_u8)
            summary.append(f"{band}\tcheesed_counts\t{counts_used}\tsum+overlay\t{png}")

        # exposure mosaic: unchanged log-scaled grayscale.
        exp_paths = [Path(r[index["cheesed_exposure_vig"]]) for r in band_rows]
        exp_acc, exp_used = mosaic(exp_paths, "sum")
        if exp_acc is not None:
            png = outdir / f"{band}_cheesed_exposure_vig_mosaic.png"
            write_gray_png(png, scaled_log_image(exp_acc))
            summary.append(f"{band}\tcheesed_exposure_vig\t{exp_used}\tsum\t{png}")

    (outdir / "cheese_qc_summary.tsv").write_text(
        "\n".join(summary) + "\n", encoding="utf-8"
    )

    # Per-frame diagnostic TSV: for each (inst, band, base) row, count
    # detected sources (from the SRCLIST or eboxlist), the FOV pixel area
    # (cheese_mask == 1 + cheese_mask == 0 inside FOV), and the fraction
    # of FOV pixels excluded by the cheese mask. This is the right level
    # of granularity to spot frames where detection failed (n_sources=0)
    # or where the mask excised an unreasonable fraction of the FOV.
    per_frame = [
        "inst\tband\tbase\tn_sources\tfov_pixels\texcluded_pixels\texcluded_frac"
    ]
    for row in rows:
        inst = row[index["inst"]]
        band = row[index["band"]]
        base = row[index["base"]]
        src_path = (
            Path(row[index["source_list"]])
            if index.get("source_list") is not None
            else None
        )
        mask_path = Path(row[index["cheese_mask"]])
        key = (inst, band, base)

        n_src = 0
        if src_path is not None and src_path.is_file():
            try:
                with fits.open(src_path, memmap=False) as hdul:
                    for hdu in hdul:
                        if isinstance(hdu, (fits.BinTableHDU, fits.TableHDU)):
                            if str(hdu.name).upper() in {
                                "SRCLIST",
                                "EBOX_LIST",
                                "BOXLIST",
                            }:
                                n_src = int(len(hdu.data))
                                break
                    else:
                        for hdu in hdul:
                            if isinstance(hdu, (fits.BinTableHDU, fits.TableHDU)):
                                n_src = int(len(hdu.data))
                                break
            except Exception:
                n_src = -1

        n_fov = 0
        n_excl = 0
        if mask_path.is_file():
            mask = fits_image(mask_path)
            fov_ref = maps_exp_vig.get(key)
            if fov_ref is not None and fov_ref.is_file():
                fov = fits_image(fov_ref)
                fov_ok = np.isfinite(fov) & (fov > 0)
            else:
                fov_ok = np.ones_like(mask, dtype=bool)
            n_fov = int(fov_ok.sum())
            n_excl = int((fov_ok & (mask <= 0)).sum())
        frac = (n_excl / n_fov) if n_fov > 0 else 0.0
        per_frame.append(
            f"{inst}\t{band}\t{base}\t{n_src}\t{n_fov}\t{n_excl}\t{frac:.4f}"
        )

    (outdir / "cheese_per_frame.tsv").write_text(
        "\n".join(per_frame) + "\n", encoding="utf-8"
    )

    # Source-overlay PNG: for each band, draw the detected emldetect
    # sources (per (inst, base)) as cyan circles on the cheesed counts
    # mosaic. Useful for spot-checking that the source detection is
    # finding what the eye sees, and that the mask radius matches the
    # PSF.
    for band in bands:
        band_rows = [row for row in rows if row[index["band"]] == band]
        counts_paths: list[Path] = []
        srclist_specs: list[tuple[Path, tuple[str, str, str]]] = []
        for row in band_rows:
            key = (row[index["inst"]], band, row[index["base"]])
            orig = maps_counts.get(key)
            counts_paths.append(
                orig
                if (orig and orig.is_file())
                else Path(row[index["cheesed_counts"]])
            )
            srclist_specs.append((Path(row[index["source_list"]]), key))

        counts_acc, _ = mosaic(counts_paths, "sum")
        if counts_acc is None:
            continue
        gray = scaled_log_image(counts_acc)  # already vertically flipped
        rgb = np.stack([gray, gray, gray], axis=-1).astype(float)
        ny, nx = counts_acc.shape

        for srclist_path, _key in srclist_specs:
            if not srclist_path.is_file():
                continue
            try:
                with fits.open(srclist_path, memmap=False) as hdul:
                    table = None
                    for hdu in hdul:
                        if (
                            isinstance(hdu, (fits.BinTableHDU, fits.TableHDU))
                            and str(hdu.name).upper() == "SRCLIST"
                        ):
                            table = hdu
                            break
                    if table is None:
                        for hdu in hdul:
                            if isinstance(hdu, (fits.BinTableHDU, fits.TableHDU)):
                                table = hdu
                                break
                    if table is None:
                        continue
                    cols = {c.upper(): c for c in table.data.names}
                    if "X_IMA" not in cols or "Y_IMA" not in cols:
                        continue
                    xs = np.asarray(table.data[cols["X_IMA"]], dtype=float) - 1.0
                    ys = np.asarray(table.data[cols["Y_IMA"]], dtype=float) - 1.0
                    if "ID_BAND" in cols:
                        keep = np.asarray(table.data[cols["ID_BAND"]], dtype=int) == 1
                        xs, ys = xs[keep], ys[keep]
                    if xs.size == 0:
                        continue
                    # Match the vertical flip in scaled_log_image.
                    ys_p = (ny - 1) - ys
                    radius = 6.0  # pixels
                    yy, xx = np.indices((ny, nx))
                    for cx, cy in zip(xs, ys_p):
                        r = np.sqrt((xx - cx) ** 2 + (yy - cy) ** 2)
                        ring = (r >= radius - 0.7) & (r <= radius + 0.7)
                        rgb[ring, 0] = 0.0
                        rgb[ring, 1] = 255.0
                        rgb[ring, 2] = 255.0
            except Exception:
                continue

        rgb_u8 = np.clip(rgb, 0.0, 255.0).astype(np.uint8)
        png = outdir / f"{band}_cheese_sources_overlay.png"
        write_rgb_png(png, rgb_u8)
        summary.append(f"{band}\tsources_overlay\t{len(band_rows)}\toverlay\t{png}")

    # Rewrite the summary TSV to include the overlay rows.
    (outdir / "cheese_qc_summary.tsv").write_text(
        "\n".join(summary) + "\n", encoding="utf-8"
    )


def event_mosaic(
    manifest_dir: str,
    outdir: str,
    detectors_arg: str = "PN",
    split_ev_arg: str = "1000",
) -> None:
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
    elif cmd == "cheese-clean-labels" and len(args) == 3:
        cheese_clean_labels(args[0], int(args[1]), int(args[2]))
    elif cmd == "maps-expr" and len(args) == 3:
        maps_expression(*args)
    elif cmd == "maps-grid" and len(args) == 4:
        maps_grid(*args)
    elif cmd == "maps-grid" and len(args) == 5:
        maps_grid(args[0], args[1], args[2], args[3], args[4])
    elif cmd == "maps-qc" and len(args) in {2, 3}:
        maps_qc(args[0], args[1], args[2] if len(args) == 3 else None)
    elif cmd == "maps-background" and len(args) == 4:
        maps_background(*args)
    elif cmd == "cheese-qc" and len(args) in {2, 3}:
        cheese_qc(args[0], args[1], args[2] if len(args) == 3 else None)
    elif cmd == "fits-rows" and len(args) == 1:
        fits_rows(*args)
    elif cmd == "fits-table-hdu" and len(args) in {1, 2}:
        fits_table_hdu(args[0], args[1] if len(args) == 2 else "SRCLIST")
    elif cmd == "fits-table-has-column" and len(args) == 3:
        fits_table_has_column(*args)
    elif cmd == "apply-image-mask" and len(args) in {3, 4}:
        apply_image_mask(
            args[0], args[1], args[2], args[3] if len(args) == 4 else "image"
        )
    elif cmd == "combine-masks" and len(args) == 3:
        combine_masks(*args)
    elif cmd == "fov-mask" and len(args) == 2:
        fov_mask(*args)
    elif cmd == "mask-subtract" and len(args) == 3:
        mask_subtract(*args)
    elif cmd == "cheese-mask" and len(args) == 8:
        cheese_make_mask(*args)
    elif cmd == "event-mosaic" and len(args) in {2, 3, 4}:
        event_mosaic(
            args[0],
            args[1],
            args[2] if len(args) >= 3 else "PN",
            args[3] if len(args) == 4 else "1000",
        )
    else:
        die(
            "Usage: tools.py shell CONFIG | rewrite-path SUM.SAS PATH | "
            "clean-band-labels CONFIG | clean-band-table CONFIG | clean-expr CONFIG DETECTOR LABEL | "
            "maps-band-table CONFIG | maps-clean-labels CONFIG MAP_LABEL | maps-expr CONFIG DETECTOR LABEL | "
            "maps-grid OUTDIR BIN_PHYS PAD_FRAC MANIFEST_DIR [DETECTORS] | "
            "maps-qc MAPS_MANIFEST OUTDIR [BACKGROUND_MANIFEST] | "
            "maps-background MAPS_MANIFEST OUTDIR FRACTION OUTSIDE_RADIUS_FRACTION | "
            "cheese-qc CHEESE_MANIFEST OUTDIR [MAPS_MANIFEST] | "
            "fits-rows FITS | fits-table-hdu FITS [PREFERRED] | fits-table-has-column FITS TABLE COLUMN | "
            "apply-image-mask IMAGE MASK OUT [image|mask] | "
            "combine-masks DETMASK REGMASK OUT | "
            "fov-mask EXPOSURE OUT | "
            "cheese-mask EMLLIST FOV OUT ID_INST ID_BAND MLMIN FLUX_E14 RADIUS_PIX | "
            "event-mosaic MANIFEST_DIR OUTDIR [DETECTORS] [SPLIT_EV]"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
