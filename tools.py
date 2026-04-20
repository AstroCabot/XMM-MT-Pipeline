#!/usr/bin/env python3
"""Small helpers for reduction_v6/pipeline.sh."""

from __future__ import annotations

import json
import re
import shlex
import sys
from pathlib import Path
from typing import Any


def die(message: str) -> None:
    raise SystemExit(message)


def load_config(path: str) -> dict[str, Any]:
    try:
        with Path(path).open("r", encoding="utf-8") as handle:
            cfg = json.load(handle)
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
        "PN": ["PN"], "EPN": ["PN"],
        "M1": ["M1"], "MOS1": ["M1"], "EMOS1": ["M1"],
        "M2": ["M2"], "MOS2": ["M2"], "EMOS2": ["M2"],
        "ALL": ["PN", "M1", "M2"], "EPIC": ["PN", "M1", "M2"],
    }
    detectors: list[str] = []
    for token in tokens:
        if token.upper() not in aliases:
            die(f"Unknown detector in config: {token}")
        for detector in aliases[token.upper()]:
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
        "CLEAN_PI_MIN": cfg.get("clean_pi_min", 200),
        "CLEAN_PI_MAX": cfg.get("clean_pi_max", 12000),
        "CLEAN_GTI_ENABLED": yesno(cfg.get("clean_gti_enabled"), True),
        "CLEAN_GTI_RATE_CUT": cfg.get("clean_gti_rate_cut", 4.80001211),
        "CLEAN_GTI_TIMEBIN": cfg.get("clean_gti_timebin", 10),
        "CLEAN_GTI_PI_MIN": cfg.get("clean_gti_pi_min", 7000),
        "CLEAN_GTI_PI_MAX": cfg.get("clean_gti_pi_max", 15000),
        "CLEAN_GTI_PATTERN_MAX": cfg.get("clean_gti_pattern_max", 0),
        "PREMERGE_IMAGE_SIZE_DEG": cfg.get("premerge_image_size_deg", 1.5),
        "IMAGE_BANDS": cfg.get("image_bands", "1000:200:500;2000:500:1000;3000:1000:2000;4000:2000:4500;5000:4500:12000"),
        "IMAGE_BIN_PHYS": cfg.get("image_bin_phys", 80),
        "IMAGE_PAD_FRAC": cfg.get("image_pad_frac", 0.08),
        "EEXPMAP_ATTREBIN": cfg.get("eexpmap_attrebin", "0.020626481"),
        "BACKGROUND_INNER_RADIUS_FRAC": cfg.get("background_inner_radius_frac", cfg.get("background_exclude_radius_frac", 0.50)),
        "BACKGROUND_OUTER_RADIUS_FRAC": cfg.get("background_outer_radius_frac", 1.00),
        "BACKGROUND_MIN_PIXELS": cfg.get("background_min_pixels", 100),
        "SOURCE_DETECT_BAND": cfg.get("source_detect_band", "4000"),
        "SOURCE_MLMIN": cfg.get("source_mlmin", 6),
        "SOURCE_LIKEMIN": cfg.get("source_likemin", 5),
        "SOURCE_BOXSIZE": cfg.get("source_boxsize", 5),
        "SOURCE_NRUNS": cfg.get("source_nruns", 3),
        "SOURCE_EMASK_THRESHOLD1": cfg.get("source_emask_threshold1", 0.30),
        "SOURCE_EMASK_THRESHOLD2": cfg.get("source_emask_threshold2", 0.50),
        "SOURCE_MASK_RADIUS_PIX": cfg.get("source_mask_radius_pix", 5),
        "SOURCE_MASK_LIKE_SCALE": cfg.get("source_mask_like_scale", 0.0),
        "SOURCE_MASK_MAX_RADIUS_PIX": cfg.get("source_mask_max_radius_pix", 15),
        "SOURCE_HOTPIXEL_FILTER": yesno(cfg.get("source_hotpixel_filter"), True),
        "STACK_TARGET_ID": cfg.get("stack_target_id", "C/2025 N1"),
        "STACK_ATTMOVE_CENTER": cfg.get("stack_attmove_center", cfg.get("stack_observer_center", "500@399")),
        "STACK_OVERLAY_CENTER": cfg.get("stack_overlay_center", "@-125989"),
        "STACK_EPHEMERIS_FILE": cfg.get("stack_ephemeris_file", ""),
        "STACK_REF_RA": cfg.get("stack_ref_ra", ""),
        "STACK_REF_DEC": cfg.get("stack_ref_dec", ""),
        "STACK_TRACK_STEP_S": cfg.get("stack_track_step_s", 30),
        "STACK_ATTMOVE_GRANULARITY": cfg.get("stack_attmove_granularity", 5),
        "STACK_ATTMOVE_MINSTABLE": cfg.get("stack_attmove_minstable", 30),
        "STACK_ATTCALC_IMAGE_SIZE_DEG": cfg.get("stack_attcalc_image_size_deg", 1.5),
        "STACK_EEXPMAP_ATTREBIN": cfg.get("stack_eexpmap_attrebin", cfg.get("eexpmap_attrebin", 4.0)),
        "LINK_ODF_CONSTITUENTS": yesno(cfg.get("link_odf_constituents"), True),
        "SKIP_ODF_LINK_PATTERNS": "|".join(cfg.get("skip_odf_link_patterns", [])),
    }
    for name, value in values.items():
        emit(name, value)


def rewrite_path(summary: str, odf_path: str) -> None:
    path = odf_path.rstrip("/") + "/"
    summary_path = Path(summary)
    lines = summary_path.read_text(encoding="utf-8", errors="surrogateescape").splitlines()
    for idx, line in enumerate(lines):
        if line.startswith("PATH "):
            lines[idx] = f"PATH {path}"
            break
    else:
        lines.append(f"PATH {path}")
    summary_path.write_text("\n".join(lines) + "\n", encoding="utf-8", errors="surrogateescape")


def relative_path(root_arg: str, path_arg: str) -> None:
    root = Path(root_arg).resolve()
    path = Path(path_arg).resolve()
    try:
        print(path.relative_to(root))
        return
    except ValueError:
        pass
    matches = sorted(root.rglob(path.name))
    if not matches:
        die(f"{path_arg} is not inside {root_arg}")
    print(matches[0].resolve().relative_to(root))


def manifest_events(manifest_arg: str, detectors_arg: str = "PN") -> list[Path]:
    paths: list[Path] = []
    seen: set[Path] = set()
    source = Path(manifest_arg)
    if source.is_file():
        manifests = [source]
    else:
        manifests = []
        for detector in normalize_detectors(detectors_arg):
            for name in (f"{detector}_clean_files.txt", f"{detector}_raw.txt"):
                for parent in (source, source / "manifest"):
                    path = parent / name
                    if path.is_file():
                        manifests.append(path)
                        break
    if not manifests:
        die(f"Missing event manifest(s): {manifest_arg}")
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


def event_hdu(hdul):
    return hdul["EVENTS"] if "EVENTS" in hdul else hdul[1]


def event_columns(path: Path):
    try:
        from astropy.io import fits
    except ImportError:
        die("FITS helpers require astropy")

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

    # Reprocessed XMM event lists can contain large placeholder SKY values.
    return np.isfinite(x) & np.isfinite(y) & (np.abs(x) < 1.0e7) & (np.abs(y) < 1.0e7)


def event_bounds(events: list[Path], extra_extent: tuple[float, float, float, float] | None = None):
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
    if extra_extent:
        xmin = min(xmin, extra_extent[0])
        xmax = max(xmax, extra_extent[1])
        ymin = min(ymin, extra_extent[2])
        ymax = max(ymax, extra_extent[3])
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
    try:
        from astropy.io import fits
    except ImportError:
        die("fits-rows requires astropy")

    path = Path(path_arg)
    if not path.is_file():
        print(-1)
        return
    with fits.open(path, memmap=True) as hdul:
        hdu = hdul["EVENTS"] if "EVENTS" in hdul else hdul[1]
        data = hdu.data
        print(0 if data is None else len(data))


def flare_qc(gti_dir_arg: str, outdir_arg: str, rate_cut_arg: str) -> None:
    import matplotlib
    import numpy as np
    from astropy.io import fits

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    paths = sorted(Path(gti_dir_arg).glob("*/*flare_lc.fits"))
    if not paths:
        return

    cut = float(rate_cut_arg)
    outdir = Path(outdir_arg)
    outdir.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(len(paths), 1, figsize=(9, max(2.4, 1.8 * len(paths))), squeeze=False, constrained_layout=True)
    rows = ["lightcurve\tbins\taccepted_bins\trejected_bins\tmedian_rate_ct_s\tp95_rate_ct_s\n"]

    for ax, path in zip(axes[:, 0], paths):
        with fits.open(path, memmap=True) as hdul:
            data = hdul["RATE"].data if "RATE" in hdul else hdul[1].data
            time = np.asarray(data["TIME"], dtype=float)
            rate = np.asarray(data["RATE"], dtype=float)

        finite = np.isfinite(time) & np.isfinite(rate)
        if not finite.any():
            continue
        x = (time - time[finite].min()) / 1000.0
        good = finite & (rate <= cut)
        bad = finite & ~good

        ax.scatter(x[good], rate[good], s=8, color="#00a6d6", label="kept")
        ax.scatter(x[bad], rate[bad], s=8, color="#d62728", label="cut")
        ax.axhline(cut, color="white", linewidth=1.0, linestyle="--")
        ax.set_facecolor("black")
        ax.tick_params(colors="white", labelsize=8)
        ax.set_ylabel("ct/s", color="white", fontsize=8)
        ax.set_title(f"{path.parent.name}/{path.name}", color="white", fontsize=8)
        rows.append(
            f"{path}\t{int(finite.sum())}\t{int(good.sum())}\t{int(bad.sum())}\t"
            f"{float(np.nanmedian(rate[finite])):.6g}\t{float(np.nanpercentile(rate[finite], 95)):.6g}\n"
        )

    axes[-1, 0].set_xlabel("ks from light-curve start", color="white", fontsize=8)
    axes[0, 0].legend(loc="upper right", fontsize=7, facecolor="black", labelcolor="white")
    fig.savefig(outdir / "flare_lightcurves.png", dpi=180, facecolor="black")
    plt.close(fig)
    (outdir / "flare_lightcurves.tsv").write_text("".join(rows), encoding="utf-8")


def exposure_grid(outdir_arg: str, bin_phys_arg: str, pad_frac_arg: str, event_args: list[str]) -> None:
    import math
    import numpy as np

    events = [Path(arg) for arg in event_args if Path(arg).is_file()]
    if not events:
        die("No clean event files found for exposure grid")

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
    if not np.isfinite([xmin, xmax, ymin, ymax]).all():
        die("Could not determine exposure grid bounds")

    bin_phys = float(bin_phys_arg)
    pad_frac = float(pad_frac_arg)
    pad_x = max(pad_frac * (xmax - xmin), bin_phys)
    pad_y = max(pad_frac * (ymax - ymin), bin_phys)
    xmin = math.floor((xmin - pad_x) / bin_phys) * bin_phys
    xmax = math.ceil((xmax + pad_x) / bin_phys) * bin_phys
    ymin = math.floor((ymin - pad_y) / bin_phys) * bin_phys
    ymax = math.ceil((ymax + pad_y) / bin_phys) * bin_phys
    nx = int(round((xmax - xmin) / bin_phys))
    ny = int(round((ymax - ymin) / bin_phys))

    outdir = Path(outdir_arg)
    outdir.mkdir(parents=True, exist_ok=True)
    env = (
        f"X_MIN_PHYS={xmin}\nX_MAX_PHYS={xmax}\n"
        f"Y_MIN_PHYS={ymin}\nY_MAX_PHYS={ymax}\n"
        f"BIN_PHYS={bin_phys}\nNX={nx}\nNY={ny}\n"
    )
    (outdir / "grid.env").write_text(env, encoding="utf-8")
    (outdir / "grid.json").write_text(
        json.dumps(
            {"x_min": xmin, "x_max": xmax, "y_min": ymin, "y_max": ymax, "bin_phys": bin_phys, "nx": nx, "ny": ny},
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )


def event_time_info(path: Path) -> tuple[float, float, str]:
    from astropy.io import fits

    with fits.open(path, memmap=True) as hdul:
        header = (hdul["EVENTS"] if "EVENTS" in hdul else hdul[1]).header
        start = header.get("TSTART")
        stop = header.get("TSTOP")
        source = header.get("EXPIDSTR") or path.stem.split("_")[-2]
    if start is None or stop is None:
        die(f"Missing TSTART/TSTOP in {path}")
    return float(start), float(stop), str(source)


def clean_groups(clean_dir_arg: str) -> dict[str, list[tuple[float, float, str, Path]]]:
    clean_dir = Path(clean_dir_arg)
    groups: dict[str, list[tuple[float, float, str, Path]]] = {}
    manifests_by_inst: dict[str, Path] = {}
    for manifest in sorted(clean_dir.glob("*_clean_files.txt")):
        inst = manifest.name.removesuffix("_clean_files.txt")
        manifests_by_inst[inst] = manifest
    for manifest in sorted((clean_dir / "manifest").glob("*_clean_files.txt")):
        inst = manifest.name.removesuffix("_clean_files.txt")
        manifests_by_inst.setdefault(inst, manifest)

    for inst, manifest in sorted(manifests_by_inst.items()):
        rows = []
        for line in manifest.read_text(encoding="utf-8").splitlines():
            path = Path(line.strip())
            if path.is_file():
                start, stop, source = event_time_info(path)
                rows.append((start, stop, source, path))
        if rows:
            groups[inst] = sorted(rows)
    return groups


def selected_clean(clean_dir_arg: str, detectors_arg: str):
    groups = clean_groups(clean_dir_arg)
    if not groups:
        die(f"No clean event manifests found in {clean_dir_arg}")
    detectors = normalize_detectors(detectors_arg)
    selected = {inst: groups[inst] for inst in detectors if groups.get(inst)}
    if not selected:
        die(f"No selected clean event manifests found in {clean_dir_arg}")
    return detectors, selected


def clean_span(selected: dict[str, list[tuple[float, float, str, Path]]]) -> tuple[float, float]:
    return (
        min(row[0] for rows in selected.values() for row in rows),
        max(row[1] for rows in selected.values() for row in rows),
    )


def odf_pointings(init_dir_arg: str, span: tuple[float, float]) -> list[tuple[str, float, float, str, str, str]]:
    from astropy.io import fits
    from astropy.time import Time

    init_dir = Path(init_dir_arg)
    epoch = Time("1998-01-01T00:00:00", scale="tt")

    def timestamp(value: Any) -> float | None:
        text = str(value).strip()
        if not text:
            return None
        return float((Time(text, scale="utc").tt - epoch).sec)

    grouped: dict[str, dict[str, Any]] = {}
    for path in sorted(init_dir.glob("*SC*ATS.FIT")):
        with fits.open(path, memmap=True) as hdul:
            data = hdul[1].data
            names = set(hdul[1].columns.names)
            required = {"VALTIME", "VALDUR", "PREQID", "TYPEID", "PTTIME"}
            if not required <= names:
                die(f"{path} is not a usable spacecraft attitude history file")
            for row in data:
                typ = str(row["TYPEID"]).strip()
                if typ not in {"P", "S"}:
                    continue
                start = timestamp(row["VALTIME"])
                stop = None if start is None else start + max(float(row["VALDUR"]), 0.0)
                if start is None or stop is None or stop <= start or stop <= span[0] or start >= span[1]:
                    continue
                preqid = str(row["PREQID"]).strip()
                if not preqid:
                    continue
                item = grouped.setdefault(preqid, {"start": start, "stop": stop, "types": set(), "pttime": ""})
                item["start"] = min(item["start"], start)
                item["stop"] = max(item["stop"], stop)
                item["types"].add(typ)
                item["pttime"] = item["pttime"] or str(row["PTTIME"]).strip()

    pointings: list[tuple[float, float, str, str, str]] = []
    for preqid, item in grouped.items():
        start = max(float(item["start"]), span[0])
        stop = min(float(item["stop"]), span[1])
        if stop - start > 1.0:
            pointings.append((start, stop, preqid, "".join(sorted(item["types"])), str(item["pttime"])))
    if not pointings:
        die(f"No ODF attitude pointings overlap selected clean events in {init_dir}")
    pointings.sort()
    return [
        (f"P{i:03d}", start, stop, preqid, types, pttime)
        for i, (start, stop, preqid, types, pttime) in enumerate(pointings, 1)
    ]


def slice_records(selected: dict[str, list[tuple[float, float, str, Path]]], pointings):
    for inst, rows in selected.items():
        for event_start, event_stop, _event_source, path in rows:
            stem = path.name.removesuffix(".fits")
            for pid, point_start, point_stop, request, _types, _pttime in pointings:
                start, stop = max(event_start, point_start), min(event_stop, point_stop)
                if stop - start > 1.0:
                    yield inst, path, f"{stem}_{pid}_{request}", pid, start, stop, request


def exposure_slices(slices_arg: str, pointings_arg: str, clean_dir_arg: str, detectors_arg: str, init_dir_arg: str) -> None:
    _detectors, selected = selected_clean(clean_dir_arg, detectors_arg)
    pointings = odf_pointings(init_dir_arg, clean_span(selected))

    Path(pointings_arg).write_text(
        "pointing\tstart\tstop\tsource\trequest\ttypes\tpttime_utc\n"
        + "".join(
            f"{pid}\t{start:.6f}\t{stop:.6f}\tSCATS\t{request}\t{types}\t{pttime}\n"
            for pid, start, stop, request, types, pttime in pointings
        ),
        encoding="utf-8",
    )

    lines = ["inst\tevent\tbase\tpointing\tstart\tstop\treference_exposure\n"] + [
        f"{inst}\t{path}\t{base}\t{pid}\t{start:.6f}\t{stop:.6f}\t{request}\n"
        for inst, path, base, pid, start, stop, request in slice_records(selected, pointings)
    ]
    Path(slices_arg).write_text("".join(lines), encoding="utf-8")


def parse_bands(value: str) -> list[tuple[str, float, float]]:
    try:
        return [(label, float(lo), float(hi)) for label, lo, hi in (entry.split(":", 2) for entry in str(value).split(";") if entry.strip())]
    except ValueError:
        die(f"Bad image_bands value: {value}")


def slice_qc(
    outdir_arg: str,
    clean_dir_arg: str,
    detectors_arg: str,
    init_dir_arg: str,
    bands_arg: str,
) -> None:
    from astropy.io import fits

    outdir = Path(outdir_arg)
    outdir.mkdir(parents=True, exist_ok=True)
    detectors, selected = selected_clean(clean_dir_arg, detectors_arg)
    pointings = odf_pointings(init_dir_arg, clean_span(selected))
    bands = parse_bands(bands_arg)
    band_columns = "\t".join(f"{label}_events" for label, _lo, _hi in bands)

    (outdir / "pre_exposure_pointings.tsv").write_text(
        "pointing\tstart\tstop\tduration_s\tsource\trequest\ttypes\tpttime_utc\n"
        + "".join(
            f"{pid}\t{start:.6f}\t{stop:.6f}\t{stop - start:.3f}\tSCATS\t{request}\t{types}\t{pttime}\n"
            for pid, start, stop, request, types, pttime in pointings
        ),
        encoding="utf-8",
    )

    slice_lines = [f"inst\tevent\tpointing\trequest\tstart\tstop\tduration_s\ttotal_events\t{band_columns}\n"]
    coverage_lines = ["inst\tevent\tspan_s\tassigned_s\tunassigned_s\n"]
    by_pointing: dict[str, dict[str, float]] = {
        pid: {
            "duration": stop - start,
            "assigned": 0.0,
            "total": 0.0,
            **{label: 0.0 for label, _lo, _hi in bands},
        }
        for pid, start, stop, _request, _types, _pttime in pointings
    }
    requests = {pid: request for pid, _start, _stop, request, _types, _pttime in pointings}

    for inst in detectors:
        for event_start, event_stop, _event_source, path in selected.get(inst, []):
            with fits.open(path, memmap=True) as hdul:
                data = event_hdu(hdul).data
                time, pi = data["TIME"], data["PI"]
                assigned = 0.0
                for _inst, _path, _base, pid, start, stop, request in slice_records({inst: [(event_start, event_stop, "", path)]}, pointings):
                    mask = (time >= start) & (time <= stop)
                    band_counts = [int((mask & (pi >= lo) & (pi <= hi)).sum()) for _label, lo, hi in bands]
                    total = int(mask.sum())
                    assigned += stop - start
                    by_pointing[pid]["assigned"] += stop - start
                    by_pointing[pid]["total"] += total
                    for (label, _lo, _hi), count in zip(bands, band_counts):
                        by_pointing[pid][label] += count
                    slice_lines.append(f"{inst}\t{path.name}\t{pid}\t{request}\t{start:.6f}\t{stop:.6f}\t{stop - start:.3f}\t{total}\t" + "\t".join(str(count) for count in band_counts) + "\n")
            coverage_lines.append(
                f"{inst}\t{path.name}\t{event_stop - event_start:.3f}\t"
                f"{assigned:.3f}\t{event_stop - event_start - assigned:.3f}\n"
            )

    by_lines = [f"pointing\trequest\tduration_s\tevent_overlap_s\ttotal_events\t{band_columns}\n"] + [
        f"{pid}\t{requests[pid]}\t{values['duration']:.3f}\t{values['assigned']:.3f}\t{int(values['total'])}\t"
        + "\t".join(str(int(values[label])) for label, _lo, _hi in bands) + "\n"
        for pid, values in by_pointing.items()
    ]
    durations = sorted(stop - start for _pid, start, stop, _request, _types, _pttime in pointings)
    median = durations[len(durations) // 2]
    summary_lines = [f"pointings {len(pointings)}", f"slices {len(slice_lines) - 1}", f"detectors {' '.join(detectors)}", f"median_pointing_s {median:.3f}", f"longest_pointing_s {max(durations):.3f}"]
    summary_lines += [f"warning {pid} has no selected-detector clean-event overlap" for pid, values in by_pointing.items() if values["assigned"] == 0]
    summary_lines += [f"warning {pid} duration is >3x median" for pid, values in by_pointing.items() if values["duration"] > 3.0 * median]
    for name, lines in {
        "pre_exposure_slices.tsv": slice_lines,
        "pre_exposure_event_coverage.tsv": coverage_lines,
        "pre_exposure_by_pointing.tsv": by_lines,
        "pre_exposure_summary.txt": [*summary_lines, ""],
    }.items():
        (outdir / name).write_text("".join(lines) if name.endswith(".tsv") else "\n".join(lines), encoding="utf-8")


def combine_maps(out_counts_arg: str, out_exposure_arg: str, out_rate_arg: str, inputs: list[str]) -> None:
    import numpy as np
    from astropy.io import fits

    if len(inputs) % 2:
        die("combine-maps expects count/exposure pairs")
    pairs = [(Path(inputs[i]), Path(inputs[i + 1])) for i in range(0, len(inputs), 2)]
    if not pairs:
        die("No maps to combine")

    counts_sum = None
    exposure_sum = None
    header = None
    for counts_path, exposure_path in pairs:
        counts, next_header = fits.getdata(counts_path, header=True)
        exposure = fits.getdata(exposure_path).astype(float)
        header = next_header if header is None else header
        counts_sum = counts.astype(float) if counts_sum is None else counts_sum + counts
        exposure_sum = exposure if exposure_sum is None else exposure_sum + exposure

    with np.errstate(divide="ignore", invalid="ignore"):
        rate = np.where(exposure_sum > 0, counts_sum / exposure_sum, np.nan)
    fits.PrimaryHDU(counts_sum, header).writeto(out_counts_arg, overwrite=True)
    fits.PrimaryHDU(exposure_sum, header).writeto(out_exposure_arg, overwrite=True)
    header["BUNIT"] = "count / s"
    fits.PrimaryHDU(rate, header).writeto(out_rate_arg, overwrite=True)


def combine_exposures(out_exposure_arg: str, inputs: list[str]) -> None:
    import numpy as np
    from astropy.io import fits

    paths = [Path(item) for item in inputs]
    if not paths:
        die("combine-exposures expects at least one exposure map")

    exposure_sum = None
    header = None
    for path in paths:
        exposure, next_header = fits.getdata(path, header=True)
        exposure = np.asarray(exposure, dtype=float)
        header = next_header if header is None else header
        exposure_sum = exposure if exposure_sum is None else exposure_sum + exposure

    fits.PrimaryHDU(exposure_sum, header).writeto(out_exposure_arg, overwrite=True)


def background_map(args: list[str]) -> None:
    import numpy as np
    from astropy.io import fits

    if len(args) != 11:
        die("background-map expects COUNTS EXPOSURE BACKGROUND NET MINPIX INNER_FRAC OUTER_FRAC SUMMARY BAND INST SOURCE")
    counts_path, exp_path, out_bkg, out_net, minpix, inner_frac, outer_frac, summary, band, inst, source = args
    counts, header = fits.getdata(counts_path, header=True)
    exposure = fits.getdata(exp_path).astype(float)
    if counts.shape != exposure.shape:
        die(f"Shape mismatch: {counts_path} vs {exp_path}")

    counts = np.asarray(counts, dtype=float)
    status = "ok"
    footprint_pixels = estimate_pixels = finite_pixels = level = 0.0
    cx = cy = detector_radius = inner_radius = outer_radius = np.nan
    background = np.zeros_like(counts, dtype=float)
    footprint = np.isfinite(exposure) & (exposure > 0)
    net = np.full_like(counts, np.nan, dtype=float)

    footprint_pixels = int(footprint.sum())
    if footprint_pixels < int(minpix):
        status = "zero_exposure" if footprint_pixels == 0 else "low_footprint"
    else:
        yy, xx = np.indices(footprint.shape)
        sy, sx = np.nonzero(footprint)
        cx = float(sx.mean())
        cy = float(sy.mean())
        footprint_r = np.hypot(sx - cx, sy - cy)
        detector_radius = float(np.percentile(footprint_r, 99.5))
        inner = max(0.0, float(inner_frac))
        outer = max(inner, float(outer_frac))
        inner_radius = inner * detector_radius
        outer_radius = outer * detector_radius
        rr = np.hypot(xx - cx, yy - cy)
        estimate = footprint & (rr >= inner_radius) & (rr <= outer_radius)
        estimate_pixels = int(estimate.sum())
        if estimate_pixels < int(minpix):
            status = "low_estimate"
        else:
            values = counts[estimate]
            values = values[np.isfinite(values)]
            finite_pixels = int(values.size)
            if finite_pixels < int(minpix):
                status = "low_counts"
            else:
                level = float(np.mean(values))
                background = np.where(footprint, level, 0.0)
                net = np.where(footprint, counts - background, np.nan)

    header["BUNIT"] = "count"
    fits.PrimaryHDU(background, header).writeto(out_bkg, overwrite=True)
    fits.PrimaryHDU(net, header).writeto(out_net, overwrite=True)

    row = (
        f"{band}\t{inst}\t{source}\t{status}\t{footprint_pixels}\t{estimate_pixels}\t{finite_pixels}\t"
        f"mean\t{level:.8g}\t{cx:.3f}\t{cy:.3f}\t"
        f"{detector_radius:.3f}\t{inner_radius:.3f}\t{outer_radius:.3f}\n"
    )
    path = Path(summary)
    if not path.exists():
        path.write_text(
            "band\tinst\tsource\tstatus\tfootprint_pixels\testimate_pixels\tfinite_estimate_pixels\tstatistic\t"
            "background_per_pixel\tcenter_x_px\tcenter_y_px\t"
            "detector_radius_px\tinner_radius_px\touter_radius_px\n",
            encoding="utf-8",
        )
    with path.open("a", encoding="utf-8") as handle:
        handle.write(row)


def combine_background(args: list[str]) -> None:
    import numpy as np
    from astropy.io import fits

    if len(args) < 4 or len(args[2:]) % 2:
        die("combine-background expects OUT_BACKGROUND OUT_NET BACKGROUND NET...")
    out_bkg, out_net = args[:2]
    pairs = [(Path(args[i]), Path(args[i + 1])) for i in range(2, len(args), 2)]
    background_sum = net_sum = valid = None
    header = None
    for bkg_path, net_path in pairs:
        background, next_header = fits.getdata(bkg_path, header=True)
        background = background.astype(float)
        net = fits.getdata(net_path).astype(float)
        finite = np.isfinite(net)
        header = next_header if header is None else header
        background_sum = background if background_sum is None else background_sum + background
        net_sum = np.where(finite, net, 0.0) if net_sum is None else net_sum + np.where(finite, net, 0.0)
        valid = finite if valid is None else valid | finite

    net = np.where(valid, net_sum, np.nan)
    header["BUNIT"] = "count"
    fits.PrimaryHDU(background_sum, header).writeto(out_bkg, overwrite=True)
    fits.PrimaryHDU(net, header).writeto(out_net, overwrite=True)


def net_rate(net_arg: str, exposure_arg: str, out_arg: str) -> None:
    import numpy as np
    from astropy.io import fits

    net, header = fits.getdata(net_arg, header=True)
    exposure = fits.getdata(exposure_arg).astype(float)
    with np.errstate(divide="ignore", invalid="ignore"):
        rate = np.where(exposure > 0, net.astype(float) / exposure, np.nan)
    header["BUNIT"] = "count / s"
    fits.PrimaryHDU(rate, header).writeto(out_arg, overwrite=True)


def sas_image(in_arg: str, out_arg: str) -> None:
    import numpy as np
    from astropy.io import fits

    data, header = fits.getdata(in_arg, header=True)
    data = np.nan_to_num(np.asarray(data, dtype=float), nan=0.0, posinf=0.0, neginf=0.0).astype(np.float32)
    fits.PrimaryHDU(data, header).writeto(out_arg, overwrite=True)


def stack_eventsets(mosaic_dir_arg: str, list_arg: str, summary_arg: str) -> None:
    import csv
    from astropy.io import fits

    rows: list[dict[str, Any]] = []
    for path in sorted(Path(mosaic_dir_arg).rglob("*")):
        name = path.name.lower()
        if not path.is_file() or "_p" not in name or "gti" in name or "tmp" in name:
            continue
        if not (name.endswith(".ds") or name.endswith(".fits")):
            continue
        try:
            with fits.open(path, memmap=True) as hdul:
                hdu = hdul["EVENTS"] if "EVENTS" in hdul else hdul[1]
                events = 0 if hdu.data is None else len(hdu.data)
                gti_rows = sum(0 if ext.data is None else len(ext.data) for ext in hdul if ext.name.upper().startswith("STDGTI"))
                header = hdu.header
        except Exception:
            continue
        if events <= 0:
            continue
        rows.append(
            {
                "path": str(path),
                "events": events,
                "gti_rows": gti_rows,
                "instrument": header.get("INSTRUME", ""),
                "obs_id": header.get("OBS_ID", ""),
                "exp_id": header.get("EXP_ID", header.get("EXPIDSTR", "")),
            }
        )

    Path(list_arg).write_text("\n".join(row["path"] for row in rows) + ("\n" if rows else ""), encoding="utf-8")
    with Path(summary_arg).open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["instrument", "obs_id", "exp_id", "events", "gti_rows", "path"], delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def source_mask(args: list[str]) -> None:
    import numpy as np
    from astropy.io import fits
    from astropy.wcs import WCS

    if len(args) not in {8, 9}:
        die("source-mask expects SRCLIST TEMPLATE KEEP_MASK SOURCE_PIXELS CATALOG BASE_RADIUS LIKE_SCALE MAX_RADIUS [MIN_ML]")
    boxlist_arg, template_arg, keep_arg, source_arg, catalog_arg, base_radius, like_scale, max_radius, *rest = args
    base_radius = float(base_radius)
    like_scale = float(like_scale)
    max_radius = float(max_radius)
    min_ml = float(rest[0]) if rest else 0.0
    template, header = fits.getdata(template_arg, header=True)
    ny, nx = template.shape
    keep = np.ones((ny, nx), dtype=np.uint8)
    source = np.zeros((ny, nx), dtype=np.uint8)
    yy, xx = np.indices((ny, nx))
    wcs = WCS(header).celestial
    rows = ["id\tx_pix\ty_pix\tradius_pix\tlikelihood\tra\tdec\n"]
    seen: set[tuple[int, int]] = set()

    with fits.open(boxlist_arg, memmap=True) as hdul:
        data = hdul["SRCLIST"].data if "SRCLIST" in hdul else hdul[1].data
        if data is not None and len(data):
            names = set(data.names)
            like_col = "LIKE" if "LIKE" in names else "DET_ML"
            for row in data:
                like = float(row[like_col]) if like_col in names else 0.0
                if like < min_ml:
                    continue
                if {"ID_INST", "ID_BAND"} <= names and (int(row["ID_INST"]) != 0 or int(row["ID_BAND"]) != 0):
                    continue
                ra = float(row["RA"]) if "RA" in names and np.isfinite(row["RA"]) else np.nan
                dec = float(row["DEC"]) if "DEC" in names and np.isfinite(row["DEC"]) else np.nan
                if np.isfinite([ra, dec]).all() and wcs.has_celestial:
                    x, y = wcs.world_to_pixel_values(ra, dec)
                    x = float(np.asarray(x))
                    y = float(np.asarray(y))
                elif {"X_IMA", "Y_IMA"} <= names:
                    x = float(row["X_IMA"]) - 1.0
                    y = float(row["Y_IMA"]) - 1.0
                else:
                    continue
                if not np.isfinite([x, y, like]).all():
                    continue
                radius = min(max_radius, base_radius + like_scale * np.sqrt(max(like, 0.0)))
                if x + radius < 0 or x - radius >= nx or y + radius < 0 or y - radius >= ny:
                    continue
                key = (round(x), round(y))
                if key in seen:
                    continue
                seen.add(key)
                inside = (xx - x) ** 2 + (yy - y) ** 2 <= radius**2
                keep[inside] = 0
                source[inside] = 1
                rows.append(f"{len(rows)}\t{x + 1.0:.3f}\t{y + 1.0:.3f}\t{radius:.3f}\t{like:.6g}\t{ra:.8g}\t{dec:.8g}\n")

    header["BUNIT"] = "mask"
    fits.PrimaryHDU(keep, header).writeto(keep_arg, overwrite=True)
    fits.PrimaryHDU(source, header).writeto(source_arg, overwrite=True)
    Path(catalog_arg).write_text("".join(rows), encoding="utf-8")


def apply_source_mask(args: list[str]) -> None:
    import numpy as np
    from astropy.io import fits

    if len(args) != 5:
        die("apply-source-mask expects KEEP_MASK IN_FITS OUT_FITS SUMMARY LABEL")
    mask_arg, in_arg, out_arg, summary_arg, label = args
    keep = fits.getdata(mask_arg).astype(bool)
    data, header = fits.getdata(in_arg, header=True)
    data = np.asarray(data, dtype=float)
    if data.shape != keep.shape:
        die(f"Mask shape mismatch: {mask_arg} vs {in_arg}")
    masked = np.where(keep, data, np.nan)
    fits.PrimaryHDU(masked, header).writeto(out_arg, overwrite=True)

    path = Path(summary_arg)
    if not path.exists():
        path.write_text("label\tinput\toutput\tmasked_pixels\tkept_pixels\n", encoding="utf-8")
    with path.open("a", encoding="utf-8") as handle:
        handle.write(f"{label}\t{in_arg}\t{out_arg}\t{int((~keep).sum())}\t{int(keep.sum())}\n")


def filter_events_mask(args: list[str]) -> None:
    import numpy as np
    from astropy.io import fits

    if len(args) != 6:
        die("filter-events-mask expects KEEP_MASK GRID_JSON IN_EVENT OUT_EVENT SUMMARY LABEL")
    mask_arg, grid_arg, in_arg, out_arg, summary_arg, label = args
    keep_image = fits.getdata(mask_arg).astype(bool)
    grid = json.loads(Path(grid_arg).read_text(encoding="utf-8"))
    x0, y0, bin_phys = float(grid["x_min"]), float(grid["y_min"]), float(grid["bin_phys"])

    with fits.open(in_arg, memmap=False) as hdul:
        evt = event_hdu(hdul)
        data = evt.data
        if data is None or len(data) == 0:
            fits.HDUList([hdu.copy() for hdu in hdul]).writeto(out_arg, overwrite=True)
            removed = kept = 0
        else:
            x = np.asarray(data["X"], dtype=float)
            y = np.asarray(data["Y"], dtype=float)
            ix = np.floor((x - x0) / bin_phys).astype(int)
            iy = np.floor((y - y0) / bin_phys).astype(int)
            inside = (ix >= 0) & (ix < keep_image.shape[1]) & (iy >= 0) & (iy < keep_image.shape[0])
            keep = np.ones(len(data), dtype=bool)
            keep[inside] = keep_image[iy[inside], ix[inside]]
            removed = int((~keep).sum())
            kept = int(keep.sum())
            out_hdus = []
            for hdu in hdul:
                if hdu is evt:
                    new = fits.BinTableHDU(data=data[keep], header=hdu.header.copy(), name=hdu.name)
                    new.header["SRCMASK"] = (True, "Sky-frame source mask applied before comet attcalc")
                    new.header["SRCREM"] = (removed, "Events removed by source mask")
                    out_hdus.append(new)
                else:
                    out_hdus.append(hdu.copy())
            Path(out_arg).parent.mkdir(parents=True, exist_ok=True)
            fits.HDUList(out_hdus).writeto(out_arg, overwrite=True)

    path = Path(summary_arg)
    if not path.exists():
        path.write_text("label\tinput\toutput\tremoved_events\tkept_events\n", encoding="utf-8")
    with path.open("a", encoding="utf-8") as handle:
        handle.write(f"{label}\t{in_arg}\t{out_arg}\t{removed}\t{kept}\n")


def source_trail_mask(args: list[str]) -> None:
    import numpy as np
    from astropy.io import fits

    if len(args) != 6:
        die("source-trail-mask expects SOURCE_CATALOG SKY_GRID TRACK_TSV TEMPLATE OUT_MASK SUMMARY")
    catalog_arg, sky_grid_arg, track_arg, template_arg, out_arg, summary_arg = args
    rows = read_tsv(catalog_arg) if Path(catalog_arg).exists() else []
    track = read_tsv(track_arg) if Path(track_arg).exists() else []
    data, header = fits.getdata(template_arg, header=True)
    keep = np.ones(np.asarray(data).shape, dtype=np.uint8)
    if not rows or not track:
        fits.PrimaryHDU(keep, header).writeto(out_arg, overwrite=True)
        Path(summary_arg).write_text("sources\ttrack_samples\tmasked_pixels\n0\t0\t0\n", encoding="utf-8")
        return

    sky_grid = json.loads(Path(sky_grid_arg).read_text(encoding="utf-8"))
    stack_grid = {
        "x_min": float(header.get("XIMIN", sky_grid["x_min"])),
        "y_min": float(header.get("YIMIN", sky_grid["y_min"])),
        "bin_phys": float(sky_grid["bin_phys"]),
    }
    if Path(template_arg).parent.parent.joinpath("grid.json").exists():
        stack_grid.update(json.loads(Path(template_arg).parent.parent.joinpath("grid.json").read_text(encoding="utf-8")))
    sky_bin = float(sky_grid["bin_phys"])
    stack_bin = float(stack_grid["bin_phys"])
    ny, nx = keep.shape
    masked = 0
    for src in rows:
        if src.get("id") == "id":
            continue
        try:
            x_sky = float(sky_grid["x_min"]) + (float(src["x_pix"]) - 0.5) * sky_bin
            y_sky = float(sky_grid["y_min"]) + (float(src["y_pix"]) - 0.5) * sky_bin
            radius = max(1.0, float(src["radius_pix"]) * sky_bin / stack_bin)
        except (KeyError, ValueError):
            continue
        r = int(np.ceil(radius))
        for tr in track:
            x = (x_sky + float(tr["dx_pix"]) * sky_bin - float(stack_grid["x_min"])) / stack_bin
            y = (y_sky + float(tr["dy_pix"]) * sky_bin - float(stack_grid["y_min"])) / stack_bin
            ix0, ix1 = max(0, int(np.floor(x - r))), min(nx, int(np.ceil(x + r + 1)))
            iy0, iy1 = max(0, int(np.floor(y - r))), min(ny, int(np.ceil(y + r + 1)))
            if ix0 >= ix1 or iy0 >= iy1:
                continue
            yy, xx = np.ogrid[iy0:iy1, ix0:ix1]
            inside = (xx - x) ** 2 + (yy - y) ** 2 <= radius**2
            before = int(keep[iy0:iy1, ix0:ix1].sum())
            keep[iy0:iy1, ix0:ix1][inside] = 0
            masked += before - int(keep[iy0:iy1, ix0:ix1].sum())
    header["BUNIT"] = "mask"
    fits.PrimaryHDU(keep, header).writeto(out_arg, overwrite=True)
    Path(summary_arg).write_text(f"sources\ttrack_samples\tmasked_pixels\n{len(rows)}\t{len(track)}\t{int((keep == 0).sum())}\n", encoding="utf-8")


def stamp_stack_map(args: list[str]) -> None:
    from astropy.io import fits

    if len(args) != 3:
        die("stamp-stack-map expects FITS REF_RA REF_DEC")
    path_arg, ref_ra, ref_dec = args
    with fits.open(path_arg, mode="update") as hdul:
        for hdu in hdul:
            if getattr(hdu, "data", None) is not None:
                hdu.header["COMFRAME"] = "COMET"
                hdu.header["REFRA"] = float(ref_ra)
                hdu.header["REFDEC"] = float(ref_dec)
        hdul.flush()


def read_tsv(path: str) -> list[dict[str, str]]:
    import csv

    with Path(path).open("r", newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def event_time_reference(event_path: str) -> tuple[float, str]:
    from astropy.io import fits

    header = fits.getheader(event_path, "EVENTS")
    mjdref = header.get("MJDREF")
    if mjdref is None:
        mjdref = float(header.get("MJDREFI", 0.0)) + float(header.get("MJDREFF", 0.0))
    return float(mjdref), str(header.get("TIMESYS", "TT")).lower()


def event_time_span(events: list[str]) -> tuple[float, float, float, str]:
    import numpy as np
    from astropy.io import fits

    tmin = np.inf
    tmax = -np.inf
    mjdref = None
    timesys = "tt"
    for event in events:
        with fits.open(event, memmap=True) as hdul:
            hdu = event_hdu(hdul)
            data = hdu.data
            if data is None or len(data) == 0:
                continue
            t = np.asarray(data["TIME"], dtype=float)
            good = np.isfinite(t)
            if not good.any():
                continue
            tmin = min(tmin, float(t[good].min()))
            tmax = max(tmax, float(t[good].max()))
            header = hdu.header
            if mjdref is None:
                value = header.get("MJDREF")
                mjdref = float(value) if value is not None else float(header.get("MJDREFI", 0.0)) + float(header.get("MJDREFF", 0.0))
                timesys = str(header.get("TIMESYS", "TT")).lower()
    if not np.isfinite([tmin, tmax]).all() or mjdref is None:
        die("Could not determine clean-event time span for stack track")
    return tmin, tmax, mjdref, timesys


def query_horizons(jd_tt: list[float], target: str, center: str, raw_arg: str) -> list[tuple[float, float, float]]:
    import csv
    import json
    import urllib.parse
    import urllib.request

    command = target if target.startswith("'") else f"'{target}'"
    results: list[str] = []
    rows: list[tuple[float, float, float]] = []
    for start in range(0, len(jd_tt), 50):
        chunk = jd_tt[start : start + 50]
        params = {
            "format": "json",
            "COMMAND": command,
            "OBJ_DATA": "NO",
            "MAKE_EPHEM": "YES",
            "EPHEM_TYPE": "OBSERVER",
            "CENTER": center,
            "TLIST": ",".join(f"{jd:.10f}" for jd in chunk),
            "TLIST_TYPE": "JD",
            "TIME_TYPE": "TT",
            "TIME_DIGITS": "FRACSEC",
            "QUANTITIES": "'1,20'",
            "ANG_FORMAT": "DEG",
            "CSV_FORMAT": "YES",
        }
        url = "https://ssd.jpl.nasa.gov/api/horizons.api?" + urllib.parse.urlencode(params)
        with urllib.request.urlopen(url, timeout=120) as response:
            payload = json.load(response)
        if payload.get("error"):
            die(payload["error"].replace("\n", "; "))
        result = payload.get("result", "")
        results.append(result)

        in_table = False
        for line in result.splitlines():
            if "$$SOE" in line:
                in_table = True
                continue
            if "$$EOE" in line:
                break
            if not in_table:
                continue
            fields = [field.strip() for field in next(csv.reader([line]))]
            nums = []
            for field in fields[1:]:
                try:
                    nums.append(float(field))
                except ValueError:
                    pass
            if len(nums) >= 2:
                rows.append((nums[0], nums[1], nums[2] if len(nums) >= 3 else 1.0))
    Path(raw_arg).write_text("\n\n".join(results), encoding="utf-8")
    if len(rows) != len(jd_tt):
        die(f"Horizons returned {len(rows)} ephemeris rows for {len(jd_tt)} requested times")
    return rows


def local_ephemeris(jd_tt: list[float], path_arg: str) -> list[tuple[float, float, float]]:
    import numpy as np

    rows = read_tsv(path_arg)
    if not rows:
        die(f"Empty ephemeris file: {path_arg}")
    if "jd_tt" in rows[0]:
        src_jd = np.array([float(row["jd_tt"]) for row in rows], dtype=float)
    elif "mjd_tt" in rows[0]:
        src_jd = np.array([float(row["mjd_tt"]) + 2400000.5 for row in rows], dtype=float)
    else:
        die("Ephemeris file needs jd_tt or mjd_tt plus ra_deg and dec_deg columns")
    ra = np.array([float(row["ra_deg"]) for row in rows], dtype=float)
    dec = np.array([float(row["dec_deg"]) for row in rows], dtype=float)
    delta = np.array([float(row.get("delta_au", 1.0)) for row in rows], dtype=float)
    return list(zip(np.interp(jd_tt, src_jd, ra), np.interp(jd_tt, src_jd, dec), np.interp(jd_tt, src_jd, delta)))


def stack_track(args: list[str]) -> None:
    import json
    import numpy as np
    from astropy.io import fits
    from astropy.time import Time
    from astropy.time import TimeDelta
    from astropy.wcs import WCS

    if len(args) != 17:
        die("stack-track expects CLEAN_DIR DETECTORS TEMPLATE ATTMOVE_TSV ATTMOVE_FITS ATTMOVE_RAW OVERLAY_TSV OVERLAY_RAW META_JSON REF_ENV TARGET ATTMOVE_CENTER OVERLAY_CENTER EPHEM_FILE REF_RA REF_DEC STEP_S")
    (
        clean_dir,
        detectors,
        template_arg,
        att_tsv_arg,
        fits_arg,
        att_raw_arg,
        overlay_tsv_arg,
        overlay_raw_arg,
        meta_arg,
        ref_env_arg,
        target,
        attmove_center,
        overlay_center,
        ephem_file,
        ref_ra_arg,
        ref_dec_arg,
        step_arg,
    ) = args
    events = [str(path) for path in manifest_events(clean_dir, detectors)]
    t0, t1, mjdref, timesys = event_time_span(events)
    step = max(float(step_arg), 1.0)
    met = np.arange(t0, t1 + 0.5 * step, step, dtype=float)
    if met[-1] < t1:
        met = np.append(met, t1)
    times = Time(mjdref, format="mjd", scale=timesys) + TimeDelta(met, format="sec")
    jd_tt = [float(value) for value in times.tt.jd]
    att_eph = local_ephemeris(jd_tt, ephem_file) if ephem_file else query_horizons(jd_tt, target, attmove_center, att_raw_arg)
    overlay_eph = att_eph if ephem_file or overlay_center == attmove_center else query_horizons(jd_tt, target, overlay_center, overlay_raw_arg)

    mid_index = int(np.argsort(jd_tt)[len(jd_tt) // 2])
    ref_ra = float(ref_ra_arg) if ref_ra_arg else float(att_eph[mid_index][0])
    ref_dec = float(ref_dec_arg) if ref_dec_arg else float(att_eph[mid_index][1])
    wcs = WCS(fits.getheader(template_arg)).celestial
    ref_x, ref_y = wcs.world_to_pixel_values(ref_ra, ref_dec)
    ref_x = float(np.asarray(ref_x))
    ref_y = float(np.asarray(ref_y))

    def write_tsv(path_arg: str, eph: list[tuple[float, float, float]], center: str) -> None:
        lines = ["time_met\tmjd_utc\tjd_tt\tra_deg\tdec_deg\tdelta_au\tref_ra_deg\tref_dec_deg\tcomet_x\tcomet_y\tref_x\tref_y\tdx_pix\tdy_pix\tcenter\n"]
        for t_met, time, jd, (ra, dec, delta) in zip(met, times, jd_tt, eph):
            x, y = wcs.world_to_pixel_values(ra, dec)
            x = float(np.asarray(x))
            y = float(np.asarray(y))
            lines.append(
                f"{t_met:.6f}\t{time.utc.mjd:.10f}\t{jd:.10f}\t"
                f"{ra:.10f}\t{dec:.10f}\t{delta:.8g}\t{ref_ra:.10f}\t{ref_dec:.10f}\t{x:.6f}\t{y:.6f}\t{ref_x:.6f}\t{ref_y:.6f}\t"
                f"{ref_x - x:.6f}\t{ref_y - y:.6f}\t{center}\n"
            )
        Path(path_arg).write_text("".join(lines), encoding="utf-8")

    write_tsv(att_tsv_arg, att_eph, attmove_center)
    write_tsv(overlay_tsv_arg, overlay_eph, overlay_center)
    cols = fits.ColDefs(
        [
            fits.Column(name="MJD", format="D", unit="d", array=np.asarray([t.utc.mjd for t in times], dtype=np.float64)),
            fits.Column(name="RA", format="E", unit="deg", array=np.asarray([row[0] for row in att_eph], dtype=np.float32)),
            fits.Column(name="DEC", format="E", unit="deg", array=np.asarray([row[1] for row in att_eph], dtype=np.float32)),
            fits.Column(name="DELTA", format="E", unit="AU", array=np.asarray([row[2] for row in att_eph], dtype=np.float32)),
        ]
    )
    hdu = fits.BinTableHDU.from_columns(cols, name="OBJTRACK")
    hdu.header["OBS_T0"] = float(t0)
    hdu.header["OBS_T1"] = float(t1)
    hdu.header["MJDREF"] = float(mjdref)
    hdu.header["TIMESYS"] = timesys.upper()
    hdu.header["REFFRAME"] = "ICRF"
    hdu.header["CENTER"] = (attmove_center, "Horizons center for SAS attmove track")
    fits.HDUList([fits.PrimaryHDU(), hdu]).writeto(fits_arg, overwrite=True)
    Path(ref_env_arg).write_text(f'COMET_REF_RA={ref_ra:.10f}\nCOMET_REF_DEC={ref_dec:.10f}\n', encoding="utf-8")
    Path(meta_arg).write_text(
        json.dumps(
            {
                "target": target,
                "attmove_center": attmove_center,
                "overlay_center": overlay_center,
                "track_samples": int(len(met)),
                "track_step_s": step,
                "timesys": timesys.upper(),
                "mjdref": mjdref,
                "tstart": t0,
                "tstop": t1,
                "reference_ra_deg": ref_ra,
                "reference_dec_deg": ref_dec,
                "reference_x_pix": ref_x,
                "reference_y_pix": ref_y,
                "ephemeris_file": ephem_file,
                "attmove_track": str(Path(att_tsv_arg).resolve()),
                "overlay_track": str(Path(overlay_tsv_arg).resolve()),
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )


def fits_png(
    fits_arg: str,
    png_arg: str,
    cmap_arg: str = "magma",
    scale_arg: str = "linear",
    grid_arg: str = "",
    manifest_arg: str = "",
    detectors_arg: str = "PN",
    colorbar_arg: str = "",
    circle_summary_arg: str = "",
    circle_band_arg: str = "",
) -> None:
    import matplotlib
    import numpy as np
    from astropy.io import fits
    from matplotlib.patches import Circle
    from matplotlib.colors import LogNorm

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    data, header = fits.getdata(fits_arg, header=True)
    data = np.asarray(data, dtype=float)
    image_extent, plot_extent = mosaic_extents(data.shape, grid_arg, manifest_arg, detectors_arg)
    ratio = (plot_extent[3] - plot_extent[2]) / max(plot_extent[1] - plot_extent[0], 1.0)

    finite = np.isfinite(data)
    positive = data[finite & (data > 0)]
    image = np.ma.masked_where(~finite | (data <= 0), data)
    cmap = plt.get_cmap(cmap_arg).copy()
    cmap.set_bad("black")
    cmap.set_under("black")
    norm = None
    vmin = vmax = None
    if scale_arg == "log" and positive.size:
        vmin = max(float(np.percentile(positive, 1.0)), np.finfo(float).tiny)
        vmax = max(float(np.percentile(positive, 99.7)), vmin * 1.01)
        norm = LogNorm(vmin=vmin, vmax=vmax)
        vmin = vmax = None
    elif positive.size:
        vmin = 0.0
        vmax = float(np.percentile(positive, 99.7))
        if vmax <= vmin:
            vmax = float(positive.max())
    if scale_arg.startswith("floor"):
        match = re.fullmatch(r"floor(?:(\d+(?:\.\d+)?))?", scale_arg)
        percentile = float(match.group(1)) if match and match.group(1) else 99.7
        image = np.ma.masked_where(~finite, np.maximum(data, 0.0))
        positive = image.compressed()
        positive = positive[positive > 0]
        vmin = 0.0
        vmax = float(np.nanpercentile(positive, percentile)) if positive.size else 1.0
        vmax = vmax if vmax > 0 else 1.0
        norm = None
    if scale_arg.startswith("signed"):
        finite_values = data[finite]
        match = re.fullmatch(r"signed(?:(\d+(?:\.\d+)?))?", scale_arg)
        percentile = float(match.group(1)) if match and match.group(1) else 99.0
        span = float(np.nanpercentile(np.abs(finite_values), percentile)) if finite_values.size else 1.0
        span = span if span > 0 else 1.0
        image = np.ma.masked_where(~finite, data)
        norm = None
        vmin = -span
        vmax = span

    fig, ax = plt.subplots(figsize=(8, max(3, 8 * ratio)), constrained_layout=True)
    im = ax.imshow(image, origin="lower", extent=image_extent, cmap=cmap, norm=norm, vmin=vmin, vmax=vmax, interpolation="nearest")
    if circle_summary_arg and Path(circle_summary_arg).exists():
        ny, nx = data.shape
        dx = (image_extent[1] - image_extent[0]) / nx
        dy = (image_extent[3] - image_extent[2]) / ny
        for row in read_tsv(circle_summary_arg):
            if circle_band_arg and row.get("band") != circle_band_arg:
                continue
            try:
                cx = float(row["center_x_px"])
                cy = float(row["center_y_px"])
                inner = float(row.get("inner_radius_px", row.get("exclude_radius_px", "nan")))
                outer = float(row.get("outer_radius_px", "nan"))
            except (KeyError, ValueError):
                continue
            if not np.isfinite([cx, cy, inner]).all() or inner <= 0:
                continue
            x = image_extent[0] + (cx + 0.5) * dx
            y = image_extent[2] + (cy + 0.5) * dy
            for radius, alpha in ((inner, 0.95), (outer, 0.75)):
                if not np.isfinite(radius) or radius <= 0:
                    continue
                r = radius * 0.5 * (abs(dx) + abs(dy))
                ax.add_patch(Circle((x, y), r, fill=False, edgecolor="white", linewidth=0.9, alpha=alpha, zorder=8))
                ax.add_patch(Circle((x, y), r, fill=False, edgecolor="cyan", linewidth=0.55, alpha=alpha, zorder=9))
    ax.set_xlim(plot_extent[0], plot_extent[1])
    ax.set_ylim(plot_extent[2], plot_extent[3])
    ax.set_aspect("equal", adjustable="box")
    ax.set_axis_off()
    if colorbar_arg in {"colorbar", "cbar"}:
        cbar = fig.colorbar(im, ax=ax, fraction=0.035, pad=0.02)
        cbar.ax.tick_params(colors="white", labelsize=7, length=2)
        cbar.outline.set_edgecolor("white")
        if header.get("BUNIT"):
            cbar.set_label(header["BUNIT"], color="white", fontsize=8)
    Path(png_arg).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(png_arg, dpi=180, facecolor="black", pad_inches=0.05)
    plt.close(fig)


def stack_png(args: list[str]) -> None:
    import matplotlib
    import numpy as np
    from astropy.io import fits
    from astropy.wcs import WCS
    from matplotlib.colors import LogNorm

    if len(args) not in {6, 9}:
        die("stack-png expects FITS PNG TRACK MODE CMAP SCALE [GRID_JSON MANIFEST DETECTORS]")
    fits_arg, png_arg, track_arg, mode, cmap_arg, scale_arg, *rest = args
    grid_arg, manifest_arg, detectors_arg = (rest + ["", "", "PN"])[:3]

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    data, header = fits.getdata(fits_arg, header=True)
    data = np.asarray(data, dtype=float)
    image_extent, plot_extent = mosaic_extents(data.shape, grid_arg, manifest_arg, detectors_arg)
    ratio = (plot_extent[3] - plot_extent[2]) / max(plot_extent[1] - plot_extent[0], 1.0)
    finite = np.isfinite(data)
    positive = data[finite & (data > 0)]
    image = np.ma.masked_where(~finite | (data <= 0), data)
    cmap = plt.get_cmap(cmap_arg).copy()
    cmap.set_bad("black")
    cmap.set_under("black")
    norm = None
    vmin = vmax = None
    if scale_arg == "log" and positive.size:
        vmin = max(float(np.nanpercentile(positive, 1.0)), np.finfo(float).tiny)
        vmax = max(float(np.nanpercentile(positive, 99.7)), vmin * 1.01)
        norm = LogNorm(vmin=vmin, vmax=vmax)
        vmin = vmax = None
    elif scale_arg.startswith("floor"):
        image = np.ma.masked_where(~finite, np.maximum(data, 0.0))
        positive = image.compressed()
        positive = positive[positive > 0]
        vmin = 0.0
        vmax = float(np.nanpercentile(positive, 99.0)) if positive.size else 1.0
    elif positive.size:
        vmin = 0.0
        vmax = float(np.nanpercentile(positive, 99.7))

    def to_plot(x: Any, y: Any) -> tuple[np.ndarray, np.ndarray]:
        ny, nx = data.shape
        return (
            image_extent[0] + (np.asarray(x, dtype=float) + 0.5) * (image_extent[1] - image_extent[0]) / nx,
            image_extent[2] + (np.asarray(y, dtype=float) + 0.5) * (image_extent[3] - image_extent[2]) / ny,
        )

    track = read_tsv(track_arg)
    fig, ax = plt.subplots(figsize=(8, max(3, 8 * ratio)), constrained_layout=True)
    im = ax.imshow(image, origin="lower", extent=image_extent, cmap=cmap, norm=norm, vmin=vmin, vmax=vmax, interpolation="nearest")
    if track:
        if mode == "trajectory":
            xs, ys = to_plot([row["comet_x"] for row in track], [row["comet_y"] for row in track])
            ax.plot(xs, ys, color="cyan", lw=1.0, alpha=0.9)
            ax.scatter(xs[0], ys[0], s=12, color="lime", marker="o")
            ax.scatter(xs[-1], ys[-1], s=18, color="red", marker="x")
        else:
            ref_ra = float(header.get("REFRA", track[0]["ref_ra_deg"]))
            ref_dec = float(header.get("REFDEC", track[0]["ref_dec_deg"]))
            if "CRPIX1" in header and "CRPIX2" in header and abs(float(header.get("CRVAL1", ref_ra)) - ref_ra) < 1e-6 and abs(float(header.get("CRVAL2", ref_dec)) - ref_dec) < 1e-6:
                x = float(header["CRPIX1"]) - 1.0
                y = float(header["CRPIX2"]) - 1.0
            else:
                try:
                    x, y = WCS(header).celestial.world_to_pixel_values(ref_ra, ref_dec)
                    x = float(np.asarray(x))
                    y = float(np.asarray(y))
                except Exception:
                    x = float(track[0]["ref_x"]) - float(header.get("STKX0", 0.0))
                    y = float(track[0]["ref_y"]) - float(header.get("STKY0", 0.0))
            if not (np.isfinite(x) and np.isfinite(y)):
                x = float(header.get("CRPIX1", data.shape[1] / 2.0)) - 1.0
                y = float(header.get("CRPIX2", data.shape[0] / 2.0)) - 1.0
            xs, ys = to_plot(x, y)
            xs = float(np.asarray(xs))
            ys = float(np.asarray(ys))
            ax.scatter([xs], [ys], s=260, facecolors="none", edgecolors="white", linewidths=1.2, zorder=8)
            ax.scatter([xs], [ys], s=170, facecolors="none", edgecolors="cyan", linewidths=1.1, zorder=9)
            ax.scatter([xs], [ys], s=95, color="white", marker="+", linewidths=1.6, zorder=10)
            ax.scatter([xs], [ys], s=50, color="cyan", marker="+", linewidths=1.3, zorder=11)
    ax.set_xlim(plot_extent[0], plot_extent[1])
    ax.set_ylim(plot_extent[2], plot_extent[3])
    ax.set_aspect("equal", adjustable="box")
    ax.set_axis_off()
    cbar = fig.colorbar(im, ax=ax, fraction=0.035, pad=0.02)
    cbar.ax.tick_params(colors="white", labelsize=7, length=2)
    cbar.outline.set_edgecolor("white")
    if header.get("BUNIT"):
        cbar.set_label(header["BUNIT"], color="white", fontsize=8)
    Path(png_arg).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(png_arg, dpi=180, facecolor="black", pad_inches=0.05)
    plt.close(fig)


def mosaic_extents(
    shape: tuple[int, int],
    grid_arg: str = "",
    manifest_arg: str = "",
    detectors_arg: str = "PN",
) -> tuple[tuple[float, float, float, float], tuple[float, float, float, float]]:
    ny, nx = shape
    image_extent = (0.0, float(nx), 0.0, float(ny))
    plot_extent = image_extent
    if grid_arg and manifest_arg:
        grid = json.loads(Path(grid_arg).read_text(encoding="utf-8"))
        image_extent = (grid["x_min"], grid["x_max"], grid["y_min"], grid["y_max"])
        _bounds, plot_extent, _nx, _ny = event_bounds(manifest_events(manifest_arg, detectors_arg), image_extent)
    return image_extent, plot_extent


def footprint_outlines(args: list[str]) -> None:
    import matplotlib
    import numpy as np
    from astropy.io import fits

    if len(args) < 5:
        die("footprint-outlines expects PNG GRID_JSON MANIFEST DETECTORS EXPOSURE...")
    png_arg, grid_arg, manifest_arg, detectors_arg, *exposure_args = args

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    first = np.asarray(fits.getdata(exposure_args[0]), dtype=float)
    ny, nx = first.shape
    image_extent, plot_extent = mosaic_extents(first.shape, grid_arg, manifest_arg, detectors_arg)
    ratio = (plot_extent[3] - plot_extent[2]) / max(plot_extent[1] - plot_extent[0], 1.0)
    fig, ax = plt.subplots(figsize=(8, max(3, 8 * ratio)), constrained_layout=True)
    colors = plt.get_cmap("magma")(np.linspace(0.05, 0.95, len(exposure_args)))
    dx = (image_extent[1] - image_extent[0]) / nx
    dy = (image_extent[3] - image_extent[2]) / ny

    for color, exposure_arg in zip(colors, exposure_args):
        data = np.asarray(fits.getdata(exposure_arg), dtype=float)
        footprint = np.isfinite(data) & (data > 0)
        if not footprint.any():
            continue
        ax.contour(footprint.astype(float), levels=[0.5], origin="lower", extent=image_extent, colors=[color], linewidths=0.8)
        sy, sx = np.nonzero(footprint)
        label = re.search(r"_(P\d{3})_", Path(exposure_arg).name)
        if label:
            ax.text(
                image_extent[0] + (sx.mean() + 0.5) * dx,
                image_extent[2] + (sy.mean() + 0.5) * dy,
                label.group(1),
                color=color,
                fontsize=5,
                ha="center",
                va="center",
            )

    ax.set_xlim(plot_extent[0], plot_extent[1])
    ax.set_ylim(plot_extent[2], plot_extent[3])
    ax.set_aspect("equal", adjustable="box")
    ax.set_axis_off()
    Path(png_arg).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(png_arg, dpi=180, facecolor="black", pad_inches=0.05)
    plt.close(fig)


def event_mosaic(manifest_dir: str, outdir: str, detectors_arg: str = "PN") -> None:
    import matplotlib
    import numpy as np
    from matplotlib.colors import LogNorm

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    detectors = normalize_detectors(detectors_arg)
    events = manifest_events(manifest_dir, " ".join(detectors))
    if not events:
        die(f"No event lists found in {manifest_dir}")

    data_bounds, extent, nx, ny = event_bounds(events)

    bands = {"soft_lt1kev": (None, 1000, "magma"), "hard_gt1kev": (1000, None, "magma")}
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)
    summary: list[str] = [
        f"detectors {' '.join(detectors)}",
        f"event_lists {len(events)}",
        f"data_bounds {data_bounds}",
        f"plot_extent {extent}",
        f"plot_pixels {nx} {ny}",
    ]

    for tag, (lo, hi, cmap) in bands.items():
        image = np.zeros((ny, nx), dtype=float)
        total = 0
        for path in events:
            hdul, x, y, pi = event_columns(path)
            good = valid_sky(x, y) & np.isfinite(pi)
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

        positive = image[image > 0]
        vmax = max(float(np.percentile(positive, 99.7)), 1.0) if positive.size else 1.0
        fig, ax = plt.subplots(figsize=(8, 8 * ny / nx), constrained_layout=True)
        ax.imshow(
            image,
            origin="lower",
            extent=extent,
            cmap=cmap,
            norm=LogNorm(vmin=1, vmax=vmax) if positive.size else None,
            interpolation="nearest",
        )
        ax.set_xlim(extent[0], extent[1])
        ax.set_ylim(extent[2], extent[3])
        ax.set_aspect("equal", adjustable="box")
        ax.set_axis_off()
        fig.savefig(out / f"{tag}_mosaic.png", dpi=180, facecolor="black", pad_inches=0.05)
        plt.close(fig)
        summary.append(f"{tag}_events {total}")

    (out / "mosaic_summary.txt").write_text("\n".join(summary) + "\n", encoding="utf-8")


def main(argv: list[str]) -> int:
    cmd, args = (argv[1], argv[2:]) if len(argv) > 1 else ("", [])
    if cmd == "shell" and len(args) == 1:
        shell(*args)
    elif cmd == "relative-path" and len(args) == 2:
        relative_path(*args)
    elif cmd == "rewrite-path" and len(args) == 2:
        rewrite_path(*args)
    elif cmd == "event-mosaic" and len(args) in {2, 3}:
        event_mosaic(args[0], args[1], args[2] if len(args) == 3 else "PN")
    elif cmd == "fits-rows" and len(args) == 1:
        fits_rows(*args)
    elif cmd == "flare-qc" and len(args) == 3:
        flare_qc(*args)
    elif cmd == "exposure-grid" and len(args) >= 4:
        exposure_grid(args[0], args[1], args[2], args[3:])
    elif cmd == "exposure-slices" and len(args) == 5:
        exposure_slices(*args)
    elif cmd == "slice-qc" and len(args) == 5:
        slice_qc(*args)
    elif cmd == "combine-maps" and len(args) >= 5:
        combine_maps(args[0], args[1], args[2], args[3:])
    elif cmd == "combine-exposures" and len(args) >= 2:
        combine_exposures(args[0], args[1:])
    elif cmd == "background-map" and len(args) == 11:
        background_map(args)
    elif cmd == "combine-background" and len(args) >= 4:
        combine_background(args)
    elif cmd == "net-rate" and len(args) == 3:
        net_rate(*args)
    elif cmd == "sas-image" and len(args) == 2:
        sas_image(*args)
    elif cmd == "stack-eventsets" and len(args) == 3:
        stack_eventsets(*args)
    elif cmd == "source-mask" and len(args) in {8, 9}:
        source_mask(args)
    elif cmd == "apply-source-mask" and len(args) == 5:
        apply_source_mask(args)
    elif cmd == "filter-events-mask" and len(args) == 6:
        filter_events_mask(args)
    elif cmd == "source-trail-mask" and len(args) == 6:
        source_trail_mask(args)
    elif cmd == "stamp-stack-map" and len(args) == 3:
        stamp_stack_map(args)
    elif cmd == "stack-track" and len(args) == 17:
        stack_track(args)
    elif cmd == "stack-png" and len(args) in {6, 9}:
        stack_png(args)
    elif cmd == "footprint-outlines" and len(args) >= 5:
        footprint_outlines(args)
    elif cmd == "fits-png" and 2 <= len(args) <= 10:
        fits_png(*args)
    else:
        die("Usage: tools.py shell CONFIG | relative-path ROOT PATH | rewrite-path SUM.SAS PATH | event-mosaic MANIFEST OUTDIR [DETECTORS] | fits-rows FITS | flare-qc GTI_DIR OUTDIR RATE_CUT | exposure-grid OUTDIR BIN PAD EVENTS... | exposure-slices SLICES_TSV POINTINGS_TSV CLEAN_DIR DETECTORS INIT_DIR | slice-qc OUTDIR CLEAN_DIR DETECTORS INIT_DIR BANDS | combine-maps OUT_COUNTS OUT_EXP OUT_RATE COUNTS EXP... | combine-exposures OUT_EXP EXP... | background-map COUNTS EXPOSURE BACKGROUND NET MINPIX INNER_FRAC OUTER_FRAC SUMMARY BAND INST SOURCE | combine-background OUT_BACKGROUND OUT_NET BACKGROUND NET... | net-rate NET EXPOSURE OUT | sas-image IN_FITS OUT_FITS | stack-eventsets MOSAIC_DIR LIST SUMMARY | source-mask SRCLIST TEMPLATE KEEP_MASK SOURCE_PIXELS CATALOG BASE_RADIUS LIKE_SCALE MAX_RADIUS [MIN_ML] | apply-source-mask KEEP_MASK IN_FITS OUT_FITS SUMMARY LABEL | filter-events-mask KEEP_MASK GRID_JSON IN_EVENT OUT_EVENT SUMMARY LABEL | source-trail-mask SOURCE_CATALOG SKY_GRID TRACK_TSV TEMPLATE OUT_MASK SUMMARY | stamp-stack-map FITS REF_RA REF_DEC | stack-track CLEAN_DIR DETECTORS TEMPLATE ATTMOVE_TSV ATTMOVE_FITS ATTMOVE_RAW OVERLAY_TSV OVERLAY_RAW META_JSON REF_ENV TARGET ATTMOVE_CENTER OVERLAY_CENTER EPHEM_FILE REF_RA REF_DEC STEP_S | stack-png FITS PNG TRACK MODE CMAP SCALE [GRID_JSON MANIFEST DETECTORS] | footprint-outlines PNG GRID_JSON MANIFEST DETECTORS EXPOSURE... | fits-png FITS PNG [CMAP] [linear|log|signed|floor] [GRID_JSON MANIFEST DETECTORS [colorbar [CIRCLE_SUMMARY BAND]]]")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
