#!/usr/bin/env python3
"""Helper commands for reduction_v6/pipeline.sh.

Keep this file deliberately small.  The shell script owns the pipeline; this
module handles JSON/config emission, FITS bookkeeping, map combination, and QC
plotting that would be painful in shell.
"""

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
        "PN": ["PN"], "EPN": ["PN"],
        "M1": ["M1"], "MOS1": ["M1"], "EMOS1": ["M1"],
        "M2": ["M2"], "MOS2": ["M2"], "EMOS2": ["M2"],
        "ALL": ["PN", "M1", "M2"], "EPIC": ["PN", "M1", "M2"],
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
        "CLEAN_PI_MIN": cfg.get("clean_pi_min", 200),
        "CLEAN_PI_MAX": cfg.get("clean_pi_max", 12000),
        "CLEAN_GTI_ENABLED": yesno(cfg.get("clean_gti_enabled"), True),
        "CLEAN_GTI_RATE_CUT": cfg.get("clean_gti_rate_cut", 4.80001211),
        "CLEAN_GTI_TIMEBIN": cfg.get("clean_gti_timebin", 10),
        "CLEAN_GTI_PI_MIN": cfg.get("clean_gti_pi_min", 7000),
        "CLEAN_GTI_PI_MAX": cfg.get("clean_gti_pi_max", 15000),
        "CLEAN_GTI_PATTERN_MAX": cfg.get("clean_gti_pattern_max", 0),
        "IMAGE_BANDS": cfg.get("image_bands", "1000:200:500;2000:500:1000;3000:1000:2000;4000:2000:4500;5000:4500:12000"),
        "IMAGE_BIN_PHYS": cfg.get("image_bin_phys", 80),
        "IMAGE_PAD_FRAC": cfg.get("image_pad_frac", 0.08),
        "EEXPMAP_ATTREBIN": cfg.get("eexpmap_attrebin", "0.020626481"),
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


def manifest_events(manifest_arg: str, detectors_arg: str = "PN") -> list[Path]:
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

    paths: list[Path] = []
    seen: set[Path] = set()
    for manifest in manifests:
        for line in manifest.read_text(encoding="utf-8").splitlines():
            path = Path(line.strip())
            if path and path not in seen:
                seen.add(path)
                paths.append(path)
    return paths


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
    from astropy.io import fits

    path = Path(path_arg)
    if not path.is_file():
        print(-1)
        return
    with fits.open(path, memmap=True) as hdul:
        hdu = event_hdu(hdul)
        print(0 if hdu.data is None else len(hdu.data))


def flare_qc(lc_dir_arg: str, outdir_arg: str, rate_cut_arg: str) -> None:
    import matplotlib
    import numpy as np
    from astropy.io import fits

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    paths = sorted(Path(lc_dir_arg).glob("*/*flare_lc.fits"))
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
    (outdir / "grid.env").write_text(
        f"X_MIN_PHYS={xmin}\nX_MAX_PHYS={xmax}\n"
        f"Y_MIN_PHYS={ymin}\nY_MAX_PHYS={ymax}\n"
        f"BIN_PHYS={bin_phys}\nNX={nx}\nNY={ny}\n",
        encoding="utf-8",
    )
    (outdir / "grid.json").write_text(
        json.dumps({"x_min": xmin, "x_max": xmax, "y_min": ymin, "y_max": ymax, "bin_phys": bin_phys, "nx": nx, "ny": ny}, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


def event_time_info(path: Path) -> tuple[float, float, str]:
    from astropy.io import fits

    with fits.open(path, memmap=True) as hdul:
        header = event_hdu(hdul).header
        start = header.get("TSTART")
        stop = header.get("TSTOP")
        source = header.get("EXPIDSTR") or path.stem.split("_")[-2]
    if start is None or stop is None:
        die(f"Missing TSTART/TSTOP in {path}")
    return float(start), float(stop), str(source)


def clean_groups(clean_dir_arg: str) -> dict[str, list[tuple[float, float, str, Path]]]:
    clean_dir = Path(clean_dir_arg)
    manifests: dict[str, Path] = {}
    for manifest in sorted(clean_dir.glob("*_clean_files.txt")):
        manifests[manifest.name.removesuffix("_clean_files.txt")] = manifest
    for manifest in sorted((clean_dir / "manifest").glob("*_clean_files.txt")):
        manifests.setdefault(manifest.name.removesuffix("_clean_files.txt"), manifest)

    groups: dict[str, list[tuple[float, float, str, Path]]] = {}
    for inst, manifest in sorted(manifests.items()):
        rows = []
        for line in manifest.read_text(encoding="utf-8").splitlines():
            path = Path(line.strip())
            if path.is_file():
                rows.append((*event_time_info(path), path))
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
        return None if not text else float((Time(text, scale="utc").tt - epoch).sec)

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
    return [(f"P{i:03d}", start, stop, preqid, types, pttime) for i, (start, stop, preqid, types, pttime) in enumerate(pointings, 1)]


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
        + "".join(f"{pid}\t{start:.6f}\t{stop:.6f}\tSCATS\t{request}\t{types}\t{pttime}\n" for pid, start, stop, request, types, pttime in pointings),
        encoding="utf-8",
    )
    Path(slices_arg).write_text(
        "inst\tevent\tbase\tpointing\tstart\tstop\treference_exposure\n"
        + "".join(f"{inst}\t{path}\t{base}\t{pid}\t{start:.6f}\t{stop:.6f}\t{request}\n" for inst, path, base, pid, start, stop, request in slice_records(selected, pointings)),
        encoding="utf-8",
    )


def parse_bands(value: str) -> list[tuple[str, float, float]]:
    try:
        return [(label, float(lo), float(hi)) for label, lo, hi in (entry.split(":", 2) for entry in str(value).split(";") if entry.strip())]
    except ValueError:
        die(f"Bad image_bands value: {value}")


def slice_qc(outdir_arg: str, clean_dir_arg: str, detectors_arg: str, init_dir_arg: str, bands_arg: str) -> None:
    from astropy.io import fits

    outdir = Path(outdir_arg)
    outdir.mkdir(parents=True, exist_ok=True)
    detectors, selected = selected_clean(clean_dir_arg, detectors_arg)
    pointings = odf_pointings(init_dir_arg, clean_span(selected))
    bands = parse_bands(bands_arg)
    band_columns = "\t".join(f"{label}_events" for label, _lo, _hi in bands)

    (outdir / "pre_exposure_pointings.tsv").write_text(
        "pointing\tstart\tstop\tduration_s\tsource\trequest\ttypes\tpttime_utc\n"
        + "".join(f"{pid}\t{start:.6f}\t{stop:.6f}\t{stop - start:.3f}\tSCATS\t{request}\t{types}\t{pttime}\n" for pid, start, stop, request, types, pttime in pointings),
        encoding="utf-8",
    )

    slice_lines = [f"inst\tevent\tpointing\trequest\tstart\tstop\tduration_s\ttotal_events\t{band_columns}\n"]
    coverage_lines = ["inst\tevent\tspan_s\tassigned_s\tunassigned_s\n"]
    by_pointing: dict[str, dict[str, float]] = {
        pid: {"duration": stop - start, "assigned": 0.0, "total": 0.0, **{label: 0.0 for label, _lo, _hi in bands}}
        for pid, start, stop, _request, _types, _pttime in pointings
    }
    requests = {pid: request for pid, _start, _stop, request, _types, _pttime in pointings}

    for inst in detectors:
        for event_start, event_stop, _source, path in selected.get(inst, []):
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
            coverage_lines.append(f"{inst}\t{path.name}\t{event_stop - event_start:.3f}\t{assigned:.3f}\t{event_stop - event_start - assigned:.3f}\n")

    by_lines = [f"pointing\trequest\tduration_s\tevent_overlap_s\ttotal_events\t{band_columns}\n"] + [
        f"{pid}\t{requests[pid]}\t{values['duration']:.3f}\t{values['assigned']:.3f}\t{int(values['total'])}\t"
        + "\t".join(str(int(values[label])) for label, _lo, _hi in bands) + "\n"
        for pid, values in by_pointing.items()
    ]
    durations = sorted(stop - start for _pid, start, stop, _request, _types, _pttime in pointings)
    median = durations[len(durations) // 2]
    summary = [f"pointings {len(pointings)}", f"slices {len(slice_lines) - 1}", f"detectors {' '.join(detectors)}", f"median_pointing_s {median:.3f}", f"longest_pointing_s {max(durations):.3f}"]
    summary += [f"warning {pid} has no selected-detector clean-event overlap" for pid, values in by_pointing.items() if values["assigned"] == 0]
    summary += [f"warning {pid} duration is >3x median" for pid, values in by_pointing.items() if values["duration"] > 3.0 * median]

    (outdir / "pre_exposure_slices.tsv").write_text("".join(slice_lines), encoding="utf-8")
    (outdir / "pre_exposure_event_coverage.tsv").write_text("".join(coverage_lines), encoding="utf-8")
    (outdir / "pre_exposure_by_pointing.tsv").write_text("".join(by_lines), encoding="utf-8")
    (outdir / "pre_exposure_summary.txt").write_text("\n".join([*summary, ""]), encoding="utf-8")


def combine_maps(out_counts_arg: str, out_exposure_arg: str, out_rate_arg: str, inputs: list[str]) -> None:
    import numpy as np
    from astropy.io import fits

    if len(inputs) % 2:
        die("combine-maps expects count/exposure pairs")
    pairs = [(Path(inputs[i]), Path(inputs[i + 1])) for i in range(0, len(inputs), 2)]
    if not pairs:
        die("No maps to combine")

    counts_sum = exposure_sum = header = None
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

    if not inputs:
        die("combine-exposures expects at least one exposure map")
    exposure_sum = header = None
    for item in inputs:
        exposure, next_header = fits.getdata(item, header=True)
        exposure = np.asarray(exposure, dtype=float)
        header = next_header if header is None else header
        exposure_sum = exposure if exposure_sum is None else exposure_sum + exposure
    fits.PrimaryHDU(exposure_sum, header).writeto(out_exposure_arg, overwrite=True)


def mosaic_extents(shape: tuple[int, int], grid_arg: str = "", manifest_arg: str = "", detectors_arg: str = "PN"):
    ny, nx = shape
    image_extent = (0.0, float(nx), 0.0, float(ny))
    plot_extent = image_extent
    if grid_arg and manifest_arg:
        grid = json.loads(Path(grid_arg).read_text(encoding="utf-8"))
        image_extent = (grid["x_min"], grid["x_max"], grid["y_min"], grid["y_max"])
        _bounds, plot_extent, _nx, _ny = event_bounds(manifest_events(manifest_arg, detectors_arg), image_extent)
    return image_extent, plot_extent


def fits_png(fits_arg: str, png_arg: str, cmap_arg: str = "magma", scale_arg: str = "linear", grid_arg: str = "", manifest_arg: str = "", detectors_arg: str = "PN") -> None:
    import matplotlib
    import numpy as np
    from astropy.io import fits
    from matplotlib.colors import LogNorm

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    data = np.asarray(fits.getdata(fits_arg), dtype=float)
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
        vmax = max(float(np.percentile(positive, 99.7)), 1.0)

    fig, ax = plt.subplots(figsize=(8, max(3, 8 * ratio)), constrained_layout=True)
    ax.imshow(image, origin="lower", extent=image_extent, cmap=cmap, norm=norm, vmin=vmin, vmax=vmax, interpolation="nearest")
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
    data_bounds, extent, nx, ny = event_bounds(events)

    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)
    summary = [
        f"detectors {' '.join(detectors)}",
        f"event_lists {len(events)}",
        f"data_bounds {data_bounds}",
        f"plot_extent {extent}",
        f"plot_pixels {nx} {ny}",
    ]
    for tag, (lo, hi) in {"soft_lt1kev": (None, 1000), "hard_gt1kev": (1000, None)}.items():
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
                hist, _xe, _ye = np.histogram2d(x[good], y[good], bins=(nx, ny), range=((extent[0], extent[1]), (extent[2], extent[3])))
                image += hist.T
                total += int(good.sum())
            hdul.close()
        positive = image[image > 0]
        vmax = max(float(np.percentile(positive, 99.7)), 1.0) if positive.size else 1.0
        fig, ax = plt.subplots(figsize=(8, 8 * ny / nx), constrained_layout=True)
        ax.imshow(image, origin="lower", extent=extent, cmap="magma", norm=LogNorm(vmin=1, vmax=vmax) if positive.size else None, interpolation="nearest")
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
    elif cmd == "rewrite-path" and len(args) == 2:
        rewrite_path(*args)
    elif cmd == "fits-rows" and len(args) == 1:
        fits_rows(*args)
    elif cmd == "event-mosaic" and len(args) in {2, 3}:
        event_mosaic(args[0], args[1], args[2] if len(args) == 3 else "PN")
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
    elif cmd == "fits-png" and 2 <= len(args) <= 7:
        fits_png(*args)
    else:
        die(
            "Usage: tools.py shell CONFIG | rewrite-path SUM.SAS PATH | fits-rows FITS | "
            "event-mosaic MANIFEST OUTDIR [DETECTORS] | flare-qc LC_DIR OUTDIR RATE_CUT | "
            "exposure-grid OUTDIR BIN PAD EVENTS... | exposure-slices SLICES POINTINGS CLEAN_DIR DETECTORS INIT_DIR | "
            "slice-qc OUTDIR CLEAN_DIR DETECTORS INIT_DIR BANDS | combine-maps OUT_COUNTS OUT_EXP OUT_RATE COUNTS EXP... | "
            "combine-exposures OUT_EXP EXP... | fits-png FITS PNG [CMAP] [SCALE] [GRID_JSON MANIFEST DETECTORS]"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
