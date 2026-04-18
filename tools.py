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
    if isinstance(text, list):
        tokens = [str(item) for item in text]
    else:
        tokens = str(text).replace(",", " ").replace(";", " ").split()

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
        try:
            expanded = aliases[token.upper()]
        except KeyError:
            die(f"Unknown detector in config: {token}")
        for detector in expanded:
            if detector not in detectors:
                detectors.append(detector)
    return detectors or ["PN"]


def emit(name: str, value: Any) -> None:
    print(f"{name}={quote(value)}")


def shell(path: str) -> None:
    cfg = load_config(path)
    emit("WORKDIR", str(cfg["workdir"]).rstrip("/"))
    emit("ODFDIR", str(cfg["odfdir"]).rstrip("/"))
    emit("SHORTLINK_NAME", cfg.get("shortlink_name", ""))
    emit("KEEP_SHORTLINK", yesno(cfg.get("keep_shortlink"), True))
    emit("SHORTLINK_MAX_PATH", cfg.get("shortlink_max_path", 60))
    emit("SAS_SETUP_SCRIPT", cfg.get("sas_setup_script", ""))
    emit("SAS_CCFPATH_CONFIG", cfg.get("sas_ccfpath", ""))
    emit("SAS_VERBOSITY_CONFIG", cfg.get("sas_verbosity", ""))
    emit("DETECTORS", " ".join(normalize_detectors(cfg.get("detectors"))))
    emit("AHF_INPUT", cfg.get("ahf_input", ""))
    emit("CLEAN_PI_MIN", cfg.get("clean_pi_min", 200))
    emit("CLEAN_PI_MAX", cfg.get("clean_pi_max", 12000))
    emit("PREMERGE_IMAGE_SIZE_DEG", cfg.get("premerge_image_size_deg", 1.5))
    emit("IMAGE_BANDS", cfg.get("image_bands", "soft:200:999;hard:1001:12000"))
    emit("IMAGE_BIN_PHYS", cfg.get("image_bin_phys", 80))
    emit("IMAGE_PAD_FRAC", cfg.get("image_pad_frac", 0.08))
    emit("EEXPMAP_ATTREBIN", cfg.get("eexpmap_attrebin", "0.020626481"))
    emit("BACKGROUND_COUNTS_QUANTILE", cfg.get("background_counts_quantile", 0.10))
    emit("BACKGROUND_EXCLUDE_RADIUS_FRAC", cfg.get("background_exclude_radius_frac", 0.50))
    emit("BACKGROUND_MIN_PIXELS", cfg.get("background_min_pixels", 100))
    emit("LINK_ODF_CONSTITUENTS", yesno(cfg.get("link_odf_constituents"), True))
    emit("SKIP_ODF_LINK_PATTERNS", "|".join(cfg.get("skip_odf_link_patterns", [])))


def rewrite_path(summary: str, odf_path: str) -> None:
    path = odf_path.rstrip("/") + "/"
    summary_path = Path(summary)
    lines = summary_path.read_text(encoding="utf-8", errors="surrogateescape").splitlines()
    found = False
    for idx, line in enumerate(lines):
        if line.startswith("PATH "):
            lines[idx] = f"PATH {path}"
            found = True
            break
    if not found:
        lines.append(f"PATH {path}")
    summary_path.write_text("\n".join(lines) + "\n", encoding="utf-8", errors="surrogateescape")


def manifest_paths(source: Path, detectors: list[str]) -> list[Path]:
    if source.is_file():
        return [source]
    manifests: list[Path] = []
    for detector in detectors:
        for name in (f"{detector}_clean_files.txt", f"{detector}_raw.txt"):
            path = source / name
            if path.is_file():
                manifests.append(path)
                break
    return manifests


def manifest_events(manifest_arg: str, detectors_arg: str = "PN") -> list[Path]:
    paths: list[Path] = []
    seen: set[Path] = set()
    source = Path(manifest_arg)
    manifests = manifest_paths(source, normalize_detectors(detectors_arg))
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


def event_columns(path: Path):
    try:
        from astropy.io import fits
    except ImportError:
        die("FITS helpers require astropy")

    hdul = fits.open(path, memmap=True)
    hdu = hdul["EVENTS"] if "EVENTS" in hdul else hdul[1]
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
    groups: dict[str, list[tuple[float, float, str, Path]]] = {}
    for manifest in sorted(Path(clean_dir_arg).glob("*_clean_files.txt")):
        inst = manifest.name.removesuffix("_clean_files.txt")
        rows = []
        for line in manifest.read_text(encoding="utf-8").splitlines():
            path = Path(line.strip())
            if path.is_file():
                start, stop, source = event_time_info(path)
                rows.append((start, stop, source, path))
        if rows:
            groups[inst] = sorted(rows)
    return groups


def odf_groups(init_dir_arg: str) -> dict[str, list[tuple[float, float, str, Path]]]:
    from astropy.io import fits
    from astropy.time import Time

    init_dir = Path(init_dir_arg)
    epoch = Time("1998-01-01T00:00:00", scale="tt")
    groups: dict[str, list[tuple[float, float, str, Path]]] = {}
    for path in sorted(init_dir.glob("*.FIT")):
        match = re.search(r"_((?:M[12])|PN)([SU]\d{3})00AUX\.FIT$", path.name)
        if not match:
            continue
        inst, source = match.groups()
        with fits.open(path, memmap=True) as hdul:
            header = hdul[1].header
            if "DATE-OBS" not in header or "DATE-END" not in header:
                continue
            start = (Time(header["DATE-OBS"], scale="utc").tt - epoch).sec
            stop = (Time(header["DATE-END"], scale="utc").tt - epoch).sec
        groups.setdefault(inst, []).append((float(start), float(stop), source, path))
    return {inst: sorted(rows) for inst, rows in groups.items()}


def exposure_slices(
    slices_arg: str,
    pointings_arg: str,
    clean_dir_arg: str,
    detectors_arg: str,
    init_dir_arg: str = "",
) -> None:
    groups = clean_groups(clean_dir_arg)
    ref_groups = (odf_groups(init_dir_arg) | groups) if init_dir_arg else groups
    if not groups:
        die(f"No clean event manifests found in {clean_dir_arg}")

    order = {"M2": 0, "M1": 1, "PN": 2}
    ref_inst = sorted(ref_groups, key=lambda item: (-len(ref_groups[item]), order.get(item, 9), item))[0]
    pointings = [(f"P{i:03d}", *row[:3]) for i, row in enumerate(ref_groups[ref_inst], 1)]

    Path(pointings_arg).write_text(
        "pointing\tstart\tstop\treference_detector\treference_exposure\n"
        + "".join(f"{pid}\t{start:.6f}\t{stop:.6f}\t{ref_inst}\t{source}\n" for pid, start, stop, source in pointings),
        encoding="utf-8",
    )

    lines = ["inst\tevent\tbase\tpointing\tstart\tstop\treference_exposure\n"]
    for inst in normalize_detectors(detectors_arg):
        for event_start, event_stop, _event_source, path in groups.get(inst, []):
            stem = path.name.removesuffix(".fits")
            for pid, point_start, point_stop, point_source in pointings:
                start = max(event_start, point_start)
                stop = min(event_stop, point_stop)
                if stop - start > 1.0:
                    base = f"{stem}_{pid}_{point_source}"
                    lines.append(f"{inst}\t{path}\t{base}\t{pid}\t{start:.6f}\t{stop:.6f}\t{point_source}\n")
    Path(slices_arg).write_text("".join(lines), encoding="utf-8")


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


def background_map(args: list[str]) -> None:
    import numpy as np
    from astropy.io import fits

    if len(args) != 11:
        die("background-map expects COUNTS EXPOSURE BACKGROUND NET COUNTS_QUANTILE MINPIX EXCLUDE_FRAC SUMMARY BAND INST SOURCE")
    counts_path, exp_path, out_bkg, out_net, quantile, minpix, exclude_frac, summary, band, inst, source = args
    counts, header = fits.getdata(counts_path, header=True)
    exposure = fits.getdata(exp_path).astype(float)
    if counts.shape != exposure.shape:
        die(f"Shape mismatch: {counts_path} vs {exp_path}")

    counts = np.asarray(counts, dtype=float)
    q = float(quantile)
    if not 0.0 <= q <= 1.0:
        die(f"Background quantile must be between 0 and 1: {q}")

    status = "ok"
    footprint_pixels = estimate_pixels = positive_pixels = level = 0.0
    cx = cy = detector_radius = exclude_radius = np.nan
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
        exclude_radius = float(exclude_frac) * detector_radius
        estimate = footprint & (np.hypot(xx - cx, yy - cy) >= exclude_radius)
        estimate_pixels = int(estimate.sum())
        if estimate_pixels < int(minpix):
            status = "low_estimate"
        else:
            values = counts[estimate]
            values = values[np.isfinite(values) & (values > 0)]
            positive_pixels = int(values.size)
            if positive_pixels < int(minpix):
                status = "low_counts"
            else:
                level = float(np.quantile(values, q))
                background = np.where(footprint, level, 0.0)
                net = np.where(footprint, counts - background, np.nan)

    header["BUNIT"] = "count"
    fits.PrimaryHDU(background, header).writeto(out_bkg, overwrite=True)
    fits.PrimaryHDU(net, header).writeto(out_net, overwrite=True)

    row = (
        f"{band}\t{inst}\t{source}\t{status}\t{footprint_pixels}\t{estimate_pixels}\t{positive_pixels}\t"
        f"{q:.4g}\t{level:.8g}\t{cx:.3f}\t{cy:.3f}\t"
        f"{detector_radius:.3f}\t{exclude_radius:.3f}\n"
    )
    path = Path(summary)
    if not path.exists():
        path.write_text(
            "band\tinst\tsource\tstatus\tfootprint_pixels\testimate_pixels\tpositive_pixels\tcounts_quantile\t"
            "background_per_pixel\tcenter_x_px\tcenter_y_px\t"
            "detector_radius_px\texclude_radius_px\n",
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


def fits_png(
    fits_arg: str,
    png_arg: str,
    cmap_arg: str = "viridis",
    scale_arg: str = "linear",
    grid_arg: str = "",
    manifest_arg: str = "",
    detectors_arg: str = "PN",
) -> None:
    import matplotlib
    import numpy as np
    from astropy.io import fits
    from matplotlib.colors import LogNorm

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    data = np.asarray(fits.getdata(fits_arg), dtype=float)
    ny, nx = data.shape
    image_extent = (0.0, float(nx), 0.0, float(ny))
    plot_extent = image_extent
    if grid_arg and manifest_arg:
        grid = json.loads(Path(grid_arg).read_text(encoding="utf-8"))
        image_extent = (grid["x_min"], grid["x_max"], grid["y_min"], grid["y_max"])
        _bounds, plot_extent, _nx, _ny = event_bounds(manifest_events(manifest_arg, detectors_arg), image_extent)

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
    if scale_arg == "signed":
        finite_values = data[finite]
        span = float(np.nanpercentile(np.abs(finite_values), 99.0)) if finite_values.size else 1.0
        span = span if span > 0 else 1.0
        image = np.ma.masked_where(~finite, data)
        norm = None
        vmin = -span
        vmax = span

    ratio = (plot_extent[3] - plot_extent[2]) / max(plot_extent[1] - plot_extent[0], 1.0)
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
    if not events:
        die(f"No event lists found in {manifest_dir}")

    data_bounds, extent, nx, ny = event_bounds(events)

    bands = {"soft_lt1kev": (None, 1000, "magma"), "hard_gt1kev": (1000, None, "viridis")}
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
    if len(argv) == 3 and argv[1] == "shell":
        shell(argv[2])
        return 0
    if len(argv) == 4 and argv[1] == "rewrite-path":
        rewrite_path(argv[2], argv[3])
        return 0
    if len(argv) in {4, 5} and argv[1] in {"event-mosaic", "repro-mosaic"}:
        event_mosaic(argv[2], argv[3], argv[4] if len(argv) == 5 else "PN")
        return 0
    if len(argv) == 3 and argv[1] == "fits-rows":
        fits_rows(argv[2])
        return 0
    if len(argv) >= 6 and argv[1] == "exposure-grid":
        exposure_grid(argv[2], argv[3], argv[4], argv[5:])
        return 0
    if len(argv) in {6, 7} and argv[1] == "exposure-slices":
        exposure_slices(argv[2], argv[3], argv[4], argv[5], argv[6] if len(argv) == 7 else "")
        return 0
    if len(argv) >= 7 and argv[1] == "combine-maps":
        combine_maps(argv[2], argv[3], argv[4], argv[5:])
        return 0
    if len(argv) == 13 and argv[1] == "background-map":
        background_map(argv[2:])
        return 0
    if len(argv) >= 6 and argv[1] == "combine-background":
        combine_background(argv[2:])
        return 0
    if len(argv) == 5 and argv[1] == "net-rate":
        net_rate(argv[2], argv[3], argv[4])
        return 0
    if 4 <= len(argv) <= 9 and argv[1] == "fits-png":
        fits_png(
            argv[2],
            argv[3],
            argv[4] if len(argv) >= 5 else "viridis",
            argv[5] if len(argv) >= 6 else "linear",
            argv[6] if len(argv) >= 7 else "",
            argv[7] if len(argv) >= 8 else "",
            argv[8] if len(argv) >= 9 else "PN",
        )
        return 0
    die("Usage: tools.py shell CONFIG | rewrite-path SUM.SAS PATH | event-mosaic MANIFEST OUTDIR [DETECTORS] | fits-rows FITS | exposure-grid OUTDIR BIN PAD EVENTS... | exposure-slices SLICES_TSV POINTINGS_TSV CLEAN_DIR DETECTORS [INIT_DIR] | combine-maps OUT_COUNTS OUT_EXP OUT_RATE COUNTS EXP... | background-map COUNTS EXPOSURE BACKGROUND NET COUNTS_QUANTILE MINPIX EXCLUDE_FRAC SUMMARY BAND INST SOURCE | combine-background OUT_BACKGROUND OUT_NET BACKGROUND NET... | net-rate NET EXPOSURE OUT | fits-png FITS PNG [CMAP] [linear|log|signed] [GRID_JSON MANIFEST DETECTORS]")


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
