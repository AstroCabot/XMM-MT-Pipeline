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


def fits_png(fits_arg: str, png_arg: str, cmap_arg: str = "viridis", scale_arg: str = "linear") -> None:
    import matplotlib
    import numpy as np
    from astropy.io import fits
    from matplotlib.colors import LogNorm

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    data = np.asarray(fits.getdata(fits_arg), dtype=float)
    finite = np.isfinite(data)
    positive = data[finite & (data > 0)]
    cmap = plt.get_cmap(cmap_arg).copy()
    cmap.set_bad("black")
    cmap.set_under("black")

    norm = None
    vmin = vmax = None
    image = np.ma.masked_where(~finite | (data <= 0), data)
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

    ny, nx = data.shape
    fig, ax = plt.subplots(figsize=(8, 8 * ny / nx), constrained_layout=True)
    ax.imshow(image, origin="lower", cmap=cmap, norm=norm, vmin=vmin, vmax=vmax, interpolation="nearest")
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

    xmin = ymin = np.inf
    xmax = ymax = -np.inf
    for path in events:
        hdul, x, y, _pi = event_columns(path)
        good = valid_sky(x, y)
        if np.any(good):
            xg = x[good]
            yg = y[good]
            bounds = (float(xg.min()), float(xg.max()), float(yg.min()), float(yg.max()))
            xmin = min(xmin, bounds[0])
            xmax = max(xmax, bounds[1])
            ymin = min(ymin, bounds[2])
            ymax = max(ymax, bounds[3])
        hdul.close()

    if not np.isfinite([xmin, xmax, ymin, ymax]).all() or xmin == xmax or ymin == ymax:
        die("Could not determine a valid mosaic extent from repro event lists")

    data_bounds = (xmin, xmax, ymin, ymax)
    pad_x = max(0.08 * (xmax - xmin), 500.0)
    pad_y = max(0.08 * (ymax - ymin), 500.0)
    extent = (xmin - pad_x, xmax + pad_x, ymin - pad_y, ymax + pad_y)
    nx = 900
    ny = max(200, int(round(nx * (extent[3] - extent[2]) / (extent[1] - extent[0]))))
    ny = min(ny, 1200)

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
    if len(argv) >= 7 and argv[1] == "combine-maps":
        combine_maps(argv[2], argv[3], argv[4], argv[5:])
        return 0
    if len(argv) in {4, 5, 6} and argv[1] == "fits-png":
        fits_png(argv[2], argv[3], argv[4] if len(argv) >= 5 else "viridis", argv[5] if len(argv) >= 6 else "linear")
        return 0
    die("Usage: tools.py shell CONFIG | rewrite-path SUM.SAS PATH | event-mosaic MANIFEST OUTDIR [DETECTORS] | fits-rows FITS | exposure-grid OUTDIR BIN PAD EVENTS... | combine-maps OUT_COUNTS OUT_EXP OUT_RATE COUNTS EXP... | fits-png FITS PNG [CMAP] [linear|log]")


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
