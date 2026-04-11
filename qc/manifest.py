#!/usr/bin/env python3
"""Generate the QC stage manifest and QC index.

Subcommands
-----------
manifest  – scan pipeline outputs and write manifest.png/md/json
index     – collect all QC PNGs into STAGE_QC.md
repro     – report per-instrument calibrated event-list candidates
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits

# Sibling imports — tolerate both package and standalone invocation.
try:
    from qc.utils import load_shell_env, read_list, fits_rows, nonempty
    from qc import plot_utils as qpu
except ImportError:
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    from utils import load_shell_env, read_list, fits_rows, nonempty
    import plot_utils as qpu


def stage_status(ok: bool, warn: bool = False) -> str:
    if ok and (not warn):
        return "ok"
    if ok and warn:
        return "warn"
    return "missing"


# ------------------------------------------------------------------
# manifest subcommand
# ------------------------------------------------------------------


def cmd_manifest(args: argparse.Namespace) -> int:
    env = load_shell_env(args.config)
    workdir = Path(env["WORKDIR"])
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    stages: dict[str, dict[str, object]] = {}

    # --- init ---
    init_ok = nonempty(workdir / "sas_setup.env") and nonempty(
        workdir / "init" / "ccf.cif"
    )
    stages["init"] = {
        "status": stage_status(init_ok),
        "sas_setup_env": str(workdir / "sas_setup.env"),
        "ccf_cif": str(workdir / "init" / "ccf.cif"),
    }

    # --- repro ---
    repro_data: dict[str, dict[str, object]] = {}
    raw_counts = []
    for inst in ("PN", "M1", "M2"):
        files = read_list(workdir / "repro" / "manifest" / f"{inst}_raw.txt")
        raw_counts.append(len(files))
        repro_data[inst] = {"n_raw_eventlists": len(files), "example": files[:3]}
    attitude_ok = nonempty(workdir / "repro" / "atthk.dat") and nonempty(
        workdir / "attitude.env"
    )
    stages["repro"] = {
        "status": stage_status(any(c > 0 for c in raw_counts) and attitude_ok),
        "attitude_ok": attitude_ok,
        "per_instrument": repro_data,
    }

    # --- clean ---
    clean_data: dict[str, dict[str, object]] = {}
    clean_counts = []
    merged_rows = []
    for inst in ("PN", "M1", "M2"):
        files = read_list(workdir / "clean" / f"{inst}_clean_files.txt")
        merged = workdir / "clean" / f"{inst}_clean_merged.fits"
        rows = fits_rows(merged)
        clean_counts.append(len(files))
        merged_rows.append(rows if rows is not None else math.nan)
        clean_data[inst] = {
            "n_clean_eventlists": len(files),
            "merged_eventlist": str(merged),
            "merged_rows": rows,
        }
    clean_warn = any(
        raw_counts[i] != clean_counts[i] and raw_counts[i] > 0 for i in range(3)
    )
    stages["clean"] = {
        "status": stage_status(any(c > 0 for c in clean_counts), warn=clean_warn),
        "per_instrument": clean_data,
    }

    # --- track ---
    track_fits = workdir / "track" / "comet_track.fits"
    seg_csv = workdir / "track" / "motion_segments.csv"
    n_segments = 0
    if seg_csv.exists():
        n_segments = max(0, len(seg_csv.read_text(encoding="utf-8").splitlines()) - 1)
    stages["track"] = {
        "status": stage_status(nonempty(track_fits) and nonempty(seg_csv)),
        "track_fits": str(track_fits),
        "motion_segments_csv": str(seg_csv),
        "n_motion_segments": n_segments,
    }

    # --- detect ---
    pseudoexps = read_list(workdir / "detect" / "stack_eventsets.txt")
    detected_rows = fits_rows(workdir / "detect" / "field_sources_all.fits")
    curated_rows = fits_rows(workdir / "detect" / "field_sources_curated.fits")
    detect_ok = (
        nonempty(workdir / "detect" / "stack_srclist.fits") and len(pseudoexps) > 0
    )
    stages["detect"] = {
        "status": stage_status(detect_ok),
        "n_pseudoexposures": len(pseudoexps),
        "n_detected_sources": detected_rows,
        "n_curated_sources": curated_rows,
    }

    # --- comet ---
    comet_ok = nonempty(workdir / "comet" / "moved_sas_setup.env")
    comet_data: dict[str, object] = {
        "status": stage_status(comet_ok),
        "moved_sas_setup_env": str(workdir / "comet" / "moved_sas_setup.env"),
    }
    for inst in ("PN", "M1", "M2"):
        rows = fits_rows(workdir / "comet" / f"{inst}_comet.fits")
        comet_data[inst] = {
            "rows": rows,
            "path": str(workdir / "comet" / f"{inst}_comet.fits"),
        }
    stages["comet"] = comet_data

    # --- contam ---
    contam_csv = workdir / "contam" / "contamination_report.csv"
    gti = workdir / "contam" / "science_gti.fits"
    contam_rows = 0
    if contam_csv.exists():
        contam_rows = max(
            0, len(contam_csv.read_text(encoding="utf-8").splitlines()) - 1
        )
    summary_json = workdir / "contam" / "contamination_summary.json"
    summary = _load_json(summary_json)
    stages["contam"] = {
        "status": stage_status(nonempty(gti)),
        "science_gti": str(gti),
        "science_policy": summary.get("science_policy"),
        "n_contam_sources": contam_rows,
    }

    # --- image ---
    image_bands = []
    for entry in env.get("IMAGE_BANDS", "broad:300:2000").split(";"):
        entry = entry.strip()
        if not entry:
            continue
        image_bands.append(entry.split(":")[0].strip())
    image_products = {}
    image_ok = False
    for band in image_bands:
        rate = workdir / "images" / band / "EPIC_rate.fits"
        exp = workdir / "images" / band / "EPIC_exp.fits"
        ok = nonempty(rate) and nonempty(exp)
        image_ok = image_ok or ok
        image_products[band] = {"rate": str(rate), "exp": str(exp), "ok": ok}
    stages["image"] = {"status": stage_status(image_ok), "bands": image_products}

    # --- lcurve ---
    lcurve_ok = False
    lcurve_warn = False
    for inst in ("PN", "M1", "M2"):
        for suffix in ("_corr_abs.fits", "_corr_relonly.fits", "_corr.fits"):
            if nonempty(workdir / "lcurve" / f"{inst}{suffix}"):
                lcurve_ok = True
                break
    stages["lcurve"] = {
        "status": stage_status(lcurve_ok, warn=lcurve_warn),
    }

    # --- spectrum ---
    spec_ok = any(
        nonempty(workdir / "spectra" / f"{inst}_src_spec.fits")
        or nonempty(workdir / "spectra" / f"{inst}_src_spec_grp.fits")
        for inst in ("PN", "M1", "M2")
    )
    stages["spectrum"] = {"status": stage_status(spec_ok)}

    # --- merge ---
    merge_ok = any(
        nonempty(workdir / "final" / f"{inst}_full_instid.fits")
        for inst in ("PN", "M1", "M2")
    )
    stages["merge"] = {"status": stage_status(merge_ok)}

    # ---- Write JSON manifest ----
    payload = {"workdir": str(workdir), "stages": stages}
    (outdir / "manifest.json").write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )

    # ---- Write Markdown manifest ----
    md_lines = [
        "# Stage manifest summary",
        "",
        f"Workdir: `{workdir}`",
        "",
        "| Stage | Status | Key check |",
        "|---|---|---|",
    ]
    for name in [
        "init",
        "repro",
        "clean",
        "track",
        "detect",
        "comet",
        "contam",
        "image",
        "lcurve",
        "spectrum",
        "merge",
    ]:
        st = stages[name]
        status = str(st.get("status", "missing"))
        if name == "repro":
            key = ", ".join(
                f"{inst}:{st['per_instrument'][inst]['n_raw_eventlists']} raw"
                for inst in ("PN", "M1", "M2")
            )
        elif name == "clean":
            key = ", ".join(
                f"{inst}:{st['per_instrument'][inst]['n_clean_eventlists']} clean"
                for inst in ("PN", "M1", "M2")
            )
        elif name == "detect":
            key = (
                f"{st.get('n_pseudoexposures', 0)} pseudoexps, "
                f"{st.get('n_detected_sources', 'n/a')} detected, "
                f"{st.get('n_curated_sources', 'n/a')} curated"
            )
        elif name == "lcurve":
            key = f"combined mode={st.get('combined_mode', 'n/a')}"
        elif name == "track":
            key = f"{st.get('n_motion_segments', 0)} motion segments"
        else:
            key = "outputs present" if status != "missing" else "missing"
        md_lines.append(f"| {name} | {status} | {key} |")
    (outdir / "manifest.md").write_text("\n".join(md_lines) + "\n", encoding="utf-8")

    # ---- Write PNG summary ----
    labels = [
        "PN raw",
        "M1 raw",
        "M2 raw",
        "PN clean",
        "M1 clean",
        "M2 clean",
        "pseudoexp",
        "curated src",
    ]
    values = [
        raw_counts[0],
        raw_counts[1],
        raw_counts[2],
        clean_counts[0],
        clean_counts[1],
        clean_counts[2],
        len(pseudoexps),
        int(curated_rows or 0),
    ]

    fig, ax = plt.subplots(figsize=(10.5, 4.6), constrained_layout=True)
    x = np.arange(len(labels))
    ax.bar(x, values)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=25, ha="right", fontsize=9)
    ax.set_ylabel("count")
    ax.set_title("Pipeline manifest overview")
    for i, v in enumerate(values):
        ax.text(
            i,
            v + max(values + [1]) * 0.02,
            str(v),
            ha="center",
            va="bottom",
            fontsize=8,
        )
    txt = [
        f"lcurve mode : {combined_mode or 'n/a'}",
        f"motion segs : {n_segments}",
        f"merged rows : PN={clean_data['PN']['merged_rows']} "
        f"M1={clean_data['M1']['merged_rows']} "
        f"M2={clean_data['M2']['merged_rows']}",
    ]
    ax.text(
        0.99,
        0.98,
        "\n".join(txt),
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=8,
        family="monospace",
        bbox=dict(boxstyle="round", facecolor="white", alpha=0.7, edgecolor="none"),
    )
    qpu.savefig(fig, outdir / "manifest.png")
    print(outdir / "manifest.png")
    return 0


# ------------------------------------------------------------------
# index subcommand
# ------------------------------------------------------------------

# Map from check keys → new semantic filenames used in the qc/ directory.
_ARTIFACT_MAP = {
    "init": ["manifest.md"],
    "repro": ["manifest.md", "manifest.png"],
    "clean": ["clean_gti.png", "manifest.md", "manifest.png"],
    "track": ["track_ephemeris.png", "qc_report.md"],
    "detect": ["detect_sources.png", "manifest.md"],
    "contam": ["contam_timeline.png", "qc_report.md"],
    "lcurve": ["lcurve.png", "region_support_PN.json", "qc_report.md"],
    "spectrum": ["spectrum.png", "qc_report.md"],
    "merge": ["manifest.md"],
}


def cmd_index(args: argparse.Namespace) -> int:
    env = load_shell_env(args.config)
    workdir = Path(env["WORKDIR"])
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    manifest = _load_json(outdir / "manifest.json")
    stages = manifest.get("stages", {})

    band = (
        env.get("QC_SPOTCHECK_BAND", "")
        or env.get("IMAGE_BANDS", "broad:300:2000").split(";")[0].split(":")[0]
    )

    # Build the per-stage artifact map dynamically for entries that depend on band.
    artifact_map = dict(_ARTIFACT_MAP)
    artifact_map["comet"] = [f"spotcheck_{band}.png", "manifest.md"]
    artifact_map["image"] = [f"image_{band}.png", "qc_report.md"]

    report_md = outdir / "qc_report.md"
    region_reports = {
        inst: outdir / f"region_support_{inst}.json"
        for inst in ("EPIC", "PN", "M1", "M2")
    }

    lines = [
        "# Stage-by-stage QC index",
        "",
        f"Workdir: `{workdir}`",
        "",
        "This file is intended to be the one place to inspect after a run.  "
        "Each section points to the main image and/or text artifact for that stage.",
        "",
        "## Summary table",
        "",
        "| Stage | Status | Main artifact(s) |",
        "|---|---|---|",
    ]

    for stage in [
        "init",
        "repro",
        "clean",
        "track",
        "detect",
        "comet",
        "contam",
        "image",
        "lcurve",
        "spectrum",
        "merge",
    ]:
        status = str(stages.get(stage, {}).get("status", "missing"))
        arts = ", ".join(f"`{a}`" for a in artifact_map.get(stage, []))
        lines.append(f"| {stage} | {status} | {arts} |")

    lines += ["", "## What to look for", ""]
    lines += [
        "### init / repro / clean",
        "Use `manifest.md` and `manifest.png` to confirm that all raw event lists "
        "were found, cleaned, and merged per instrument.",
        "",
        "### track",
        "Inspect `track_ephemeris.png`.  The track should be smooth, and the "
        "motion-segmentation curve should not show pathological jumps.",
        "",
        "### detect",
        "Inspect `detect_sources.png`.  The visible portion of the comet track "
        "should lie over the detect mosaics.",
        "",
        "### comet",
        f"Inspect the spot-check montage `spotcheck_{band}.png`.  In a good comet "
        "frame the comet stays compact while stationary field sources trail.",
        "",
        "### contam",
        "Inspect `contam_timeline.png` and the contamination section in `qc_report.md`.",
        "",
        "### image",
        f"Inspect `image_{band}.png`.  The counts, exposure, rate, and masked-rate "
        "panels should all be finite in the comet-centered footprint.",
        "",
        "### lcurve",
        "Inspect `lcurve.png`, `lcurve/lc_summary.tsv`, and the "
        "`region_support_*.json` files.",
        "",
        "### spectrum",
        "Inspect `spectrum.png`.  The combined source spectrum should sit above "
        "the scaled background over at least part of the configured band.",
        "",
    ]

    lines.append("## Region-support quick notes")
    lines.append("")
    for inst, path in region_reports.items():
        if not path.exists():
            continue
        try:
            data = json.loads(path.read_text(encoding="utf-8"))
        except Exception:
            continue
        src_med = data.get("source", {}).get("median")
        bkg_med = data.get("background", {}).get("median")
        src_zero = data.get("source", {}).get("zero_fraction")
        bkg_zero = data.get("background", {}).get("zero_fraction")
        lines.append(
            f"- **{inst}**: source median exp={src_med}, background median exp={bkg_med}, "
            f"source zero-fraction={src_zero}, background zero-fraction={bkg_zero}"
        )
    lines.append("")

    if report_md.exists():
        lines.append("## Primary QC report")
        lines.append("")
        lines.append(
            "The companion `qc_report.md` contains the detailed stage-level "
            "measurements produced by `checks.py`."
        )
        lines.append("")

    (outdir / "STAGE_QC.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(outdir / "STAGE_QC.md")
    return 0


# ------------------------------------------------------------------
# repro subcommand
# ------------------------------------------------------------------


def _find_candidates(repro: Path, patterns: list[str]) -> list[str]:
    out: list[str] = []
    for patt in patterns:
        out.extend(str(p) for p in sorted(repro.rglob(patt)) if p.is_file())
    seen: set[str] = set()
    uniq: list[str] = []
    for p in out:
        if p not in seen:
            uniq.append(p)
            seen.add(p)
    return uniq


def cmd_repro(args: argparse.Namespace) -> int:
    env = load_shell_env(args.config)
    workdir = Path(env["WORKDIR"])
    repro = workdir / "repro"
    data = {
        "PN": _find_candidates(
            repro,
            [
                "*EPN*ImagingEvts.ds",
                "*EPN*ImagingEvts*.FIT*",
                "*PIEVLI*.FIT*",
                "*EPN*EVLI*.FIT*",
            ],
        ),
        "M1": _find_candidates(
            repro,
            [
                "*EMOS1*ImagingEvts.ds",
                "*EMOS1*ImagingEvts*.FIT*",
                "*M1EVLI*.FIT*",
                "*EMOS1*EVLI*.FIT*",
            ],
        ),
        "M2": _find_candidates(
            repro,
            [
                "*EMOS2*ImagingEvts.ds",
                "*EMOS2*ImagingEvts*.FIT*",
                "*M2EVLI*.FIT*",
                "*EMOS2*EVLI*.FIT*",
            ],
        ),
    }
    print(
        json.dumps(
            {k: {"count": len(v), "files": v} for k, v in data.items()}, indent=2
        )
    )
    return 0


# ------------------------------------------------------------------
# helpers
# ------------------------------------------------------------------


def _load_json(path: Path) -> dict:
    if not path.exists():
        return {}
    return json.loads(path.read_text(encoding="utf-8"))


# ------------------------------------------------------------------
# CLI
# ------------------------------------------------------------------


def main() -> int:
    ap = argparse.ArgumentParser(
        description="Manifest, QC index, and repro diagnostics for the XMM comet package."
    )
    sub = ap.add_subparsers(dest="command")
    sub.required = True

    p_manifest = sub.add_parser(
        "manifest", help="Write stage-level manifest JSON/MD/PNG"
    )
    p_manifest.add_argument("--config", required=True)
    p_manifest.add_argument("--outdir", required=True)

    p_index = sub.add_parser("index", help="Build stage-by-stage QC navigation index")
    p_index.add_argument("--config", required=True)
    p_index.add_argument("--outdir", required=True)

    p_repro = sub.add_parser(
        "repro", help="Report per-instrument calibrated event-list candidates"
    )
    p_repro.add_argument("--config", required=True)

    args = ap.parse_args()
    if args.command == "manifest":
        return cmd_manifest(args)
    elif args.command == "index":
        return cmd_index(args)
    elif args.command == "repro":
        return cmd_repro(args)
    return 1


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(1)
