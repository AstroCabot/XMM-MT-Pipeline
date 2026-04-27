#!/usr/bin/env python3
"""Helper toolkit for reduction_v10/pipeline.sh.

Subcommands:
  shell CONFIG                                Emit shell variables for pipeline.sh.
  rewrite-sum-path SUM.SAS NEW_PATH           Rewrite the PATH line of an ODF SUM.SAS.
  repro-manifest REPRODIR DETECTOR OUTPATH    Write a raw event-list manifest.
  frames-split REPRODIR FRAMESDIR ATTHK DETECTORS THRESHOLD_AMIN MIN_DURATION_S
                                              Split multi-pointing repro events.
  clean-band-table CONFIG                     TSV of per-band/per-instrument filter expressions.
  clean-manifest EVENTS_DIR OUTPATH           List *_clean.fits paths into a manifest file.
  fits-rows PATH                              Print EVENTS row count (-1 if missing).
  qc-init INITDIR OUTDIR [SAS_SETUP_ENV]      Build init QC products.
  qc-repro REPRODIR DETECTORS OUTDIR SPLIT_EV [ATTITUDE_ENV]
                                              Build repro QC + loE/hiE mosaics.
  qc-frames FRAMESDIR DETECTORS OUTDIR SPLIT_EV
                                              Build frames QC + mosaics + timeline.
  qc-clean CLEANDIR FRAMESDIR DETECTORS SPLIT_EV OUTDIR
                                              loE/hiE QC mosaics + flare lightcurves.
  cut-run CLEANDIR FRAMESDIR OUTDIR CONFIG    Build cut PI-bands with sigma-clip masking.
  qc-cut CUTDIR FRAMESDIR DETECTORS OUTDIR    Cut QC: per-(det,band) mosaics + frame boxes.
  maps-counts CUTDIR MAPSDIR LOGDIR CONFIG    evselect per-frame counts images.
  maps-exposure CUTDIR MAPSDIR ATTHK LOGDIR CONFIG
                                              eexpmap per-frame vignetted exposure maps.
  maps-background MAPSDIR CONFIG              a + b*E fit per frame on off-source pixels.
  maps-corrected MAPSDIR CONFIG               (counts - bkg) / exp_vig per frame.
  qc-maps MAPSDIR FRAMESDIR DETECTORS OUTDIR  Per-(det,band) nobkg-corrected + corrected mosaics.
  build-track FRAMES_MANIFEST TARGET ID_TYPE OBSERVER STEP TRACK_INPUT OUT_FITS OUT_ENV OUT_JSON
                                              Build comet ephemeris over the obs window.
  qc-track TRACK_FITS CUTDIR FRAMESDIR DETECTORS OUTDIR
                                              Cut mosaics with comet trajectory overlay.
"""

from __future__ import annotations

import json
import math
import shlex
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable

import numpy as np
from astropy.io import fits
from astropy.time import Time, TimeDelta
from astropy.wcs import WCS

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.patches import Rectangle

# PN -> EPN, M1 -> EMOS1, M2 -> EMOS2.  proc_task is the SAS reprocessing
# task; events_glob matches that detector's ImagingEvts files.
INSTRUMENTS: dict[str, dict[str, str]] = {
    "PN": {"sas": "EPN",   "proc_task": "epproc", "events_glob": "*EPN*ImagingEvts.ds"},
    "M1": {"sas": "EMOS1", "proc_task": "emproc", "events_glob": "*EMOS1*ImagingEvts.ds"},
    "M2": {"sas": "EMOS2", "proc_task": "emproc", "events_glob": "*EMOS2*ImagingEvts.ds"},
}

DETECTOR_ALIASES = {
    "PN": "PN", "EPN": "PN",
    "M1": "M1", "MOS1": "M1", "EMOS1": "M1",
    "M2": "M2", "MOS2": "M2", "EMOS2": "M2",
}


def die(msg: str) -> None:
    raise SystemExit(msg)


def _family(detector: str) -> str:
    """PN events use the 'pn' filter spec; MOS1/MOS2 share the 'mos' spec."""
    return "pn" if DETECTOR_ALIASES[detector.upper()] == "PN" else "mos"


def _clean_bands(cfg: dict[str, Any]) -> list[dict[str, Any]]:
    """Validate and return cfg['clean_bands'].

    Each band needs: label, description, pi_ranges, and per-family (pn, mos)
    {pattern_min, pattern_max, flag, extra} blocks. Validated at load so any
    typo fails fast before SAS runs the multi-hour evselect/attcalc loop.
    """
    bands = cfg.get("clean_bands")
    if not isinstance(bands, list) or not bands:
        die("config missing non-empty clean_bands list")
    seen: set[str] = set()
    for band in bands:
        label = str(band.get("label", "")).strip()
        if not label or label in seen:
            die(f"clean_bands needs unique non-empty labels (bad: {label!r})")
        seen.add(label)
        ranges = band.get("pi_ranges")
        if not isinstance(ranges, list) or not ranges:
            die(f"band {label}: pi_ranges must be a non-empty list")
        for lo_hi in ranges:
            if not (isinstance(lo_hi, list) and len(lo_hi) == 2):
                die(f"band {label}: bad pi_ranges entry {lo_hi}")
        for fam in ("pn", "mos"):
            spec = band.get(fam)
            if not isinstance(spec, dict):
                die(f"band {label}: missing {fam} filter block")
            if not isinstance(spec.get("extra", []), list):
                die(f"band {label}: {fam}.extra must be a list")
    return bands


def _band_expr(band: dict[str, Any], family: str) -> str:
    """Build the SAS evselect expression for one (band, instrument-family)."""
    spec = band[family]
    pmin = int(spec.get("pattern_min", 0))
    pmax = int(spec.get("pattern_max", pmin))
    pi_terms = "||".join(f"PI in [{int(lo)}:{int(hi)}]" for lo, hi in band["pi_ranges"])
    terms = [f"({pi_terms})",
             (f"(PATTERN=={pmin})" if pmin == pmax
              else f"((PATTERN>={pmin})&&(PATTERN<={pmax}))")]
    flag = str(spec.get("flag", "")).strip()
    if flag:
        terms.append(flag if flag.startswith("#") else f"({flag})")
    for extra in spec.get("extra", []):
        if str(extra).strip():
            terms.append(f"({str(extra).strip()})")
    return "&&".join(terms)


def _cut_bands(cfg: dict[str, Any]) -> list[dict[str, Any]]:
    """Validate cfg['cut_bands']: each needs label, pi_min, pi_max."""
    bands = cfg.get("cut_bands")
    if not isinstance(bands, list) or not bands:
        die("config missing non-empty cut_bands list")
    seen: set[str] = set()
    for b in bands:
        for key in ("label", "pi_min", "pi_max"):
            if key not in b:
                die(f"cut_bands entry missing {key}: {b}")
        label = str(b["label"]).strip()
        if not label or label in seen:
            die(f"cut_bands needs unique non-empty labels (bad: {label!r})")
        seen.add(label)
        if int(b["pi_min"]) > int(b["pi_max"]):
            die(f"cut_bands {label}: pi_min > pi_max")
    return bands


def _flare_lc_expr(cfg: dict[str, Any], detector: str) -> str:
    """Lightcurve filter expression used for the per-frame flare GTI."""
    pmax = int(cfg.get("clean_gti_pattern_max", 0))
    pi_lo = int(cfg.get("clean_gti_pi_min", 7000))
    pi_hi = int(cfg.get("clean_gti_pi_max", 15000))
    qual = "(FLAG==0)" if _family(detector) == "pn" else "#XMMEA_EM"
    return f"(PATTERN<={pmax})&&(PI in [{pi_lo}:{pi_hi}])&&{qual}"


def load_config(path: str) -> dict[str, Any]:
    cfg = json.loads(Path(path).read_text(encoding="utf-8"))
    if not isinstance(cfg, dict):
        die("config must be a JSON object")
    for key in ("workdir", "odfdir"):
        if not cfg.get(key):
            die(f"config missing required key: {key}")
    return cfg


def normalize_detectors(value: Any) -> list[str]:
    tokens = value if isinstance(value, list) else str(value or "").replace(",", " ").split()
    out: list[str] = []
    for token in tokens:
        key = str(token).strip().upper()
        if key not in DETECTOR_ALIASES:
            die(f"unknown detector: {token}")
        det = DETECTOR_ALIASES[key]
        if det not in out:
            out.append(det)
    return out or ["PN"]


def shell(config_path: str) -> None:
    cfg = load_config(config_path)
    detectors = normalize_detectors(cfg.get("detectors", "PN"))
    bands = _clean_bands(cfg)
    band_labels = [str(b["label"]) for b in bands]
    cut_bands = _cut_bands(cfg)
    cut_labels = [str(b["label"]) for b in cut_bands]
    rate_cuts = cfg.get("clean_gti_rate_cut", {})
    if not isinstance(rate_cuts, dict):
        die("clean_gti_rate_cut must be a {detector: rate} mapping")
    for det in detectors:
        if det not in rate_cuts:
            die(f"clean_gti_rate_cut: no rate for detector {det}")

    values: dict[str, Any] = {
        "WORKDIR": str(cfg["workdir"]).rstrip("/"),
        "ODFDIR": str(cfg["odfdir"]).rstrip("/"),
        "SAS_SETUP_SCRIPT": cfg.get("sas_setup_script", ""),
        "SAS_CCFPATH_CONFIG": cfg.get("sas_ccfpath", ""),
        "SAS_VERBOSITY_CONFIG": cfg.get("sas_verbosity", ""),
        "DETECTORS": " ".join(detectors),
        "PROC_TASKS": " ".join(f"{d}:{INSTRUMENTS[d]['proc_task']}" for d in detectors),
        "EVENT_GLOBS": " ".join(f"{d}:{INSTRUMENTS[d]['events_glob']}" for d in detectors),
        "QC_SPLIT_EV": cfg.get("qc_split_ev", 1000),
        "FRAMES_THRESHOLD_AMIN": cfg.get("frames_threshold_amin", 1.0),
        "FRAMES_MIN_DURATION_S": cfg.get("frames_min_duration_s", 1000.0),
        "CLEAN_BAND_LABELS": " ".join(band_labels),
        "CLEAN_GTI_TIMEBIN": cfg.get("clean_gti_timebin", 10),
        "CUT_BAND_LABELS":   " ".join(cut_labels),
        "TARGET_ID":      cfg.get("target_id", ""),
        "TARGET_ID_TYPE": cfg.get("target_id_type", "smallbody"),
        "TRACK_OBSERVER": cfg.get("track_observer", "500@-125989"),
        "TRACK_STEP":     cfg.get("track_step", "1m"),
        "TRACK_INPUT":    cfg.get("track_input", ""),
    }
    for det in detectors:
        values[f"CLEAN_RATE_{det}"]   = rate_cuts[det]
        values[f"FLARE_EXPR_{det}"]   = _flare_lc_expr(cfg, det)
        for band in bands:
            values[f"CLEAN_EXPR_{det}_{band['label']}"] = _band_expr(band, _family(det))

    for name, val in values.items():
        print(f"{name}={shlex.quote(str(val))}")


def rewrite_sum_path(sum_path: str, new_dir: str) -> None:
    """Rewrite the PATH line in an ODF SUM.SAS file.

    odfingest writes an absolute PATH; when we adopt an existing init/
    directory the recorded path is from a previous workdir and SAS will
    refuse to find the ODF constituents.
    """
    target = Path(sum_path)
    new_path = str(new_dir).rstrip("/") + "/"
    text = target.read_text(encoding="utf-8", errors="surrogateescape")
    lines = text.splitlines()
    for i, line in enumerate(lines):
        if line.startswith("PATH "):
            lines[i] = f"PATH {new_path}"
            break
    else:
        lines.append(f"PATH {new_path}")
    target.write_text("\n".join(lines) + "\n", encoding="utf-8", errors="surrogateescape")


def repro_manifest(repro_dir: str, detector: str, out_path: str) -> None:
    """Find this detector's ImagingEvts files under repro_dir and list them."""
    detector = DETECTOR_ALIASES[detector.strip().upper()]
    glob = INSTRUMENTS[detector]["events_glob"]
    root = Path(repro_dir)
    if not root.is_dir():
        die(f"repro dir does not exist: {root}")
    paths = sorted(p.resolve() for p in root.glob(glob) if p.is_file())
    out = Path(out_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text("\n".join(str(p) for p in paths) + ("\n" if paths else ""),
                   encoding="utf-8")


FRAMES_TSV_HEADER = (
    "inst\tsource_base\tbase\tpointing\tsplit\tevent\t"
    "n_rows\tt_start\tt_stop\tduration_s\tra\tdec\tn_samples"
)


def _cluster_pointings(atthk: Path, t_start: float, t_stop: float,
                       threshold_amin: float, min_duration_s: float) -> list[dict]:
    """Cluster ATTHK rows into stable pointings inside [t_start, t_stop].

    A new cluster starts whenever the running centroid drifts more than
    threshold_amin arcminutes from the new sample. Clusters shorter than
    min_duration_s are dropped (slew transitions); after that, near-equal
    centroids are coalesced so a target re-visited later in the observation
    is recognised as the same pointing.
    """
    with fits.open(atthk, memmap=True) as hdul:
        data = hdul["ATTHK"].data
        t   = np.asarray(data["TIME"],   dtype=float)
        ra  = np.asarray(data["AHFRA"],  dtype=float)
        dec = np.asarray(data["AHFDEC"], dtype=float)
    keep = (t >= t_start) & (t <= t_stop)
    t, ra, dec = t[keep], ra[keep], dec[keep]
    if t.size == 0:
        return []

    def sep_amin(r1, d1, r2, d2):
        return float(np.hypot((r1 - r2) * np.cos(np.radians(d2)), d1 - d2) * 60.0)

    raw: list[dict] = []
    cur: dict | None = None
    for ti, ri, di in zip(t, ra, dec):
        if cur and sep_amin(ri, di,
                            float(np.mean(cur["ras"])),
                            float(np.mean(cur["decs"]))) <= threshold_amin:
            cur["tmax"] = float(ti)
            cur["ras"].append(float(ri)); cur["decs"].append(float(di))
        else:
            if cur: raw.append(cur)
            cur = {"tmin": float(ti), "tmax": float(ti),
                   "ras": [float(ri)], "decs": [float(di)]}
    if cur: raw.append(cur)

    stable = [c for c in raw if (c["tmax"] - c["tmin"]) >= min_duration_s]
    if not stable:
        return []

    # Coalesce stable clusters whose centroids fall within threshold.
    used = [False] * len(stable)
    pointings: list[dict] = []
    for i, ci in enumerate(stable):
        if used[i]:
            continue
        used[i] = True
        group = [ci]
        cra, cdec = float(np.mean(ci["ras"])), float(np.mean(ci["decs"]))
        for j in range(i + 1, len(stable)):
            if used[j]:
                continue
            cj = stable[j]
            if sep_amin(float(np.mean(cj["ras"])), float(np.mean(cj["decs"])),
                        cra, cdec) <= threshold_amin:
                used[j] = True
                group.append(cj)
        ras  = [r for c in group for r in c["ras"]]
        decs = [d for c in group for d in c["decs"]]
        segs = sorted([(c["tmin"], c["tmax"]) for c in group])
        pointings.append({
            "ra": float(sum(ras) / len(ras)),
            "dec": float(sum(decs) / len(decs)),
            "t_start": segs[0][0],
            "t_stop":  segs[-1][1],
            "duration": sum(b - a for a, b in segs),
            "n_samples": len(ras),
        })
    pointings.sort(key=lambda p: p["t_start"])
    for k, p in enumerate(pointings):
        p["id"] = f"P{k:02d}"
    return pointings


def _slice_event_file(src: Path, dst: Path, t1: float, t2: float) -> int:
    """Write src restricted to TIME in [t1, t2]; return EVENTS rows kept.

    GTI-style HDUs are intersected with [t1, t2]; TSTART/TSTOP keys clip to
    the window so eexpmap and friends see the right exposure window.
    """
    dst.parent.mkdir(parents=True, exist_ok=True)
    with fits.open(src, memmap=True) as hdul:
        out_hdus: list = []
        kept = 0
        for hdu in hdul:
            name = (hdu.name or "").upper()
            if name == "EVENTS" and hdu.data is not None:
                m = (hdu.data["TIME"] >= t1) & (hdu.data["TIME"] <= t2)
                kept = int(m.sum())
                new = fits.BinTableHDU(data=hdu.data[m],
                                       header=hdu.header.copy(), name="EVENTS")
                new.header["TSTART"] = float(t1); new.header["TSTOP"] = float(t2)
                out_hdus.append(new)
            elif (name.startswith("STDGTI") or name.startswith("GTI")) and hdu.data is not None and len(hdu.data):
                cols = hdu.data.columns.names
                sc = "START" if "START" in cols else cols[0]
                ec = "STOP"  if "STOP"  in cols else cols[1]
                starts = np.maximum(hdu.data[sc].astype(float), t1)
                stops  = np.minimum(hdu.data[ec].astype(float), t2)
                keep = stops > starts
                if keep.any():
                    sub = hdu.data[keep].copy()
                    sub[sc] = starts[keep]; sub[ec] = stops[keep]
                    new = fits.BinTableHDU(data=sub, header=hdu.header.copy(),
                                           name=hdu.name)
                else:
                    new = fits.BinTableHDU(data=np.zeros(0, dtype=hdu.data.dtype),
                                           header=hdu.header.copy(), name=hdu.name)
                new.header["TSTART"] = float(t1); new.header["TSTOP"] = float(t2)
                out_hdus.append(new)
            else:
                out_hdus.append(hdu.copy())
        fits.HDUList(out_hdus).writeto(dst, overwrite=True)
    return kept


def _frame_base(base: str, pid: str) -> str:
    """Insert a pointing id before the trailing _ImagingEvts in a base name."""
    suffix = "_ImagingEvts"
    return f"{base[:-len(suffix)]}_{pid}{suffix}" if base.endswith(suffix) \
        else f"{base}_{pid}"


def _frame_row(det: str, src_base: str, frame_base: str, split: str,
               event_path: Path, n_rows: int, t_start: float, t_stop: float,
               pointing: dict | None) -> str:
    """Build one TSV row for frames.tsv.

    pointing is None for single-pointing pass-throughs that didn't go through
    the attitude clusterer; otherwise it carries id, ra, dec, n_samples and
    duration. Same column ordering as FRAMES_TSV_HEADER.
    """
    duration = pointing["duration"] if pointing else max(t_stop - t_start, 0.0)
    return "\t".join([
        det, src_base, frame_base,
        pointing["id"] if pointing else "",
        split, str(event_path),
        str(n_rows), f"{t_start:.3f}", f"{t_stop:.3f}", f"{duration:.3f}",
        f"{pointing['ra']:.6f}"  if pointing else "",
        f"{pointing['dec']:.6f}" if pointing else "",
        str(pointing["n_samples"]) if pointing else "",
    ])


def frames_split(repro_dir: str, frames_dir: str, atthk_path: str,
                 detectors_arg: str, threshold_amin_arg: str,
                 min_duration_s_arg: str) -> None:
    """For each detector, split multi-pointing repro events by attitude.

    Single-pointing event lists pass through (their absolute path is written
    into <DET>_frames.txt). Multi-pointing event lists are time-sliced into
    one file per stable pointing under frames_dir/events/<DET>/<source_base>/.
    """
    repro = Path(repro_dir); frames = Path(frames_dir); atthk = Path(atthk_path)
    if not atthk.is_file():
        die(f"atthk file not found: {atthk}")
    detectors = normalize_detectors(detectors_arg)
    threshold_amin = float(threshold_amin_arg)
    min_duration_s = float(min_duration_s_arg)

    (frames / "manifest").mkdir(parents=True, exist_ok=True)
    summary_rows = [FRAMES_TSV_HEADER]

    for det in detectors:
        raw_manifest = repro / "manifest" / f"{det}_raw.txt"
        events = _read_manifest(raw_manifest)
        if not events:
            die(f"empty raw manifest: {raw_manifest}")
        out_paths: list[Path] = []

        for src in events:
            if not src.is_file():
                die(f"missing event file: {src}")
            base = src.stem
            with fits.open(src, memmap=True) as hdul:
                evt = hdul["EVENTS"].data
                n_rows = 0 if evt is None else int(len(evt))
                if n_rows:
                    times = evt["TIME"].astype(float)
                    t_lo, t_hi = float(times.min()), float(times.max())
                else:
                    t_lo = t_hi = 0.0
            pointings = (_cluster_pointings(atthk, t_lo, t_hi,
                                            threshold_amin, min_duration_s)
                         if n_rows else [])

            if len(pointings) <= 1:
                p = pointings[0] if pointings else None
                out_paths.append(src.resolve())
                summary_rows.append(_frame_row(
                    det, base, base, "no", src.resolve(), n_rows, t_lo, t_hi, p,
                ))
                continue

            for p in pointings:
                fb = _frame_base(base, p["id"])
                dst = frames / "events" / det / base / f"{fb}{src.suffix or '.fits'}"
                kept = _slice_event_file(src, dst, p["t_start"], p["t_stop"])
                if kept > 0:
                    out_paths.append(dst.resolve())
                summary_rows.append(_frame_row(
                    det, base, fb, "yes" if kept > 0 else "zero",
                    dst.resolve(), kept, p["t_start"], p["t_stop"], p,
                ))

        (frames / "manifest" / f"{det}_frames.txt").write_text(
            "\n".join(str(p) for p in out_paths) + ("\n" if out_paths else ""),
            encoding="utf-8",
        )

    (frames / "frames.tsv").write_text("\n".join(summary_rows) + "\n",
                                       encoding="utf-8")


def _file_kind(name: str) -> str:
    n = name.upper()
    if n.endswith("SUM.SAS"):    return "SUM.SAS"
    if n.endswith("SUM.ASC"):    return "SUM.ASC"
    if n == "CCF.CIF":           return "ccf.cif"
    if n.startswith("MANIFEST"): return "MANIFEST"
    for ext in (".FIT", ".FTZ", ".FITS", ".DS", ".ASC"):
        if n.endswith(ext):      return ext.lstrip(".").lower()
    return "other"


def _write_listing(out_path: Path, paths: Iterable[Path]) -> None:
    out_path.write_text("\n".join(str(p) for p in paths) + "\n", encoding="utf-8")


def qc_init(init_dir: str, out_dir: str, sas_env: str | None = None) -> None:
    init = Path(init_dir)
    out = Path(out_dir); out.mkdir(parents=True, exist_ok=True)

    files = sorted(p for p in init.iterdir() if p.is_file() or p.is_symlink())
    _write_listing(out / "files.txt", files)

    counts: dict[str, int] = {}
    for p in files:
        counts[_file_kind(p.name)] = counts.get(_file_kind(p.name), 0) + 1
    (out / "file_type_counts.txt").write_text(
        "".join(f"{k}\t{counts[k]}\n" for k in sorted(counts)), encoding="utf-8"
    )

    sums  = sorted(init.glob("*SUM.SAS"))
    has_ccf = (init / "ccf.cif").exists()
    status = [
        f"ccf.cif\t{'ok' if has_ccf else 'missing'}",
        f"SUM.SAS\t{'ok' if sums else 'missing'}\t{sums[0].name if sums else ''}",
        f"odf_links\t{sum(1 for p in files if p.is_symlink())}",
    ]
    if sas_env:
        status.append(f"sas_setup.env\t{'ok' if Path(sas_env).is_file() else 'missing'}")
    (out / "status.txt").write_text("\n".join(status) + "\n", encoding="utf-8")


def _stage_qc(stage_dir: Path, manifest_subdir: str, manifest_suffix: str,
              detectors: list[str], out: Path) -> dict[str, list[Path]]:
    """Shared QC scaffolding: file listing, per-detector counts, status TSVs.

    Returns the per-detector resolved-event-list dict so callers can build
    mosaics and stage-specific summaries on top.
    """
    files = sorted(p for p in stage_dir.rglob("*") if p.is_file())
    _write_listing(out / "files.txt", files)

    manifests: dict[str, list[Path]] = {}
    counts: list[str] = []
    status: list[str] = []
    for det in detectors:
        manifest = stage_dir / manifest_subdir / f"{det}_{manifest_suffix}.txt"
        evts = _read_manifest(manifest)
        manifests[det] = evts
        ok = bool(evts) and all(p.is_file() for p in evts)
        counts.append(f"{det}\t{len(evts)}")
        status.append(f"{det}_{manifest_suffix}.txt\t{'ok' if ok else 'missing_or_stale'}")
    (out / "manifest_counts.txt").write_text("\n".join(counts) + "\n", encoding="utf-8")
    (out / "status.txt").write_text("\n".join(status) + "\n", encoding="utf-8")
    return manifests


def qc_repro(repro_dir: str, detectors_arg: str, out_dir: str,
             split_ev_arg: str, attitude_env: str | None = None) -> None:
    repro = Path(repro_dir)
    out = Path(out_dir); out.mkdir(parents=True, exist_ok=True)
    detectors = normalize_detectors(detectors_arg)

    manifests = _stage_qc(repro, "manifest", "raw", detectors, out)

    extras = [f"atthk.dat\t{'ok' if (repro / 'atthk.dat').exists() else 'missing'}"]
    if attitude_env:
        extras.append(
            f"attitude.env\t{'ok' if Path(attitude_env).is_file() else 'missing'}"
        )
    with (out / "status.txt").open("a", encoding="utf-8") as fh:
        fh.write("\n".join(extras) + "\n")

    _event_mosaics(manifests, out, float(split_ev_arg))


def clean_band_table(config_path: str) -> None:
    """Emit one row per band with both per-instrument expressions."""
    cfg = load_config(config_path)
    print("label\tdescription\tpn_expression\tmos_expression")
    for band in _clean_bands(cfg):
        print("\t".join([
            str(band["label"]), str(band.get("description", "")),
            _band_expr(band, "pn"), _band_expr(band, "mos"),
        ]))


def clean_manifest(events_dir: str, out_path: str) -> None:
    """Write the per-(detector, band) manifest of cleaned event files."""
    root = Path(events_dir)
    paths = sorted(p.resolve() for p in root.glob("*_clean.fits") if p.is_file())
    out = Path(out_path); out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text("\n".join(str(p) for p in paths) + ("\n" if paths else ""),
                   encoding="utf-8")


def fits_rows(path_arg: str) -> None:
    """Print EVENTS row count, or -1 if the file is missing."""
    path = Path(path_arg)
    if not path.is_file():
        print(-1); return
    with fits.open(path, memmap=True) as hdul:
        hdu = hdul["EVENTS"] if "EVENTS" in hdul else hdul[1]
        print(0 if hdu.data is None else len(hdu.data))


def qc_frames(frames_dir: str, detectors_arg: str, out_dir: str,
              split_ev_arg: str) -> None:
    frames = Path(frames_dir)
    out = Path(out_dir); out.mkdir(parents=True, exist_ok=True)
    detectors = normalize_detectors(detectors_arg)

    manifests = _stage_qc(frames, "manifest", "frames", detectors, out)

    summary = frames / "frames.tsv"
    if summary.is_file():
        (out / "frames.tsv").write_bytes(summary.read_bytes())
        _frames_timeline(summary, out)

    _event_mosaics(manifests, out, float(split_ev_arg),
                   box_events_by_det=manifests)


def _read_manifest(path: Path) -> list[Path]:
    if not path.is_file():
        return []
    return [Path(line.strip()) for line in path.read_text(encoding="utf-8").splitlines()
            if line.strip()]


def _frames_timeline(summary_tsv: Path, out_dir: Path) -> None:
    """Render frames.tsv as a per-detector Gantt-chart timeline.

    Each row is a (detector, source_base) lane; bars span [t_start, t_stop]
    coloured by pointing id. Bars labelled 'no' (single-pointing pass-through)
    are drawn in grey to distinguish them from genuine splits.
    """
    rows = summary_tsv.read_text(encoding="utf-8").strip().splitlines()
    if len(rows) < 2:
        return
    header = rows[0].split("\t")
    idx = {name: i for i, name in enumerate(header)}
    records = [r.split("\t") for r in rows[1:]]

    lanes: list[tuple[str, str]] = []
    for r in records:
        key = (r[idx["inst"]], r[idx["source_base"]])
        if key not in lanes:
            lanes.append(key)

    pids = sorted({r[idx["pointing"]] for r in records if r[idx["pointing"]]})
    cmap = plt.get_cmap("tab20")
    pid_color = {pid: cmap(i % cmap.N) for i, pid in enumerate(pids)}

    t0 = min(float(r[idx["t_start"]]) for r in records if r[idx["t_start"]])
    fig_h = max(2.5, 0.45 * len(lanes) + 1.5)
    with plt.style.context("dark_background"):
        fig, ax = plt.subplots(figsize=(11, fig_h), dpi=120)
        for r in records:
            try:
                t_a = float(r[idx["t_start"]]) - t0
                t_b = float(r[idx["t_stop"]])  - t0
            except ValueError:
                continue
            if t_b <= t_a:
                continue
            lane = lanes.index((r[idx["inst"]], r[idx["source_base"]]))
            split = r[idx["split"]]
            pid = r[idx["pointing"]]
            color = pid_color.get(pid, (0.6, 0.6, 0.6, 1.0)) if split == "yes" \
                    else (0.45, 0.45, 0.45, 1.0)
            ax.barh(lane, t_b - t_a, left=t_a, height=0.65,
                    color=color, edgecolor="white", linewidth=0.3)
            if pid:
                ax.text(t_a + 0.5 * (t_b - t_a), lane, pid,
                        ha="center", va="center", fontsize=7, color="white")

        ax.set_yticks(range(len(lanes)))
        ax.set_yticklabels([f"{i} {b}" for i, b in lanes], fontsize=8)
        ax.invert_yaxis()
        ax.set_xlabel(f"time since t0 = {t0:.0f} [s]")
        ax.set_title(f"frames timeline  ({len(records)} entries, "
                     f"{sum(1 for r in records if r[idx['split']] == 'yes')} split)")
        ax.grid(axis="x", linestyle=":", alpha=0.3, color="white")
        fig.tight_layout()
        fig.savefig(out_dir / "frames_timeline.png")
        plt.close(fig)


def _mosaic_grid(events: list[Path]
                 ) -> tuple[tuple[float, float, float, float], int, int] | None:
    """Compute the same X/Y extent + grid size used by the QC mosaics.

    Returns (extent, nx, ny) or None when there are no usable events.
    Sharing this with the QC renderer keeps the cut sigma-clip's notion of
    a "5-sigma pixel" identical to the pixel the user inspects in qc-cut.
    """
    xs_all: list[np.ndarray] = []
    ys_all: list[np.ndarray] = []
    for p in events:
        if not p.is_file():
            continue
        x, y, _ = _load_event_xypi(p)
        if x.size:
            xs_all.append(x); ys_all.append(y)
    if not xs_all:
        return None
    x = np.concatenate(xs_all); y = np.concatenate(ys_all)
    xlo, xhi = float(x.min()), float(x.max())
    ylo, yhi = float(y.min()), float(y.max())
    pad_x = max(0.05 * (xhi - xlo), 250.0)
    pad_y = max(0.05 * (yhi - ylo), 250.0)
    extent = (xlo - pad_x, xhi + pad_x, ylo - pad_y, yhi + pad_y)
    nx = 900
    ny = max(200, min(1200,
        int(round(nx * (extent[3] - extent[2]) / (extent[1] - extent[0])))))
    return extent, nx, ny


def cut_run(clean_dir: str, frames_dir: str, out_dir: str,
            config_path: str) -> None:
    """For each (detector, cut band): concat the cleaned per-PPS-band events,
    apply the cut PI window, sigma-clip events sitting in pixels brighter
    than median + N*sigma_MAD of the combined per-detector mosaic, and write
    per-frame cut FITS preserving the cleaned events' header (so REFXCRVL
    etc. stay valid for downstream WCS).
    """
    cfg = load_config(config_path)
    cut_bands  = _cut_bands(cfg)
    sigma_thr  = float(cfg.get("cut_sigma_clip", 5.0))
    detectors  = normalize_detectors(cfg.get("detectors", "PN"))
    clean = Path(clean_dir); frames = Path(frames_dir); out = Path(out_dir)
    (out / "manifest").mkdir(parents=True, exist_ok=True)

    summary = ["det\tband\tpi_min\tpi_max\tn_in\tn_out\tn_clipped\tn_bad_pix\tn_frames"]
    (out / "cut_bands.tsv").write_text(
        "label\tpi_min\tpi_max\n"
        + "\n".join(f"{b['label']}\t{int(b['pi_min'])}\t{int(b['pi_max'])}"
                    for b in cut_bands) + "\n",
        encoding="utf-8",
    )
    for det in detectors:
        frame_paths = _read_manifest(frames / "manifest" / f"{det}_frames.txt")
        if not frame_paths:
            continue
        clean_band_dirs = [d for d in sorted((clean / "events" / det).iterdir())
                           if d.is_dir()]
        for band in cut_bands:
            row = _cut_one_band(det, frame_paths, band,
                                clean_band_dirs, sigma_thr, out)
            if row:
                summary.append(row)
    (out / "cut_summary.tsv").write_text("\n".join(summary) + "\n", encoding="utf-8")


def _cut_one_band(det: str, frame_paths: list[Path], band: dict[str, Any],
                  clean_band_dirs: list[Path], sigma_thr: float,
                  out_dir: Path) -> str | None:
    label = str(band["label"])
    pi_lo, pi_hi = int(band["pi_min"]), int(band["pi_max"])

    per_frame: list[tuple[Path, np.ndarray, Path]] = []
    for frame in frame_paths:
        base = frame.stem
        files = [d / f"{base}_{d.name}_clean.fits" for d in clean_band_dirs]
        files = [f for f in files if f.is_file()]
        chunks: list[np.ndarray] = []
        src: Path | None = None
        for f in files:
            with fits.open(f, memmap=True) as hdul:
                e = hdul["EVENTS"]
                if e.data is None or len(e.data) == 0:
                    if src is None:
                        src = f
                    continue
                chunks.append(np.asarray(e.data))
                if src is None:
                    src = f
        if not chunks or src is None:
            continue
        events = np.concatenate(chunks)
        events = events[(events["PI"] >= pi_lo) & (events["PI"] <= pi_hi)]
        if len(events):
            per_frame.append((frame, events, src))
    if not per_frame:
        return None

    # Combined mosaic, same grid as qc-cut so the user sees exactly which
    # pixels the clip is acting on.
    all_x = np.concatenate([np.asarray(ev["X"], float) for _, ev, _ in per_frame])
    all_y = np.concatenate([np.asarray(ev["Y"], float) for _, ev, _ in per_frame])
    fin = (np.isfinite(all_x) & np.isfinite(all_y)
           & (np.abs(all_x) < 1e7) & (np.abs(all_y) < 1e7))
    xs, ys = all_x[fin], all_y[fin]
    if xs.size == 0:
        return None
    xlo, xhi = float(xs.min()), float(xs.max())
    ylo, yhi = float(ys.min()), float(ys.max())
    pad_x = max(0.05 * (xhi - xlo), 250.0)
    pad_y = max(0.05 * (yhi - ylo), 250.0)
    extent = (xlo - pad_x, xhi + pad_x, ylo - pad_y, yhi + pad_y)
    nx = 900
    ny = max(200, min(1200,
        int(round(nx * (extent[3] - extent[2]) / (extent[1] - extent[0])))))
    img, _xe, _ye = np.histogram2d(
        xs, ys, bins=(nx, ny),
        range=((extent[0], extent[1]), (extent[2], extent[3])),
    )
    img = img.T  # (ny, nx)

    # Poisson σ ≈ sqrt(λ) on the typical pixel: take λ = median of the
    # positive pixels (robust to bright sources/hot pixels that we are
    # actually trying to clip), then mask everything above λ + N*sqrt(λ).
    # For deep PN data this is essentially Gaussian σ; for sparse MOS data
    # MAD would collapse to zero whereas this stays meaningful.
    positive = img[img > 0]
    lam = float(np.median(positive)) if positive.size else 0.0
    threshold = lam + sigma_thr * np.sqrt(max(lam, 1.0))
    bad_mask  = img > threshold
    n_bad_pix = int(bad_mask.sum())

    bin_x = (extent[1] - extent[0]) / nx
    bin_y = (extent[3] - extent[2]) / ny
    events_dir = out_dir / "events" / det / label
    events_dir.mkdir(parents=True, exist_ok=True)
    out_paths: list[Path] = []
    n_in = n_out = n_clipped = 0
    for frame, events, src in per_frame:
        n_in += len(events)
        x = np.asarray(events["X"], float); y = np.asarray(events["Y"], float)
        bx = np.floor((x - extent[0]) / bin_x).astype(int)
        by = np.floor((y - extent[2]) / bin_y).astype(int)
        in_bounds = (bx >= 0) & (bx < nx) & (by >= 0) & (by < ny)
        keep = np.ones(len(events), dtype=bool)
        if n_bad_pix:
            in_bad = np.zeros(len(events), dtype=bool)
            in_bad[in_bounds] = bad_mask[by[in_bounds], bx[in_bounds]]
            keep &= ~in_bad
            n_clipped += int(in_bad.sum())
        kept = events[keep]
        n_out += len(kept)
        out_path = events_dir / f"{frame.stem}_{label}_cut.fits"
        out_paths.append(out_path.resolve())
        # Preserve the source clean file's full HDU structure (PrimaryHDU,
        # OFFSETS, EXPOSU*, STDGTI*, GTI* etc.) so SAS tasks like evselect
        # and eexpmap that index by the dataset summary keep working;
        # only the EVENTS table data gets replaced. We must re-stamp the
        # column-WCS keys (TCTYP/TCRVL/TCRPX/TCDLT/TCUNI) afterwards
        # because astropy regenerates the EVENTS column header from the
        # new data dtype on write and would otherwise drop them.
        _COL_WCS_PREFIXES = ("TCTYP", "TCRVL", "TCRPX", "TCDLT", "TCUNI")
        with fits.open(src, memmap=False) as src_hdul:
            new_hdulist = fits.HDUList([h.copy() for h in src_hdul])
        evt_hdu = next(h for h in new_hdulist if h.name == "EVENTS")
        col_wcs = {k: evt_hdu.header[k] for k in list(evt_hdu.header)
                   if any(k.startswith(p) for p in _COL_WCS_PREFIXES)}
        evt_hdu.data = kept
        for k, v in col_wcs.items():
            evt_hdu.header[k] = v
        new_hdulist.writeto(out_path, overwrite=True)

    manifest = out_dir / "manifest" / f"{det}_{label}_cut.txt"
    manifest.write_text("\n".join(str(p) for p in out_paths) + "\n",
                        encoding="utf-8")

    return (f"{det}\t{label}\t{pi_lo}\t{pi_hi}\t{n_in}\t{n_out}\t"
            f"{n_clipped}\t{n_bad_pix}\t{len(per_frame)}")


def _band_mosaics(events_by_det_band: dict[str, dict[str, list[Path]]],
                  out_dir: Path, *,
                  box_events_by_det: dict[str, list[Path]] | None = None,
                  track_xy_by_det: dict[str, tuple[np.ndarray, np.ndarray]] | None = None
                  ) -> None:
    """Render one mosaic per (detector, pre-binned band) — no PI re-split.

    Each detector gets a shared extent across its bands so soft and hard
    PNGs align side-by-side; per-band stretches are independent.
    """
    summary = []
    for det, by_band in events_by_det_band.items():
        all_events = [p for evs in by_band.values() for p in evs if p.is_file()]
        grid = _mosaic_grid(all_events)
        if grid is None:
            continue
        extent, nx, ny = grid
        boxes = (_boxes_from_events(box_events_by_det[det])
                 if box_events_by_det and det in box_events_by_det else [])
        track_xy = track_xy_by_det.get(det) if track_xy_by_det else None
        summary.append(
            f"{det}\textent={extent}\tpixels={nx}x{ny}\tboxes={len(boxes)}"
        )
        for label, events in by_band.items():
            events = [p for p in events if p.is_file()]
            if not events:
                continue
            xs, ys, _ = zip(*(_load_event_xypi(p) for p in events))
            img, total = _hist2d(list(xs), list(ys),
                                 [np.ones_like(x, dtype=bool) for x in xs],
                                 nx, ny, extent)
            title = (f"{det}  {label}  {total:,} events"
                     + (f"  frames={len(boxes)}" if boxes else "")
                     + ("  +track" if track_xy is not None else ""))
            _save_mosaic(img, extent, boxes,
                         out_dir / f"{det}_{label}_mosaic.png", title,
                         track_xy=track_xy)
            summary.append(f"{det}\t{label}_events={total}")
    (out_dir / "mosaic_summary.txt").write_text(
        "\n".join(summary) + "\n", encoding="utf-8",
    )


def _cut_events_by_det_band(cut_dir: Path, detectors: list[str]
                            ) -> dict[str, dict[str, list[Path]]]:
    """Build {det: {band: [event_paths]}} from on-disk cut manifests."""
    out: dict[str, dict[str, list[Path]]] = {}
    for det in detectors:
        out[det] = {}
        for m in sorted((cut_dir / "manifest").glob(f"{det}_*_cut.txt")):
            label = m.stem.removeprefix(f"{det}_").removesuffix("_cut")
            out[det][label] = _read_manifest(m)
    return out


def qc_cut(cut_dir: str, frames_dir: str, detectors_arg: str,
           out_dir: str) -> None:
    """Cut QC: file listing, per-(det, band) counts, status, mosaics with
    upstream per-frame box overlay. The cut bands are pre-segmented by PI
    so each (det, band) gets its own mosaic — no PI re-split."""
    cut = Path(cut_dir); frames = Path(frames_dir)
    out = Path(out_dir); out.mkdir(parents=True, exist_ok=True)
    detectors = normalize_detectors(detectors_arg)

    files = sorted(p for p in cut.rglob("*") if p.is_file())
    _write_listing(out / "files.txt", files)

    events_by_det_band = _cut_events_by_det_band(cut, detectors)
    counts: list[str] = []
    status: list[str] = []
    for det, by_band in events_by_det_band.items():
        for label, evts in by_band.items():
            ok = bool(evts) and all(p.is_file() for p in evts)
            counts.append(f"{det}\t{label}\t{len(evts)}")
            status.append(f"{det}_{label}_cut.txt\t"
                          f"{'ok' if ok else 'missing_or_stale'}")
    (out / "manifest_counts.txt").write_text(
        "\n".join(counts) + "\n", encoding="utf-8")
    (out / "status.txt").write_text(
        "\n".join(status) + "\n", encoding="utf-8")
    for tsv in ("cut_bands.tsv", "cut_summary.tsv"):
        src = cut / tsv
        if src.is_file():
            (out / tsv).write_bytes(src.read_bytes())

    box_events = {det: _read_manifest(frames / "manifest" / f"{det}_frames.txt")
                  for det in detectors}
    _band_mosaics(events_by_det_band, out, box_events_by_det=box_events)


def _run_sas(tag: str, log_dir: Path, cmd: list[str]) -> None:
    """Run a SAS command, tee output to a per-task log."""
    log_dir.mkdir(parents=True, exist_ok=True)
    log_path = log_dir / f"{tag}.log"
    print(f"[{datetime.now(timezone.utc).strftime('%FT%TZ')}] {tag}: "
          + " ".join(cmd), flush=True)
    with log_path.open("w", encoding="utf-8") as logf:
        result = subprocess.run(cmd, stdout=logf, stderr=subprocess.STDOUT)
    if result.returncode != 0:
        die(f"{tag} failed (rc={result.returncode}); see {log_path}")


def _frame_grid(event_path: Path, bin_arcsec: float, pad_pix: int
                ) -> dict[str, float | int] | None:
    """Per-frame physical-coordinate image grid from the event list X/Y bounds.

    The bin size is bin_arcsec / abs(REFXCDLT in arcsec); xmin/xmax/ymin/ymax
    are aligned to that bin and padded by pad_pix on each side so eexpmap
    has a clean margin.
    """
    with fits.open(event_path, memmap=True) as hdul:
        hdu = hdul["EVENTS"]
        h = hdu.header
        d = hdu.data
    if d is None or len(d) == 0:
        return None
    x = np.asarray(d["X"], dtype=float)
    y = np.asarray(d["Y"], dtype=float)
    good = np.isfinite(x) & np.isfinite(y) & (np.abs(x) < 1e7) & (np.abs(y) < 1e7)
    if not good.any():
        return None
    pix_arcsec = abs(float(h["REFXCDLT"])) * 3600.0
    bin_phys = bin_arcsec / pix_arcsec
    xmin = (math.floor(float(x[good].min()) / bin_phys) - pad_pix) * bin_phys
    xmax = (math.ceil (float(x[good].max()) / bin_phys) + pad_pix) * bin_phys
    ymin = (math.floor(float(y[good].min()) / bin_phys) - pad_pix) * bin_phys
    ymax = (math.ceil (float(y[good].max()) / bin_phys) + pad_pix) * bin_phys
    nx = int(round((xmax - xmin) / bin_phys))
    ny = int(round((ymax - ymin) / bin_phys))
    return {"xmin": xmin, "xmax": xmax, "ymin": ymin, "ymax": ymax,
            "bin": bin_phys, "nx": nx, "ny": ny}


def _frame_base_from_cut(evt_path: Path, label: str) -> str:
    s = evt_path.stem
    suf = f"_{label}_cut"
    return s[:-len(suf)] if s.endswith(suf) else s


def maps_counts(cut_dir: str, maps_dir: str, log_dir: str,
                config_path: str) -> None:
    """Bin each cut event list into a counts image on a per-frame grid."""
    cfg = load_config(config_path)
    cut_bands = _cut_bands(cfg)
    detectors = normalize_detectors(cfg.get("detectors", "PN"))
    bin_arcsec = float(cfg.get("map_bin_arcsec", 4.0))
    pad_pix = int(cfg.get("map_pad_pix", 1))

    cut = Path(cut_dir); maps = Path(maps_dir); log = Path(log_dir)
    grid_rows = ["det\tband\tbase\tximagemin\tximagemax\tyimagemin"
                 "\tyimagemax\tbin\tnx\tny"]
    for det in detectors:
        for band in cut_bands:
            label = str(band["label"])
            events = _read_manifest(cut / "manifest" / f"{det}_{label}_cut.txt")
            band_dir = maps / det / label
            band_dir.mkdir(parents=True, exist_ok=True)
            for evt in events:
                base = _frame_base_from_cut(evt, label)
                grid = _frame_grid(evt, bin_arcsec, pad_pix)
                if grid is None:
                    continue
                grid_rows.append(
                    f"{det}\t{label}\t{base}\t{grid['xmin']:g}\t{grid['xmax']:g}"
                    f"\t{grid['ymin']:g}\t{grid['ymax']:g}\t{grid['bin']:g}"
                    f"\t{grid['nx']}\t{grid['ny']}"
                )
                outp = band_dir / f"{base}_{label}_counts.fits"
                if outp.is_file():
                    continue
                _run_sas(f"maps_counts_{det}_{label}_{base}", log, [
                    "evselect",
                    f"table={evt}:EVENTS",
                    f"imageset={outp}",
                    "withimageset=yes",
                    "xcolumn=X", "ycolumn=Y",
                    "imagebinning=binSize",
                    f"ximagebinsize={grid['bin']:g}",
                    f"yimagebinsize={grid['bin']:g}",
                    "withxranges=yes",
                    f"ximagemin={grid['xmin']:g}",
                    f"ximagemax={grid['xmax']:g}",
                    "withyranges=yes",
                    f"yimagemin={grid['ymin']:g}",
                    f"yimagemax={grid['ymax']:g}",
                    "ignorelegallimits=yes", "writedss=yes",
                ])
    (maps / "maps_grid.tsv").write_text("\n".join(grid_rows) + "\n",
                                        encoding="utf-8")
    _build_maps_manifest(maps, detectors, cut_bands)


def maps_exposure(cut_dir: str, maps_dir: str, atthk_path: str,
                  log_dir: str, config_path: str) -> None:
    """Build per-frame vignetted exposure maps via eexpmap."""
    cfg = load_config(config_path)
    cut_bands = _cut_bands(cfg)
    detectors = normalize_detectors(cfg.get("detectors", "PN"))
    attrebin = str(cfg.get("map_eexpmap_attrebin", "0.020626481"))

    cut = Path(cut_dir); maps = Path(maps_dir); log = Path(log_dir)
    if not Path(atthk_path).is_file():
        die(f"atthk file missing: {atthk_path}")

    for det in detectors:
        for band in cut_bands:
            label = str(band["label"])
            pi_lo, pi_hi = int(band["pi_min"]), int(band["pi_max"])
            band_dir = maps / det / label
            if not band_dir.is_dir():
                continue
            for counts in sorted(band_dir.glob(f"*_{label}_counts.fits")):
                base = counts.stem.removesuffix(f"_{label}_counts")
                evt = cut / "events" / det / label / f"{base}_{label}_cut.fits"
                if not evt.is_file():
                    continue
                outp = band_dir / f"{base}_{label}_exp_vig.fits"
                if outp.is_file():
                    continue
                _run_sas(f"maps_exposure_{det}_{label}_{base}", log, [
                    "eexpmap",
                    f"imageset={counts}",
                    f"attitudeset={atthk_path}",
                    f"eventset={evt}",
                    f"expimageset={outp}",
                    f"pimin={pi_lo}", f"pimax={pi_hi}",
                    "withvignetting=yes",
                    f"attrebin={attrebin}",
                ])
    _build_maps_manifest(maps, detectors, cut_bands)


def _read_first2d(path: Path) -> tuple[np.ndarray, fits.Header]:
    with fits.open(path, memmap=True) as hdul:
        for hdu in hdul:
            d = getattr(hdu, "data", None)
            if d is not None and np.asarray(d).ndim >= 2:
                return np.asarray(d, dtype=float), hdu.header.copy()
    die(f"no 2D image in {path}")


def _solve_bkg(c: np.ndarray, e: np.ndarray, mode: str) -> tuple[float, float, float, str]:
    """Return (a, b, rmse, fit_mode) given off-source counts/exposure samples."""
    if c.size == 0:
        return 0.0, 0.0, 0.0, "empty"
    if mode == "flat":
        a = float(c.mean())
        return a, 0.0, float(np.sqrt(np.mean((c - a) ** 2))), "flat"
    have_span = c.size > 1 and (e.max() - e.min()) > 0
    if have_span:
        A = np.column_stack([np.ones(c.size), e])
        coeffs, *_ = np.linalg.lstsq(A, c, rcond=None)
        a_u, b_u = float(coeffs[0]), float(coeffs[1])
    else:
        a_u, b_u = float(c.mean()), 0.0
    if mode == "none":
        rmse = float(np.sqrt(np.mean((c - (a_u + b_u * e)) ** 2)))
        return a_u, b_u, rmse, "unconstrained"
    # nonneg: pick best of {unconstrained-if-pos, intercept_zero, slope_zero, zero}
    cands: list[tuple[str, float, float, float]] = []
    if a_u >= 0 and b_u >= 0:
        rmse = float(np.sqrt(np.mean((c - (a_u + b_u * e)) ** 2)))
        cands.append(("unconstrained", a_u, b_u, rmse))
    denom = float(np.dot(e, e))
    b_only = max(float(np.dot(e, c)) / denom, 0.0) if denom > 0 else 0.0
    cands.append(("intercept_zero", 0.0, b_only,
                  float(np.sqrt(np.mean((c - b_only * e) ** 2)))))
    a_only = max(float(c.mean()), 0.0)
    cands.append(("slope_zero", a_only, 0.0,
                  float(np.sqrt(np.mean((c - a_only) ** 2)))))
    cands.append(("zero", 0.0, 0.0,
                  float(np.sqrt(np.mean(c ** 2)))))
    name, a, b, rmse = min(cands, key=lambda r: r[3])
    return a, b, rmse, name


def maps_background(maps_dir: str, config_path: str) -> None:
    """Per-frame a + b·E background fit on off-source sample pixels."""
    cfg = load_config(config_path)
    cut_bands = _cut_bands(cfg)
    detectors = normalize_detectors(cfg.get("detectors", "PN"))
    mode = str(cfg.get("map_bkg_constrain", "flat")).strip().lower()
    if mode not in {"none", "nonneg", "flat"}:
        die(f"map_bkg_constrain must be none/nonneg/flat (got {mode!r})")
    outside_frac = float(cfg.get("map_bkg_outside_radius_fraction", 0.5))

    maps = Path(maps_dir)
    summary = ["det\tband\tbase\tfit_a\tfit_b\trmse\tn_sample\t"
               "center_x\tcenter_y\tradius_pix\tfit_mode"]

    for det in detectors:
        for band in cut_bands:
            label = str(band["label"])
            band_dir = maps / det / label
            if not band_dir.is_dir():
                continue
            for counts_path in sorted(band_dir.glob(f"*_{label}_counts.fits")):
                base = counts_path.stem.removesuffix(f"_{label}_counts")
                exp_path = band_dir / f"{base}_{label}_exp_vig.fits"
                if not exp_path.is_file():
                    continue
                bkg_path = band_dir / f"{base}_{label}_bkg.fits"
                row = _fit_one_bkg(counts_path, exp_path, bkg_path,
                                   mode, outside_frac)
                summary.append(f"{det}\t{label}\t{base}\t{row}")
    (maps / "bkg_summary.tsv").write_text("\n".join(summary) + "\n",
                                          encoding="utf-8")
    _build_maps_manifest(maps, detectors, cut_bands)


def _fit_one_bkg(counts_path: Path, exp_path: Path, bkg_path: Path,
                 mode: str, outside_frac: float) -> str:
    counts_img, c_header = _read_first2d(counts_path)
    exp_img, _ = _read_first2d(exp_path)
    if counts_img.shape != exp_img.shape:
        die(f"counts/exposure shape mismatch: {counts_path} vs {exp_path}")
    footprint = np.isfinite(exp_img) & (exp_img > 0)
    ny, nx = counts_img.shape
    if footprint.any():
        ys, xs = np.where(footprint)
        cy = float(ys.mean()); cx = float(xs.mean())
        half = 0.5 * max(ys.max() - ys.min(), xs.max() - xs.min())
        radius = outside_frac * half
        yy, xx = np.ogrid[:ny, :nx]
        outside_circle = (yy - cy) ** 2 + (xx - cx) ** 2 >= radius * radius
        sample = footprint & outside_circle & np.isfinite(counts_img)
        if not sample.any():
            sample = footprint & np.isfinite(counts_img)
    else:
        cy = (ny - 1) / 2.0; cx = (nx - 1) / 2.0
        radius = 0.0
        sample = np.zeros_like(footprint, dtype=bool)
    n_sample = int(sample.sum())
    c = counts_img[sample].astype(float) if n_sample else np.zeros(0)
    e = exp_img[sample].astype(float) if n_sample else np.zeros(0)
    a, b, rmse, fit_mode = _solve_bkg(c, e, mode)
    bg = (a + b * exp_img).astype(np.float32)
    bg = np.where(footprint, bg, 0.0).astype(np.float32)

    new_header = fits.Header()
    for key in c_header:
        if key in {"SIMPLE", "BITPIX", "NAXIS", "NAXIS1", "NAXIS2", "EXTEND"}:
            continue
        try:
            new_header[key] = c_header[key]
        except (ValueError, KeyError):
            continue
    out_hdu = fits.PrimaryHDU(data=bg, header=new_header)
    out_hdu.header["BGMODEL"] = ("a+bE", "background model")
    out_hdu.header["BGMODE"]  = (fit_mode, "fit constraint mode")
    out_hdu.header["BGA"]     = (a, "fit intercept")
    out_hdu.header["BGB"]     = (b, "fit slope vs exp_vig")
    out_hdu.header["BGRMSE"]  = (rmse, "sample RMS residual")
    out_hdu.header["BG_CX"]   = (cx, "exclusion center x [pix]")
    out_hdu.header["BG_CY"]   = (cy, "exclusion center y [pix]")
    out_hdu.header["BG_RAD"]  = (radius, "exclusion radius [pix]")
    out_hdu.header["BG_NSAM"] = (n_sample, "n sample pixels")
    out_hdu.writeto(bkg_path, overwrite=True)
    return (f"{a:.6g}\t{b:.6g}\t{rmse:.6g}\t{n_sample}"
            f"\t{cx:.2f}\t{cy:.2f}\t{radius:.2f}\t{fit_mode}")


def maps_corrected(maps_dir: str, config_path: str) -> None:
    """(counts - bkg) / exp_vig per frame; not clipped."""
    cfg = load_config(config_path)
    cut_bands = _cut_bands(cfg)
    detectors = normalize_detectors(cfg.get("detectors", "PN"))
    maps = Path(maps_dir)
    for det in detectors:
        for band in cut_bands:
            label = str(band["label"])
            band_dir = maps / det / label
            if not band_dir.is_dir():
                continue
            for counts in sorted(band_dir.glob(f"*_{label}_counts.fits")):
                base = counts.stem.removesuffix(f"_{label}_counts")
                exp = band_dir / f"{base}_{label}_exp_vig.fits"
                bkg = band_dir / f"{base}_{label}_bkg.fits"
                if not (exp.is_file() and bkg.is_file()):
                    continue
                out = band_dir / f"{base}_{label}_corrected.fits"
                _make_corrected(counts, exp, bkg, out)
    _build_maps_manifest(maps, detectors, cut_bands)


def _make_corrected(counts_path: Path, exp_path: Path, bkg_path: Path,
                    out_path: Path) -> None:
    c, c_header = _read_first2d(counts_path)
    e, _ = _read_first2d(exp_path)
    b, _ = _read_first2d(bkg_path)
    with np.errstate(divide="ignore", invalid="ignore"):
        rate = (c - b) / e
    rate = np.where(np.isfinite(rate) & (e > 0), rate, 0.0).astype(np.float32)
    new_header = fits.Header()
    for key in c_header:
        if key in {"SIMPLE", "BITPIX", "NAXIS", "NAXIS1", "NAXIS2", "EXTEND"}:
            continue
        try:
            new_header[key] = c_header[key]
        except (ValueError, KeyError):
            continue
    fits.PrimaryHDU(data=rate, header=new_header).writeto(out_path, overwrite=True)


def _build_maps_manifest(maps_dir: Path, detectors: list[str],
                         cut_bands: list[dict[str, Any]]) -> None:
    rows = ["det\tband\tbase\tcounts\texp_vig\tbkg\tcorrected"]
    for det in detectors:
        for band in cut_bands:
            label = str(band["label"])
            band_dir = maps_dir / det / label
            if not band_dir.is_dir():
                continue
            for counts in sorted(band_dir.glob(f"*_{label}_counts.fits")):
                base = counts.stem.removesuffix(f"_{label}_counts")
                paths = {
                    "counts":    counts,
                    "exp_vig":   band_dir / f"{base}_{label}_exp_vig.fits",
                    "bkg":       band_dir / f"{base}_{label}_bkg.fits",
                    "corrected": band_dir / f"{base}_{label}_corrected.fits",
                }
                row = [det, label, base]
                for k in ("counts", "exp_vig", "bkg", "corrected"):
                    p = paths[k]
                    row.append(str(p.resolve()) if p.is_file() else "")
                rows.append("\t".join(row))
    (maps_dir / "maps_manifest.tsv").write_text(
        "\n".join(rows) + "\n", encoding="utf-8")


def _read_maps_grid(maps_dir: Path) -> dict[tuple[str, str, str], dict[str, float]]:
    p = maps_dir / "maps_grid.tsv"
    if not p.is_file():
        return {}
    lines = p.read_text(encoding="utf-8").splitlines()
    if not lines:
        return {}
    head = lines[0].split("\t")
    idx = {n: i for i, n in enumerate(head)}
    out: dict[tuple[str, str, str], dict[str, float]] = {}
    for line in lines[1:]:
        parts = line.split("\t")
        if len(parts) < len(head):
            continue
        out[(parts[idx["det"]], parts[idx["band"]], parts[idx["base"]])] = {
            "xmin": float(parts[idx["ximagemin"]]),
            "ymin": float(parts[idx["yimagemin"]]),
            "bin":  float(parts[idx["bin"]]),
            "nx":   int(parts[idx["nx"]]),
            "ny":   int(parts[idx["ny"]]),
        }
    return out


def qc_maps(maps_dir: str, frames_dir: str, detectors_arg: str,
            out_dir: str) -> None:
    """Per-(detector, band) mosaics of the nobkg-corrected (counts/exp_vig)
    and the corrected ((counts - bkg)/exp_vig) rate maps. Rates are co-added
    across frames as sum(num)/sum(exp) on a viewing grid built from the
    union of the per-frame physical extents.
    """
    maps = Path(maps_dir); frames = Path(frames_dir)
    out = Path(out_dir); out.mkdir(parents=True, exist_ok=True)
    detectors = normalize_detectors(detectors_arg)

    files = sorted(p for p in maps.rglob("*") if p.is_file())
    _write_listing(out / "files.txt", files)

    for tsv in ("maps_grid.tsv", "maps_manifest.tsv", "bkg_summary.tsv"):
        src = maps / tsv
        if src.is_file():
            (out / tsv).write_bytes(src.read_bytes())

    grids = _read_maps_grid(maps)
    rows = _read_maps_manifest_rows(maps)

    summary = ["det\tband\tnframes\tnx\tny\tnobkg_corrected\tcorrected"]
    box_events = {det: _read_manifest(frames / "manifest" / f"{det}_frames.txt")
                  for det in detectors}

    bands: list[str] = []
    for r in rows:
        if r["band"] not in bands:
            bands.append(r["band"])

    for det in detectors:
        for band in bands:
            band_rows = [r for r in rows if r["det"] == det and r["band"] == band]
            band_rows = [r for r in band_rows
                         if r["counts"] and r["exp_vig"] and r["bkg"]]
            if not band_rows:
                continue
            keys = [(det, band, r["base"]) for r in band_rows]
            keys = [k for k in keys if k in grids]
            if not keys:
                continue
            bin_phys = grids[keys[0]]["bin"]
            xmin = min(grids[k]["xmin"] for k in keys)
            ymin = min(grids[k]["ymin"] for k in keys)
            xmax = max(grids[k]["xmin"] + grids[k]["nx"] * grids[k]["bin"]
                       for k in keys)
            ymax = max(grids[k]["ymin"] + grids[k]["ny"] * grids[k]["bin"]
                       for k in keys)
            nx = int(round((xmax - xmin) / bin_phys))
            ny = int(round((ymax - ymin) / bin_phys))
            sum_c = np.zeros((ny, nx), dtype=float)
            sum_e = np.zeros((ny, nx), dtype=float)
            sum_b = np.zeros((ny, nx), dtype=float)
            for r in band_rows:
                key = (det, band, r["base"])
                if key not in grids:
                    continue
                g = grids[key]
                ix = int(round((g["xmin"] - xmin) / bin_phys))
                iy = int(round((g["ymin"] - ymin) / bin_phys))
                cimg, _ = _read_first2d(Path(r["counts"]))
                eimg, _ = _read_first2d(Path(r["exp_vig"]))
                bimg, _ = _read_first2d(Path(r["bkg"]))
                fy, fx = cimg.shape
                sum_c[iy:iy+fy, ix:ix+fx] += cimg
                sum_e[iy:iy+fy, ix:ix+fx] += eimg
                sum_b[iy:iy+fy, ix:ix+fx] += bimg
            with np.errstate(divide="ignore", invalid="ignore"):
                nobkg = np.where(sum_e > 0, sum_c / sum_e, 0.0)
                corr  = np.where(sum_e > 0, (sum_c - sum_b) / sum_e, 0.0)
            extent = (xmin, xmin + nx * bin_phys,
                      ymin, ymin + ny * bin_phys)
            boxes = _boxes_from_events(box_events.get(det, []))

            for tag, img in (("nobkg-corrected", nobkg), ("corrected", corr)):
                positive = img[img > 0]
                if positive.size:
                    vmin = float(np.percentile(positive, 5))
                    vmax = float(np.percentile(positive, 99.5))
                else:
                    vmin, vmax = 1.0e-6, 1.0e-3
                if not (vmax > vmin > 0):
                    vmax = max(vmax, vmin * 10.0, 1e-6)
                title = (f"{det}  {band}  {tag}  "
                         f"frames={len(band_rows)}  bin={bin_phys:g} phys")
                _save_mosaic(img, extent, boxes,
                             out / f"{det}_{band}_{tag}_mosaic.png", title,
                             vmin=vmin, vmax=vmax)
            summary.append(
                f"{det}\t{band}\t{len(band_rows)}\t{nx}\t{ny}"
                f"\t{det}_{band}_nobkg-corrected_mosaic.png"
                f"\t{det}_{band}_corrected_mosaic.png"
            )
    (out / "qc_maps_summary.tsv").write_text("\n".join(summary) + "\n",
                                             encoding="utf-8")


def _read_maps_manifest_rows(maps_dir: Path) -> list[dict[str, str]]:
    p = maps_dir / "maps_manifest.tsv"
    if not p.is_file():
        return []
    lines = p.read_text(encoding="utf-8").splitlines()
    if not lines:
        return []
    head = lines[0].split("\t")
    out: list[dict[str, str]] = []
    for line in lines[1:]:
        parts = line.split("\t")
        if len(parts) < len(head):
            continue
        out.append({n: parts[i] for i, n in enumerate(head)})
    return out


def _observation_window(event_paths: list[Path]) -> tuple[float, float, float, str]:
    """Read TSTART/TSTOP/MJDREF/TIMESYS from event headers and union them."""
    tmin, tmax = float("inf"), float("-inf")
    mjdref: float | None = None
    timesys: str | None = None
    for p in event_paths:
        with fits.open(p, memmap=True) as hdul:
            h = hdul["EVENTS"].header
        tmin = min(tmin, float(h["TSTART"]))
        tmax = max(tmax, float(h["TSTOP"]))
        mref = float(h.get("MJDREF",
                           float(h.get("MJDREFI", 0)) + float(h.get("MJDREFF", 0))))
        tsys = str(h.get("TIMESYS", "TT")).strip().upper()
        if mjdref is None:
            mjdref, timesys = mref, tsys
        elif abs(mref - mjdref) > 1e-9 or tsys != timesys:
            die(f"event time-system mismatch in {p}")
    if not (np.isfinite(tmin) and np.isfinite(tmax)) or mjdref is None:
        die("could not derive observation window from event headers")
    return tmin, tmax, mjdref, timesys  # type: ignore[return-value]


def _xmm_seconds_to_isot(seconds: float, mjdref: float, timesys: str) -> str:
    return (Time(mjdref, format="mjd", scale=timesys.lower())
            + TimeDelta(seconds, format="sec")).utc.isot


def _step_seconds(step: str) -> float:
    """Parse a Horizons-style step ('30s', '1m', '1h', '1d') into seconds."""
    s = step.strip().lower()
    for unit, mul in (("d", 86400.0), ("h", 3600.0), ("m", 60.0), ("s", 1.0)):
        if s.endswith(unit):
            return float(s[: -len(unit)].strip()) * mul
    return float(s)


def _query_horizons(target_id: str, id_type: str, observer: str,
                    start_utc: str, stop_utc: str, step: str
                    ) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Query JPL Horizons for (MJD, RA, DEC, Δ) over the window.

    observer is a Horizons location code; for XMM-Newton observations use
    "500@-125989" so the apparent positions include XMM's orbital parallax
    (geocentric "500@399" is offset by up to ~0.3° from a moving target at
    XMM's apogee).
    """
    from astroquery.jplhorizons import Horizons  # heavy import, only here

    obj = Horizons(id=target_id, id_type=id_type, location=observer,
                   epochs={"start": start_utc, "stop": stop_utc, "step": step.strip()})
    eph = obj.ephemerides(quantities="2,20", extra_precision=True)

    def _col(names: list[str]) -> str:
        for want in names:
            for col in eph.colnames:
                if col.lower() == want.lower():
                    return col
        die(f"Horizons output missing one of {names}; got {eph.colnames}")

    jd  = np.asarray(eph[_col(["datetime_jd", "JD"])], dtype=float)
    ra  = np.asarray(eph[_col(["RA_app", "RA_ICRF_app", "RA"])], dtype=float)
    dec = np.asarray(eph[_col(["DEC_app", "DEC_ICRF_app", "DEC"])], dtype=float)
    delta = np.asarray(eph[_col(["delta"])], dtype=float)
    return jd - 2400000.5, ra, dec, delta


def _load_track_input(path_arg: str
                      ) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Load (MJD, RA, DEC, DELTA) from a FITS table or a CSV with named columns."""
    p = Path(path_arg)
    if not p.is_file():
        die(f"track_input not found: {p}")
    if p.suffix.lower() in {".csv", ".txt"}:
        data = np.genfromtxt(p, delimiter=",", names=True, dtype=None, encoding=None)
        names = data.dtype.names or ()
        def _get(*want):
            for w in want:
                for n in names:
                    if n.lower() == w.lower():
                        return n
            die(f"{p}: missing column {want}")
        if "mjd" in (n.lower() for n in names):
            mjd = np.asarray(data[_get("MJD")], dtype=float)
        else:
            mjd = np.asarray(Time(np.asarray(data[_get("time_iso", "utc_iso", "isot")]),
                                  format="isot", scale="utc").mjd, dtype=float)
        return (mjd,
                np.asarray(data[_get("ra_deg", "ra")], dtype=float),
                np.asarray(data[_get("dec_deg", "dec")], dtype=float),
                np.asarray(data[_get("delta_au", "delta")], dtype=float))
    with fits.open(p, memmap=False) as hdul:
        for hdu in hdul[1:]:
            if getattr(hdu, "data", None) is not None and len(hdu.data):
                t = hdu.data
                cols = {c.upper(): c for c in t.columns.names}
                for need in ("MJD", "RA", "DEC", "DELTA"):
                    if need not in cols:
                        die(f"{p}: track FITS missing {need}")
                return (np.asarray(t[cols["MJD"]], dtype=float),
                        np.asarray(t[cols["RA"]],  dtype=float),
                        np.asarray(t[cols["DEC"]], dtype=float),
                        np.asarray(t[cols["DELTA"]], dtype=float))
    die(f"{p}: no usable track table")


def build_track(frames_manifest_glob: str, target_id: str, id_type: str,
                observer: str, step: str, track_input: str, out_fits: str,
                out_env: str, out_json: str) -> None:
    """Resolve the observation window, query Horizons (or load a local file),
    sort and dedupe by MJD, then write track.fits + track.env + track.json.
    """
    paths = sorted({Path(line.strip()).resolve()
                    for glob_pat in frames_manifest_glob.split(",")
                    for line in Path(glob_pat).read_text(encoding="utf-8").splitlines()
                    if line.strip()})
    if not paths:
        die("no frame events to define observation window")
    tmin, tmax, mjdref, timesys = _observation_window(paths)
    start_utc = _xmm_seconds_to_isot(tmin, mjdref, timesys)
    stop_utc  = _xmm_seconds_to_isot(tmax, mjdref, timesys)

    # Horizons aligns sample times to step boundaries, so its first/last
    # samples can fall inside [start_utc, stop_utc] by up to one step.
    # Pad the query window by one step on each side so the resulting track
    # still brackets every event time.
    step_sec = _step_seconds(step)
    pad = TimeDelta(step_sec, format="sec")
    if track_input.strip():
        mjd, ra, dec, delta = _load_track_input(track_input.strip())
        source = f"file:{Path(track_input).name}"
    else:
        if not target_id.strip():
            die("config target_id is empty and no track_input provided")
        mjd, ra, dec, delta = _query_horizons(
            target_id.strip(), (id_type or "smallbody").strip(),
            observer.strip() or "500@-125989",
            (Time(start_utc, format="isot", scale="utc") - pad).utc.isot,
            (Time(stop_utc,  format="isot", scale="utc") + pad).utc.isot,
            step,
        )
        source = f"Horizons:{target_id.strip()}@{observer}"

    order = np.argsort(mjd)
    mjd, ra, dec, delta = (a[order] for a in (mjd, ra, dec, delta))
    keep = np.concatenate([[True], np.diff(mjd) > 0]) if len(mjd) > 1 \
           else np.ones(len(mjd), dtype=bool)
    mjd, ra, dec, delta = (a[keep] for a in (mjd, ra, dec, delta))
    obs_mjd_lo = float(Time(start_utc, format="isot", scale="utc").mjd)
    obs_mjd_hi = float(Time(stop_utc,  format="isot", scale="utc").mjd)
    if len(mjd) == 0 or mjd[0] > obs_mjd_lo or mjd[-1] < obs_mjd_hi:
        die("track does not cover the full observation window")

    mid = len(mjd) // 2
    ref_ra, ref_dec = float(ra[mid]), float(dec[mid])

    out = Path(out_fits); out.parent.mkdir(parents=True, exist_ok=True)
    hdu = fits.BinTableHDU.from_columns([
        fits.Column(name="MJD",   format="D", array=mjd),
        fits.Column(name="RA",    format="D", array=ra),
        fits.Column(name="DEC",   format="D", array=dec),
        fits.Column(name="DELTA", format="D", array=delta),
    ], name="OBJTRACK")
    hdu.header["ORIGINTRK"] = (source[:68], "Track source")
    hdu.header["MJDREF"]  = (mjdref, "Event-file MJD reference")
    hdu.header["TIMESYS"] = (timesys, "Event-file time scale")
    hdu.header["TSTART"]  = (float(tmin), "Observation start [s]")
    hdu.header["TSTOP"]   = (float(tmax), "Observation stop [s]")
    fits.HDUList([fits.PrimaryHDU(), hdu]).writeto(out, overwrite=True)

    Path(out_env).write_text(
        f'export TRACK_FILE={shlex.quote(str(out.resolve()))}\n'
        f'export COMET_REF_RA="{ref_ra:.10f}"\n'
        f'export COMET_REF_DEC="{ref_dec:.10f}"\n',
        encoding="utf-8",
    )
    Path(out_json).write_text(json.dumps({
        "source": source,
        "target_id": target_id.strip(),
        "observer": observer.strip(),
        "track_input": track_input.strip(),
        "step": step,
        "n_rows": int(len(mjd)),
        "observation_tstart_s": float(tmin),
        "observation_tstop_s":  float(tmax),
        "observation_start_utc": start_utc,
        "observation_stop_utc":  stop_utc,
        "track_mjd_min": float(mjd[0]),
        "track_mjd_max": float(mjd[-1]),
        "ref_ra_deg":  ref_ra,
        "ref_dec_deg": ref_dec,
    }, indent=2) + "\n", encoding="utf-8")


def _read_track(track_fits: str
                ) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    with fits.open(track_fits, memmap=False) as hdul:
        t = hdul["OBJTRACK"].data
        return (np.asarray(t["MJD"], dtype=float),
                np.asarray(t["RA"],  dtype=float),
                np.asarray(t["DEC"], dtype=float),
                np.asarray(t["DELTA"], dtype=float))


def _event_xy_wcs(event_path: Path) -> WCS:
    """Build the celestial WCS from XMM event REFX*/REFY* keywords."""
    with fits.open(event_path, memmap=True) as hdul:
        h = hdul["EVENTS"].header
    w = WCS(naxis=2)
    w.wcs.crval = [float(h["REFXCRVL"]), float(h["REFYCRVL"])]
    w.wcs.crpix = [float(h["REFXCRPX"]), float(h["REFYCRPX"])]
    w.wcs.cdelt = [float(h["REFXCDLT"]), float(h["REFYCDLT"])]
    w.wcs.ctype = [str(h["REFXCTYP"]).strip(), str(h["REFYCTYP"]).strip()]
    return w


def qc_track(track_fits: str, cut_dir: str, frames_dir: str,
             detectors_arg: str, out_dir: str) -> None:
    """Overlay the comet trajectory on the cut per-(detector, band) mosaics.

    All cut events of one observation share the same celestial reference
    (inherited from the cleaned events' REFXCRVL/REFYCRVL/...), so we
    project the track once per detector and overlay the resulting X(t), Y(t)
    curve on every cut band PNG.
    """
    cut = Path(cut_dir); frames = Path(frames_dir)
    out = Path(out_dir); out.mkdir(parents=True, exist_ok=True)
    detectors = normalize_detectors(detectors_arg)

    mjd, ra, dec, _delta = _read_track(track_fits)
    (out / "track_window.txt").write_text(
        f"track_mjd_min\t{mjd[0]:.6f}\ntrack_mjd_max\t{mjd[-1]:.6f}\n"
        f"track_rows\t{len(mjd)}\n", encoding="utf-8",
    )

    events_by_det_band = _cut_events_by_det_band(cut, detectors)
    box_events = {det: _read_manifest(frames / "manifest" / f"{det}_frames.txt")
                  for det in detectors}
    track_xy_by_det: dict[str, tuple[np.ndarray, np.ndarray]] = {}
    for det, by_band in events_by_det_band.items():
        event_path = next(
            (p for evs in by_band.values() for p in evs if p.is_file()), None,
        )
        if event_path is None:
            continue
        track_xy_by_det[det] = _event_xy_wcs(event_path).wcs_world2pix(ra, dec, 1)

    _band_mosaics(events_by_det_band, out,
                  box_events_by_det=box_events,
                  track_xy_by_det=track_xy_by_det)


def _clean_events_by_det(clean_dir: Path,
                         detectors: list[str]) -> dict[str, list[Path]]:
    """Combine all per-band cleaned event manifests into one list per detector."""
    out: dict[str, list[Path]] = {}
    for det in detectors:
        events: list[Path] = []
        for m in sorted((clean_dir / "manifest").glob(f"{det}_*_clean.txt")):
            events.extend(_read_manifest(m))
        out[det] = events
    return out


def qc_clean(clean_dir: str, frames_dir: str, detectors_arg: str,
             split_ev_arg: str, out_dir: str) -> None:
    """Clean QC: file listing, per-(det, band) counts, status, loE/hiE
    mosaics with per-frame boxes, and per-detector flare lightcurve grid.
    """
    clean = Path(clean_dir)
    frames = Path(frames_dir)
    out = Path(out_dir); out.mkdir(parents=True, exist_ok=True)
    detectors = normalize_detectors(detectors_arg)

    files = sorted(p for p in clean.rglob("*") if p.is_file())
    _write_listing(out / "files.txt", files)

    # Discover band labels from on-disk manifests so QC stays independent
    # of config drift between runs.
    labels: list[str] = []
    for det in detectors:
        for m in sorted((clean / "manifest").glob(f"{det}_*_clean.txt")):
            label = m.stem.removeprefix(f"{det}_").removesuffix("_clean")
            if label and label not in labels:
                labels.append(label)

    counts: list[str] = []
    status: list[str] = []
    for det in detectors:
        for label in labels:
            manifest = clean / "manifest" / f"{det}_{label}_clean.txt"
            evts = _read_manifest(manifest)
            ok = bool(evts) and all(p.is_file() for p in evts)
            counts.append(f"{det}\t{label}\t{len(evts)}")
            status.append(f"{det}_{label}_clean.txt\t"
                          f"{'ok' if ok else 'missing_or_stale'}")
    (out / "manifest_counts.txt").write_text(
        "\n".join(counts) + "\n", encoding="utf-8")
    (out / "status.txt").write_text(
        "\n".join(status) + "\n", encoding="utf-8")

    for tsv in ("clean_band_filters.tsv", "flare_gti_summary.tsv"):
        src = clean / tsv
        if src.is_file():
            (out / tsv).write_bytes(src.read_bytes())

    events_by_det = _clean_events_by_det(clean, detectors)
    box_events    = {det: _read_manifest(frames / "manifest" / f"{det}_frames.txt")
                     for det in detectors}
    _event_mosaics(events_by_det, out, float(split_ev_arg),
                   box_events_by_det=box_events)
    summary = clean / "flare_gti_summary.tsv"
    if summary.is_file():
        _flare_lc_grid(summary, detectors, out)


def _flare_lc_grid(summary_tsv: Path, detectors: list[str],
                   out_dir: Path) -> None:
    """Per-detector grid of flare lightcurves with rate-cut threshold and GTI overlay.

    Each panel is one frame: bin-centred rates as a step plot, the
    per-instrument rate cut as a dashed red line, and accepted-time
    intervals (from the STDGTI table) shaded green. Times are referenced
    to each frame's own t0 so panels share a comparable axis.
    """
    rows = summary_tsv.read_text(encoding="utf-8").strip().splitlines()
    if len(rows) < 2:
        return
    header = rows[0].split("\t")
    idx = {name: i for i, name in enumerate(header)}
    by_det: dict[str, list[list[str]]] = {d: [] for d in detectors}
    for r in rows[1:]:
        parts = r.split("\t")
        if parts[idx["inst"]] in by_det:
            by_det[parts[idx["inst"]]].append(parts)

    for det, recs in by_det.items():
        if not recs:
            continue
        ncols = 4
        nrows = (len(recs) + ncols - 1) // ncols
        with plt.style.context("dark_background"):
            fig, axes = plt.subplots(nrows, ncols,
                                     figsize=(ncols * 3.4, nrows * 2.0),
                                     dpi=120, squeeze=False)
            for ax in axes.flat:
                ax.set_visible(False)
            for ax, rec in zip(axes.flat, recs):
                rate_cut = float(rec[idx["rate_cut"]])
                lc_path = Path(rec[idx["lightcurve"]])
                gti_path = Path(rec[idx["gti"]])
                if not lc_path.is_file():
                    continue
                with fits.open(lc_path, memmap=True) as hdul:
                    rate = hdul["RATE"].data
                    t_abs = np.asarray(rate["TIME"], dtype=float)
                    r = np.asarray(rate["RATE"], dtype=float)
                if t_abs.size == 0:
                    continue
                ax.set_visible(True)
                t0 = float(t_abs[0])
                t = t_abs - t0
                ax.step(t, r, where="mid", color="white", linewidth=0.6)
                ax.axhline(rate_cut, color="#ff5555", linestyle="--",
                           linewidth=0.7)
                if gti_path.is_file():
                    with fits.open(gti_path, memmap=True) as hdul:
                        gti_hdu = hdul["STDGTI"] if "STDGTI" in hdul else hdul[1]
                        cols = gti_hdu.data.columns.names
                        starts = np.asarray(gti_hdu.data[cols[0]], dtype=float) - t0
                        stops  = np.asarray(gti_hdu.data[cols[1]], dtype=float) - t0
                    for s, e in zip(starts, stops):
                        ax.axvspan(s, e, color="#33aa33", alpha=0.2, lw=0)
                ax.set_title(rec[idx["frame"]], fontsize=6)
                ax.tick_params(labelsize=6)
            fig.suptitle(f"{det} flare lightcurves  "
                         f"(red --- = rate cut, green = accepted GTI)",
                         fontsize=10)
            fig.tight_layout()
            fig.savefig(out_dir / f"{det}_flare_lcs.png")
            plt.close(fig)


def _load_event_xypi(path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Read X, Y, PI from an XMM event file and drop invalid sky positions."""
    with fits.open(path, memmap=True) as hdul:
        hdu = hdul["EVENTS"] if "EVENTS" in hdul else hdul[1]
        x  = np.asarray(hdu.data["X"],  dtype=float)
        y  = np.asarray(hdu.data["Y"],  dtype=float)
        pi = np.asarray(hdu.data["PI"], dtype=float)
    good = np.isfinite(x) & np.isfinite(y) & (np.abs(x) < 1e7) & (np.abs(y) < 1e7)
    return x[good], y[good], pi[good]


def _prepare_mosaic(events: list[Path], draw_boxes: bool):
    """Read events, return (xs, ys, pis, extent, nx, ny, boxes) or None.

    Used by both the repro/frames loE/hiE renderer and the clean per-band
    renderer so all stages share the same extent/grid/box logic.
    """
    events = [p for p in events if p.is_file()]
    if not events:
        return None
    xs, ys, pis = zip(*(_load_event_xypi(p) for p in events))
    if not any(arr.size for arr in xs):
        return None
    x_all = np.concatenate(xs); y_all = np.concatenate(ys)
    xlo, xhi = float(x_all.min()), float(x_all.max())
    ylo, yhi = float(y_all.min()), float(y_all.max())
    pad_x = max(0.05 * (xhi - xlo), 250.0)
    pad_y = max(0.05 * (yhi - ylo), 250.0)
    extent = (xlo - pad_x, xhi + pad_x, ylo - pad_y, yhi + pad_y)
    nx = 900
    ny = max(200, min(1200,
        int(round(nx * (extent[3] - extent[2]) / (extent[1] - extent[0])))))
    boxes = ([(float(x.min()), float(x.max()), float(y.min()), float(y.max()))
              for x, y in zip(xs, ys) if x.size > 0]
             if draw_boxes else [])
    return list(xs), list(ys), list(pis), extent, nx, ny, boxes


def _hist2d(xs: list[np.ndarray], ys: list[np.ndarray],
            masks: list[np.ndarray], nx: int, ny: int,
            extent: tuple[float, float, float, float]
            ) -> tuple[np.ndarray, int]:
    img = np.zeros((ny, nx), dtype=float)
    total = 0
    for x, y, m in zip(xs, ys, masks):
        if not m.any():
            continue
        h, _xe, _ye = np.histogram2d(
            x[m], y[m], bins=(nx, ny),
            range=((extent[0], extent[1]), (extent[2], extent[3])),
        )
        img += h.T
        total += int(m.sum())
    return img, total


def _save_mosaic(img: np.ndarray, extent: tuple[float, float, float, float],
                 boxes: list[tuple[float, float, float, float]],
                 out_path: Path, title: str,
                 track_xy: tuple[np.ndarray, np.ndarray] | None = None,
                 *, vmin: float | None = None, vmax: float | None = None) -> None:
    ny, nx = img.shape
    with plt.style.context("dark_background"):
        fig, ax = plt.subplots(figsize=(9, 9 * ny / nx + 0.6), dpi=120)
        if vmin is None:
            vmin = max(img.max() * 1e-4, 0.5)
        if vmax is None:
            vmax = max(img.max(), 1.0)
        norm = LogNorm(vmin=vmin, vmax=vmax)
        im = ax.imshow(img, origin="lower", cmap="gray", norm=norm,
                       extent=extent, aspect="equal", interpolation="nearest")
        fig.colorbar(im, ax=ax, label="counts/bin")
        for x0, x1, y0, y1 in boxes:
            ax.add_patch(Rectangle(
                (x0, y0), x1 - x0, y1 - y0,
                fill=False, edgecolor="#ffd24a", linewidth=0.4, alpha=0.5,
            ))
        if track_xy is not None:
            tx, ty = track_xy
            ax.plot(tx, ty, color="#33ddff", linewidth=1.0, alpha=0.9, zorder=5)
            ax.scatter([tx[0]],  [ty[0]],  marker="o", s=24,
                       facecolor="#33ff66", edgecolor="black", zorder=6,
                       label="track start")
            ax.scatter([tx[-1]], [ty[-1]], marker="X", s=28,
                       facecolor="#ff8833", edgecolor="black", zorder=6,
                       label="track end")
            ax.legend(loc="upper right", fontsize=7, framealpha=0.5)
        ax.set_xlabel("X [det]"); ax.set_ylabel("Y [det]")
        ax.set_title(title)
        fig.tight_layout()
        fig.savefig(out_path)
        plt.close(fig)


def _boxes_from_events(events: list[Path]) -> list[tuple[float, float, float, float]]:
    """One sky-X/Y bounding box per event list (used as per-frame outlines)."""
    out: list[tuple[float, float, float, float]] = []
    for p in events:
        if not p.is_file():
            continue
        x, y, _ = _load_event_xypi(p)
        if x.size:
            out.append((float(x.min()), float(x.max()),
                        float(y.min()), float(y.max())))
    return out


def _event_mosaics(events_by_det: dict[str, list[Path]], out_dir: Path,
                   split_ev: float, *,
                   box_events_by_det: dict[str, list[Path]] | None = None,
                   track_xy_by_det: dict[str, tuple[np.ndarray, np.ndarray]] | None = None
                   ) -> None:
    """Render one loE and one hiE mosaic per detector (broad PI split).

    Each detector gets its own extent and stretch (their FOVs and sensitivities
    differ enough that a shared one is misleading).

    box_events_by_det optionally provides the event lists used to compute
    per-frame bounding boxes; the natural source from clean and downstream is
    the upstream frames manifest, so each frame contributes one outline even
    after band/cleaning splits the events into multiple files.

    track_xy_by_det optionally provides per-detector projected-ephemeris
    (xs, ys) arrays drawn over each mosaic.
    """
    summary: list[str] = [f"qc_split_ev\t{split_ev:g}"]
    for det, events in events_by_det.items():
        prepared = _prepare_mosaic(events, draw_boxes=False)
        if prepared is None:
            continue
        xs, ys, pis, extent, nx, ny, _ = prepared
        boxes = (_boxes_from_events(box_events_by_det[det])
                 if box_events_by_det and det in box_events_by_det else [])
        track_xy = track_xy_by_det.get(det) if track_xy_by_det else None
        summary.append(
            f"{det}\tevent_lists={len(xs)}\textent={extent}\t"
            f"pixels={nx}x{ny}\tboxes={len(boxes)}"
        )
        # loE/hiE here are the *broad* PI-split QC bands -- distinct from the
        # cut stage's "soft"/"hard" science bands (200-1000 / 1001-4000 eV).
        for tag, lo, hi in (("loE", None, split_ev), ("hiE", split_ev, None)):
            masks = []
            for pi in pis:
                k = np.ones_like(pi, dtype=bool)
                if lo is not None: k &= pi > lo
                if hi is not None: k &= pi < hi
                masks.append(k)
            img, total = _hist2d(xs, ys, masks, nx, ny, extent)
            band = (f"PI<{split_ev:g}" if tag == "loE" else f"PI>{split_ev:g}")
            title = (f"{det}  {tag} ({band})  {total:,} events"
                     + (f"  frames={len(boxes)}" if boxes else ""))
            _save_mosaic(img, extent, boxes,
                         out_dir / f"{det}_{tag}_mosaic.png", title,
                         track_xy=track_xy)
            summary.append(f"{det}\t{tag}_events={total}")
    (out_dir / "mosaic_summary.txt").write_text(
        "\n".join(summary) + "\n", encoding="utf-8",
    )


def main(argv: list[str]) -> int:
    if len(argv) < 2:
        die(__doc__)
    cmd, args = argv[1], argv[2:]
    if cmd == "shell" and len(args) == 1:
        shell(*args)
    elif cmd == "rewrite-sum-path" and len(args) == 2:
        rewrite_sum_path(*args)
    elif cmd == "repro-manifest" and len(args) == 3:
        repro_manifest(*args)
    elif cmd == "frames-split" and len(args) == 6:
        frames_split(*args)
    elif cmd == "clean-band-table" and len(args) == 1:
        clean_band_table(*args)
    elif cmd == "clean-manifest" and len(args) == 2:
        clean_manifest(*args)
    elif cmd == "fits-rows" and len(args) == 1:
        fits_rows(*args)
    elif cmd == "qc-init" and len(args) in {2, 3}:
        qc_init(args[0], args[1], args[2] if len(args) == 3 else None)
    elif cmd == "qc-repro" and len(args) in {4, 5}:
        qc_repro(args[0], args[1], args[2], args[3],
                 args[4] if len(args) == 5 else None)
    elif cmd == "qc-frames" and len(args) == 4:
        qc_frames(*args)
    elif cmd == "qc-clean" and len(args) == 5:
        qc_clean(*args)
    elif cmd == "cut-run" and len(args) == 4:
        cut_run(*args)
    elif cmd == "qc-cut" and len(args) == 4:
        qc_cut(*args)
    elif cmd == "build-track" and len(args) == 9:
        build_track(*args)
    elif cmd == "qc-track" and len(args) == 5:
        qc_track(*args)
    elif cmd == "maps-counts" and len(args) == 4:
        maps_counts(*args)
    elif cmd == "maps-exposure" and len(args) == 5:
        maps_exposure(*args)
    elif cmd == "maps-background" and len(args) == 2:
        maps_background(*args)
    elif cmd == "maps-corrected" and len(args) == 2:
        maps_corrected(*args)
    elif cmd == "qc-maps" and len(args) == 4:
        qc_maps(*args)
    else:
        die(__doc__)
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
