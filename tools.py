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
  cut-run CLEANDIR FRAMESDIR OUTDIR LOGDIR CONFIG [DETECTORS]
                                              Build cut PI-bands; emanom CCD drop on MOS.
  qc-cut CUTDIR FRAMESDIR DETECTORS OUTDIR    Cut QC: per-(det,band) mosaics + frame boxes.
  maps-counts CUTDIR MAPSDIR FRAMESDIR LOGDIR CONFIG [DETECTORS]
                                              evselect per-frame counts images on a single
                                              per-(det,frame) grid (broadband-derived).
  maps-exposure CUTDIR MAPSDIR ATTHK LOGDIR CONFIG
                                              eexpmap per-frame vignetted exposure maps.
  maps-background MAPSDIR CONFIG              a + b*E fit per frame on off-source pixels.
  maps-corrected MAPSDIR CONFIG               (counts - bkg) / exp_vig per frame.
  qc-maps MAPSDIR FRAMESDIR DETECTORS OUTDIR  Per-(det,band) nobkg-corrected + corrected mosaics.
  cheese-detect MAPSDIR CUTDIR CHEESEDIR LOGDIR CONFIG [DETECTORS]
                                              Detect sources via SAS eboxdetect per frame on
                                              cheese_detection_detectors (default PN); write
                                              sources.tsv + mask.fits.
  cheese-mask MAPSDIR CUTDIR CHEESEDIR CONFIG [DETECTORS]
                                              Apply cheese mask to events + counts/exp/bkg per frame.
  qc-cheese CHEESEDIR FRAMESDIR DETECTORS OUTDIR
                                              Per-(det,band) cheesed mosaics with sources circled.
  stack-events CHEESEDIR TRACK_FITS TRACK_ENV STACKDIR CONFIG [DETECTORS]
                                              Per-event time-resolved shift to comet co-moving frame.
  stack-coadd CHEESEDIR STACKDIR TRACK_FITS TRACK_ENV CONFIG [DETECTORS]
                                              Per-(det,band) stacked counts/exp/bkg/corrected.
  stack-merge STACKDIR CONFIG [DETECTORS]
                                              Pure-Python merge of per-detector stacks onto a
                                              common (union) grid; bilinear resample and sum;
                                              writes merged_<band>_{counts,exp_vig,bkg,corrected}.
  qc-stack STACKDIR DETECTORS OUTDIR          Per-(det,band) and merged stacked mosaics.
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
from scipy.ndimage import correlate as _ndi_correlate
from astropy.io import fits
from astropy.time import Time, TimeDelta
from astropy.wcs import WCS

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, Normalize
from matplotlib.patches import Circle, Rectangle

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


# Structural FITS keys astropy regenerates from the data dtype on write —
# never copy these from a source header verbatim.
_STRUCT_KEYS = frozenset({"SIMPLE", "BITPIX", "NAXIS", "NAXIS1", "NAXIS2", "EXTEND"})

# Column-WCS prefixes (TCTYP6 etc. on EVENTS) astropy drops when the
# table data is reassigned; preserved by _write_filtered_events.
_COL_WCS_PREFIXES = ("TCTYP", "TCRVL", "TCRPX", "TCDLT", "TCUNI")


def _resolve_detectors(cfg: dict[str, Any], detectors_arg: str) -> list[str]:
    """Pick detector list, honouring an optional --detectors override."""
    return (normalize_detectors(detectors_arg) if detectors_arg.strip()
            else normalize_detectors(cfg.get("detectors", "PN")))


def _image_header(src_header: fits.Header) -> fits.Header:
    """Copy a FITS header, dropping the structural keys astropy will rewrite."""
    new = fits.Header()
    for key in src_header:
        if key in _STRUCT_KEYS:
            continue
        try:
            new[key] = src_header[key]
        except (ValueError, KeyError):
            continue
    return new


def _update_cut_dss(out_path: Path, pi_lo: int, pi_hi: int,
                    pattern_hi: int) -> None:
    """Rewrite the EVENTS-extension DSS PI and PATTERN values in-place.

    Cut event files inherit their HDUs (and DSS) from the first clean-band
    sub-file concatenated in (clean-band 1000), which describes only that
    sub-band's filter — wrong for the wider cut band. Downstream SAS tools
    that read DSS (arfgen, rmfgen, evselect with DSS-aware expressions)
    would otherwise see e.g. PI=201:500 on a cut file that actually spans
    201..1000 (soft) or 1001..4000 (hard).
    """
    with fits.open(out_path, mode="update") as hdul:
        evt = next(h for h in hdul if h.name == "EVENTS")
        hdr = evt.header
        for k in list(hdr):
            if not (k.startswith("DSTYP") and k[5:].isdigit()):
                continue
            idx = k[5:]
            v = str(hdr[k]).strip().upper()
            if v == "PI":
                hdr[f"DSVAL{idx}"] = f"{pi_lo}:{pi_hi}"
            elif v == "PATTERN":
                hdr[f"DSVAL{idx}"] = f"0:{pattern_hi}"
        hdr["BNDDSSFX"] = ("yes",
                           "PI/PATTERN DSS rewritten to cut-band range")


def _write_filtered_events(src_path: Path, kept_rows: np.ndarray,
                           out_path: Path) -> None:
    """Write src's full HDU structure to out_path with EVENTS data swapped
    for kept_rows; re-stamp column-WCS keys (TCTYP6 etc.) that astropy
    drops on data reassignment. Used by both cut and cheese for events."""
    with fits.open(src_path, memmap=False) as src_hdul:
        new_hdulist = fits.HDUList([h.copy() for h in src_hdul])
    evt = next(h for h in new_hdulist if h.name == "EVENTS")
    col_wcs = {k: evt.header[k] for k in list(evt.header)
               if any(k.startswith(p) for p in _COL_WCS_PREFIXES)}
    evt.data = kept_rows
    for k, v in col_wcs.items():
        evt.header[k] = v
    out_path.parent.mkdir(parents=True, exist_ok=True)
    new_hdulist.writeto(out_path, overwrite=True)


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


_SCIENCE_FILTERS = {"thin1", "thin2", "medium", "thick"}


def _frame_filter(clean_event_path: Path) -> str:
    """Read the FILTER keyword from the cleaned event file's primary header."""
    with fits.open(clean_event_path, memmap=True) as hdul:
        return str(hdul[0].header.get("FILTER", "")).strip()


def _frame_ontime(clean_event_path: Path) -> float:
    """Read ONTIME from the EVENTS extension; fall back to TSTOP-TSTART."""
    with fits.open(clean_event_path, memmap=True) as hdul:
        h = hdul["EVENTS"].header
        ontime = h.get("ONTIME")
        if ontime is None:
            t0, t1 = h.get("TSTART"), h.get("TSTOP")
            if t0 is not None and t1 is not None:
                return float(t1) - float(t0)
            return 0.0
        return float(ontime)


def _emanom_bad_ccds(frame_event_path: Path, log_dir: Path,
                     drop_codes: set[str]) -> tuple[list[int], str]:
    """Run emanom on a MOS frame event list, return (bad_ccd_numbers, ccd_states).

    Input must be a broad-PI event list (the post-frames file is right —
    a narrow band-1000 cleaned event has too few corner counts to determine
    state). emanom writes per-CCD ANOMFL<n> for CCDs 2-7 (the outer ring;
    CCD1 is the central science chip and is assumed always good).

    drop_codes is a subset of {'G','I','B','O','U'}: 'B' (bad), 'O' (off)
    are the unambiguous reject states; 'U' (undetermined) and 'I'
    (intermediate) are typically kept by default and can be added later if
    a band needs to be conservative.
    """
    if not frame_event_path.is_file():
        die(f"emanom: frame event file missing: {frame_event_path}")
    tag = f"cut_emanom_{frame_event_path.stem}"
    corner_file = log_dir / f"{tag}_corner.fits"
    _run_sas(tag, log_dir, [
        "emanom",
        f"eventfile={frame_event_path}",
        f"cornerfile={corner_file}",
        "keepcorner=no",
        "writekeys=yes",
        "writelog=no",
    ])
    with fits.open(frame_event_path, memmap=True) as hdul:
        hdr = hdul[0].header
        states = "G" + "".join(
            str(hdr.get(f"ANOMFL{ccd}", "U")).strip()[:1].upper()
            for ccd in range(2, 8)
        )
    bad = [i + 1 for i, c in enumerate(states) if c.upper() in drop_codes]
    return bad, states


def _merge_tsv_keep_other(out_path: Path, header: str, key_indices: tuple[int, ...],
                          new_rows: list[str], dets_in_run: set[str]) -> None:
    """Write header + new_rows + any preserved rows from existing file.

    A row is preserved if its first column (the detector) is NOT in
    dets_in_run — this is what supports `--detectors M1,M2` rerunning a
    stage without trampling PN's already-recorded summary rows.
    """
    preserved: list[str] = []
    if out_path.is_file():
        existing = out_path.read_text(encoding="utf-8").splitlines()
        if existing and existing[0] == header:
            for line in existing[1:]:
                parts = line.split("\t")
                if parts and parts[0] not in dets_in_run:
                    preserved.append(line)
    out_path.write_text(
        header + "\n" + "\n".join(preserved + new_rows) + "\n",
        encoding="utf-8",
    )


def cut_run(clean_dir: str, frames_dir: str, out_dir: str, log_dir: str,
            config_path: str, detectors_arg: str = "") -> None:
    """For each (detector, cut band): concat the cleaned per-PPS-band events,
    apply the cut PI window, drop events from MOS CCDs flagged anomalous by
    emanom, sigma-clip events sitting in pixels brighter than median +
    N*sqrt(median) of the combined per-detector mosaic, and write per-frame
    cut FITS preserving the cleaned events' header (so REFXCRVL etc. stay
    valid for downstream WCS).
    """
    cfg = load_config(config_path)
    cut_bands  = _cut_bands(cfg)
    sigma_thr  = float(cfg.get("cut_sigma_clip", 5.0))
    drop_codes = set(str(s).upper() for s in cfg.get("cut_emanom_drop", ["B", "O", "U"]))
    min_ontime = float(cfg.get("cut_min_ontime", 0.0))
    detectors = _resolve_detectors(cfg, detectors_arg)
    clean = Path(clean_dir); frames = Path(frames_dir)
    out = Path(out_dir); log = Path(log_dir)
    (out / "manifest").mkdir(parents=True, exist_ok=True)

    (out / "cut_bands.tsv").write_text(
        "label\tpi_min\tpi_max\n"
        + "\n".join(f"{b['label']}\t{int(b['pi_min'])}\t{int(b['pi_max'])}"
                    for b in cut_bands) + "\n",
        encoding="utf-8",
    )

    summary_rows: list[str] = []
    emanom_rows:  list[str] = []
    summary_header = ("det\tband\tpi_min\tpi_max\tn_in\tn_out\t"
                      "n_clipped\tn_bad_pix\tn_emanom_dropped\tn_frames")
    emanom_header  = "det\tframe\tfilter\tdecision\tccd_states\tbad_ccds"

    for det in detectors:
        frame_paths = _read_manifest(frames / "manifest" / f"{det}_frames.txt")
        if not frame_paths:
            continue
        clean_band_dirs = [d for d in sorted((clean / "events" / det).iterdir())
                           if d.is_dir()]

        # Per-frame filter check + (MOS-only) emanom anomaly state.
        # Frames with non-science FILTER (Closed, CalClosed, CalThick, ...)
        # are dropped entirely from cut so they don't contaminate the maps.
        kept_frames: list[Path] = []
        bad_per_frame: dict[str, list[int]] = {}
        for frame in frame_paths:
            base = frame.stem
            if not frame.is_file():
                continue
            filt = _frame_filter(frame)
            if filt.lower() not in _SCIENCE_FILTERS:
                emanom_rows.append(
                    f"{det}\t{base}\t{filt}\tskipped_non_science\t\t-")
                continue
            ontime = _frame_ontime(frame)
            if ontime < min_ontime:
                emanom_rows.append(
                    f"{det}\t{base}\t{filt}\t"
                    f"skipped_short_ontime[{ontime:.2f}s]\t\t-")
                continue
            kept_frames.append(frame)
            states = ""
            bad: list[int] = []
            if _family(det) == "mos":
                bad, states = _emanom_bad_ccds(frame, log, drop_codes)
                bad_per_frame[base] = bad
            emanom_rows.append(
                f"{det}\t{base}\t{filt}\tkept_science\t{states}\t"
                + (",".join(str(c) for c in bad) if bad else "-")
            )

        for band in cut_bands:
            row = _cut_one_band(det, kept_frames, band, clean_band_dirs,
                                sigma_thr, bad_per_frame, out)
            if row:
                summary_rows.append(row)

    _merge_tsv_keep_other(out / "cut_summary.tsv", summary_header,
                          (0,), summary_rows, set(detectors))
    _merge_tsv_keep_other(out / "cut_emanom.tsv", emanom_header,
                          (0,), emanom_rows, set(detectors))


def _cut_one_band(det: str, frame_paths: list[Path], band: dict[str, Any],
                  clean_band_dirs: list[Path], sigma_thr: float,
                  bad_per_frame: dict[str, list[int]],
                  out_dir: Path) -> str | None:
    label = str(band["label"])
    pi_lo, pi_hi = int(band["pi_min"]), int(band["pi_max"])

    per_frame: list[tuple[Path, np.ndarray, Path]] = []
    n_emanom_dropped = 0
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
        bad_ccds = bad_per_frame.get(base, [])
        if bad_ccds and "CCDNR" in events.dtype.names:
            n_before = len(events)
            events = events[~np.isin(events["CCDNR"], bad_ccds)]
            n_emanom_dropped += n_before - len(events)
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
        _write_filtered_events(src, kept, out_path)
        _update_cut_dss(out_path, pi_lo, pi_hi,
                        4 if _family(det) == "pn" else 12)

    manifest = out_dir / "manifest" / f"{det}_{label}_cut.txt"
    manifest.write_text("\n".join(str(p) for p in out_paths) + "\n",
                        encoding="utf-8")

    return (f"{det}\t{label}\t{pi_lo}\t{pi_hi}\t{n_in}\t{n_out}\t"
            f"{n_clipped}\t{n_bad_pix}\t{n_emanom_dropped}\t{len(per_frame)}")


def _band_mosaics(events_by_det_band: dict[str, dict[str, list[Path]]],
                  out_dir: Path, *,
                  box_events_by_det: dict[str, list[Path]] | None = None,
                  track_xy_by_det: dict[str, tuple[np.ndarray, np.ndarray]] | None = None,
                  circles: list[tuple[float, float, float]] | None = None,
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
                         track_xy=track_xy, circles=circles)
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


def _coadd_band_maps(maps_dir: Path, detectors: list[str], band: str,
                     grids: dict[tuple[str, str, str], dict[str, float]],
                     excluded: set[str] | None = None
                     ) -> tuple[np.ndarray, np.ndarray, float, float, float, int, int]:
    """Co-add per-frame counts and exp_vig across detectors for one band onto
    a sky-fixed viewing grid (aligned to the per-frame physical grids' bin).

    Returns (sum_counts, sum_exp, xmin, ymin, bin, nx, ny).
    """
    rows = _read_maps_manifest_rows(maps_dir)
    band_rows = [r for r in rows
                 if r["band"] == band and r["det"] in detectors
                 and r["counts"] and r["exp_vig"]]
    if excluded:
        band_rows = [r for r in band_rows if r["base"] not in excluded]
    keys = [(r["det"], r["band"], r["base"]) for r in band_rows
            if (r["det"], r["band"], r["base"]) in grids]
    if not keys:
        die(f"cheese: no maps grid info for band={band}; run maps-counts first")
    bin_phys = grids[keys[0]]["bin"]
    xmin = min(grids[k]["xmin"] for k in keys)
    ymin = min(grids[k]["ymin"] for k in keys)

    frame_imgs: list[tuple[dict, np.ndarray, np.ndarray]] = []
    xmax = ymax = float("-inf")
    for r in band_rows:
        key = (r["det"], r["band"], r["base"])
        if key not in grids:
            continue
        g = grids[key]
        c, _ = _read_first2d(Path(r["counts"]))
        e, _ = _read_first2d(Path(r["exp_vig"]))
        fy, fx = c.shape
        xmax = max(xmax, g["xmin"] + fx * bin_phys)
        ymax = max(ymax, g["ymin"] + fy * bin_phys)
        frame_imgs.append((g, c, e))
    nx = int(round((xmax - xmin) / bin_phys))
    ny = int(round((ymax - ymin) / bin_phys))
    sum_c = np.zeros((ny, nx), dtype=float)
    sum_e = np.zeros((ny, nx), dtype=float)
    for g, c, e in frame_imgs:
        ix = int(round((g["xmin"] - xmin) / bin_phys))
        iy = int(round((g["ymin"] - ymin) / bin_phys))
        fy, fx = c.shape
        sum_c[iy:iy+fy, ix:ix+fx] += c
        sum_e[iy:iy+fy, ix:ix+fx] += e
    return sum_c, sum_e, xmin, ymin, bin_phys, nx, ny


def _read_eboxdetect_sources(srclist_path: Path) -> list[dict[str, float]]:
    """Read SAS eboxdetect SRCLIST. Returns list of {ra, dec, like, rate, flux}."""
    if not srclist_path.is_file():
        return []
    with fits.open(srclist_path, memmap=False) as hdul:
        tab = next((h for h in hdul if isinstance(h, fits.BinTableHDU)
                    and getattr(h, "data", None) is not None
                    and len(h.data) > 0), None)
        if tab is None:
            return []
        names = {n.upper(): n for n in tab.columns.names}
        ra_col   = names.get("RA")
        dec_col  = names.get("DEC")
        if ra_col is None or dec_col is None:
            return []
        like_col = names.get("LIKE") or names.get("ML_RATE")
        rate_col = names.get("RATE") or names.get("ML_RATE")
        flux_col = names.get("FLUX")
        out: list[dict[str, float]] = []
        for row in tab.data:
            out.append({
                "ra":   float(row[ra_col]),
                "dec":  float(row[dec_col]),
                "like": float(row[like_col]) if like_col else 0.0,
                "rate": float(row[rate_col]) if rate_col else 0.0,
                "flux": float(row[flux_col]) if flux_col else 0.0,
            })
        return out


def _dedupe_sources(sources: list[dict[str, float]],
                    radius_arcsec: float) -> list[dict[str, float]]:
    """Drop sources within radius_arcsec of an earlier (more significant) one."""
    rad_deg = radius_arcsec / 3600.0
    kept: list[dict[str, float]] = []
    for s in sorted(sources, key=lambda r: -r["like"]):
        cd = math.cos(math.radians(s["dec"]))
        is_dup = False
        for k in kept:
            d_ra  = (s["ra"]  - k["ra"]) * cd
            d_dec = s["dec"] - k["dec"]
            if d_ra * d_ra + d_dec * d_dec < rad_deg * rad_deg:
                is_dup = True
                break
        if not is_dup:
            kept.append(s)
    return kept


def _coadd_for_detection(det_rows: list[dict[str, str]],
                         grids: dict[tuple[str, str, str], dict[str, float]],
                         out_counts: Path, out_exp: Path) -> None:
    """Co-add per-frame counts and exp_vig for one (detector, band) onto a
    single sky-fixed image with SAS-friendly WCS keys (CRPIX/LTV adjusted
    for the new grid origin), so eboxdetect can run on the result.
    """
    keys = [(r["det"], r["band"], r["base"]) for r in det_rows
            if (r["det"], r["band"], r["base"]) in grids]
    bin_phys = grids[keys[0]]["bin"]
    xmin = min(grids[k]["xmin"] for k in keys)
    ymin = min(grids[k]["ymin"] for k in keys)

    sample_header: fits.Header | None = None
    sample_g: dict[str, float] | None = None
    frame_imgs: list[tuple[dict[str, float], np.ndarray, np.ndarray]] = []
    xmax = ymax = float("-inf")
    for r in det_rows:
        key = (r["det"], r["band"], r["base"])
        if key not in grids:
            continue
        g = grids[key]
        c, c_hdr = _read_first2d(Path(r["counts"]))
        e, _ = _read_first2d(Path(r["exp_vig"]))
        if sample_header is None:
            sample_header, sample_g = c_hdr, g
        fy, fx = c.shape
        xmax = max(xmax, g["xmin"] + fx * bin_phys)
        ymax = max(ymax, g["ymin"] + fy * bin_phys)
        frame_imgs.append((g, c, e))
    nx = int(round((xmax - xmin) / bin_phys))
    ny = int(round((ymax - ymin) / bin_phys))
    sum_c = np.zeros((ny, nx), dtype=np.float32)
    sum_e = np.zeros((ny, nx), dtype=np.float32)
    for g, c, e in frame_imgs:
        ix = int(round((g["xmin"] - xmin) / bin_phys))
        iy = int(round((g["ymin"] - ymin) / bin_phys))
        fy, fx = c.shape
        sum_c[iy:iy + fy, ix:ix + fx] += c
        sum_e[iy:iy + fy, ix:ix + fx] += e

    # WCS for the viewing grid: keep the celestial reference (CRVAL/CDELT/
    # CTYPE) from the sample frame; recompute CRPIX and LTV for the new
    # grid origin so SAS image-pixel↔physical-coord conversions stay
    # consistent.
    assert sample_header is not None and sample_g is not None
    sample_crpix1 = float(sample_header["CRPIX1"])
    sample_crpix2 = float(sample_header["CRPIX2"])
    refxcrpx = (sample_crpix1 - 0.5) * bin_phys + sample_g["xmin"]
    refycrpx = (sample_crpix2 - 0.5) * bin_phys + sample_g["ymin"]
    new_hdr = _image_header(sample_header)
    new_hdr["CRPIX1"] = (refxcrpx - xmin) / bin_phys + 0.5
    new_hdr["CRPIX2"] = (refycrpx - ymin) / bin_phys + 0.5
    new_hdr["LTV1"] = 0.5 - xmin / bin_phys
    new_hdr["LTV2"] = 0.5 - ymin / bin_phys

    fits.PrimaryHDU(data=sum_c, header=new_hdr.copy()).writeto(out_counts, overwrite=True)
    fits.PrimaryHDU(data=sum_e, header=new_hdr.copy()).writeto(out_exp,    overwrite=True)


def cheese_detect(maps_dir: str, cut_dir: str, cheese_dir: str,
                  log_dir: str, config_path: str,
                  detectors_arg: str = "") -> None:
    """Stage 1: run SAS eboxdetect on every per-frame
    (counts, exp_vig, bkg) image of each detection-detector (default PN
    only), union the source lists, and write sources.tsv + sky-fixed
    mask.fits on the same physical grid as the maps. The mask is then
    applied to all detectors at stage 2 — PN is the most sensitive
    detector and dominates source recall.
    """
    cfg = load_config(config_path)
    detection_dets = normalize_detectors(
        cfg.get("cheese_detection_detectors", ["PN"]))
    band = str(cfg.get("cheese_band", "hard"))
    radius_arcsec = float(cfg.get("cheese_radius_arcsec", 25.0))
    likemin = float(cfg.get("cheese_eboxdetect_likemin", 8))
    boxsize = int(cfg.get("cheese_eboxdetect_boxsize", 5))
    nruns   = int(cfg.get("cheese_eboxdetect_nruns", 3))
    excluded = set(str(s) for s in cfg.get("excluded_frames", []))

    cut_bands_cfg = _cut_bands(cfg)
    band_info = next((b for b in cut_bands_cfg if b["label"] == band), None)
    if band_info is None:
        die(f"cheese_band {band!r} is not in cut_bands")
    pi_min = int(band_info["pi_min"])
    pi_max = int(band_info["pi_max"])

    maps = Path(maps_dir); cheese = Path(cheese_dir); log = Path(log_dir)
    cheese.mkdir(parents=True, exist_ok=True)
    src_dir = cheese / "srclists"; src_dir.mkdir(parents=True, exist_ok=True)

    # Per-frame eboxdetect on each detection detector's hard-band products.
    rows = _read_maps_manifest_rows(maps)
    # Per-detector co-added counts/exp_vig (sky-fixed grid, full SAS WCS),
    # then ONE eboxdetect call per detection-detector. This gives the depth
    # of the whole observation in a single image so likemin=8 picks out
    # real sources without per-frame noise contamination.
    grids = _read_maps_grid(maps)
    coadd_dir = cheese / "coadd"; coadd_dir.mkdir(parents=True, exist_ok=True)
    raw_sources: list[dict[str, float]] = []
    for det in detection_dets:
        det_rows = [r for r in rows
                    if r["det"] == det and r["band"] == band
                    and r["counts"] and r["exp_vig"]]
        if excluded:
            before = len(det_rows)
            det_rows = [r for r in det_rows if r["base"] not in excluded]
            dropped = before - len(det_rows)
            if dropped:
                print(f"  cheese-detect: {det} excluding {dropped} "
                      f"frame(s) per cfg.excluded_frames", flush=True)
        det_rows = [
            r for r in det_rows
            if (lambda img: bool((img > 0).any()))(
                _read_first2d(Path(r["exp_vig"]))[0])
        ]
        if not det_rows:
            continue
        co_counts = coadd_dir / f"{det}_{band}_counts.fits"
        co_exp    = coadd_dir / f"{det}_{band}_exp_vig.fits"
        _coadd_for_detection(det_rows, grids, co_counts, co_exp)
        srclist = src_dir / f"{det}_{band}_eboxlist.fits"
        if not srclist.is_file():
            _run_sas(f"cheese_eboxdetect_{det}_{band}", log, [
                "eboxdetect",
                f"imagesets={co_counts}",
                f"expimagesets={co_exp}",
                "withexpimage=yes",
                "usemap=no",
                f"likemin={likemin}",
                f"boxsize={boxsize}",
                f"nruns={nruns}",
                f"pimin={pi_min}",
                f"pimax={pi_max}",
                f"boxlistset={srclist}",
            ])
        raw_sources.extend(_read_eboxdetect_sources(srclist))

    # Mask-grid: union of every detector's frame footprints (so we cover
    # the whole observation, not just the detection-detector's FOV).
    apply_dets = _resolve_detectors(cfg, detectors_arg)
    _, _, xmin, ymin, bin_phys, nx, ny = _coadd_band_maps(
        maps, apply_dets, band, grids, excluded=excluded)

    sample_counts = next(
        (Path(r["counts"]) for r in rows
         if r["counts"] and Path(r["counts"]).is_file()), None,
    )
    if sample_counts is None:
        die("cheese: no maps counts files on disk")
    with fits.open(sample_counts, memmap=True) as h:
        sample_header = h[0].header.copy()
    pixel_arcsec = abs(float(sample_header["CDELT1"])) * 3600.0
    radius_pix   = radius_arcsec / pixel_arcsec
    radius_phys  = radius_pix * bin_phys

    # Dedupe by RA/DEC within the mask radius.
    deduped = _dedupe_sources(raw_sources, radius_arcsec)

    cut_event = next((f for f in Path(cut_dir).rglob("*_cut.fits")
                      if f.is_file()), None)
    if cut_event is None:
        die("cheese: no cut events to anchor RA/DEC -> physical X/Y WCS")
    wcs = _event_xy_wcs(cut_event)

    yy, xx = np.indices((ny, nx))
    mask = np.ones((ny, nx), dtype=np.uint8)
    src_lines = ["det\tx_pix\ty_pix\tx_phys\ty_phys\tra\tdec\t"
                 "radius_pix\tradius_phys\tlike\trate\tflux"]
    for s in deduped:
        x_phys, y_phys = (float(v) for v in
                          wcs.wcs_world2pix(s["ra"], s["dec"], 1))
        x_pix = (x_phys - xmin) / bin_phys - 0.5
        y_pix = (y_phys - ymin) / bin_phys - 0.5
        if 0 <= x_pix < nx and 0 <= y_pix < ny:
            mask[(yy - y_pix) ** 2 + (xx - x_pix) ** 2
                 <= radius_pix * radius_pix] = 0
        src_lines.append(
            f"combined\t{int(round(x_pix))}\t{int(round(y_pix))}\t"
            f"{x_phys:.2f}\t{y_phys:.2f}\t"
            f"{s['ra']:.6f}\t{s['dec']:.6f}\t"
            f"{radius_pix:.3f}\t{radius_phys:.2f}\t"
            f"{s['like']:.2f}\t{s['rate']:.6g}\t{s['flux']:.6g}"
        )
    (cheese / "sources.tsv").write_text(
        "\n".join(src_lines) + "\n", encoding="utf-8")

    mask_hdu = fits.PrimaryHDU(data=mask, header=_image_header(sample_header))
    mask_hdu.header["MASKXMIN"] = (float(xmin), "physical-X grid origin")
    mask_hdu.header["MASKYMIN"] = (float(ymin), "physical-Y grid origin")
    mask_hdu.header["MASKBIN"]  = (float(bin_phys), "physical units per pixel")
    mask_hdu.header["MASKNX"]   = (int(nx), "image pixels in X")
    mask_hdu.header["MASKNY"]   = (int(ny), "image pixels in Y")
    mask_hdu.header["NSOURCES"] = (len(deduped), "number of detected sources")
    mask_hdu.header["LIKEMIN"]  = (float(likemin), "eboxdetect LIKE threshold")
    mask_hdu.header["RADPIX"]   = (float(radius_pix), "mask radius in pixels")
    mask_hdu.writeto(cheese / "mask.fits", overwrite=True)

    n_excl = int((mask == 0).sum())
    n_total = int(mask.size)
    (cheese / "cheese_summary.tsv").write_text(
        "metric\tvalue\n"
        f"detection_band\t{band}\n"
        f"detection_detectors\t{','.join(detection_dets)}\n"
        f"apply_detectors\t{','.join(apply_dets)}\n"
        f"likemin\t{likemin}\n"
        f"boxsize\t{boxsize}\n"
        f"nruns\t{nruns}\n"
        f"mask_radius_arcsec\t{radius_arcsec}\n"
        f"mask_radius_pix\t{radius_pix:.3f}\n"
        f"n_eboxdetect_runs\t{len(detection_dets)}\n"
        f"n_raw_sources\t{len(raw_sources)}\n"
        f"n_unique_sources\t{len(deduped)}\n"
        f"n_excluded_pixels\t{n_excl}\n"
        f"n_total_pixels\t{n_total}\n"
        f"excluded_fraction\t{n_excl / n_total:.4f}\n",
        encoding="utf-8",
    )


def _read_cheese_mask(cheese_dir: Path
                      ) -> tuple[np.ndarray, float, float, float]:
    mask_path = cheese_dir / "mask.fits"
    if not mask_path.is_file():
        die(f"cheese mask missing: {mask_path}")
    with fits.open(mask_path, memmap=False) as h:
        hdr = h[0].header
        return (np.asarray(h[0].data, dtype=np.uint8),
                float(hdr["MASKXMIN"]), float(hdr["MASKYMIN"]),
                float(hdr["MASKBIN"]))


def _read_sources(cheese_dir: Path) -> list[tuple[float, float, float]]:
    p = cheese_dir / "sources.tsv"
    if not p.is_file():
        return []
    lines = p.read_text(encoding="utf-8").splitlines()
    if len(lines) < 2:
        return []
    head = lines[0].split("\t")
    idx = {n: i for i, n in enumerate(head)}
    out: list[tuple[float, float, float]] = []
    for line in lines[1:]:
        parts = line.split("\t")
        if len(parts) < len(head):
            continue
        out.append((float(parts[idx["x_phys"]]),
                    float(parts[idx["y_phys"]]),
                    float(parts[idx["radius_phys"]])))
    return out


def cheese_mask(maps_dir: str, cut_dir: str, cheese_dir: str,
                config_path: str, detectors_arg: str = "") -> None:
    """Stage 2 of cheese: drop event rows inside any source circle and
    zero the corresponding pixels in per-frame counts/exp_vig/bkg images
    (corrected is recomputed from the masked counts and exp_vig)."""
    cfg = load_config(config_path)
    detectors = _resolve_detectors(cfg, detectors_arg)
    cut_bands = _cut_bands(cfg)
    excluded = set(str(s) for s in cfg.get("excluded_frames", []))

    maps = Path(maps_dir); cut = Path(cut_dir); cheese = Path(cheese_dir)
    sources = _read_sources(cheese)
    mask, mx, my, mbin = _read_cheese_mask(cheese)

    summary_rows: list[str] = []
    for det in detectors:
        for band in cut_bands:
            label = str(band["label"])
            evt_paths = _read_manifest(
                cut / "manifest" / f"{det}_{label}_cut.txt")
            if excluded:
                before = len(evt_paths)
                evt_paths = [
                    p for p in evt_paths
                    if p.stem.removesuffix(f"_{label}_cut") not in excluded
                ]
                dropped = before - len(evt_paths)
                if dropped:
                    print(f"  cheese-mask: {det} {label} excluding {dropped} "
                          f"frame(s) per cfg.excluded_frames", flush=True)
            n_in = n_kept = n_frames = 0
            out_paths: list[Path] = []
            for src_path in evt_paths:
                if not src_path.is_file():
                    continue
                with fits.open(src_path, memmap=True) as h:
                    data = np.asarray(h["EVENTS"].data)
                n_in += len(data)
                if len(data) and sources:
                    xs = np.asarray(data["X"], dtype=float)
                    ys = np.asarray(data["Y"], dtype=float)
                    keep = np.ones(len(data), dtype=bool)
                    for sx, sy, sr in sources:
                        keep &= (xs - sx) ** 2 + (ys - sy) ** 2 > sr * sr
                    kept = data[keep]
                else:
                    kept = data
                n_kept += len(kept)
                n_frames += 1
                out_path = (cheese / "events" / det / label
                            / f"{src_path.stem}_cheesed.fits")
                _write_filtered_events(src_path, kept, out_path)
                out_paths.append(out_path.resolve())
            (cheese / "manifest").mkdir(parents=True, exist_ok=True)
            (cheese / "manifest" / f"{det}_{label}_cheesed.txt").write_text(
                "\n".join(str(p) for p in out_paths)
                + ("\n" if out_paths else ""),
                encoding="utf-8",
            )
            summary_rows.append(
                f"{det}\t{label}\t{n_in}\t{n_kept}\t{n_in - n_kept}\t{n_frames}"
            )

    grids = _read_maps_grid(maps)
    map_summary: list[str] = []
    for r in _read_maps_manifest_rows(maps):
        det = r["det"]
        if det not in detectors:
            continue
        if r["base"] in excluded:
            continue
        if not (r["counts"] and r["exp_vig"] and r["bkg"]):
            continue
        key = (det, r["band"], r["base"])
        if key not in grids:
            continue
        g = grids[key]
        c, c_hdr = _read_first2d(Path(r["counts"]))
        e, _ = _read_first2d(Path(r["exp_vig"]))
        b, _ = _read_first2d(Path(r["bkg"]))
        ix = int(round((g["xmin"] - mx) / mbin))
        iy = int(round((g["ymin"] - my) / mbin))
        fy, fx = c.shape
        slab = mask[iy:iy + fy, ix:ix + fx].astype(float)
        if slab.shape != c.shape:
            full = np.ones_like(c)
            sy = min(slab.shape[0], fy); sx = min(slab.shape[1], fx)
            full[:sy, :sx] = slab[:sy, :sx]
            slab = full
        out_dir = cheese / "maps" / det / r["band"]
        out_dir.mkdir(parents=True, exist_ok=True)
        new_hdr = _image_header(c_hdr)
        new_hdr["CHEESE"] = ("yes", "point-source mask applied")
        c_out = (c * slab).astype(np.float32)
        e_out = (e * slab).astype(np.float32)
        b_out = (b * slab).astype(np.float32)
        with np.errstate(divide="ignore", invalid="ignore"):
            corr_out = np.where(e_out > 0, (c_out - b_out) / e_out, 0.0
                                ).astype(np.float32)
        for arr, suffix in ((c_out,    "counts"),
                            (e_out,    "exp_vig"),
                            (b_out,    "bkg"),
                            (corr_out, "corrected")):
            out = out_dir / f"{r['base']}_{r['band']}_{suffix}_cheesed.fits"
            fits.PrimaryHDU(data=arr, header=new_hdr.copy()).writeto(
                out, overwrite=True)
        map_summary.append(
            f"{det}\t{r['band']}\t{r['base']}\t"
            f"{out_dir / (r['base'] + '_' + r['band'] + '_counts_cheesed.fits')}\t"
            f"{out_dir / (r['base'] + '_' + r['band'] + '_exp_vig_cheesed.fits')}\t"
            f"{out_dir / (r['base'] + '_' + r['band'] + '_bkg_cheesed.fits')}\t"
            f"{out_dir / (r['base'] + '_' + r['band'] + '_corrected_cheesed.fits')}"
        )

    _merge_tsv_keep_other(
        cheese / "events_summary.tsv",
        "det\tband\tn_in\tn_kept\tn_dropped\tn_frames",
        (0,), summary_rows, set(detectors))
    _merge_tsv_keep_other(
        cheese / "maps_manifest.tsv",
        "det\tband\tbase\tcounts\texp_vig\tbkg\tcorrected",
        (0,), map_summary, set(detectors))


def qc_cheese(cheese_dir: str, frames_dir: str, detectors_arg: str,
              out_dir: str, config_path: str = "") -> None:
    """Per-(det, band) mosaics of the *cut* events (pre-cheese, so the
    underlying source structure is visible) with translucent red circles
    overplotted to show what gets masked out at the cheese stage."""
    cheese = Path(cheese_dir); frames = Path(frames_dir)
    out = Path(out_dir); out.mkdir(parents=True, exist_ok=True)
    detectors = normalize_detectors(detectors_arg)
    excluded: set[str] = set()
    if config_path:
        try:
            excluded = set(str(s) for s in
                           load_config(config_path).get("excluded_frames", []))
        except Exception:
            pass

    files = sorted(p for p in cheese.rglob("*") if p.is_file())
    _write_listing(out / "files.txt", files)
    for tsv in ("sources.tsv", "cheese_summary.tsv",
                "events_summary.tsv", "maps_manifest.tsv"):
        src = cheese / tsv
        if src.is_file():
            (out / tsv).write_bytes(src.read_bytes())

    cut_dir = cheese.parent / "cut"
    events_by_det_band = _cut_events_by_det_band(cut_dir, detectors)
    if excluded:
        for det, by_band in events_by_det_band.items():
            for label, paths in by_band.items():
                before = len(paths)
                by_band[label] = [
                    p for p in paths
                    if p.stem.removesuffix(f"_{label}_cut") not in excluded
                ]
                dropped = before - len(by_band[label])
                if dropped:
                    print(f"  qc-cheese: {det} {label} excluding {dropped} "
                          f"frame(s) per cfg.excluded_frames", flush=True)
    box_events = {det: _read_manifest(frames / "manifest" / f"{det}_frames.txt")
                  for det in detectors}
    circles = _read_sources(cheese)
    _band_mosaics(events_by_det_band, out, box_events_by_det=box_events,
                  circles=circles)


# ============================================================
# Stack stage
# ============================================================

def _read_track_env(env_path: Path) -> tuple[float, float]:
    """Parse output/track/track.env for COMET_REF_RA / COMET_REF_DEC."""
    ra = dec = None
    for line in Path(env_path).read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if line.startswith("export COMET_REF_RA="):
            ra = float(line.split("=", 1)[1].strip().strip('"'))
        elif line.startswith("export COMET_REF_DEC="):
            dec = float(line.split("=", 1)[1].strip().strip('"'))
    if ra is None or dec is None:
        die(f"track env missing COMET_REF_RA / COMET_REF_DEC: {env_path}")
    return ra, dec


def _xmm_seconds_to_mjd_utc(seconds: np.ndarray, mjdref: float,
                            timesys: str) -> np.ndarray:
    """Convert XMM mission-elapsed-time TIME (in TIMESYS scale) to MJD UTC."""
    t = (Time(mjdref, format="mjd", scale=timesys.lower())
         + TimeDelta(np.asarray(seconds, dtype=float), format="sec"))
    return np.asarray(t.utc.mjd, dtype=float)


def _interpolate_track(mjd_target: np.ndarray, mjd_track: np.ndarray,
                       ra_track: np.ndarray, dec_track: np.ndarray
                       ) -> tuple[np.ndarray, np.ndarray]:
    """Interpolate track RA/DEC at the requested MJDs. Track must be sorted
    ascending; RA is unwrapped before interp to handle 0/360 boundary."""
    ra_unwrapped = np.unwrap(np.deg2rad(ra_track))
    ra_at_t  = np.rad2deg(np.interp(mjd_target, mjd_track, ra_unwrapped)) % 360.0
    dec_at_t = np.interp(mjd_target, mjd_track, dec_track)
    return ra_at_t, dec_at_t


def _shift_events(src_path: Path, dst_path: Path,
                  mjd_track: np.ndarray, ra_track: np.ndarray,
                  dec_track: np.ndarray, ref_ra: float, ref_dec: float,
                  wcs: WCS) -> tuple[int, float, float]:
    """Apply a per-event time-resolved shift so the comet stays at the
    track's reference RA/DEC. Returns (n_events, mid_dx_phys, mid_dy_phys).
    """
    with fits.open(src_path, memmap=True) as h:
        evt_hdr = h["EVENTS"].header
        data = np.asarray(h["EVENTS"].data)
        mjdref = float(evt_hdr.get("MJDREF",
                       float(evt_hdr.get("MJDREFI", 0))
                       + float(evt_hdr.get("MJDREFF", 0))))
        timesys = str(evt_hdr.get("TIMESYS", "TT")).strip().upper()
    n = len(data)
    if n == 0:
        _write_filtered_events(src_path, data, dst_path)
        return 0, 0.0, 0.0

    times = np.asarray(data["TIME"], dtype=float)
    mjd_utc = _xmm_seconds_to_mjd_utc(times, mjdref, timesys)
    ra_t, dec_t = _interpolate_track(mjd_utc, mjd_track, ra_track, dec_track)
    x_c, y_c = wcs.wcs_world2pix(ra_t, dec_t, 1)
    x_ref, y_ref = (float(v) for v in wcs.wcs_world2pix(ref_ra, ref_dec, 1))

    new_data = data.copy()
    new_data["X"] = data["X"] - x_c + x_ref
    new_data["Y"] = data["Y"] - y_c + y_ref
    _write_filtered_events(src_path, new_data, dst_path)

    # Mid-window shift used by the per-frame map shift in stack-coadd.
    mid_idx = n // 2
    return n, float(x_c[mid_idx] - x_ref), float(y_c[mid_idx] - y_ref)


def stack_events(cheese_dir: str, track_fits: str, track_env: str,
                 stack_dir: str, config_path: str,
                 detectors_arg: str = "") -> None:
    """For each cheesed cut event file, shift X/Y per event so the comet
    stays at the track reference RA/DEC. Output goes to
    stack/events/<DET>/<BAND>/<frame>_<band>_stack.fits with the original
    HDU structure preserved.
    """
    cfg = load_config(config_path)
    cut_bands = _cut_bands(cfg)
    detectors = _resolve_detectors(cfg, detectors_arg)
    excluded = set(str(s) for s in cfg.get("excluded_frames", []))

    cheese = Path(cheese_dir); stack = Path(stack_dir)
    stack.mkdir(parents=True, exist_ok=True)
    (stack / "events").mkdir(parents=True, exist_ok=True)
    (stack / "manifest").mkdir(parents=True, exist_ok=True)

    mjd_track, ra_track, dec_track, _ = _read_track(track_fits)
    ref_ra, ref_dec = _read_track_env(Path(track_env))

    sample_event = next(
        (p for det in detectors for label in (b["label"] for b in cut_bands)
         for p in _read_manifest(cheese / "manifest"
                                 / f"{det}_{label}_cheesed.txt")
         if p.is_file()),
        None,
    )
    if sample_event is None:
        die("stack-events: no cheesed events on disk; run cheese first")
    wcs = _event_xy_wcs(sample_event)

    summary_rows = ["det\tband\tn_events\tn_frames\tmid_dx_phys\tmid_dy_phys"]
    for det in detectors:
        for band in cut_bands:
            label = str(band["label"])
            srcs = _read_manifest(
                cheese / "manifest" / f"{det}_{label}_cheesed.txt")
            if excluded:
                before = len(srcs)
                srcs = [
                    p for p in srcs
                    if p.stem.removesuffix("_cheesed").removesuffix(
                        f"_{label}_cut") not in excluded
                ]
                dropped = before - len(srcs)
                if dropped:
                    print(f"  stack-events: {det} {label} excluding {dropped} "
                          f"frame(s) per cfg.excluded_frames", flush=True)
            out_dir = stack / "events" / det / label
            out_dir.mkdir(parents=True, exist_ok=True)
            paths: list[Path] = []
            tot_events = 0
            sum_dx = sum_dy = 0.0
            for src_path in srcs:
                if not src_path.is_file():
                    continue
                base = src_path.stem
                if base.endswith("_cheesed"):
                    base = base[: -len("_cheesed")]
                if base.endswith(f"_{label}_cut"):
                    base = base[: -len(f"_{label}_cut")]
                dst = out_dir / f"{base}_{label}_stack.fits"
                n, dx, dy = _shift_events(src_path, dst,
                                          mjd_track, ra_track, dec_track,
                                          ref_ra, ref_dec, wcs)
                paths.append(dst.resolve())
                tot_events += n
                sum_dx += dx; sum_dy += dy
            (stack / "manifest" / f"{det}_{label}_stack.txt").write_text(
                "\n".join(str(p) for p in paths)
                + ("\n" if paths else ""), encoding="utf-8",
            )
            n_frames = len(paths)
            summary_rows.append(
                f"{det}\t{label}\t{tot_events}\t{n_frames}\t"
                f"{(sum_dx / max(n_frames, 1)):.2f}\t"
                f"{(sum_dy / max(n_frames, 1)):.2f}"
            )
    _merge_tsv_keep_other(
        stack / "events_summary.tsv",
        "det\tband\tn_events\tn_frames\tmid_dx_phys\tmid_dy_phys",
        (0,), summary_rows[1:], set(detectors),
    )


def _stack_grid_extent(per_frame_offsets: list[tuple[dict[str, float], float, float]],
                       bin_phys: float
                       ) -> tuple[float, float, int, int]:
    """Given a list of (frame-grid-dict, dx_phys, dy_phys) offsets to apply
    in the stack frame, return (xmin, ymin, nx, ny) of the union extent."""
    xmin = min(g["xmin"] - dx for g, dx, _ in per_frame_offsets)
    ymin = min(g["ymin"] - dy for g, _, dy in per_frame_offsets)
    xmax = max(g["xmin"] - dx + g["nx"] * bin_phys
               for g, dx, _ in per_frame_offsets)
    ymax = max(g["ymin"] - dy + g["ny"] * bin_phys
               for g, _, dy in per_frame_offsets)
    nx = int(round((xmax - xmin) / bin_phys))
    ny = int(round((ymax - ymin) / bin_phys))
    return xmin, ymin, nx, ny


def _frame_midtime_mjd_utc(event_path: Path) -> float:
    """Return the MJD UTC at the midpoint of an event file's [TSTART, TSTOP]."""
    with fits.open(event_path, memmap=True) as h:
        hdr = h["EVENTS"].header
        tstart = float(hdr["TSTART"])
        tstop = float(hdr["TSTOP"])
        mjdref = float(hdr.get("MJDREF",
                               float(hdr.get("MJDREFI", 0))
                               + float(hdr.get("MJDREFF", 0))))
        timesys = str(hdr.get("TIMESYS", "TT")).strip().upper()
    return float(_xmm_seconds_to_mjd_utc(
        np.array([0.5 * (tstart + tstop)]), mjdref, timesys)[0])


def _trail_kernel(evt_path: Path, mid_mjd: float,
                  mjd_track: np.ndarray, ra_track: np.ndarray,
                  dec_track: np.ndarray, wcs: WCS, bin_phys: float,
                  step_sec: float = 30.0
                  ) -> np.ndarray:
    """Time-weighted distribution of the comet's pixel-offset relative to the
    frame midtime, in stack-grid pixels (one pixel = bin_phys phys-px).
    Sums to 1, odd shape, centered. Sampled across the union of STDGTI*
    intervals so GTI gaps within the frame get zero weight."""
    with fits.open(evt_path, memmap=True) as h:
        evt_hdr = h["EVENTS"].header
        mjdref = float(evt_hdr.get("MJDREF",
                       float(evt_hdr.get("MJDREFI", 0))
                       + float(evt_hdr.get("MJDREFF", 0))))
        timesys = str(evt_hdr.get("TIMESYS", "TT")).strip().upper()
        ivs: list[tuple[float, float]] = []
        for hdu in h[1:]:
            name = (hdu.name or "").upper()
            if name.startswith("STDGTI") and hdu.data is not None:
                for s, e in zip(hdu.data["START"], hdu.data["STOP"]):
                    if e > s:
                        ivs.append((float(s), float(e)))
    if not ivs:
        return np.array([[1.0]], dtype=np.float32)

    samples = []
    for s0, s1 in ivs:
        n = max(2, int(np.ceil((s1 - s0) / step_sec)))
        samples.append(np.linspace(s0, s1, n))
    times = np.concatenate(samples)
    mjd_utc = _xmm_seconds_to_mjd_utc(times, mjdref, timesys)
    ra_t, dec_t = _interpolate_track(mjd_utc, mjd_track, ra_track, dec_track)
    x_t, y_t = wcs.wcs_world2pix(ra_t, dec_t, 1)

    ra_m, dec_m = _interpolate_track(np.array([mid_mjd]),
                                     mjd_track, ra_track, dec_track)
    x_m, y_m = (float(v) for v in wcs.wcs_world2pix(ra_m[0], dec_m[0], 1))

    dx = (np.asarray(x_t, dtype=float) - x_m) / bin_phys
    dy = (np.asarray(y_t, dtype=float) - y_m) / bin_phys
    half_x = int(np.ceil(max(abs(float(dx.min())),
                             abs(float(dx.max())))))
    half_y = int(np.ceil(max(abs(float(dy.min())),
                             abs(float(dy.max())))))
    nx_k, ny_k = 2 * half_x + 1, 2 * half_y + 1
    kernel, _, _ = np.histogram2d(
        dx, dy, bins=(nx_k, ny_k),
        range=((-half_x - 0.5, half_x + 0.5),
               (-half_y - 0.5, half_y + 0.5)))
    kernel = kernel.T  # (ny, nx)
    s = float(kernel.sum())
    if s <= 0:
        return np.array([[1.0]], dtype=np.float32)
    return (kernel / s).astype(np.float32)


def stack_coadd(cheese_dir: str, stack_dir: str, track_fits: str,
                track_env: str, config_path: str,
                detectors_arg: str = "") -> None:
    """Build the per-(det, band) stacked counts/exp_vig/bkg/corrected maps
    in the comet co-moving frame.

    counts come from the per-event-shifted events written by stack-events
    (binned via histogram2d). exp_vig and bkg come from the cheesed maps
    integer-pixel-shifted by the comet position at each frame's midtime
    (good enough at 4 arcsec/pix; comet drifts ~1 pixel per frame).
    """
    cfg = load_config(config_path)
    cut_bands = _cut_bands(cfg)
    detectors = _resolve_detectors(cfg, detectors_arg)
    quantile = float(cfg.get("exp_floor_quantile", 0.01))
    min_exp_max = float(cfg.get("stack_min_exp_max", 0.0))
    excluded = set(str(s) for s in cfg.get("excluded_frames", []))

    cheese = Path(cheese_dir); stack = Path(stack_dir)
    grids = _read_maps_grid(cheese.parent / "maps")

    mjd_track, ra_track, dec_track, _ = _read_track(track_fits)
    ref_ra, ref_dec = _read_track_env(Path(track_env))

    sample_event = next(
        (p for det in detectors for label in (b["label"] for b in cut_bands)
         for p in _read_manifest(stack / "manifest"
                                 / f"{det}_{label}_stack.txt")
         if p.is_file()),
        None,
    )
    if sample_event is None:
        die("stack-coadd: no shifted events; run stack-events first")
    wcs = _event_xy_wcs(sample_event)
    x_ref, y_ref = (float(v) for v in wcs.wcs_world2pix(ref_ra, ref_dec, 1))

    summary = ["det\tband\tnx\tny\tn_events\ttotal_exposure\tbkg_total"]
    for det in detectors:
        for band in cut_bands:
            label = str(band["label"])
            stack_paths = _read_manifest(
                stack / "manifest" / f"{det}_{label}_stack.txt")
            stack_paths = [p for p in stack_paths if p.is_file()]
            if excluded:
                before = len(stack_paths)
                stack_paths = [
                    p for p in stack_paths
                    if p.stem.removesuffix(f"_{label}_stack") not in excluded
                ]
                dropped = before - len(stack_paths)
                if dropped:
                    print(f"  stack-coadd: {det} {label} excluding {dropped} "
                          f"frame(s) per cfg.excluded_frames", flush=True)
            if not stack_paths:
                continue

            # Drop short-exposure / slew-class frames whose bkg model would
            # paint a near-empty footprint over a large area, leaving a
            # negative halo in the comoving stack. The threshold is on the
            # cheesed exp_vig peak — that's the per-pixel exposure most
            # representative of how much we should trust this frame.
            cheese_maps_dir = cheese / "maps" / det / label
            if min_exp_max > 0:
                kept: list[Path] = []
                for p in stack_paths:
                    base = p.stem.removesuffix(f"_{label}_stack")
                    cexp = (cheese_maps_dir
                            / f"{base}_{label}_exp_vig_cheesed.fits")
                    if cexp.is_file():
                        with fits.open(cexp) as h:
                            e_max = float(np.asarray(h[0].data).max())
                        if e_max < min_exp_max:
                            print(
                                f"  stack-coadd skip {det} {label} {base}: "
                                f"exp_max={e_max:.0f}s < "
                                f"stack_min_exp_max={min_exp_max:.0f}s",
                                flush=True)
                            continue
                    kept.append(p)
                stack_paths = kept
                if not stack_paths:
                    continue

            # Per-frame mid-time shift (for the maps integer shift AND to
            # size the stack grid via the shifted footprints) plus a
            # within-frame trail kernel (for smearing exp_vig/bkg along the
            # comet path so they match the events' per-event-shift geometry).
            offsets: list[tuple[dict[str, float], float, float, Path,
                                np.ndarray]] = []
            max_kx = max_ky = 0
            for evt_path in stack_paths:
                base = evt_path.stem.removesuffix(f"_{label}_stack")
                key = (det, label, base)
                if key not in grids:
                    continue
                mid_mjd = _frame_midtime_mjd_utc(evt_path)
                ra_m, dec_m = _interpolate_track(np.array([mid_mjd]),
                                                 mjd_track, ra_track,
                                                 dec_track)
                xc, yc = (float(v) for v in
                          wcs.wcs_world2pix(ra_m[0], dec_m[0], 1))
                kernel = _trail_kernel(evt_path, mid_mjd, mjd_track,
                                       ra_track, dec_track, wcs,
                                       grids[key]["bin"])
                max_kx = max(max_kx, (kernel.shape[1] - 1) // 2)
                max_ky = max(max_ky, (kernel.shape[0] - 1) // 2)
                offsets.append((grids[key], xc - x_ref, yc - y_ref,
                                evt_path, kernel))
            if not offsets:
                continue

            bin_phys = offsets[0][0]["bin"]
            xmin, ymin, nx, ny = _stack_grid_extent(
                [(g, dx, dy) for g, dx, dy, _, _ in offsets], bin_phys)
            xmin -= max_kx * bin_phys
            ymin -= max_ky * bin_phys
            nx += 2 * max_kx
            ny += 2 * max_ky

            sum_c = np.zeros((ny, nx), dtype=np.float32)
            sum_e = np.zeros((ny, nx), dtype=np.float32)
            sum_b = np.zeros((ny, nx), dtype=np.float32)
            cheese_maps = cheese / "maps" / det / label
            n_events_total = 0

            for g, dx, dy, evt_path, kernel in offsets:
                ix = int(round((g["xmin"] - dx - xmin) / bin_phys))
                iy = int(round((g["ymin"] - dy - ymin) / bin_phys))

                with fits.open(evt_path, memmap=True) as h:
                    data = np.asarray(h["EVENTS"].data)
                if len(data):
                    xs = np.asarray(data["X"], dtype=float)
                    ys = np.asarray(data["Y"], dtype=float)
                    good = (np.isfinite(xs) & np.isfinite(ys)
                            & (np.abs(xs) < 1e7) & (np.abs(ys) < 1e7))
                    h2, _xe, _ye = np.histogram2d(
                        xs[good], ys[good], bins=(nx, ny),
                        range=((xmin, xmin + nx * bin_phys),
                               (ymin, ymin + ny * bin_phys)),
                    )
                    sum_c += h2.T.astype(np.float32)
                    n_events_total += int(good.sum())

                # Smear exp_vig/bkg along the within-frame comet path so
                # chip gaps and cheese disks match the events' geometry.
                kh, kw = kernel.shape
                py, px = (kh - 1) // 2, (kw - 1) // 2

                base = evt_path.stem.removesuffix(f"_{label}_stack")
                exp_path = cheese_maps / f"{base}_{label}_exp_vig_cheesed.fits"
                bkg_path = cheese_maps / f"{base}_{label}_bkg_cheesed.fits"
                # evselect rounds nx/ny slightly differently from a naive
                # (xmax-xmin)/bin, so clip the slab to the actual image
                # shape and to the stack-grid bounds.
                for arr_in, accum in ((exp_path, sum_e), (bkg_path, sum_b)):
                    if not arr_in.is_file():
                        continue
                    img, _ = _read_first2d(arr_in)
                    if kh > 1 or kw > 1:
                        img = _ndi_correlate(
                            np.pad(img.astype(np.float32),
                                   ((py, py), (px, px)),
                                   mode="constant", constant_values=0.0),
                            kernel, mode="constant", cval=0.0)
                        ix_eff, iy_eff = ix - px, iy - py
                    else:
                        img = img.astype(np.float32)
                        ix_eff, iy_eff = ix, iy
                    fy, fx = img.shape
                    sx0, sy0 = max(0, ix_eff), max(0, iy_eff)
                    sx1 = min(nx, ix_eff + fx); sy1 = min(ny, iy_eff + fy)
                    if sx1 <= sx0 or sy1 <= sy0:
                        continue
                    accum[sy0:sy1, sx0:sx1] += img[
                        sy0 - iy_eff : sy0 - iy_eff + (sy1 - sy0),
                        sx0 - ix_eff : sx0 - ix_eff + (sx1 - sx0),
                    ].astype(np.float32)

            corrected = _safe_rate(sum_c - sum_b, sum_e, quantile)

            # Sample header from the first frame's counts to inherit
            # celestial WCS keys (CRVAL/CDELT/CTYPE) — those are universal,
            # only CRPIX/LTV need to shift to the stack grid origin.
            first_g, _, _, first_evt, _ = offsets[0]
            first_base = first_evt.stem.removesuffix(f"_{label}_stack")
            sample_counts = (cheese.parent / "maps" / det / label
                             / f"{first_base}_{label}_counts.fits")
            with fits.open(sample_counts, memmap=True) as h:
                sample_header = h[0].header.copy()
            new_hdr = _image_header(sample_header)
            refxcrpx = (float(sample_header["CRPIX1"]) - 0.5) * bin_phys + first_g["xmin"]
            refycrpx = (float(sample_header["CRPIX2"]) - 0.5) * bin_phys + first_g["ymin"]
            new_hdr["CRPIX1"] = (refxcrpx - xmin) / bin_phys + 0.5
            new_hdr["CRPIX2"] = (refycrpx - ymin) / bin_phys + 0.5
            new_hdr["LTV1"] = 0.5 - xmin / bin_phys
            new_hdr["LTV2"] = 0.5 - ymin / bin_phys
            new_hdr["STACKREF_RA"]  = (ref_ra,  "comet reference RA [deg]")
            new_hdr["STACKREF_DEC"] = (ref_dec, "comet reference DEC [deg]")

            for arr, name in ((sum_c,    "counts"),
                              (sum_e,    "exp_vig"),
                              (sum_b,    "bkg"),
                              (corrected, "corrected")):
                out = stack / f"{det}_{label}_{name}.fits"
                fits.PrimaryHDU(data=arr, header=new_hdr.copy()).writeto(
                    out, overwrite=True)

            summary.append(
                f"{det}\t{label}\t{nx}\t{ny}\t{n_events_total}\t"
                f"{float(sum_e.sum()):.6g}\t{float(sum_b.sum()):.6g}"
            )

    _merge_tsv_keep_other(
        stack / "stack_summary.tsv",
        "det\tband\tnx\tny\tn_events\ttotal_exposure\tbkg_total",
        (0,), summary[1:], set(detectors),
    )


def stack_merge(stack_dir: str, config_path: str,
                detectors_arg: str = "") -> None:
    """Pure-Python merge of per-detector stacks: a common grid is built as
    the union of per-detector footprints (all share CRVAL/CDELT, so the
    pixel offset between stacks is just CRPIX_d - CRPIX_ref); each
    detector's counts/exp_vig/bkg is bilinearly resampled onto that grid
    and summed. exp_vig is pre-multiplied per detector by its relative
    effective area (cfg.stack_merge_weights, default {PN:1.0, M1:0.4,
    M2:0.4} matching eimagecombine) so the resulting rate map is in
    PN-equivalent count rate, uniform across PN chip gaps where MOS
    covers."""
    from math import ceil, floor
    from scipy.ndimage import map_coordinates
    cfg = load_config(config_path)
    cut_bands = _cut_bands(cfg)
    detectors = _resolve_detectors(cfg, detectors_arg)
    weights = {str(k).upper(): float(v) for k, v in
               cfg.get("stack_merge_weights",
                       {"PN": 1.0, "M1": 0.4, "M2": 0.4}).items()}
    quantile = float(cfg.get("exp_floor_quantile", 0.01))
    stack = Path(stack_dir)

    summary = ["band\ttype\tn_pixels_positive\ttotal\tmin\tmax"]
    for band in cut_bands:
        label = str(band["label"])
        per_det = []
        for d in detectors:
            cp = stack / f"{d}_{label}_counts.fits"
            ep = stack / f"{d}_{label}_exp_vig.fits"
            bp = stack / f"{d}_{label}_bkg.fits"
            if not all(p.is_file() for p in (cp, ep, bp)):
                continue
            with fits.open(cp) as h:
                hdr = h[0].header.copy()
                shape = h[0].data.shape
            per_det.append({"det": d, "hdr": hdr, "shape": shape,
                            "paths": (cp, ep, bp)})
        if not per_det:
            continue

        # Reference = first detector's grid; the offset between any other
        # detector's grid and this is CRPIX_d - CRPIX_ref (real-valued, so
        # use bilinear interp on the slab).
        ref_hdr = per_det[0]["hdr"]
        rcx, rcy = float(ref_hdr["CRPIX1"]), float(ref_hdr["CRPIX2"])
        x_lo = y_lo = float("inf"); x_hi = y_hi = float("-inf")
        for d in per_det:
            cx = float(d["hdr"]["CRPIX1"]); cy = float(d["hdr"]["CRPIX2"])
            ny_d, nx_d = d["shape"]
            sx, sy = rcx - cx, rcy - cy
            d["shift"] = (sx, sy)
            x_lo = min(x_lo, 1 + sx); x_hi = max(x_hi, nx_d + sx)
            y_lo = min(y_lo, 1 + sy); y_hi = max(y_hi, ny_d + sy)

        x0 = floor(x_lo); y0 = floor(y_lo)
        nx = int(ceil(x_hi) - x0) + 1
        ny = int(ceil(y_hi) - y0) + 1

        sum_c = np.zeros((ny, nx), dtype=np.float32)
        sum_e = np.zeros((ny, nx), dtype=np.float32)
        sum_b = np.zeros((ny, nx), dtype=np.float32)
        i_grid, j_grid = np.indices((ny, nx))

        for d in per_det:
            sx, sy = d["shift"]
            w = float(weights.get(d["det"].upper(), 1.0))
            # Common pixel (j+1, i+1) (1-indexed) maps to ref pixel
            # (j+1+x0-1, i+1+y0-1) = (j+x0, i+y0). Detector d's pixel for
            # the same celestial position is shifted back by (sx, sy):
            #   d_pix_1based = (j+x0-sx, i+y0-sy)
            # 0-indexed coords for map_coordinates:
            j_in_d = j_grid + x0 - sx - 1.0
            i_in_d = i_grid + y0 - sy - 1.0
            # Counts and bkg unweighted; exp_vig pre-multiplied by w_d so
            # the merged rate is PN-equivalent (eimagecombine recipe).
            scales = (1.0, w, 1.0)
            for path, accum, scale in zip(d["paths"],
                                          (sum_c, sum_e, sum_b), scales):
                with fits.open(path) as h:
                    arr = h[0].data.astype(np.float32)
                if scale != 1.0:
                    arr = arr * np.float32(scale)
                accum += map_coordinates(
                    arr, [i_in_d, j_in_d], order=1,
                    mode="constant", cval=0.0).astype(np.float32)

        corr = _safe_rate(sum_c - sum_b, sum_e, quantile)

        new_hdr = ref_hdr.copy()
        new_hdr["CRPIX1"] = rcx - x0 + 1
        new_hdr["CRPIX2"] = rcy - y0 + 1
        for k in ("LTV1", "LTV2"):
            if k in new_hdr:
                del new_hdr[k]

        for arr, name in ((sum_c, "counts"), (sum_e, "exp_vig"),
                          (sum_b, "bkg"), (corr, "corrected")):
            out = stack / f"merged_{label}_{name}.fits"
            fits.PrimaryHDU(data=arr, header=new_hdr.copy()).writeto(
                out, overwrite=True)
            pos = arr[arr > 0]
            summary.append(
                f"{label}\t{name}\t{int(pos.size)}\t{float(arr.sum()):.6g}\t"
                f"{float(pos.min()) if pos.size else 0:.6g}\t"
                f"{float(pos.max()) if pos.size else 0:.6g}"
            )

    (stack / "stack_merge_summary.tsv").write_text(
        "\n".join(summary) + "\n", encoding="utf-8")


def qc_stack(stack_dir: str, detectors_arg: str, out_dir: str,
             config_path: str = "") -> None:
    """Per-(det, band) stacked-frame mosaics of counts, exp_vig, bkg, and
    corrected, plus the merged products if present. No frame boxes
    (stacked products lose per-frame structure) and no comet overlay
    (comet sits at the reference by construction).

    Corrected maps are rendered with a 0-centered linear `vanimo`
    diverging colormap; a Gaussian-smoothed companion is also written for
    each detector + the merged image."""
    from scipy.ndimage import gaussian_filter
    stack = Path(stack_dir)
    out = Path(out_dir); out.mkdir(parents=True, exist_ok=True)
    detectors = normalize_detectors(detectors_arg)

    sigma_px = 10.0
    if config_path:
        try:
            sigma_px = float(load_config(config_path).get(
                "qc_stack_smooth_sigma_px", 10.0))
        except Exception:
            pass

    files = sorted(p for p in stack.rglob("*") if p.is_file()
                   if "/tmp/" not in str(p))
    _write_listing(out / "files.txt", files)
    for tsv in ("events_summary.tsv", "stack_summary.tsv",
                "stack_merge_summary.tsv"):
        src = stack / tsv
        if src.is_file():
            (out / tsv).write_bytes(src.read_bytes())

    # "merged" is treated as a virtual detector if outputs exist.
    det_tags = list(detectors)
    if any((stack / f"merged_{b}_counts.fits").is_file()
           for b in ("soft", "hard")):
        det_tags.append("merged")

    def _div_norm(arr: np.ndarray) -> tuple[Normalize, float]:
        """0-centered linear norm with extent = max(|p1|, |p99|).
        Computed on the supplied (unsmoothed) array so both the regular
        and smoothed mosaics share the same color range."""
        finite = arr[np.isfinite(arr)]
        if finite.size:
            p1 = float(np.percentile(finite, 1))
            p99 = float(np.percentile(finite, 99))
            v = max(abs(p1), abs(p99))
        else:
            v = 0.0
        if v == 0.0:
            v = 1e-6
        return Normalize(vmin=-v, vmax=v), v

    summary = ["det\tband\ttype\tn_pixels_positive\tmin\tmax"]
    for det in det_tags:
        for band in ("soft", "hard"):
            for tag in ("counts", "exp_vig", "bkg", "corrected"):
                p = stack / f"{det}_{band}_{tag}.fits"
                if not p.is_file():
                    continue
                img, hdr = _read_first2d(p)
                ny, nx = img.shape
                extent = (0.0, float(nx), 0.0, float(ny))

                if tag == "corrected":
                    norm, v = _div_norm(img)
                    title = (f"{det}  {band}  corrected [stacked]  "
                             f"{nx}x{ny}  ±{v:.3g}")
                    _save_mosaic(img, extent, [],
                                 out / f"{det}_{band}_corrected_mosaic.png",
                                 title, norm=norm, cmap="vanimo",
                                 colorbar_label="rate")
                    img_s = gaussian_filter(img.astype(np.float32),
                                            sigma=sigma_px).astype(np.float32)
                    title_s = (f"{det}  {band}  corrected smoothed "
                               f"σ={sigma_px:g}px [stacked]  "
                               f"{nx}x{ny}  ±{v:.3g}")
                    _save_mosaic(
                        img_s, extent, [],
                        out / f"{det}_{band}_corrected_smoothed_mosaic.png",
                        title_s, norm=norm, cmap="vanimo",
                        colorbar_label="rate")
                    fits.PrimaryHDU(data=img_s, header=hdr.copy()).writeto(
                        out / f"{det}_{band}_corrected_smoothed.fits",
                        overwrite=True)
                    summary.append(
                        f"{det}\t{band}\tcorrected\t-\t"
                        f"{float(img.min()):.6g}\t{float(img.max()):.6g}")
                else:
                    positive = img[img > 0]
                    vmin = (float(np.percentile(positive, 5))
                            if positive.size else 1e-6)
                    vmax = (float(np.percentile(positive, 99.5))
                            if positive.size else 1.0)
                    if not (vmax > vmin > 0):
                        vmax = max(vmax, vmin * 10.0, 1e-6)
                    title = (f"{det}  {band}  {tag} [stacked]  "
                             f"{nx}x{ny}  "
                             f"{int(positive.size)} positive px")
                    _save_mosaic(img, extent, [],
                                 out / f"{det}_{band}_{tag}_mosaic.png",
                                 title, vmin=vmin, vmax=vmax)
                    summary.append(
                        f"{det}\t{band}\t{tag}\t{positive.size}\t"
                        f"{float(positive.min()) if positive.size else 0:.6g}\t"
                        f"{float(positive.max()) if positive.size else 0:.6g}"
                    )
    (out / "qc_summary.tsv").write_text("\n".join(summary) + "\n",
                                         encoding="utf-8")


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
    """Per-frame physical-coordinate image grid derived from the supplied
    event list's X/Y bounds. Pass the **broadband** frame event list (post-
    attitude-cluster, pre-band-cut) so the grid covers the full detector
    footprint at this attitude — band-cut event lists are sparse and would
    truncate the grid to where photons happened to land.

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


def maps_counts(cut_dir: str, maps_dir: str, frames_dir: str, log_dir: str,
                config_path: str, detectors_arg: str = "") -> None:
    """Bin each cut event list into a counts image on a per-frame grid.

    The grid is derived once per (detector, frame) from the broadband
    post-attitude-cluster frame event list (frames/manifest/<DET>_frames.txt),
    so the same xmin/xmax/ymin/ymax/nx/ny are used for every cut band of
    that frame. This guarantees that exposed pixels with zero events in a
    sparse band (e.g. MOS soft) are still represented in the grid, instead
    of being truncated by the band's event extent.
    """
    cfg = load_config(config_path)
    cut_bands = _cut_bands(cfg)
    detectors = _resolve_detectors(cfg, detectors_arg)
    bin_arcsec = float(cfg.get("map_bin_arcsec", 4.0))
    pad_pix = int(cfg.get("map_pad_pix", 1))
    excluded = set(str(s) for s in cfg.get("excluded_frames", []))

    cut = Path(cut_dir); maps = Path(maps_dir); log = Path(log_dir)
    frames = Path(frames_dir)
    grid_header = ("det\tband\tbase\tximagemin\tximagemax\tyimagemin"
                   "\tyimagemax\tbin\tnx\tny")
    grid_rows: list[str] = []
    for det in detectors:
        # One grid per (det, frame) from the broadband frame event list.
        frame_paths = _read_manifest(
            frames / "manifest" / f"{det}_frames.txt")
        grid_by_base: dict[str, dict[str, float | int]] = {}
        for fp in frame_paths:
            if not fp.is_file():
                continue
            g = _frame_grid(fp, bin_arcsec, pad_pix)
            if g is not None:
                grid_by_base[fp.stem] = g

        for band in cut_bands:
            label = str(band["label"])
            events = _read_manifest(cut / "manifest" / f"{det}_{label}_cut.txt")
            band_dir = maps / det / label
            band_dir.mkdir(parents=True, exist_ok=True)
            for evt in events:
                base = _frame_base_from_cut(evt, label)
                if base in excluded:
                    print(f"  maps-counts skip excluded: {base}", flush=True)
                    continue
                grid = grid_by_base.get(base)
                if grid is None:
                    # Fallback: derive from the cut events (only happens if
                    # the frame manifest doesn't list this base).
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
    _merge_tsv_keep_other(maps / "maps_grid.tsv", grid_header,
                          (0,), grid_rows, set(detectors))
    _build_maps_manifest(maps, detectors, cut_bands)


def maps_exposure(cut_dir: str, maps_dir: str, atthk_path: str,
                  log_dir: str, config_path: str,
                  detectors_arg: str = "") -> None:
    """Build per-frame exposure maps via eexpmap. Two are produced per
    frame and band: `*_exp_vig.fits` (vignetting on, the standard exposure
    used for rate calculations) and `*_exp_unvig.fits` (vignetting off,
    used by the two-component bkg fit to capture the unvignetted
    instrumental/QPB component separately from the vignetted sky bkg).
    """
    cfg = load_config(config_path)
    cut_bands = _cut_bands(cfg)
    detectors = _resolve_detectors(cfg, detectors_arg)
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
                for tag, vig in (("exp_vig", "yes"), ("exp_unvig", "no")):
                    outp = band_dir / f"{base}_{label}_{tag}.fits"
                    if outp.is_file():
                        continue
                    _run_sas(
                        f"maps_exposure_{tag}_{det}_{label}_{base}", log, [
                            "eexpmap",
                            f"imageset={counts}",
                            f"attitudeset={atthk_path}",
                            f"eventset={evt}",
                            f"expimageset={outp}",
                            f"pimin={pi_lo}", f"pimax={pi_hi}",
                            f"withvignetting={vig}",
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


def _solve_bkg2(c: np.ndarray, e_unvig: np.ndarray, e_vig: np.ndarray,
                mode: str) -> tuple[float, float, float, str]:
    """Two-component bkg fit: c ≈ q*e_unvig + s*e_vig.

    `q` is the unvignetted (instrumental / QPB-like) per-second-per-pixel
    coefficient; `s` is the vignetted (sky / SWCX-like) coefficient. Modes:
        "none"   – ordinary least squares (q, s can be any sign)
        "nonneg" – constrained to q, s ≥ 0 via NNLS (default-ish for
                   physical bkg components)
    Returns (q, s, rmse, fit_mode_label).
    """
    if c.size < 2:
        return 0.0, 0.0, 0.0, "empty"
    A = np.column_stack([e_unvig, e_vig])
    if mode == "nonneg":
        from scipy.optimize import nnls
        sol, _ = nnls(A, c)
        q, s = float(sol[0]), float(sol[1])
        label = "nonneg"
    else:
        sol, *_ = np.linalg.lstsq(A, c, rcond=None)
        q, s = float(sol[0]), float(sol[1])
        label = "lstsq"
    pred = q * e_unvig + s * e_vig
    rmse = float(np.sqrt(np.mean((c - pred) ** 2)))
    return q, s, rmse, label


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


def maps_background(maps_dir: str, config_path: str,
                    detectors_arg: str = "") -> None:
    """Per-frame two-component bkg fit: counts ≈ q·E_unvig + s·E_vig.

    `q` represents the unvignetted instrumental/QPB component (constant
    per-pixel-per-time, since particles aren't focused by the mirror); `s`
    represents the vignetted sky / SWCX component (scales with effective-
    area×time = exp_vig). Sample pixels are taken outside an off-source
    radius (rough comet/source exclusion), with low-exposure pixels
    excluded, and 5σ Poisson outliers iteratively dropped.
    """
    cfg = load_config(config_path)
    cut_bands = _cut_bands(cfg)
    detectors = _resolve_detectors(cfg, detectors_arg)
    mode = str(cfg.get("map_bkg_constrain", "nonneg")).strip().lower()
    if mode not in {"none", "nonneg"}:
        die(f"map_bkg_constrain must be none/nonneg (got {mode!r})")
    outside_frac = float(cfg.get("map_bkg_outside_radius_fraction", 0.5))
    scale_cfg = cfg.get("map_bkg_scale", 1.0)
    sigma_clip = float(cfg.get("map_bkg_sigma_clip", 5.0))
    low_exp_q = float(cfg.get("map_bkg_low_exp_quantile", 0.05))

    maps = Path(maps_dir)
    summary_header = ("det\tband\tbase\tfit_q\tfit_s\trmse\tn_sample\t"
                      "n_clipped\tcenter_x\tcenter_y\tradius_pix\tfit_mode")
    summary_rows: list[str] = []

    for det in detectors:
        for band in cut_bands:
            label = str(band["label"])
            band_dir = maps / det / label
            if not band_dir.is_dir():
                continue
            scale = _bkg_scale_for(scale_cfg, det, label)
            for counts_path in sorted(band_dir.glob(f"*_{label}_counts.fits")):
                base = counts_path.stem.removesuffix(f"_{label}_counts")
                vig_path   = band_dir / f"{base}_{label}_exp_vig.fits"
                unvig_path = band_dir / f"{base}_{label}_exp_unvig.fits"
                if not (vig_path.is_file() and unvig_path.is_file()):
                    continue
                bkg_path = band_dir / f"{base}_{label}_bkg.fits"
                row = _fit_one_bkg(counts_path, vig_path, unvig_path,
                                   bkg_path, mode, outside_frac, scale,
                                   sigma_clip, low_exp_q)
                summary_rows.append(f"{det}\t{label}\t{base}\t{row}")
    _merge_tsv_keep_other(maps / "bkg_summary.tsv", summary_header,
                          (0,), summary_rows, set(detectors))
    _build_maps_manifest(maps, detectors, cut_bands)


def _bkg_scale_for(cfg_val: Any, det: str, label: str) -> float:
    """Resolve cfg.map_bkg_scale to a scalar for one (detector, band).
    Accepts a number (applies everywhere) or a dict {DET: scalar} or
    {DET: {BAND: scalar}}. Falls back to 1.0 if no match."""
    if isinstance(cfg_val, (int, float)):
        return float(cfg_val)
    if isinstance(cfg_val, dict):
        per_det = (cfg_val.get(det)
                   or cfg_val.get(det.upper())
                   or cfg_val.get(det.lower()))
        if isinstance(per_det, (int, float)):
            return float(per_det)
        if isinstance(per_det, dict):
            v = (per_det.get(label)
                 or per_det.get(label.upper())
                 or per_det.get(label.lower()))
            if isinstance(v, (int, float)):
                return float(v)
    return 1.0


def _fit_one_bkg(counts_path: Path, vig_path: Path, unvig_path: Path,
                 bkg_path: Path, mode: str, outside_frac: float,
                 scale: float = 1.0, sigma_clip: float = 5.0,
                 low_exp_q: float = 0.05) -> str:
    """Two-component bkg fit and bkg image writer.

    Reads counts, vignetted exposure (exp_vig) and unvignetted exposure
    (exp_unvig). Solves c ≈ q·E_unvig + s·E_vig on a sample of valid
    off-source pixels (footprint ∩ outside-radius ∩ finite ∩ exp_vig in
    the upper (1-low_exp_q) of positive exposures). Iteratively drops
    >sigma_clip·σ Poisson outliers from the fit. Then writes the bkg
    image as scale·(q·E_unvig + s·E_vig), zeroed outside the footprint.
    No post-fit taper — the spatial profile comes from the model itself.
    """
    counts_img, c_header = _read_first2d(counts_path)
    exp_vig, _ = _read_first2d(vig_path)
    exp_unvig, _ = _read_first2d(unvig_path)
    if not (counts_img.shape == exp_vig.shape == exp_unvig.shape):
        die(f"shape mismatch among counts/vig/unvig: {counts_path}")

    footprint = np.isfinite(exp_vig) & (exp_vig > 0)
    ny, nx = counts_img.shape

    # Source/comet exclusion via outside-radius circle (radius derived
    # from the footprint extent).
    if footprint.any():
        ys, xs = np.where(footprint)
        cy = float(ys.mean()); cx = float(xs.mean())
        half = 0.5 * max(ys.max() - ys.min(), xs.max() - xs.min())
        radius = outside_frac * half
        yy, xx = np.ogrid[:ny, :nx]
        outside_circle = (yy - cy) ** 2 + (xx - cx) ** 2 >= radius * radius
        sample = (footprint & outside_circle
                  & np.isfinite(counts_img) & np.isfinite(exp_unvig))
        if not sample.any():
            sample = (footprint & np.isfinite(counts_img)
                      & np.isfinite(exp_unvig))
    else:
        cy = (ny - 1) / 2.0; cx = (nx - 1) / 2.0
        radius = 0.0
        sample = np.zeros_like(footprint, dtype=bool)

    # Drop the bottom low_exp_q fraction of positive exposures from the
    # fit — those vignetting-tail pixels are individually unreliable.
    if low_exp_q > 0 and footprint.any():
        thresh = float(np.quantile(exp_vig[footprint], low_exp_q))
        sample &= exp_vig >= thresh

    n_sample = int(sample.sum())
    c_all = counts_img[sample].astype(float) if n_sample else np.zeros(0)
    eu_all = exp_unvig[sample].astype(float) if n_sample else np.zeros(0)
    ev_all = exp_vig[sample].astype(float)   if n_sample else np.zeros(0)

    # Initial fit; iteratively drop >sigma_clip·σ Poisson outliers.
    keep = np.ones(n_sample, dtype=bool)
    q, s, rmse, fit_mode = _solve_bkg2(c_all, eu_all, ev_all, mode)
    if sigma_clip > 0 and n_sample > 1:
        for _ in range(3):
            predicted = q * eu_all + s * ev_all
            sigma = np.sqrt(np.maximum(predicted, 1.0))
            new_keep = (c_all - predicted) <= sigma_clip * sigma
            if np.array_equal(new_keep, keep):
                break
            keep = new_keep
            if int(keep.sum()) < 2:
                break
            q, s, rmse, fit_mode = _solve_bkg2(
                c_all[keep], eu_all[keep], ev_all[keep], mode)
    n_used = int(keep.sum()) if n_sample else 0
    n_clipped = n_sample - n_used

    # Build bkg image: q·E_unvig + s·E_vig, zeroed outside footprint.
    bg = (scale * (q * exp_unvig + s * exp_vig)).astype(np.float32)
    bg = np.where(footprint, bg, 0.0).astype(np.float32)

    out_hdu = fits.PrimaryHDU(data=bg, header=_image_header(c_header))
    out_hdu.header["BGMODEL"]  = ("q*Eu+s*Ev", "two-component bkg model")
    out_hdu.header["BGMODE"]   = (fit_mode, "fit constraint mode")
    out_hdu.header["BGQ"]      = (q, "unvig (instrumental) coefficient")
    out_hdu.header["BGS"]      = (s, "vig (sky) coefficient")
    out_hdu.header["BGSCALE"]  = (scale, "post-fit multiplicative scale")
    out_hdu.header["BGSIGCL"]  = (sigma_clip, "sample sigma-clip threshold")
    out_hdu.header["BGNCLIP"]  = (n_clipped, "n pixels clipped from fit")
    out_hdu.header["BGLOWQ"]   = (low_exp_q, "low-exposure mask quantile")
    out_hdu.header["BGRMSE"]   = (rmse, "sample RMS residual")
    out_hdu.header["BG_CX"]    = (cx, "exclusion center x [pix]")
    out_hdu.header["BG_CY"]    = (cy, "exclusion center y [pix]")
    out_hdu.header["BG_RAD"]   = (radius, "exclusion radius [pix]")
    out_hdu.header["BG_NSAM"]  = (n_used, "n pixels used in final fit")
    out_hdu.writeto(bkg_path, overwrite=True)
    return (f"{q:.6g}\t{s:.6g}\t{rmse:.6g}\t{n_used}\t{n_clipped}"
            f"\t{cx:.2f}\t{cy:.2f}\t{radius:.2f}\t{fit_mode}")


def maps_corrected(maps_dir: str, config_path: str,
                   detectors_arg: str = "") -> None:
    """(counts - bkg) / exp_vig per frame; not clipped."""
    cfg = load_config(config_path)
    cut_bands = _cut_bands(cfg)
    detectors = _resolve_detectors(cfg, detectors_arg)
    quantile = float(cfg.get("exp_floor_quantile", 0.01))
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
                _make_corrected(counts, exp, bkg, out, quantile)
    _build_maps_manifest(maps, detectors, cut_bands)


def _safe_rate(num: np.ndarray, exp: np.ndarray,
               quantile: float = 0.01) -> np.ndarray:
    """num / exp with the exposure map floored at its `quantile`-th
    percentile of positive values (in-memory; the exposure FITS itself is
    untouched). Pixels where exp == 0 stay at 0 (truly outside FOV); the
    floor only affects the deep-vignetting tail where (num/exp) would
    otherwise blow up to enormous magnitudes."""
    pos = exp[exp > 0]
    if pos.size == 0:
        return np.zeros_like(num, dtype=np.float32)
    floor = float(np.quantile(pos, quantile))
    exp_safe = np.where(exp > 0, np.maximum(exp, floor), 0.0)
    with np.errstate(divide="ignore", invalid="ignore"):
        rate = np.where(exp_safe > 0, num / exp_safe, 0.0)
    return rate.astype(np.float32)


def _make_corrected(counts_path: Path, exp_path: Path, bkg_path: Path,
                    out_path: Path, quantile: float = 0.01) -> None:
    c, c_header = _read_first2d(counts_path)
    e, _ = _read_first2d(exp_path)
    b, _ = _read_first2d(bkg_path)
    rate = _safe_rate(c - b, e, quantile)
    fits.PrimaryHDU(data=rate, header=_image_header(c_header)).writeto(
        out_path, overwrite=True)


def _build_maps_manifest(maps_dir: Path, _detectors_unused: list[str],
                         cut_bands: list[dict[str, Any]]) -> None:
    """Scan every <DET>/<BAND>/ subdir on disk so partial reruns (e.g.
    --detectors M1,M2) leave PN rows intact in the manifest."""
    rows = ["det\tband\tbase\tcounts\texp_vig\tbkg\tcorrected"]
    band_labels = [str(b["label"]) for b in cut_bands]
    for det_dir in sorted(p for p in maps_dir.iterdir() if p.is_dir()):
        det = det_dir.name
        if det not in DETECTOR_ALIASES:
            continue
        for label in band_labels:
            band_dir = det_dir / label
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
            out_dir: str, config_path: str = "") -> None:
    """Per-(detector, band) mosaics of the nobkg-corrected (counts/exp_vig)
    and the corrected ((counts - bkg)/exp_vig) rate maps. Rates are co-added
    across frames as sum(num)/sum(exp) on a viewing grid built from the
    union of the per-frame physical extents.
    """
    maps = Path(maps_dir); frames = Path(frames_dir)
    out = Path(out_dir); out.mkdir(parents=True, exist_ok=True)
    detectors = normalize_detectors(detectors_arg)
    quantile = 0.01
    excluded: set[str] = set()
    if config_path:
        try:
            cfg = load_config(config_path)
            quantile = float(cfg.get("exp_floor_quantile", 0.01))
            excluded = set(str(s) for s in cfg.get("excluded_frames", []))
        except Exception:
            pass

    files = sorted(p for p in maps.rglob("*") if p.is_file())
    _write_listing(out / "files.txt", files)

    for tsv in ("maps_grid.tsv", "maps_manifest.tsv", "bkg_summary.tsv"):
        src = maps / tsv
        if src.is_file():
            (out / tsv).write_bytes(src.read_bytes())

    grids = _read_maps_grid(maps)
    rows = _read_maps_manifest_rows(maps)

    summary = ["det\tband\tnframes\tnx\tny\tcounts\texp_vig\tbkg"
               "\tnobkg_corrected\tcorrected"]
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
            if excluded:
                before = len(band_rows)
                band_rows = [r for r in band_rows if r["base"] not in excluded]
                dropped = before - len(band_rows)
                if dropped:
                    print(f"  qc-maps: {det} {band} excluding {dropped} "
                          f"frame(s) per cfg.excluded_frames", flush=True)
            if not band_rows:
                continue
            keys = [(det, band, r["base"]) for r in band_rows]
            keys = [k for k in keys if k in grids]
            if not keys:
                continue
            bin_phys = grids[keys[0]]["bin"]
            # Read actual per-frame image shapes (evselect rounds nx/ny
            # slightly differently from a naive (xmax-xmin)/bin), and use
            # those to size the viewing grid so np.sum can broadcast.
            frame_imgs: list[tuple[dict[str, str], np.ndarray, np.ndarray, np.ndarray]] = []
            xmin = min(grids[k]["xmin"] for k in keys)
            ymin = min(grids[k]["ymin"] for k in keys)
            xmax = ymax = float("-inf")
            for r in band_rows:
                key = (det, band, r["base"])
                if key not in grids:
                    continue
                g = grids[key]
                cimg, _ = _read_first2d(Path(r["counts"]))
                eimg, _ = _read_first2d(Path(r["exp_vig"]))
                bimg, _ = _read_first2d(Path(r["bkg"]))
                fy, fx = cimg.shape
                xmax = max(xmax, g["xmin"] + fx * bin_phys)
                ymax = max(ymax, g["ymin"] + fy * bin_phys)
                frame_imgs.append((r, cimg, eimg, bimg))
            nx = int(round((xmax - xmin) / bin_phys))
            ny = int(round((ymax - ymin) / bin_phys))
            sum_c = np.zeros((ny, nx), dtype=float)
            sum_e = np.zeros((ny, nx), dtype=float)
            sum_b = np.zeros((ny, nx), dtype=float)
            for r, cimg, eimg, bimg in frame_imgs:
                g = grids[(det, band, r["base"])]
                ix = int(round((g["xmin"] - xmin) / bin_phys))
                iy = int(round((g["ymin"] - ymin) / bin_phys))
                fy, fx = cimg.shape
                sum_c[iy:iy+fy, ix:ix+fx] += cimg
                sum_e[iy:iy+fy, ix:ix+fx] += eimg
                sum_b[iy:iy+fy, ix:ix+fx] += bimg
            nobkg = _safe_rate(sum_c,           sum_e, quantile)
            corr  = _safe_rate(sum_c - sum_b,   sum_e, quantile)
            extent = (xmin, xmin + nx * bin_phys,
                      ymin, ymin + ny * bin_phys)
            boxes = _boxes_from_events(box_events.get(det, []))

            for tag, img in (("counts", sum_c),
                             ("exp_vig", sum_e),
                             ("bkg", sum_b),
                             ("nobkg-corrected", nobkg),
                             ("corrected", corr)):
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

            # Diverging-colormap render of the (full bkg-subtracted)
            # corrected map: 0-centered linear vanimo, ±max(|p1|,|p99|).
            finite = corr[np.isfinite(corr)]
            if finite.size:
                p1 = float(np.percentile(finite, 1))
                p99 = float(np.percentile(finite, 99))
                v = max(abs(p1), abs(p99), 1e-6)
            else:
                v = 1e-6
            norm = Normalize(vmin=-v, vmax=v)
            title_v = (f"{det}  {band}  corrected (vanimo)  "
                       f"frames={len(band_rows)}  ±{v:.3g}")
            _save_mosaic(corr, extent, boxes,
                         out / f"{det}_{band}_corrected_vanimo_mosaic.png",
                         title_v, norm=norm, cmap="vanimo",
                         colorbar_label="rate")
            summary.append(
                f"{det}\t{band}\t{len(band_rows)}\t{nx}\t{ny}"
                f"\t{det}_{band}_counts_mosaic.png"
                f"\t{det}_{band}_exp_vig_mosaic.png"
                f"\t{det}_{band}_bkg_mosaic.png"
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
                 *, vmin: float | None = None, vmax: float | None = None,
                 circles: list[tuple[float, float, float]] | None = None,
                 cmap: str = "gray",
                 norm: "matplotlib.colors.Normalize | None" = None,
                 colorbar_label: str = "counts/bin") -> None:
    ny, nx = img.shape
    with plt.style.context("dark_background"):
        fig, ax = plt.subplots(figsize=(9, 9 * ny / nx + 0.6), dpi=120)
        if norm is None:
            if vmin is None:
                vmin = max(img.max() * 1e-4, 0.5)
            if vmax is None:
                vmax = max(img.max(), 1.0)
            norm = LogNorm(vmin=vmin, vmax=vmax)
        im = ax.imshow(img, origin="lower", cmap=cmap, norm=norm,
                       extent=extent, aspect="equal", interpolation="nearest")
        fig.colorbar(im, ax=ax, label=colorbar_label)
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
        if circles:
            for cx, cy, cr in circles:
                ax.add_patch(Circle(
                    (cx, cy), cr,
                    facecolor=(1.0, 0.33, 0.33, 0.35),
                    edgecolor="black", linewidth=0.4, zorder=7,
                ))
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
    elif cmd == "cut-run" and len(args) in {5, 6}:
        cut_run(args[0], args[1], args[2], args[3], args[4],
                args[5] if len(args) == 6 else "")
    elif cmd == "qc-cut" and len(args) == 4:
        qc_cut(*args)
    elif cmd == "cheese-detect" and len(args) in {5, 6}:
        cheese_detect(args[0], args[1], args[2], args[3], args[4],
                      args[5] if len(args) == 6 else "")
    elif cmd == "cheese-mask" and len(args) in {4, 5}:
        cheese_mask(args[0], args[1], args[2], args[3],
                    args[4] if len(args) == 5 else "")
    elif cmd == "qc-cheese" and len(args) in {4, 5}:
        qc_cheese(args[0], args[1], args[2], args[3],
                  args[4] if len(args) == 5 else "")
    elif cmd == "stack-events" and len(args) in {5, 6}:
        stack_events(args[0], args[1], args[2], args[3], args[4],
                     args[5] if len(args) == 6 else "")
    elif cmd == "stack-coadd" and len(args) in {5, 6}:
        stack_coadd(args[0], args[1], args[2], args[3], args[4],
                    args[5] if len(args) == 6 else "")
    elif cmd == "stack-merge" and len(args) in {2, 3}:
        stack_merge(args[0], args[1],
                    args[2] if len(args) == 3 else "")
    elif cmd == "qc-stack" and len(args) in {3, 4}:
        qc_stack(args[0], args[1], args[2],
                 args[3] if len(args) == 4 else "")
    elif cmd == "build-track" and len(args) == 9:
        build_track(*args)
    elif cmd == "qc-track" and len(args) == 5:
        qc_track(*args)
    elif cmd == "maps-counts" and len(args) in {5, 6}:
        maps_counts(args[0], args[1], args[2], args[3], args[4],
                    args[5] if len(args) == 6 else "")
    elif cmd == "maps-exposure" and len(args) in {5, 6}:
        maps_exposure(args[0], args[1], args[2], args[3], args[4],
                      args[5] if len(args) == 6 else "")
    elif cmd == "maps-background" and len(args) in {2, 3}:
        maps_background(args[0], args[1],
                        args[2] if len(args) == 3 else "")
    elif cmd == "maps-corrected" and len(args) in {2, 3}:
        maps_corrected(args[0], args[1],
                       args[2] if len(args) == 3 else "")
    elif cmd == "qc-maps" and len(args) in {4, 5}:
        qc_maps(args[0], args[1], args[2], args[3],
                args[4] if len(args) == 5 else "")
    else:
        die(__doc__)
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
