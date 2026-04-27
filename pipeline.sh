#!/usr/bin/env bash
# XMM-Newton diffuse moving-target pipeline (v10).
#
# Stages:
#   init   cifbuild + odfingest -> ccf.cif, *SUM.SAS, sas_setup.env
#   repro  epproc/emproc per detector + atthkgen -> event lists, atthk.dat,
#          per-detector raw manifests, attitude.env
#   frames split multi-pointing repro events into per-pointing frames using
#          the attitude-history clustering driven by atthk.dat
#   clean  per-frame flare GTI from a high-energy lightcurve, then per-band
#          PI/PATTERN/FLAG event filters (PPS standard) + attcalc to sky
#   cut    PI segmentation into soft (200-1000 eV) and hard (1001-4000 eV)
#          plus a sigma clip that drops events in 5+sigma mosaic pixels
#   track  query JPL Horizons (or load track_input) for the comet ephemeris
#          over the observation window; products feed the later stack stage
#   maps   per-frame counts/exposure_vig/bkg/corrected images on a per-frame
#          physical grid (4 arcsec/pix default); separable substages
#          maps-counts, maps-exposure, maps-background, maps-corrected
#
# Both stages adopt existing on-disk products: if the expected output files
# exist they are accepted, only book-keeping (PATH rewrite, manifests, env
# files) is regenerated, and SAS is not re-run unless --force is given.
#
# QC stages mirror the production stages (qc-init, qc-repro). qc reruns all.
set -Eeuo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONFIG="$SCRIPT_DIR/config.json"
PYTHON="${PYTHON:-python3}"
STAGE=""
FORCE=0
DETECTORS_OVERRIDE=""

usage() {
  cat <<'EOF'
Usage: ./pipeline.sh STAGE [--config FILE] [--force]

Stages:
  init        Build SAS calibration index and ODF summary.
  repro       Run epproc/emproc per detector and atthkgen.
  frames      Split multi-pointing repro events by attitude.
  clean       Per-frame flare GTI + PPS-standard band filters + attcalc.
  cut         PI segmentation (soft/hard) + sigma clip on mosaic outliers.
  track       Query Horizons (or load track_input) for the comet ephemeris.
  maps             Run all four maps substages.
  maps-counts      evselect per-frame counts on per-frame physical grid.
  maps-exposure    eexpmap per-frame vignetted exposure maps.
  maps-background  Per-frame a + b*E fit on off-source pixels.
  maps-corrected   Per-frame (counts - bkg) / exp_vig.
  all         Run init, repro, frames, clean, cut, track, maps.
  qc-init     QC products for init.
  qc-repro    QC products for repro (listing, status, loE/hiE mosaics).
  qc-frames   QC products for frames (listing, status, mosaics, timeline).
  qc-clean    QC products for clean (loE/hiE mosaics, flare lightcurves).
  qc-cut      QC products for cut (per-(det,band) mosaics with frame boxes).
  qc-track    Cut mosaics with comet trajectory overlay.
  qc-maps     Per-(det, band) nobkg-corrected and corrected mosaics.
  qc          Run all qc-* stages.

Options:
  --config FILE        JSON config (default: ./config.json).
  --force              Rebuild stage products even if outputs already exist.
  --detectors LIST     Restrict the run to a comma-separated subset of
                       detectors (e.g. M1,M2). Existing products for
                       other detectors are preserved.
EOF
}

while (($#)); do
  case "$1" in
    init|repro|frames|clean|cut|track|maps|maps-counts|maps-exposure|maps-background|maps-corrected|all|qc-init|qc-repro|qc-frames|qc-clean|qc-cut|qc-track|qc-maps|qc) STAGE="$1"; shift ;;
    --config|-c) CONFIG="$2"; shift 2 ;;
    --force|-f)  FORCE=1; shift ;;
    --detectors|-d) DETECTORS_OVERRIDE="$2"; shift 2 ;;
    -h|--help)   usage; exit 0 ;;
    *) echo "Unknown argument: $1" >&2; usage >&2; exit 2 ;;
  esac
done
[[ -n "$STAGE" ]] || { usage >&2; exit 2; }

[[ -f "$CONFIG" ]]              || { echo "Missing config: $CONFIG" >&2; exit 1; }
[[ -f "$SCRIPT_DIR/tools.py" ]] || { echo "Missing tools.py" >&2; exit 1; }
command -v "$PYTHON" >/dev/null  || { echo "Missing python: $PYTHON" >&2; exit 1; }

eval "$("$PYTHON" "$SCRIPT_DIR/tools.py" shell "$CONFIG")"
if [[ -n "$DETECTORS_OVERRIDE" ]]; then
  DETECTORS="$(echo "$DETECTORS_OVERRIDE" | tr ',' ' ')"
fi

INITDIR="$WORKDIR/init"
REPRODIR="$WORKDIR/repro"
FRAMESDIR="$WORKDIR/frames"
CLEANDIR="$WORKDIR/clean"
CUTDIR="$WORKDIR/cut"
TRACKDIR="$WORKDIR/track"
MAPSDIR="$WORKDIR/maps"
LOGDIR="$WORKDIR/logs"
QCDIR="$WORKDIR/qc"
SAS_ENV="$WORKDIR/sas_setup.env"
ATT_ENV="$WORKDIR/attitude.env"
TRACK_ENV="$TRACKDIR/track.env"
TRACK_FITS="$TRACKDIR/track.fits"
TRACK_JSON="$TRACKDIR/track.json"
mkdir -p "$INITDIR" "$REPRODIR/manifest" "$FRAMESDIR/manifest" \
         "$CLEANDIR/manifest" "$CUTDIR/manifest" "$TRACKDIR" \
         "$MAPSDIR" "$LOGDIR"

trap 'echo "[ERR] line $LINENO: $BASH_COMMAND" >&2' ERR

# ---------- helpers ----------

log_run() {
  # log_run TAG CMD ARGS...   stream stdout/stderr to log file and tee.
  local tag="$1"; shift
  local logfile="$LOGDIR/${tag}.log"
  echo "[$(date -u +%FT%TZ)] $tag: $*" | tee -a "$LOGDIR/pipeline.log"
  ( "$@" ) 2>&1 | tee "$logfile"
  return "${PIPESTATUS[0]}"
}

source_sas() {
  # Source SAS once. After this cifbuild/epproc/etc are on PATH.
  if command -v sasversion >/dev/null && [[ -n "${SAS_DIR:-}" ]]; then return 0; fi
  [[ -n "$SAS_SETUP_SCRIPT" && -f "$SAS_SETUP_SCRIPT" ]] \
    || { echo "sas_setup_script not found: $SAS_SETUP_SCRIPT" >&2; exit 1; }
  set +u; # shellcheck disable=SC1090
  source "$SAS_SETUP_SCRIPT"
  set -u
}

first_sumsas() { find "$1" -maxdepth 1 -name '*SUM.SAS' -type f | sort | head -n 1; }

write_sas_setup_env() {
  local sumsas="$1"
  {
    echo "# reduction_v10/pipeline.sh"
    echo "export SAS_CCF=$(printf '%q' "$INITDIR/ccf.cif")"
    echo "export SAS_ODF=$(printf '%q' "$sumsas")"
    [[ -n "$SAS_CCFPATH_CONFIG"  ]] && echo "export SAS_CCFPATH=$(printf '%q' "$SAS_CCFPATH_CONFIG")"
    [[ -n "$SAS_VERBOSITY_CONFIG" ]] && echo "export SAS_VERBOSITY=$(printf '%q' "$SAS_VERBOSITY_CONFIG")"
  } > "$SAS_ENV"
}

write_attitude_env() {
  local ahf
  ahf="$(find "$ODFDIR" -maxdepth 1 -type f -name '*ATS.FIT*' | sort | head -n 1)"
  [[ -n "$ahf" ]] || { echo "No *ATS.FIT* attitude file in $ODFDIR" >&2; exit 1; }
  {
    echo "ATTHK_FILE=$(printf '%q' "$REPRODIR/atthk.dat")"
    echo "AHF_FILE=$(printf '%q'  "$ahf")"
  } > "$ATT_ENV"
}

link_odf_into_init() {
  # Symlink every ODF constituent (excluding the gpg/tarball wrappers) into
  # init/. cifbuild/odfingest expect their inputs to live next to SUM.SAS.
  local f
  for f in "$ODFDIR"/*; do
    [[ -f "$f" ]] || continue
    case "${f##*/}" in *.gz.gpg|*.tar.gz|*.tar) continue ;; esac
    ln -sf "$f" "$INITDIR/$(basename "$f")"
  done
}

# ---------- stage: init ----------

init_have_products() {
  local sumsas
  sumsas="$(first_sumsas "$INITDIR")"
  [[ -s "$INITDIR/ccf.cif" && -n "$sumsas" && -s "$sumsas" ]]
}

stage_init() {
  if [[ "$FORCE" == "1" ]]; then
    rm -f "$INITDIR/ccf.cif" "$INITDIR"/*SUM.SAS "$INITDIR"/*SUM.ASC "$SAS_ENV"
  fi

  link_odf_into_init

  if ! init_have_products; then
    [[ -d "$ODFDIR" ]] || { echo "ODFDIR missing: $ODFDIR" >&2; exit 1; }
    source_sas
    export SAS_ODF="$ODFDIR"
    export SAS_CCF="$INITDIR/ccf.cif"
    [[ -n "$SAS_CCFPATH_CONFIG"  ]] && export SAS_CCFPATH="$SAS_CCFPATH_CONFIG"
    [[ -n "$SAS_VERBOSITY_CONFIG" ]] && export SAS_VERBOSITY="$SAS_VERBOSITY_CONFIG"
    pushd "$INITDIR" >/dev/null
    log_run init_cifbuild  cifbuild
    export SAS_CCF="$INITDIR/ccf.cif"
    log_run init_odfingest odfingest odfdir="$ODFDIR" outdir="$INITDIR"
    popd >/dev/null
  fi

  local sumsas
  sumsas="$(first_sumsas "$INITDIR")"
  [[ -n "$sumsas" ]] || { echo "no *SUM.SAS in $INITDIR" >&2; exit 1; }

  # Always rewrite PATH (cheap; fixes paths after directory moves).
  "$PYTHON" "$SCRIPT_DIR/tools.py" rewrite-sum-path "$sumsas" "$INITDIR"
  write_sas_setup_env "$sumsas"
  echo "init ok: $(basename "$sumsas")"
}

# ---------- stage: repro ----------

manifest_paths_exist() {
  local manifest="$1" line
  [[ -s "$manifest" ]] || return 1
  while IFS= read -r line; do
    [[ -z "$line" ]] && continue
    [[ -s "$line" ]] || return 1
  done < "$manifest"
}

repro_have_products() {
  local det manifest
  for det in $DETECTORS; do
    manifest="$REPRODIR/manifest/${det}_raw.txt"
    "$PYTHON" "$SCRIPT_DIR/tools.py" repro-manifest "$REPRODIR" "$det" "$manifest"
    manifest_paths_exist "$manifest" || return 1
  done
  [[ -s "$REPRODIR/atthk.dat" ]]
}

stage_repro() {
  init_have_products || stage_init

  if [[ "$FORCE" == "1" ]]; then
    rm -f "$REPRODIR/atthk.dat" "$REPRODIR/manifest/"*_raw.txt "$ATT_ENV"
    rm -f "$REPRODIR"/*ImagingEvts.ds "$REPRODIR"/*Badpixels.ds "$REPRODIR"/*AttHk.ds
  fi

  if [[ "$FORCE" != "1" ]] && repro_have_products; then
    write_attitude_env
    echo "repro ok: adopted existing products"
    return 0
  fi

  source_sas
  # shellcheck disable=SC1090
  source "$SAS_ENV"

  pushd "$REPRODIR" >/dev/null
  local entry det task
  for entry in $PROC_TASKS; do
    det="${entry%%:*}"; task="${entry##*:}"
    "$PYTHON" "$SCRIPT_DIR/tools.py" repro-manifest "$REPRODIR" "$det" \
      "$REPRODIR/manifest/${det}_raw.txt"
    if [[ "$FORCE" == "1" ]] || ! manifest_paths_exist "$REPRODIR/manifest/${det}_raw.txt"; then
      log_run "repro_${task}" "$task"
      "$PYTHON" "$SCRIPT_DIR/tools.py" repro-manifest "$REPRODIR" "$det" \
        "$REPRODIR/manifest/${det}_raw.txt"
    fi
    manifest_paths_exist "$REPRODIR/manifest/${det}_raw.txt" \
      || { echo "no $det events after $task" >&2; exit 1; }
  done

  if [[ "$FORCE" == "1" || ! -s "$REPRODIR/atthk.dat" ]]; then
    log_run repro_atthkgen atthkgen atthkset="$REPRODIR/atthk.dat"
  fi
  popd >/dev/null

  write_attitude_env
  echo "repro ok"
}

# ---------- stage: frames ----------

frames_have_products() {
  local det
  [[ -s "$FRAMESDIR/frames.tsv" ]] || return 1
  for det in $DETECTORS; do
    manifest_paths_exist "$FRAMESDIR/manifest/${det}_frames.txt" || return 1
  done
}

stage_frames() {
  repro_have_products || stage_repro

  if [[ "$FORCE" == "1" ]]; then
    rm -rf "$FRAMESDIR/events" "$FRAMESDIR/manifest" "$FRAMESDIR/frames.tsv"
    mkdir -p "$FRAMESDIR/manifest"
  elif frames_have_products; then
    echo "frames ok: adopted existing products"
    return 0
  fi

  [[ -s "$ATT_ENV" ]] || { echo "missing $ATT_ENV (run repro first)" >&2; exit 1; }
  # shellcheck disable=SC1090
  source "$ATT_ENV"
  [[ -s "$ATTHK_FILE" ]] || { echo "atthk file missing: $ATTHK_FILE" >&2; exit 1; }

  log_run frames_split \
    "$PYTHON" "$SCRIPT_DIR/tools.py" frames-split \
      "$REPRODIR" "$FRAMESDIR" "$ATTHK_FILE" \
      "$DETECTORS" "$FRAMES_THRESHOLD_AMIN" "$FRAMES_MIN_DURATION_S"
  echo "frames ok"
}

# ---------- stage: clean ----------

clean_have_products() {
  [[ -s "$CLEANDIR/clean_band_filters.tsv" ]] || return 1
  [[ -s "$CLEANDIR/flare_gti_summary.tsv" ]] || return 1
  local det label
  for det in $DETECTORS; do
    for label in $CLEAN_BAND_LABELS; do
      manifest_paths_exist "$CLEANDIR/manifest/${det}_${label}_clean.txt" || return 1
    done
  done
}

stage_clean() {
  frames_have_products || stage_frames

  if [[ "$FORCE" == "1" ]]; then
    rm -rf "$CLEANDIR"
    mkdir -p "$CLEANDIR/manifest"
  elif clean_have_products; then
    echo "clean ok: adopted existing products"
    return 0
  fi

  source_sas
  # shellcheck disable=SC1090
  source "$SAS_ENV"
  # shellcheck disable=SC1090
  source "$ATT_ENV"
  [[ -s "$ATTHK_FILE" ]] || { echo "atthk missing: $ATTHK_FILE" >&2; exit 1; }

  "$PYTHON" "$SCRIPT_DIR/tools.py" clean-band-table "$CONFIG" \
    > "$CLEANDIR/clean_band_filters.tsv"
  printf 'inst\tframe\trate_cut\ttimebin\tlightcurve\tgti\n' \
    > "$CLEANDIR/flare_gti_summary.tsv"

  local det rate_var rate flare_var flare frame base lc gti expr_var expr clean rows
  for det in $DETECTORS; do
    rate_var="CLEAN_RATE_${det}";  rate="${!rate_var}"
    flare_var="FLARE_EXPR_${det}"; flare="${!flare_var}"
    mkdir -p "$CLEANDIR/lightcurves/$det" "$CLEANDIR/gti/$det"

    while IFS= read -r frame; do
      [[ -n "$frame" ]] || continue
      base="$(basename "${frame%.*}")"
      lc="$CLEANDIR/lightcurves/$det/${base}_flare_lc.fits"
      gti="$CLEANDIR/gti/$det/${base}_flare_gti.fits"

      log_run "clean_lc_${det}_${base}" \
        evselect table="${frame}:EVENTS" \
          withrateset=yes rateset="$lc" \
          maketimecolumn=yes timecolumn=TIME \
          timebinsize="$CLEAN_GTI_TIMEBIN" \
          makeratecolumn=yes expression="$flare"
      log_run "clean_gti_${det}_${base}" \
        tabgtigen table="${lc}:RATE" gtiset="$gti" \
          expression="RATE <= $rate"
      printf '%s\t%s\t%s\t%s\t%s\t%s\n' \
        "$det" "$base" "$rate" "$CLEAN_GTI_TIMEBIN" "$lc" "$gti" \
        >> "$CLEANDIR/flare_gti_summary.tsv"

      for label in $CLEAN_BAND_LABELS; do
        mkdir -p "$CLEANDIR/events/$det/$label"
        clean="$CLEANDIR/events/$det/$label/${base}_${label}_clean.fits"
        expr_var="CLEAN_EXPR_${det}_${label}"; expr="${!expr_var}"
        log_run "clean_evt_${det}_${base}_${label}" \
          evselect table="${frame}:EVENTS" \
            withfilteredset=yes filteredset="$clean" \
            destruct=yes keepfilteroutput=yes updateexposure=yes \
            writedss=yes \
            expression="(${expr})&&gti(${gti},TIME)"
        rows="$("$PYTHON" "$SCRIPT_DIR/tools.py" fits-rows "$clean")"
        if [[ "$rows" -gt 0 ]]; then
          log_run "clean_attcalc_${det}_${base}_${label}" \
            attcalc eventset="$clean" \
              attitudelabel=ahf refpointlabel=pnt \
              withatthkset=yes atthkset="$ATTHK_FILE"
        fi
      done
    done < "$FRAMESDIR/manifest/${det}_frames.txt"

    for label in $CLEAN_BAND_LABELS; do
      "$PYTHON" "$SCRIPT_DIR/tools.py" clean-manifest \
        "$CLEANDIR/events/$det/$label" \
        "$CLEANDIR/manifest/${det}_${label}_clean.txt"
    done
  done

  echo "clean ok"
}

# ---------- stage: cut ----------

cut_have_products() {
  [[ -s "$CUTDIR/cut_summary.tsv" ]] || return 1
  local det label
  for det in $DETECTORS; do
    for label in $CUT_BAND_LABELS; do
      manifest_paths_exist "$CUTDIR/manifest/${det}_${label}_cut.txt" || return 1
    done
  done
}

stage_cut() {
  clean_have_products || stage_clean

  if [[ "$FORCE" == "1" ]]; then
    # When restricting to a subset of detectors, only nuke that subset's
    # cut events + manifests so others survive the rerun.
    local det label
    for det in $DETECTORS; do
      rm -rf "$CUTDIR/events/$det"
      for label in $CUT_BAND_LABELS; do
        rm -f "$CUTDIR/manifest/${det}_${label}_cut.txt"
      done
    done
  elif cut_have_products; then
    echo "cut ok: adopted existing products"
    return 0
  fi

  source_sas
  # shellcheck disable=SC1090
  source "$SAS_ENV"

  log_run cut_run \
    "$PYTHON" "$SCRIPT_DIR/tools.py" cut-run \
      "$CLEANDIR" "$FRAMESDIR" "$CUTDIR" "$LOGDIR" "$CONFIG" "$DETECTORS"
  echo "cut ok"
}

# ---------- stage: track ----------

track_have_products() {
  [[ -s "$TRACK_FITS" && -s "$TRACK_ENV" && -s "$TRACK_JSON" ]]
}

stage_track() {
  cut_have_products || stage_cut

  if [[ "$FORCE" == "1" ]]; then
    rm -f "$TRACK_FITS" "$TRACK_ENV" "$TRACK_JSON"
  elif track_have_products; then
    echo "track ok: adopted existing products"
    return 0
  fi

  local manifests=""
  for det in $DETECTORS; do
    manifests+="${manifests:+,}$FRAMESDIR/manifest/${det}_frames.txt"
  done

  log_run track_build \
    "$PYTHON" "$SCRIPT_DIR/tools.py" build-track \
      "$manifests" "$TARGET_ID" "$TARGET_ID_TYPE" "$TRACK_OBSERVER" \
      "$TRACK_STEP" "$TRACK_INPUT" \
      "$TRACK_FITS" "$TRACK_ENV" "$TRACK_JSON"
  echo "track ok"
}

# ---------- stage: maps ----------

_force_rm_per_det() {
  # Remove only the in-scope detectors' maps products with a given suffix.
  local suffix="$1" det
  for det in $DETECTORS; do
    rm -f "$MAPSDIR/$det"/*/"*${suffix}"
  done
}

stage_maps_counts() {
  cut_have_products || stage_cut
  source_sas
  # shellcheck disable=SC1090
  source "$SAS_ENV"
  if [[ "$FORCE" == "1" ]]; then
    _force_rm_per_det "_counts.fits"
  fi
  log_run maps_counts \
    "$PYTHON" "$SCRIPT_DIR/tools.py" maps-counts \
      "$CUTDIR" "$MAPSDIR" "$LOGDIR" "$CONFIG" "$DETECTORS"
  echo "maps-counts ok"
}

stage_maps_exposure() {
  [[ -s "$MAPSDIR/maps_grid.tsv" ]] || stage_maps_counts
  source_sas
  # shellcheck disable=SC1090
  source "$SAS_ENV"
  # shellcheck disable=SC1090
  source "$ATT_ENV"
  if [[ "$FORCE" == "1" ]]; then
    _force_rm_per_det "_exp_vig.fits"
  fi
  log_run maps_exposure \
    "$PYTHON" "$SCRIPT_DIR/tools.py" maps-exposure \
      "$CUTDIR" "$MAPSDIR" "$ATTHK_FILE" "$LOGDIR" "$CONFIG" "$DETECTORS"
  echo "maps-exposure ok"
}

stage_maps_background() {
  [[ -s "$MAPSDIR/maps_manifest.tsv" ]] || stage_maps_exposure
  if [[ "$FORCE" == "1" ]]; then
    _force_rm_per_det "_bkg.fits"
  fi
  log_run maps_background \
    "$PYTHON" "$SCRIPT_DIR/tools.py" maps-background \
      "$MAPSDIR" "$CONFIG" "$DETECTORS"
  echo "maps-background ok"
}

stage_maps_corrected() {
  [[ -s "$MAPSDIR/bkg_summary.tsv" ]] || stage_maps_background
  if [[ "$FORCE" == "1" ]]; then
    _force_rm_per_det "_corrected.fits"
  fi
  log_run maps_corrected \
    "$PYTHON" "$SCRIPT_DIR/tools.py" maps-corrected \
      "$MAPSDIR" "$CONFIG" "$DETECTORS"
  echo "maps-corrected ok"
}

stage_maps() {
  stage_maps_counts
  stage_maps_exposure
  stage_maps_background
  stage_maps_corrected
  echo "maps ok"
}

# ---------- QC stages ----------

stage_qc_init() {
  init_have_products || { echo "no init products to QC" >&2; exit 1; }
  local out="$QCDIR/init"
  rm -rf "$out"
  "$PYTHON" "$SCRIPT_DIR/tools.py" qc-init "$INITDIR" "$out" "$SAS_ENV"
  echo "qc-init wrote: $out"
}

stage_qc_repro() {
  [[ -d "$REPRODIR" ]] || { echo "missing $REPRODIR" >&2; exit 1; }
  local out="$QCDIR/repro"
  rm -rf "$out"
  "$PYTHON" "$SCRIPT_DIR/tools.py" qc-repro \
    "$REPRODIR" "$DETECTORS" "$out" "$QC_SPLIT_EV" "$ATT_ENV"
  echo "qc-repro wrote: $out"
}

stage_qc_frames() {
  [[ -d "$FRAMESDIR" ]] || { echo "missing $FRAMESDIR" >&2; exit 1; }
  local out="$QCDIR/frames"
  rm -rf "$out"
  "$PYTHON" "$SCRIPT_DIR/tools.py" qc-frames \
    "$FRAMESDIR" "$DETECTORS" "$out" "$QC_SPLIT_EV"
  echo "qc-frames wrote: $out"
}

stage_qc_clean() {
  [[ -d "$CLEANDIR" ]] || { echo "missing $CLEANDIR" >&2; exit 1; }
  [[ -d "$FRAMESDIR" ]] || { echo "missing $FRAMESDIR" >&2; exit 1; }
  local out="$QCDIR/clean"
  rm -rf "$out"
  "$PYTHON" "$SCRIPT_DIR/tools.py" qc-clean \
    "$CLEANDIR" "$FRAMESDIR" "$DETECTORS" "$QC_SPLIT_EV" "$out"
  echo "qc-clean wrote: $out"
}

stage_qc_cut() {
  [[ -d "$CUTDIR" ]]    || { echo "missing $CUTDIR" >&2; exit 1; }
  [[ -d "$FRAMESDIR" ]] || { echo "missing $FRAMESDIR" >&2; exit 1; }
  local out="$QCDIR/cut"
  rm -rf "$out"
  "$PYTHON" "$SCRIPT_DIR/tools.py" qc-cut \
    "$CUTDIR" "$FRAMESDIR" "$DETECTORS" "$out"
  echo "qc-cut wrote: $out"
}

stage_qc_track() {
  [[ -s "$TRACK_FITS" ]] || { echo "missing $TRACK_FITS (run track first)" >&2; exit 1; }
  [[ -d "$CUTDIR" ]]     || { echo "missing $CUTDIR" >&2; exit 1; }
  [[ -d "$FRAMESDIR" ]]  || { echo "missing $FRAMESDIR" >&2; exit 1; }
  local out="$QCDIR/track"
  rm -rf "$out"
  "$PYTHON" "$SCRIPT_DIR/tools.py" qc-track \
    "$TRACK_FITS" "$CUTDIR" "$FRAMESDIR" "$DETECTORS" "$out"
  echo "qc-track wrote: $out"
}

stage_qc_maps() {
  [[ -d "$MAPSDIR" ]]   || { echo "missing $MAPSDIR" >&2; exit 1; }
  [[ -d "$FRAMESDIR" ]] || { echo "missing $FRAMESDIR" >&2; exit 1; }
  local out="$QCDIR/maps"
  rm -rf "$out"
  "$PYTHON" "$SCRIPT_DIR/tools.py" qc-maps \
    "$MAPSDIR" "$FRAMESDIR" "$DETECTORS" "$out"
  echo "qc-maps wrote: $out"
}

# ---------- dispatch ----------

case "$STAGE" in
  init)             stage_init ;;
  repro)            stage_repro ;;
  frames)           stage_frames ;;
  clean)            stage_clean ;;
  cut)              stage_cut ;;
  track)            stage_track ;;
  maps)             stage_maps ;;
  maps-counts)      stage_maps_counts ;;
  maps-exposure)    stage_maps_exposure ;;
  maps-background)  stage_maps_background ;;
  maps-corrected)   stage_maps_corrected ;;
  all)              stage_init; stage_repro; stage_frames; stage_clean
                    stage_cut;  stage_track; stage_maps ;;
  qc-init)          stage_qc_init ;;
  qc-repro)         stage_qc_repro ;;
  qc-frames)        stage_qc_frames ;;
  qc-clean)         stage_qc_clean ;;
  qc-cut)           stage_qc_cut ;;
  qc-track)         stage_qc_track ;;
  qc-maps)          stage_qc_maps ;;
  qc)               stage_qc_init; stage_qc_repro; stage_qc_frames
                    stage_qc_clean; stage_qc_cut; stage_qc_track
                    stage_qc_maps ;;
esac
