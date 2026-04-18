#!/usr/bin/env bash
set -euo pipefail

# ==============================================================================
# Command-line interface
# ==============================================================================

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONFIG="$SCRIPT_DIR/config.json"
PYTHON="${PYTHON:-python3}"
STAGE="init"

usage() {
  cat <<'EOF'
Usage: ./qc.sh [init|repro|clean|exposure|all] [--config FILE]

Stages:
  init   Summarize initial SAS/ODF setup products.
  repro  Summarize reprocessed event lists and attitude products.
  clean  Summarize cleaned event lists and make PN soft/hard mosaics.
  exposure
         Summarize EPIC exposure products and render quick-look maps.
  all    Run QC for every implemented pipeline stage.
EOF
}

while (($#)); do
  case "$1" in
    init|repro|clean|exposure|all) STAGE="$1"; shift ;;
    --config|-c) CONFIG="$2"; shift 2 ;;
    --help|-h) usage; exit 0 ;;
    *) echo "Unknown argument: $1" >&2; usage >&2; exit 2 ;;
  esac
done

# ==============================================================================
# Configuration
# ==============================================================================

need_cmd() { command -v "$1" >/dev/null 2>&1 || { echo "Missing command: $1" >&2; exit 1; }; }
need_file() { [[ -f "$1" ]] || { echo "Missing file: $1" >&2; exit 1; }; }

need_cmd "$PYTHON"
need_file "$SCRIPT_DIR/tools.py"
need_file "$CONFIG"
eval "$("$PYTHON" "$SCRIPT_DIR/tools.py" shell "$CONFIG")"
[[ -n "${DETECTORS:-}" ]] || { echo "No detectors selected" >&2; exit 1; }

# ==============================================================================
# Stage QC: init
# ==============================================================================

file_type() {
  local name="$1"
  case "$name" in
    *SUM.SAS) echo "SUM.SAS" ;;
    ccf.cif) echo "ccf.cif" ;;
    *.FIT) echo "FIT" ;;
    *.ASC) echo "ASC" ;;
    MANIFEST.*) echo "MANIFEST" ;;
    *.*) echo "${name##*.}" ;;
    *) echo "other" ;;
  esac
}

qc_init() {
  local initdir="$WORKDIR/init"
  local qcdir="$WORKDIR/qc/init"
  local files="$qcdir/output_files.txt"
  local counts="$qcdir/file_type_counts.txt"

  [[ -d "$initdir" ]] || { echo "Missing init output directory: $initdir" >&2; exit 1; }
  mkdir -p "$qcdir"

  {
    [[ -f "$WORKDIR/sas_setup.env" ]] && echo "$WORKDIR/sas_setup.env"
    [[ -d "$WORKDIR/logs" ]] && find "$WORKDIR/logs" -maxdepth 1 -type f -name 'init_*.log' | sort
    find "$initdir" -maxdepth 1 \( -type f -o -type l \) | sort
  } > "$files"

  while IFS= read -r path; do
    file_type "$(basename "$path")"
  done < "$files" | sort | uniq -c | awk '{ printf "%s %s\n", $2, $1 }' > "$counts"

  echo "Wrote $files"
  echo "Wrote $counts"
}

# ==============================================================================
# Stage QC: repro
# ==============================================================================

count_lines() {
  [[ -s "$1" ]] && wc -l < "$1" || echo 0
}

qc_repro() {
  local reprodir="$WORKDIR/repro"
  local qcdir="$WORKDIR/qc/repro"
  local files="$qcdir/output_files.txt"
  local counts="$qcdir/manifest_counts.txt"
  local status="$qcdir/status.txt"

  [[ -d "$reprodir" ]] || { echo "Missing repro output directory: $reprodir" >&2; exit 1; }
  mkdir -p "$qcdir"

  {
    [[ -f "$WORKDIR/attitude.env" ]] && echo "$WORKDIR/attitude.env"
    [[ -d "$WORKDIR/logs" ]] && find "$WORKDIR/logs" -maxdepth 1 -type f -name 'repro_*.log' | sort
    find "$reprodir" -maxdepth 2 \( -type f -o -type l \) | sort
  } > "$files"

  {
    for inst in $DETECTORS; do
      echo "$inst $(count_lines "$reprodir/manifest/${inst}_raw.txt")"
    done
  } > "$counts"

  {
    [[ -s "$reprodir/atthk.dat" ]] && echo "atthk.dat ok" || echo "atthk.dat missing"
    [[ -s "$WORKDIR/attitude.env" ]] && echo "attitude.env ok" || echo "attitude.env missing"
    [[ "$DETECTORS" == *PN* ]] && { [[ -s "$WORKDIR/logs/repro_epproc.log" ]] && echo "repro_epproc.log ok" || echo "repro_epproc.log missing"; }
    [[ "$DETECTORS" == *M1* || "$DETECTORS" == *M2* ]] && { [[ -s "$WORKDIR/logs/repro_emproc.log" ]] && echo "repro_emproc.log ok" || echo "repro_emproc.log missing"; }
    [[ -s "$WORKDIR/logs/repro_atthkgen.log" ]] && echo "repro_atthkgen.log ok" || echo "repro_atthkgen.log missing"
  } > "$status"

  "$PYTHON" "$SCRIPT_DIR/tools.py" event-mosaic "$reprodir/manifest" "$qcdir" "$DETECTORS"

  echo "Wrote $files"
  echo "Wrote $counts"
  echo "Wrote $status"
  echo "Wrote $qcdir/soft_lt1kev_mosaic.png"
  echo "Wrote $qcdir/hard_gt1kev_mosaic.png"
}

# ==============================================================================
# Stage QC: clean
# ==============================================================================

qc_clean() {
  local cleandir="$WORKDIR/clean"
  local qcdir="$WORKDIR/qc/clean"
  local files="$qcdir/output_files.txt"
  local counts="$qcdir/manifest_counts.txt"
  local status="$qcdir/status.txt"

  [[ -d "$cleandir" ]] || { echo "Missing clean output directory: $cleandir" >&2; exit 1; }
  mkdir -p "$qcdir"

  {
    find "$cleandir" -maxdepth 2 \( -type f -o -type l \) | sort
    [[ -d "$WORKDIR/logs" ]] && find "$WORKDIR/logs" -maxdepth 1 -type f -name 'clean_*.log' | sort
  } > "$files"

  {
    for inst in $DETECTORS; do
      echo "$inst $(count_lines "$cleandir/${inst}_clean_files.txt")"
    done
  } > "$counts"

  {
    for inst in $DETECTORS; do
      [[ -s "$cleandir/${inst}_clean_merged.fits" ]] && echo "${inst}_clean_merged.fits ok" || echo "${inst}_clean_merged.fits missing"
    done
  } > "$status"

  "$PYTHON" "$SCRIPT_DIR/tools.py" event-mosaic "$cleandir" "$qcdir" "$DETECTORS"

  echo "Wrote $files"
  echo "Wrote $counts"
  echo "Wrote $status"
  echo "Wrote $qcdir/soft_lt1kev_mosaic.png"
  echo "Wrote $qcdir/hard_gt1kev_mosaic.png"
}

# ==============================================================================
# Stage QC: exposure
# ==============================================================================

image_band_labels() {
  local band label
  IFS=';' read -r -a bands <<< "$IMAGE_BANDS"
  for band in "${bands[@]}"; do
    [[ -n "$band" ]] || continue
    IFS=':' read -r label _ <<< "$band"
    [[ -n "$label" ]] && echo "$label"
  done
}

qc_exposure() {
  local expdir="$WORKDIR/exposure"
  local qcdir="$WORKDIR/qc/exposure"
  local files="$qcdir/output_files.txt"
  local status="$qcdir/status.txt"
  local label

  [[ -d "$expdir" ]] || { echo "Missing exposure output directory: $expdir" >&2; exit 1; }
  mkdir -p "$qcdir"

  {
    find "$expdir" -maxdepth 2 \( -type f -o -type l \) | sort
    [[ -d "$WORKDIR/logs" ]] && find "$WORKDIR/logs" -maxdepth 1 -type f -name 'exposure_*.log' | sort
  } > "$files"

  {
    [[ -s "$expdir/grid.env" ]] && echo "grid.env ok" || echo "grid.env missing"
    for label in $(image_band_labels); do
      for inst in $DETECTORS; do
        [[ -s "$expdir/$label/${inst}_counts.fits" ]] && echo "$label/${inst}_counts.fits ok" || echo "$label/${inst}_counts.fits missing"
        [[ -s "$expdir/$label/${inst}_exposure.fits" ]] && echo "$label/${inst}_exposure.fits ok" || echo "$label/${inst}_exposure.fits missing"
      done
      for product in EPIC_counts EPIC_exposure EPIC_rate; do
        [[ -s "$expdir/$label/${product}.fits" ]] && echo "$label/${product}.fits ok" || echo "$label/${product}.fits missing"
      done
    done
  } > "$status"

  for label in $(image_band_labels); do
    [[ -s "$expdir/$label/EPIC_exposure.fits" ]] && \
      "$PYTHON" "$SCRIPT_DIR/tools.py" fits-png "$expdir/$label/EPIC_exposure.fits" "$qcdir/${label}_exposure.png" viridis linear
    [[ -s "$expdir/$label/EPIC_rate.fits" ]] && \
      "$PYTHON" "$SCRIPT_DIR/tools.py" fits-png "$expdir/$label/EPIC_rate.fits" "$qcdir/${label}_rate.png" magma log
  done

  echo "Wrote $files"
  echo "Wrote $status"
}

# ==============================================================================
# Entrypoint
# ==============================================================================

case "$STAGE" in
  init) qc_init ;;
  repro) qc_repro ;;
  clean) qc_clean ;;
  exposure) qc_exposure ;;
  all) qc_init; qc_repro; qc_clean; qc_exposure ;;
  *) echo "Unknown QC stage: $STAGE" >&2; exit 2 ;;
esac
