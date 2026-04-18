#!/usr/bin/env bash
set -euo pipefail

# ==============================================================================
# Command-line interface
# ==============================================================================

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONFIG="$SCRIPT_DIR/config.json"
PYTHON="${PYTHON:-python3}"
STAGE="init"
FORCE=0
PRINT_ENV=0
NO_SHORTLINK=0

usage() {
  cat <<'EOF'
Usage: ./pipeline.sh [init|repro|clean|exposure|background|all|env] [--config FILE] [--force] [--print-env] [--no-shortlink]

Stages:
  init         Prepare the SAS ODF state used by all later stages:
                 output/init/ccf.cif
                 output/init/*SUM.SAS
                 output/init/<raw ODF symlinks>
                 output/sas_setup.env

  repro        Run epproc/emproc and create attitude products.

  clean        Apply standard EPIC event-quality and energy cuts.

  exposure     Build soft/hard EPIC counts, exposure, and rate maps.

  background   Build per-pointing quantile background maps.

  all          Run every implemented science-production stage in order.

  env          Print output/sas_setup.env after init has run.

Options:
  --config     JSON config file. Default: ./config.json
  --force      Rebuild init products.
  --print-env  Print the saved SAS environment after init.
  --no-shortlink
               Use the configured workdir directly instead of a /tmp shortlink.
EOF
}

while (($#)); do
  case "$1" in
    init|repro|clean|exposure|background|all|env) STAGE="$1"; shift ;;
    --config|-c) CONFIG="$2"; shift 2 ;;
    --force|-f) FORCE=1; shift ;;
    --print-env) PRINT_ENV=1; shift ;;
    --no-shortlink) NO_SHORTLINK=1; shift ;;
    --help|-h) usage; exit 0 ;;
    *) echo "Unknown argument: $1" >&2; usage >&2; exit 2 ;;
  esac
done

# ==============================================================================
# Small utilities
# ==============================================================================

need_cmd() { command -v "$1" >/dev/null 2>&1 || { echo "Missing command: $1" >&2; exit 1; }; }
need_file() { [[ -f "$1" ]] || { echo "Missing file: $1" >&2; exit 1; }; }
is_yes() { [[ "${1,,}" == "1" || "${1,,}" == "y" || "${1,,}" == "yes" || "${1,,}" == "true" ]]; }

# ==============================================================================
# Configuration and output layout
# ==============================================================================

need_cmd "$PYTHON"
need_file "$SCRIPT_DIR/tools.py"
need_file "$CONFIG"
eval "$("$PYTHON" "$SCRIPT_DIR/tools.py" shell "$CONFIG")"
[[ -n "${DETECTORS:-}" ]] || { echo "No detectors selected" >&2; exit 1; }

REAL_WORKDIR="$WORKDIR"
mkdir -p "$REAL_WORKDIR"

if [[ "$NO_SHORTLINK" != "1" && "${#REAL_WORKDIR}" -gt "$SHORTLINK_MAX_PATH" ]]; then
  SHORTLINK="/tmp/_xmm_${SHORTLINK_NAME:-v6_$$}"
  rm -f "$SHORTLINK"
  ln -sfn "$REAL_WORKDIR" "$SHORTLINK"
  WORKDIR="$SHORTLINK"
  is_yes "$KEEP_SHORTLINK" || trap 'rm -f "$SHORTLINK"' EXIT
fi

INITDIR="$WORKDIR/init"
REAL_INITDIR="$REAL_WORKDIR/init"
REPRODIR="$WORKDIR/repro"
REAL_REPRODIR="$REAL_WORKDIR/repro"
CLEANDIR="$WORKDIR/clean"
EXPOSUREDIR="$WORKDIR/exposure"
BACKGROUNDDIR="$WORKDIR/background"
LOGDIR="$WORKDIR/logs"
mkdir -p "$INITDIR" "$REPRODIR" "$CLEANDIR" "$EXPOSUREDIR" "$BACKGROUNDDIR" "$LOGDIR"

# ==============================================================================
# Logging and SAS environment setup
# ==============================================================================

runlog() {
  local name="$1"
  shift
  echo "[$(date -u +%FT%TZ)] $name" | tee -a "$LOGDIR/pipeline.log"
  "$@" 2>&1 | tee "$LOGDIR/${name}.log"
}

source_sas() {
  if command -v sasversion >/dev/null 2>&1 && [[ -n "${SAS_DIR:-}" && -n "${SAS_PATH:-}" ]]; then
    return 0
  fi

  local setup_script="${SAS_SETUP_SCRIPT:-}"
  if [[ -z "$setup_script" && -n "${SAS_DIR:-}" && -f "$SAS_DIR/setsas.sh" ]]; then
    setup_script="$SAS_DIR/setsas.sh"
  fi
  [[ -z "$setup_script" ]] && return 0
  need_file "$setup_script"
  set +e +u
  # shellcheck source=/dev/null
  source "$setup_script"
  local status=$?
  set -euo pipefail
  if [[ "$status" -ne 0 ]]; then
    echo "Failed to source SAS setup script: $setup_script" >&2
    echo "Initialize HEASOFT/SAS first, or update sas_setup_script in config.json." >&2
    exit "$status"
  fi
}

set_sas_init_env() {
  export SAS_ODF="$ODFDIR"
  export SAS_CCF="$INITDIR/ccf.cif"
  [[ -n "$SAS_CCFPATH_CONFIG" ]] && export SAS_CCFPATH="$SAS_CCFPATH_CONFIG"
  [[ -n "$SAS_VERBOSITY_CONFIG" ]] && export SAS_VERBOSITY="$SAS_VERBOSITY_CONFIG"
  source_sas

  if [[ -z "${SAS_CCFPATH:-}" ]]; then
    echo "Warning: SAS_CCFPATH is unset; cifbuild may fail without a visible CCF repository." >&2
  fi
}

set_sas_repro_env() {
  need_file "$REAL_WORKDIR/sas_setup.env"
  # shellcheck source=/dev/null
  source "$REAL_WORKDIR/sas_setup.env"
  source_sas
  export SAS_CCF="$INITDIR/ccf.cif"
  local matches sum
  matches="$(find "$INITDIR" -maxdepth 1 -type f -name '*SUM.SAS' | sort)"
  sum="${matches%%$'\n'*}"
  [[ -n "$sum" ]] && export SAS_ODF="$sum"
  need_file "${SAS_CCF:-}"
  need_file "${SAS_ODF:-}"
}

set_sas_clean_env() {
  set_sas_repro_env
  need_file "$REAL_WORKDIR/attitude.env"
  # shellcheck source=/dev/null
  source "$REAL_WORKDIR/attitude.env"
  need_file "${ATTHKGEN_FILE:-}"
}

# ==============================================================================
# Stage: init
#
# Purpose:
#   Create the initial SAS products expected by downstream reduction stages.
#
# Outputs:
#   output/init/ccf.cif
#   output/init/*SUM.SAS
#   output/init/<symlinks to raw ODF constituents>
#   output/sas_setup.env
# ==============================================================================

mirror_odf_into_init() {
  is_yes "$LINK_ODF_CONSTITUENTS" || return 0

  local f base skip pat
  IFS='|' read -r -a SKIP_PATTERNS <<< "$SKIP_ODF_LINK_PATTERNS"
  for f in "$ODFDIR"/*; do
    [[ -f "$f" ]] || continue
    base="$(basename "$f")"
    skip=0
    for pat in "${SKIP_PATTERNS[@]}"; do
      [[ "$base" == $pat ]] && { skip=1; break; }
    done
    [[ "$skip" == "1" ]] && continue
    [[ -e "$INITDIR/$base" ]] || ln -sf "$f" "$INITDIR/$base"
  done
}

write_sas_setup() {
  local sum_real="$1"
  {
    echo "# Generated by reduction_v6/pipeline.sh"
    printf 'export SAS_CCF=%q\n' "$REAL_INITDIR/ccf.cif"
    printf 'export SAS_ODF=%q\n' "$sum_real"
    [[ -n "${SAS_CCFPATH:-}" ]] && printf 'export SAS_CCFPATH=%q\n' "$SAS_CCFPATH"
    [[ -n "${SAS_VERBOSITY:-}" ]] && printf 'export SAS_VERBOSITY=%q\n' "$SAS_VERBOSITY"
  } > "$REAL_WORKDIR/sas_setup.env"
}

print_env() {
  need_file "$REAL_WORKDIR/sas_setup.env"
  cat "$REAL_WORKDIR/sas_setup.env"
}

stage_init() {
  if [[ "$FORCE" != "1" && -s "$REAL_WORKDIR/sas_setup.env" ]]; then
    # shellcheck source=/dev/null
    source "$REAL_WORKDIR/sas_setup.env"
    if [[ -s "${SAS_CCF:-}" && -s "${SAS_ODF:-}" ]]; then
      echo "Skipping init; found $REAL_WORKDIR/sas_setup.env"
      [[ "$PRINT_ENV" == "1" ]] && print_env
      return 0
    fi
    echo "Found stale sas_setup.env; rebuilding init products"
  fi

  [[ -d "$ODFDIR" ]] || { echo "ODF directory does not exist: $ODFDIR" >&2; exit 1; }

  if [[ "$FORCE" == "1" ]]; then
    rm -f "$INITDIR/ccf.cif" "$INITDIR"/*SUM.SAS "$REAL_WORKDIR/sas_setup.env"
  fi

  set_sas_init_env
  need_cmd cifbuild
  need_cmd odfingest

  pushd "$INITDIR" >/dev/null
  runlog init_cifbuild cifbuild
  export SAS_CCF="$INITDIR/ccf.cif"
  runlog init_odfingest odfingest odfdir="$ODFDIR" outdir="$INITDIR"
  popd >/dev/null

  local sum real_sum
  local matches
  matches="$(find "$INITDIR" -maxdepth 1 -type f -name '*SUM.SAS' | sort)"
  sum="${matches%%$'\n'*}"
  [[ -n "$sum" ]] || { echo "odfingest did not create a *SUM.SAS file in $INITDIR" >&2; exit 1; }

  real_sum="$REAL_INITDIR/$(basename "$sum")"
  mirror_odf_into_init
  "$PYTHON" "$SCRIPT_DIR/tools.py" rewrite-path "$sum" "$REAL_INITDIR/"
  write_sas_setup "$real_sum"

  echo "Wrote $REAL_WORKDIR/sas_setup.env"
  [[ "$PRINT_ENV" == "1" ]] && print_env
}

# ==============================================================================
# Stage: repro
#
# Purpose:
#   Reprocess EPIC ODF data into calibrated event lists and attitude products.
#
# Outputs:
#   output/repro/*ImagingEvts.ds
#   output/repro/manifest/PN_raw.txt, M1_raw.txt, M2_raw.txt
#   output/repro/atthk.dat
#   output/attitude.env
# ==============================================================================

write_manifest() {
  local outfile="$1"
  local root="$2"
  shift 2

  mkdir -p "$(dirname "$outfile")"
  shopt -s globstar nullglob
  for pattern in "$@"; do
    for f in "$root"/**/$pattern; do
      [[ -f "$f" ]] && readlink -f "$f"
    done
  done | awk '!seen[$0]++' > "$outfile"
}

choose_original_ahf() {
  local pattern found
  for pattern in '*SC*ATS.FIT*' '*SC*ATS.FT*' '*ATS.FIT*' '*ATS.FT*' '*ATTTSR*.FIT*' '*AHF*.FIT*'; do
    found="$(find "$ODFDIR" -type f -name "$pattern" | sort)"
    found="${found%%$'\n'*}"
    [[ -n "$found" ]] && { readlink -f "$found"; return 0; }
  done
  return 1
}

stage_repro() {
  set_sas_repro_env
  need_cmd atthkgen

  mkdir -p "$REPRODIR/manifest"
  pushd "$REPRODIR" >/dev/null

  local want_pn=0 want_mos=0 mos_missing=0 inst
  for inst in $DETECTORS; do
    [[ "$inst" == "PN" ]] && want_pn=1
    if [[ "$inst" == "M1" || "$inst" == "M2" ]]; then
      want_mos=1
      [[ ! -s "$REPRODIR/manifest/${inst}_raw.txt" ]] && mos_missing=1
    fi
  done

  if [[ "$want_pn" == "1" && ( "$FORCE" == "1" || ! -s "$REPRODIR/manifest/PN_raw.txt" ) ]]; then
    need_cmd epproc
    runlog repro_epproc epproc
  fi
  if [[ "$want_mos" == "1" && ( "$FORCE" == "1" || "$mos_missing" == "1" ) ]]; then
    need_cmd emproc
    runlog repro_emproc emproc
  fi

  for inst in $DETECTORS; do
    case "$inst" in
      PN) write_manifest "$REPRODIR/manifest/PN_raw.txt" "$REPRODIR" '*EPN*ImagingEvts.ds' '*EPN*ImagingEvts*.FIT*' '*PIEVLI*.FIT*' '*EPN*EVLI*.FIT*' ;;
      M1) write_manifest "$REPRODIR/manifest/M1_raw.txt" "$REPRODIR" '*EMOS1*ImagingEvts.ds' '*EMOS1*ImagingEvts*.FIT*' '*M1EVLI*.FIT*' '*EMOS1*EVLI*.FIT*' ;;
      M2) write_manifest "$REPRODIR/manifest/M2_raw.txt" "$REPRODIR" '*EMOS2*ImagingEvts.ds' '*EMOS2*ImagingEvts*.FIT*' '*M2EVLI*.FIT*' '*EMOS2*EVLI*.FIT*' ;;
    esac
  done

  local have_events=0
  for inst in $DETECTORS; do
    [[ -s "$REPRODIR/manifest/${inst}_raw.txt" ]] && have_events=1
  done
  [[ "$have_events" == "1" ]] || { echo "No selected EPIC imaging event lists found after repro" >&2; exit 1; }

  if [[ "$FORCE" == "1" || ! -s "$REPRODIR/atthk.dat" ]]; then
    runlog repro_atthkgen atthkgen atthkset="$REPRODIR/atthk.dat"
  fi

  local ahf="$AHF_INPUT"
  [[ -n "$ahf" ]] || ahf="$(choose_original_ahf)" || {
    echo "Could not locate an original spacecraft attitude file in $ODFDIR" >&2
    exit 1
  }
  [[ -f "$ahf" ]] || { echo "AHF/ATS file not found: $ahf" >&2; exit 1; }

  {
    printf 'ATTHKGEN_FILE=%q\n' "$REAL_REPRODIR/atthk.dat"
    printf 'ATTHK_FILE=%q\n' "$REAL_REPRODIR/atthk.dat"
    printf 'AHF_FILE=%q\n' "$ahf"
  } > "$REAL_WORKDIR/attitude.env"

  popd >/dev/null
  echo "Wrote $REAL_WORKDIR/attitude.env"
}

# ==============================================================================
# Stage: clean
#
# Purpose:
#   Apply standard EPIC event-quality and PI cuts, without espfilt.
#
# Outputs:
#   output/clean/<inst>/*.clean.fits
#   output/clean/<inst>_clean_files.txt
#   output/clean/<inst>_clean_merged.fits
# ==============================================================================

clean_expr() {
  local inst="$1"
  if [[ "$inst" == "PN" ]]; then
    echo "(FLAG==0)&&(PATTERN<=4)&&(PI in [$CLEAN_PI_MIN:$CLEAN_PI_MAX])"
  else
    echo "#XMMEA_EM&&(PATTERN<=12)&&(PI in [$CLEAN_PI_MIN:$CLEAN_PI_MAX])"
  fi
}

merge_pairwise() {
  local outfile="$1"
  shift
  local -a files=("$@")
  [[ ${#files[@]} -gt 0 ]] || return 1
  if [[ ${#files[@]} -eq 1 ]]; then
    cp -f "${files[0]}" "$outfile"
    return 0
  fi

  local tmp="${files[0]}" next outtmp i=1
  local -a made=()
  while [[ $i -lt ${#files[@]} ]]; do
    next="${files[$i]}"
    outtmp="$outfile.tmp${i}.fits"
    merge set1="$tmp" set2="$next" outset="$outtmp" imagesize="$PREMERGE_IMAGE_SIZE_DEG"
    made+=("$outtmp")
    tmp="$outtmp"
    i=$((i + 1))
  done
  mv -f "$tmp" "$outfile"
  for tmp in "${made[@]}"; do
    [[ "$tmp" == "$outfile" ]] || rm -f "$tmp"
  done
}

clean_one_event() {
  local inst="$1"
  local evt="$2"
  local outfile="$3"
  local expr
  expr="$(clean_expr "$inst")"

  if [[ "$FORCE" != "1" && -s "$outfile" ]]; then
    return 0
  fi

  runlog "clean_${inst}_$(basename "${outfile%.fits}")" \
    evselect table="${evt}:EVENTS" \
      withfilteredset=yes filteredset="$outfile" \
      destruct=yes keepfilteroutput=yes updateexposure=yes writedss=yes \
      expression="$expr"

  if [[ "$("$PYTHON" "$SCRIPT_DIR/tools.py" fits-rows "$outfile")" -le 0 ]]; then
    echo "No events left after cleaning: $outfile" >&2
    rm -f "$outfile"
    return 0
  fi

  attcalc eventset="$outfile" attitudelabel=ahf refpointlabel=pnt withatthkset=yes atthkset="$ATTHKGEN_FILE"
}

write_clean_manifest() {
  local inst="$1"
  find "$CLEANDIR/$inst" -maxdepth 1 -type f -name '*.clean.fits' | sort > "$CLEANDIR/${inst}_clean_files.txt"
}

stage_clean() {
  set_sas_clean_env
  need_cmd evselect
  need_cmd attcalc
  need_cmd merge

  local inst rawlist evt base outfile merged
  for inst in $DETECTORS; do
    rawlist="$REPRODIR/manifest/${inst}_raw.txt"
    [[ -s "$rawlist" ]] || continue
    mkdir -p "$CLEANDIR/$inst"

    while IFS= read -r evt; do
      [[ -n "$evt" ]] || continue
      base="$(basename "$evt")"
      base="${base%.*}"
      outfile="$CLEANDIR/$inst/${base}.clean.fits"
      clean_one_event "$inst" "$evt" "$outfile"
    done < "$rawlist"

    write_clean_manifest "$inst"
    mapfile -t clean_files < "$CLEANDIR/${inst}_clean_files.txt"
    [[ ${#clean_files[@]} -gt 0 ]] || continue

    merged="$CLEANDIR/${inst}_clean_merged.fits"
    if [[ "$FORCE" == "1" || ! -s "$merged" ]]; then
      merge_pairwise "$merged" "${clean_files[@]}"
      attcalc eventset="$merged" attitudelabel=ahf refpointlabel=pnt withatthkset=yes atthkset="$ATTHKGEN_FILE"
    fi
  done
}

# ==============================================================================
# Stage: exposure
#
# Purpose:
#   Build matched EPIC counts, exposure, and count-rate maps from cleaned events.
#
# Outputs:
#   output/exposure/grid.env, grid.json
#   output/exposure/pointings.tsv, slices.tsv
#   output/exposure/<band>/<inst>/<event>_<pointing>_counts.fits
#   output/exposure/<band>/<inst>/<event>_<pointing>_exposure.fits
#   output/exposure/<band>/<inst>_counts.fits
#   output/exposure/<band>/<inst>_exposure.fits
#   output/exposure/<band>/EPIC_counts.fits
#   output/exposure/<band>/EPIC_exposure.fits
#   output/exposure/<band>/EPIC_rate.fits
# ==============================================================================

image_band_lines() {
  local band label pimin pimax
  IFS=';' read -r -a bands <<< "$IMAGE_BANDS"
  for band in "${bands[@]}"; do
    [[ -n "$band" ]] || continue
    IFS=':' read -r label pimin pimax <<< "$band"
    if [[ -z "${label:-}" || -z "${pimin:-}" || -z "${pimax:-}" ]]; then
      echo "Bad image_bands entry: $band" >&2
      exit 1
    fi
    printf '%s %s %s\n' "$label" "$pimin" "$pimax"
  done
}

stage_exposure() {
  set_sas_clean_env
  need_cmd evselect
  need_cmd eexpmap

  local inst evt label pimin pimax banddir instdir base point start stop ref
  local counts exposure slices pointings
  local -a events=()
  for inst in $DETECTORS; do
    manifest="$CLEANDIR/${inst}_clean_files.txt"
    [[ -s "$manifest" ]] || continue
    while IFS= read -r evt; do
      [[ -n "$evt" ]] && events+=("$evt")
    done < "$manifest"
  done
  [[ ${#events[@]} -gt 0 ]] || { echo "No clean event files found in $CLEANDIR" >&2; exit 1; }

  if [[ "$FORCE" == "1" || ! -s "$EXPOSUREDIR/grid.env" ]]; then
    "$PYTHON" "$SCRIPT_DIR/tools.py" exposure-grid "$EXPOSUREDIR" "$IMAGE_BIN_PHYS" "$IMAGE_PAD_FRAC" "${events[@]}"
  fi
  slices="$EXPOSUREDIR/slices.tsv"
  pointings="$EXPOSUREDIR/pointings.tsv"
  "$PYTHON" "$SCRIPT_DIR/tools.py" exposure-slices "$slices" "$pointings" "$CLEANDIR" "$DETECTORS" "$INITDIR"
  # shellcheck source=/dev/null
  source "$EXPOSUREDIR/grid.env"
  find "$EXPOSUREDIR" -type f -name '*_events.fits' -delete

  while read -r label pimin pimax; do
    [[ -n "$label" ]] || continue
    banddir="$EXPOSUREDIR/$label"
    mkdir -p "$banddir"
    local -a epic=()

    for inst in $DETECTORS; do
      instdir="$banddir/$inst"
      mkdir -p "$instdir"
      local -a det=()

      while IFS=$'\t' read -r slice_inst evt base point start stop ref; do
        [[ "$slice_inst" == "$inst" ]] || continue
        [[ "$slice_inst" == "inst" ]] && continue
        counts="$instdir/${base}_counts.fits"
        exposure="$instdir/${base}_exposure.fits"

        if [[ "$FORCE" == "1" || ! -s "$counts" ]]; then
          runlog "exposure_${label}_${inst}_${base}_counts" \
            evselect table="${evt}:EVENTS" \
              withimageset=yes imageset="$counts" \
              xcolumn=X ycolumn=Y imagebinning=binSize \
              ximagebinsize="$BIN_PHYS" yimagebinsize="$BIN_PHYS" \
              withxranges=yes ximagemin="$X_MIN_PHYS" ximagemax="$X_MAX_PHYS" \
              withyranges=yes yimagemin="$Y_MIN_PHYS" yimagemax="$Y_MAX_PHYS" \
              ignorelegallimits=yes \
              writedss=yes \
              expression="(TIME>=$start)&&(TIME<=$stop)&&(PI in [$pimin:$pimax])"
        fi

        if [[ "$FORCE" == "1" || ! -s "$exposure" ]]; then
          runlog "exposure_${label}_${inst}_${base}_eexpmap" \
            eexpmap imageset="$counts" \
              attitudeset="$ATTHKGEN_FILE" \
              eventset="$evt" \
              expimageset="$exposure" \
              pimin="$pimin" pimax="$pimax" \
              withvignetting=yes \
              attrebin="$EEXPMAP_ATTREBIN"
        fi
        det+=("$counts" "$exposure")
      done < "$slices"

      [[ ${#det[@]} -gt 0 ]] || continue
      runlog "exposure_${label}_${inst}_combine" \
        "$PYTHON" "$SCRIPT_DIR/tools.py" combine-maps \
          "$banddir/${inst}_counts.fits" \
          "$banddir/${inst}_exposure.fits" \
          "$banddir/${inst}_rate.fits" \
          "${det[@]}"
      epic+=("$banddir/${inst}_counts.fits" "$banddir/${inst}_exposure.fits")
    done

    [[ ${#epic[@]} -gt 0 ]] || continue
    runlog "exposure_${label}_combine" \
      "$PYTHON" "$SCRIPT_DIR/tools.py" combine-maps \
        "$banddir/EPIC_counts.fits" \
        "$banddir/EPIC_exposure.fits" \
        "$banddir/EPIC_rate.fits" \
        "${epic[@]}"
  done < <(image_band_lines)
}

# ==============================================================================
# Stage: background
#
# Purpose:
#   Build constant per-pointing quantile background maps in counts space.
#
# Outputs:
#   output/background/<band>/<inst>/<event>_<pointing>_{background,net_counts}.fits
#   output/background/<band>/<inst>_{background,net_counts}.fits
#   output/background/<band>/EPIC_{background,net_counts}.fits
#   output/background/<band>/background_summary.tsv
# ==============================================================================

stage_background() {
  need_file "$EXPOSUREDIR/grid.json"
  need_file "$EXPOSUREDIR/slices.tsv"

  local label pimin pimax banddir summary inst slice_inst evt base point start stop ref instdir
  local counts exposure background net
  while read -r label pimin pimax; do
    [[ -n "$label" ]] || continue
    banddir="$BACKGROUNDDIR/$label"
    summary="$banddir/background_summary.tsv"
    mkdir -p "$banddir"
    rm -f "$summary"
    find "$banddir" -type f \( -name '*_background.fits' -o -name '*_net_counts.fits' \) -delete
    local -a epic=()

    for inst in $DETECTORS; do
      instdir="$banddir/$inst"
      mkdir -p "$instdir"
      local -a det=()

      while IFS=$'\t' read -r slice_inst evt base point start stop ref; do
        [[ "$slice_inst" == "$inst" ]] || continue
        background="$instdir/${base}_background.fits"
        net="$instdir/${base}_net_counts.fits"
        counts="$EXPOSUREDIR/$label/$inst/${base}_counts.fits"
        exposure="$EXPOSUREDIR/$label/$inst/${base}_exposure.fits"
        [[ -s "$counts" && -s "$exposure" ]] || continue

        "$PYTHON" "$SCRIPT_DIR/tools.py" background-map \
          "$counts" "$exposure" "$background" "$net" \
          "$BACKGROUND_COUNTS_QUANTILE" "$BACKGROUND_MIN_PIXELS" \
          "$BACKGROUND_EXCLUDE_RADIUS_FRAC" \
          "$summary" "$label" "$inst" "$base"
        det+=("$background" "$net")
      done < "$EXPOSUREDIR/slices.tsv"

      [[ ${#det[@]} -gt 0 ]] || continue
      "$PYTHON" "$SCRIPT_DIR/tools.py" combine-background \
        "$banddir/${inst}_background.fits" \
        "$banddir/${inst}_net_counts.fits" \
        "${det[@]}"
      epic+=("$banddir/${inst}_background.fits" "$banddir/${inst}_net_counts.fits")
    done

    [[ ${#epic[@]} -gt 0 ]] || continue
    "$PYTHON" "$SCRIPT_DIR/tools.py" combine-background \
      "$banddir/EPIC_background.fits" \
      "$banddir/EPIC_net_counts.fits" \
      "${epic[@]}"
  done < <(image_band_lines)
}

# ==============================================================================
# Entrypoint
# ==============================================================================

case "$STAGE" in
  init) stage_init ;;
  repro) stage_repro ;;
  clean) stage_clean ;;
  exposure) stage_exposure ;;
  background) stage_background ;;
  all) stage_init; stage_repro; stage_clean; stage_exposure; stage_background ;;
  env) print_env ;;
  *) echo "Unknown stage: $STAGE" >&2; exit 2 ;;
esac
