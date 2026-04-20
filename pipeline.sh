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
RUN_QC=0

usage() {
  cat <<'EOF'
Usage: ./pipeline.sh [stage] [--config FILE] [--force] [--qc] [--print-env] [--no-shortlink]

Stages:
  init        Prepare SAS ODF state.
  repro       Run epproc/emproc and create attitude products.
  clean       Apply EPIC quality, PI, and fixed high-energy flare-GTI cuts.
  exposure    Build counts, exposure, and rate maps.
  all         Run all implemented production stages.
  qc-init     QC init outputs.
  qc-repro    QC repro outputs.
  qc-clean    QC clean outputs.
  qc-exposure QC exposure outputs.
  qc          Run all QC stages.
  env         Print output/sas_setup.env.

Options:
  --config FILE   JSON config file. Default: ./config.json
  --force         Rebuild stage products.
  --qc            Run matching QC after the requested production stage.
  --print-env     Print saved SAS environment after init.
  --no-shortlink  Use configured workdir instead of a /tmp shortlink.
EOF
}

while (($#)); do
  case "$1" in
    init|repro|clean|exposure|all|qc-init|qc-repro|qc-clean|qc-exposure|qc|env) STAGE="$1"; shift ;;
    --config|-c) CONFIG="$2"; shift 2 ;;
    --force|-f) FORCE=1; shift ;;
    --qc) RUN_QC=1; shift ;;
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
GTIDIR="$CLEANDIR/gti"
LIGHTCURVEDIR="$CLEANDIR/lightcurves"
GTI_SUMMARY="$CLEANDIR/flare_gti_summary.tsv"
EXPOSUREDIR="$WORKDIR/exposure"
LOGDIR="$WORKDIR/logs"
mkdir -p "$INITDIR" "$REPRODIR" "$CLEANDIR" "$EXPOSUREDIR" "$LOGDIR"

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
# Create ccf.cif, *SUM.SAS, ODF symlinks, and sas_setup.env.
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
# Reprocess EPIC data and write raw manifests, atthk.dat, and attitude.env.
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
# Apply EPIC event-quality, PI, and optional fixed high-energy flare-GTI cuts.
# ==============================================================================

clean_expr() {
  local inst="$1"
  if [[ "$inst" == "PN" ]]; then
    echo "(FLAG==0)&&(PATTERN<=4)&&(PI in [$CLEAN_PI_MIN:$CLEAN_PI_MAX])"
  else
    echo "#XMMEA_EM&&(PATTERN<=12)&&(PI in [$CLEAN_PI_MIN:$CLEAN_PI_MAX])"
  fi
}

clean_gti_enabled() {
  is_yes "${CLEAN_GTI_ENABLED:-yes}"
}

flare_lc_expr() {
  local inst="$1" qual
  if [[ "$inst" == "PN" ]]; then
    qual="FLAG==0"
  else
    qual="#XMMEA_EM"
  fi
  echo "(PATTERN<=$CLEAN_GTI_PATTERN_MAX)&&(PI in [$CLEAN_GTI_PI_MIN:$CLEAN_GTI_PI_MAX])&&($qual)"
}

build_flare_gti() {
  local inst="$1" evt="$2" base="$3"
  local lc gti expr
  mkdir -p "$LIGHTCURVEDIR/$inst" "$GTIDIR/$inst"

  lc="$LIGHTCURVEDIR/$inst/${base}_flare_lc.fits"
  gti="$GTIDIR/$inst/${base}_flare_gti.fits"
  if [[ "$FORCE" == "1" || ! -s "$gti" ]]; then
    expr="$(flare_lc_expr "$inst")"
    runlog "lc_${inst}_${base}" \
      evselect table="${evt}:EVENTS" withrateset=yes rateset="$lc" \
        maketimecolumn=yes timecolumn=TIME timebinsize="$CLEAN_GTI_TIMEBIN" \
        makeratecolumn=yes expression="$expr"
    runlog "gti_${inst}_${base}" \
      tabgtigen table="${lc}:RATE" gtiset="$gti" expression="RATE <= $CLEAN_GTI_RATE_CUT"
  fi

  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$inst" "$base" "$CLEAN_GTI_RATE_CUT" "$CLEAN_GTI_TIMEBIN" "$CLEAN_GTI_PI_MIN" "$CLEAN_GTI_PI_MAX" "$lc" "$gti" >> "$GTI_SUMMARY"
  CLEAN_GTI_RESULT="$gti"
}

clean_one_event() {
  local inst="$1"
  local evt="$2"
  local outfile="$3"
  local expr base
  expr="$(clean_expr "$inst")"

  if [[ "$FORCE" != "1" && -s "$outfile" ]]; then
    return 0
  fi

  base="$(basename "$outfile")"
  base="${base%.fits}"
  base="${base%_clean}"
  base="${base%.clean}"
  if clean_gti_enabled; then
    CLEAN_GTI_RESULT=""
    build_flare_gti "$inst" "$evt" "$base"
    expr="($expr)&&gti($CLEAN_GTI_RESULT,TIME)"
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
  mkdir -p "$CLEANDIR/manifest"
  find "$CLEANDIR/events/$inst" -maxdepth 1 -type f -name '*_clean.fits' | sort > "$CLEANDIR/manifest/${inst}_clean_files.txt"
}

stage_clean() {
  if [[ "$FORCE" != "1" ]] && clean_stage_complete; then
    echo "Skipping clean; found existing clean products in $REAL_WORKDIR/clean"
    return 0
  fi

  set_sas_clean_env
  need_cmd evselect
  need_cmd attcalc
  clean_gti_enabled && need_cmd tabgtigen

  if clean_gti_enabled; then
    mkdir -p "$GTIDIR"
    printf 'inst\tevent\trate_cut_ct_s\ttimebin_s\tpi_min\tpi_max\tlightcurve\tgti\n' > "$GTI_SUMMARY"
  fi

  local inst rawlist evt base outfile
  for inst in $DETECTORS; do
    rawlist="$REPRODIR/manifest/${inst}_raw.txt"
    [[ -s "$rawlist" ]] || continue
    mkdir -p "$CLEANDIR/events/$inst" "$CLEANDIR/manifest"

    while IFS= read -r evt; do
      [[ -n "$evt" ]] || continue
      base="$(basename "$evt")"
      base="${base%.*}"
      outfile="$CLEANDIR/events/$inst/${base}_clean.fits"
      clean_one_event "$inst" "$evt" "$outfile"
    done < "$rawlist"

    write_clean_manifest "$inst"
  done
}

# ==============================================================================
# Stage: exposure
#
# Build ODF-pointing counts, exposure, and rate maps from cleaned events.
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

band_count_expr() {
  local label="$1" pimin="$2" pimax="$3" start="$4" stop="$5"
  local time_expr="(TIME>=$start)&&(TIME<=$stop)"
  case "$label" in
    1000)
      echo "$time_expr&&(RAWY>=13)&&(PATTERN==0)&&(PI in [201:500])&&((FLAG & 0x2fb002c)==0)"
      ;;
    2000)
      echo "$time_expr&&(RAWY>=13)&&(PATTERN<=4)&&(PI in [501:1000])&&((FLAG & 0x2fb002c)==0)"
      ;;
    3000)
      echo "$time_expr&&(RAWY>=13)&&(PATTERN<=4)&&(PI in [1001:2000])&&((FLAG & 0x2fb0024)==0)"
      ;;
    4000)
      echo "$time_expr&&(RAWY>=13)&&(PATTERN<=4)&&(PI in [2001:4500])&&((FLAG & 0x2fb0024)==0)"
      ;;
    5000)
      echo "$time_expr&&(RAWY>=13)&&(PATTERN<=4)&&((PI in [4501:7800])||(PI in [8201:12000]))&&((FLAG & 0x2fb0024)==0)"
      ;;
    *)
      echo "$time_expr&&(PI in [$pimin:$pimax])"
      ;;
  esac
}

band_exposure_ref_expr() {
  local label="$1" start="$2" stop="$3"
  local time_expr="(TIME>=$start)&&(TIME<=$stop)"
  case "$label" in
    1000|2000)
      echo "$time_expr&&(RAWY>=13)&&((FLAG & 0x2fb002c)==0)"
      ;;
    3000|4000|5000)
      echo "$time_expr&&(RAWY>=13)&&((FLAG & 0x2fb0024)==0)"
      ;;
    *)
      echo "$time_expr"
      ;;
  esac
}

clean_manifest_for_detector() {
  local inst="$1"
  local candidate
  for candidate in "$CLEANDIR/${inst}_clean_files.txt" "$CLEANDIR/manifest/${inst}_clean_files.txt"; do
    if [[ -s "$candidate" ]]; then
      echo "$candidate"
      return 0
    fi
  done
  return 1
}

clean_stage_complete() {
  local inst manifest evt have=0
  for inst in $DETECTORS; do
    manifest="$(clean_manifest_for_detector "$inst" || true)"
    [[ -n "$manifest" ]] || return 1
    while IFS= read -r evt; do
      [[ -n "$evt" ]] || continue
      [[ -s "$evt" ]] || return 1
      have=1
    done < "$manifest"
  done
  [[ "$have" == "1" ]]
}

stage_exposure() {
  set_sas_clean_env
  need_cmd evselect
  need_cmd eexpmap

  local inst evt manifest label pimin pimax banddir instdir base point start stop ref
  local counts exposure novig_exposure refimage slices pointings expr refexpr
  local -a events=()
  for inst in $DETECTORS; do
    manifest="$(clean_manifest_for_detector "$inst" || true)"
    [[ -n "$manifest" ]] || continue
    while IFS= read -r evt; do
      [[ -n "$evt" ]] && events+=("$evt")
    done < "$manifest"
  done
  [[ ${#events[@]} -gt 0 ]] || { echo "No clean event files found in $CLEANDIR" >&2; exit 1; }

  if [[ "$FORCE" == "1" ]]; then
    find "$EXPOSUREDIR" -type f \
      \( -name '*_counts.fits' -o -name '*_exposure.fits' -o -name '*_novig_exposure.fits' -o -name '*_rate.fits' -o -name '*_exposure_ref.fits' \) -delete
  fi

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
    local -a epic_novig=()

    for inst in $DETECTORS; do
      instdir="$banddir/$inst"
      mkdir -p "$instdir"
      local -a det=()
      local -a det_novig=()

      while IFS=$'\t' read -r slice_inst evt base point start stop ref; do
        [[ "$slice_inst" == "$inst" ]] || continue
        [[ "$slice_inst" == "inst" ]] && continue
        counts="$instdir/${base}_counts.fits"
        refimage="$instdir/${base}_exposure_ref.fits"
        exposure="$instdir/${base}_exposure.fits"
        novig_exposure="$instdir/${base}_novig_exposure.fits"
        expr="$(band_count_expr "$label" "$pimin" "$pimax" "$start" "$stop")"
        refexpr="$(band_exposure_ref_expr "$label" "$start" "$stop")"

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
              expression="$expr"
        fi

        if [[ "$FORCE" == "1" || ! -s "$refimage" ]]; then
          runlog "exposure_${label}_${inst}_${base}_refimage" \
            evselect table="${evt}:EVENTS" \
              withimageset=yes imageset="$refimage" \
              xcolumn=X ycolumn=Y imagebinning=binSize \
              ximagebinsize="$BIN_PHYS" yimagebinsize="$BIN_PHYS" \
              withxranges=yes ximagemin="$X_MIN_PHYS" ximagemax="$X_MAX_PHYS" \
              withyranges=yes yimagemin="$Y_MIN_PHYS" yimagemax="$Y_MAX_PHYS" \
              ignorelegallimits=yes \
              writedss=yes \
              expression="$refexpr"
        fi

        if [[ "$FORCE" == "1" || ! -s "$exposure" ]]; then
          runlog "exposure_${label}_${inst}_${base}_eexpmap" \
            eexpmap imageset="$refimage" \
              attitudeset="$ATTHKGEN_FILE" \
              eventset="$evt" \
              expimageset="$exposure" \
              pimin="$pimin" pimax="$pimax" \
              withvignetting=yes \
              attrebin="$EEXPMAP_ATTREBIN"
        fi
        if [[ "$FORCE" == "1" || ! -s "$novig_exposure" ]]; then
          runlog "exposure_${label}_${inst}_${base}_eexpmap_novig" \
            eexpmap imageset="$refimage" \
              attitudeset="$ATTHKGEN_FILE" \
              eventset="$evt" \
              expimageset="$novig_exposure" \
              pimin="$pimin" pimax="$pimax" \
              withvignetting=no \
              attrebin="$EEXPMAP_ATTREBIN"
        fi
        det+=("$counts" "$exposure")
        det_novig+=("$novig_exposure")
      done < "$slices"

      [[ ${#det[@]} -gt 0 ]] || continue
      runlog "exposure_${label}_${inst}_combine" \
        "$PYTHON" "$SCRIPT_DIR/tools.py" combine-maps \
          "$banddir/${inst}_counts.fits" \
          "$banddir/${inst}_exposure.fits" \
          "$banddir/${inst}_rate.fits" \
          "${det[@]}"
      runlog "exposure_${label}_${inst}_combine_novig" \
        "$PYTHON" "$SCRIPT_DIR/tools.py" combine-exposures \
          "$banddir/${inst}_novig_exposure.fits" \
          "${det_novig[@]}"
      epic+=("$banddir/${inst}_counts.fits" "$banddir/${inst}_exposure.fits")
      epic_novig+=("$banddir/${inst}_novig_exposure.fits")
    done

    [[ ${#epic[@]} -gt 0 ]] || continue
    runlog "exposure_${label}_combine" \
      "$PYTHON" "$SCRIPT_DIR/tools.py" combine-maps \
        "$banddir/EPIC_counts.fits" \
        "$banddir/EPIC_exposure.fits" \
        "$banddir/EPIC_rate.fits" \
        "${epic[@]}"
    runlog "exposure_${label}_combine_novig" \
      "$PYTHON" "$SCRIPT_DIR/tools.py" combine-exposures \
        "$banddir/EPIC_novig_exposure.fits" \
        "${epic_novig[@]}"
  done < <(image_band_lines)
}

# ==============================================================================
# QC helpers
# ==============================================================================

file_type() {
  local name="$1"
  case "$name" in
    *SUM.SAS) echo "SUM.SAS" ;;
    ccf.cif) echo "ccf.cif" ;;
    *.FIT|*.FITZ|*.fits|*.fits.gz) echo "FITS" ;;
    *.ASC) echo "ASC" ;;
    MANIFEST.*) echo "MANIFEST" ;;
    *.*) echo "${name##*.}" ;;
    *) echo "other" ;;
  esac
}

count_lines() {
  [[ -s "$1" ]] && wc -l < "$1" || echo 0
}

image_band_labels() {
  local band label
  IFS=';' read -r -a bands <<< "$IMAGE_BANDS"
  for band in "${bands[@]}"; do
    [[ -n "$band" ]] || continue
    IFS=':' read -r label _ <<< "$band"
    [[ -n "$label" ]] && echo "$label"
  done
}

single_detector() {
  set -- $DETECTORS
  [[ $# -eq 1 ]] && echo "$1"
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

qc_clean() {
  local cleandir="$WORKDIR/clean"
  local qcdir="$WORKDIR/qc/clean"
  local files="$qcdir/output_files.txt"
  local counts="$qcdir/manifest_counts.txt"
  local status="$qcdir/status.txt"
  local gti_summary="$qcdir/flare_gti_summary.tsv"
  local inst manifest lc_root event_count manifest_count

  [[ -d "$cleandir" ]] || { echo "Missing clean output directory: $cleandir" >&2; exit 1; }
  mkdir -p "$qcdir"

  {
    find "$cleandir" -maxdepth 5 \( -type f -o -type l \) | sort
    [[ -d "$WORKDIR/logs" ]] && find "$WORKDIR/logs" -maxdepth 1 -type f -name 'clean_*.log' | sort
  } > "$files"

  {
    for inst in $DETECTORS; do
      manifest="$(clean_manifest_for_detector "$inst" || true)"
      echo "$inst $(count_lines "$manifest")"
    done
  } > "$counts"

  {
    for inst in $DETECTORS; do
      manifest="$(clean_manifest_for_detector "$inst" || true)"
      if [[ -n "$manifest" ]]; then
        manifest_count="$(count_lines "$manifest")"
        event_count="$(find "$cleandir/events/$inst" -maxdepth 1 -type f -name '*_clean.fits' 2>/dev/null | wc -l)"
        echo "$(basename "$manifest") ok"
        echo "${inst}_manifest_events $manifest_count"
        echo "${inst}_clean_event_files_on_disk $event_count"
      else
        echo "${inst}_clean_files.txt missing"
      fi
    done
  } > "$status"

  "$PYTHON" "$SCRIPT_DIR/tools.py" event-mosaic "$cleandir" "$qcdir" "$DETECTORS"
  "$PYTHON" "$SCRIPT_DIR/tools.py" slice-qc "$qcdir" "$cleandir" "$DETECTORS" "$WORKDIR/init" "$IMAGE_BANDS"
  [[ -s "$cleandir/flare_gti_summary.tsv" ]] && cp -f "$cleandir/flare_gti_summary.tsv" "$gti_summary"
  [[ -s "$cleandir/flare_summary.tsv" ]] && cp -f "$cleandir/flare_summary.tsv" "$qcdir/flare_summary.tsv"
  lc_root="$cleandir/gti"
  [[ -d "$cleandir/lightcurves" ]] && lc_root="$cleandir/lightcurves"
  "$PYTHON" "$SCRIPT_DIR/tools.py" flare-qc "$lc_root" "$qcdir" "$CLEAN_GTI_RATE_CUT"

  echo "Wrote $files"
  echo "Wrote $counts"
  echo "Wrote $status"
  [[ -s "$gti_summary" ]] && echo "Wrote $gti_summary"
  [[ -s "$qcdir/flare_summary.tsv" ]] && echo "Wrote $qcdir/flare_summary.tsv"
  [[ -s "$qcdir/flare_lightcurves.png" ]] && echo "Wrote $qcdir/flare_lightcurves.png"
  [[ -s "$qcdir/flare_lightcurves.tsv" ]] && echo "Wrote $qcdir/flare_lightcurves.tsv"
  echo "Wrote $qcdir/soft_lt1kev_mosaic.png"
  echo "Wrote $qcdir/hard_gt1kev_mosaic.png"
  echo "Wrote $qcdir/pre_exposure_*.tsv"
  echo "Wrote $qcdir/pre_exposure_summary.txt"
}

qc_exposure() {
  local expdir="$WORKDIR/exposure"
  local qcdir="$WORKDIR/qc/exposure"
  local files="$qcdir/output_files.txt"
  local status="$qcdir/status.txt"
  local label detector counts_map exposure_map novig_map rate_map

  [[ -d "$expdir" ]] || { echo "Missing exposure output directory: $expdir" >&2; exit 1; }
  mkdir -p "$qcdir"
  detector="$(single_detector)"

  {
    find "$expdir" -maxdepth 4 \( -type f -o -type l \) | sort
    [[ -d "$WORKDIR/logs" ]] && find "$WORKDIR/logs" -maxdepth 1 -type f -name 'exposure_*.log' | sort
  } > "$files"

  {
    [[ -s "$expdir/grid.env" ]] && echo "grid.env ok" || echo "grid.env missing"
    [[ -s "$expdir/pointings.tsv" ]] && echo "pointings $(($(wc -l < "$expdir/pointings.tsv") - 1))" || echo "pointings.tsv missing"
    [[ -s "$expdir/slices.tsv" ]] && echo "slices $(($(wc -l < "$expdir/slices.tsv") - 1))" || echo "slices.tsv missing"
    for label in $(image_band_labels); do
      for inst in $DETECTORS; do
        [[ -s "$expdir/$label/${inst}_counts.fits" ]] && echo "$label/${inst}_counts.fits ok" || echo "$label/${inst}_counts.fits missing"
        [[ -s "$expdir/$label/${inst}_exposure.fits" ]] && echo "$label/${inst}_exposure.fits ok" || echo "$label/${inst}_exposure.fits missing"
        [[ -s "$expdir/$label/${inst}_novig_exposure.fits" ]] && echo "$label/${inst}_novig_exposure.fits ok" || echo "$label/${inst}_novig_exposure.fits missing"
      done
      for product in EPIC_counts EPIC_exposure EPIC_novig_exposure EPIC_rate; do
        [[ -s "$expdir/$label/${product}.fits" ]] && echo "$label/${product}.fits ok" || echo "$label/${product}.fits missing"
      done
    done
  } > "$status"

  rm -f "$qcdir"/*.png
  for label in $(image_band_labels); do
    counts_map="$expdir/$label/EPIC_counts.fits"
    exposure_map="$expdir/$label/EPIC_exposure.fits"
    novig_map="$expdir/$label/EPIC_novig_exposure.fits"
    rate_map="$expdir/$label/EPIC_rate.fits"
    if [[ -n "$detector" ]]; then
      counts_map="$expdir/$label/${detector}_counts.fits"
      exposure_map="$expdir/$label/${detector}_exposure.fits"
      novig_map="$expdir/$label/${detector}_novig_exposure.fits"
      rate_map="$expdir/$label/${detector}_rate.fits"
    fi
    [[ -s "$counts_map" ]] && \
      "$PYTHON" "$SCRIPT_DIR/tools.py" fits-png "$counts_map" "$qcdir/${label}_counts_mosaic.png" magma log "$expdir/grid.json" "$WORKDIR/clean" "$DETECTORS"
    [[ -s "$exposure_map" ]] && \
      "$PYTHON" "$SCRIPT_DIR/tools.py" fits-png "$exposure_map" "$qcdir/${label}_exposure_mosaic.png" magma linear "$expdir/grid.json" "$WORKDIR/clean" "$DETECTORS"
    [[ -s "$novig_map" ]] && \
      "$PYTHON" "$SCRIPT_DIR/tools.py" fits-png "$novig_map" "$qcdir/${label}_novig_exposure_mosaic.png" magma linear "$expdir/grid.json" "$WORKDIR/clean" "$DETECTORS"
    [[ -s "$rate_map" ]] && \
      "$PYTHON" "$SCRIPT_DIR/tools.py" fits-png "$rate_map" "$qcdir/${label}_rate_mosaic.png" magma log "$expdir/grid.json" "$WORKDIR/clean" "$DETECTORS"
  done

  echo "Wrote $files"
  echo "Wrote $status"
  echo "Wrote $qcdir/*_mosaic.png"
}

run_qc_stage() {
  case "$1" in
    init) qc_init ;;
    repro) qc_repro ;;
    clean) qc_clean ;;
    exposure) qc_exposure ;;
    all|qc) qc_init; qc_repro; qc_clean; qc_exposure ;;
    *) echo "No QC stage for: $1" >&2; exit 2 ;;
  esac
}

maybe_run_qc() {
  [[ "$RUN_QC" == "1" ]] || return 0
  run_qc_stage "$STAGE"
}

# ==============================================================================
# Entrypoint
# ==============================================================================

case "$STAGE" in
  init) stage_init; maybe_run_qc ;;
  repro) stage_repro; maybe_run_qc ;;
  clean) stage_clean; maybe_run_qc ;;
  exposure) stage_exposure; maybe_run_qc ;;
  all) stage_init; stage_repro; stage_clean; stage_exposure; maybe_run_qc ;;
  qc-init) qc_init ;;
  qc-repro) qc_repro ;;
  qc-clean) qc_clean ;;
  qc-exposure) qc_exposure ;;
  qc) run_qc_stage all ;;
  env) print_env ;;
  *) echo "Unknown stage: $STAGE" >&2; exit 2 ;;
esac
