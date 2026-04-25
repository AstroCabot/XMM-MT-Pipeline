#!/usr/bin/env bash
set -Eeuo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONFIG="$SCRIPT_DIR/config.json"
PYTHON="${PYTHON:-python3}"
STAGE="init"
FORCE=0
PRINT_ENV=0
NO_SHORTLINK=0
RUN_QC=0
MAPS_RERUN_FROM=""
CHEESE_RERUN_FROM=""

usage() {
  cat <<'EOF'
Usage: ./pipeline.sh [stage] [--config FILE] [--force] [--qc] [--print-env] [--no-shortlink]

Stages:
  init       Prepare SAS ODF state.
  repro      Run epproc/emproc and write raw EPIC manifests.
  clean      Apply flare GTIs and PPS-style band event filters.
  maps       Make SAS count, exposure, and hard-band source products.
  cheese     Make hard-band point-source masks and masked map images.
  all        Run init, repro, clean, maps, and cheese.
  qc-init    QC init outputs.
  qc-repro   QC repro outputs.
  qc-clean   QC clean outputs.
  qc-maps    QC maps outputs.
  qc-cheese  QC cheese outputs.
  qc         Run all implemented QC stages.
  env        Print output/sas_setup.env.

Options:
  --config FILE   JSON config file. Default: ./config.json
  --force         Rebuild stage products.
  --qc            Run matching QC after a production stage.
  --print-env     Print saved SAS environment after init.
  --no-shortlink  Use configured workdir instead of a /tmp shortlink.
  --maps-from STEP
                  Rebuild maps products from STEP downstream. STEP is one of:
                  grid, counts, exposure, mask, sources.
  --cheese-from STEP
                  Rebuild cheese products from STEP downstream. STEP is one of:
                  regions, mask, images.
EOF
}

while (($#)); do
  case "$1" in
    init|repro|clean|maps|cheese|all|qc-init|qc-repro|qc-clean|qc-maps|qc-cheese|qc|env) STAGE="$1"; shift ;;
    --config|-c) CONFIG="$2"; shift 2 ;;
    --force|-f) FORCE=1; shift ;;
    --qc) RUN_QC=1; shift ;;
    --print-env) PRINT_ENV=1; shift ;;
    --no-shortlink) NO_SHORTLINK=1; shift ;;
    --maps-from) MAPS_RERUN_FROM="$2"; shift 2 ;;
    --cheese-from) CHEESE_RERUN_FROM="$2"; shift 2 ;;
    --help|-h) usage; exit 0 ;;
    *) echo "Unknown argument: $1" >&2; usage >&2; exit 2 ;;
  esac
done

need_cmd() { command -v "$1" >/dev/null 2>&1 || { echo "Missing command: $1" >&2; exit 1; }; }
need_file() { [[ -f "$1" ]] || { echo "Missing file: $1" >&2; exit 1; }; }
is_yes() { [[ "${1,,}" == "1" || "${1,,}" == "y" || "${1,,}" == "yes" || "${1,,}" == "true" ]]; }

need_cmd "$PYTHON"
need_file "$SCRIPT_DIR/tools.py"
need_file "$CONFIG"
eval "$("$PYTHON" "$SCRIPT_DIR/tools.py" shell "$CONFIG")"
[[ -n "${DETECTORS:-}" ]] || { echo "No detectors selected" >&2; exit 1; }

REAL_WORKDIR="$WORKDIR"
mkdir -p "$REAL_WORKDIR"

if [[ "$NO_SHORTLINK" != "1" && "${#REAL_WORKDIR}" -gt "$SHORTLINK_MAX_PATH" ]]; then
  SHORTLINK="/tmp/_xmm_${SHORTLINK_NAME:-v6_$$}"
  if [[ -e "$SHORTLINK" && ! -L "$SHORTLINK" ]]; then
    echo "Shortlink path exists and is not a symlink: $SHORTLINK" >&2
    exit 1
  fi
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
REAL_CLEANDIR="$REAL_WORKDIR/clean"
GTIDIR="$CLEANDIR/gti"
LIGHTCURVEDIR="$CLEANDIR/lightcurves"
GTI_SUMMARY="$CLEANDIR/flare_gti_summary.tsv"
MAPSDIR="$WORKDIR/maps"
REAL_MAPSDIR="$REAL_WORKDIR/maps"
CHEESEDIR="$WORKDIR/cheese"
REAL_CHEESEDIR="$REAL_WORKDIR/cheese"
LOGDIR="$WORKDIR/logs"
mkdir -p "$INITDIR" "$REPRODIR" "$CLEANDIR" "$MAPSDIR" "$CHEESEDIR" "$LOGDIR"

on_error() {
  local status=$? line="${BASH_LINENO[0]:-?}" cmd="${BASH_COMMAND:-?}"
  set +e
  echo "[$(date -u +%FT%TZ)] ERROR status=$status line=$line command=$cmd" | tee -a "$LOGDIR/pipeline.log" >&2
  return "$status"
}
trap on_error ERR

write_if_changed() {
  local path="$1" tmp
  tmp="${path}.tmp.$$"
  mkdir -p "$(dirname "$path")"
  cat > "$tmp"
  if [[ -e "$path" ]] && cmp -s "$tmp" "$path"; then
    rm -f "$tmp"
  else
    mv -f "$tmp" "$path"
  fi
}

runlog() {
  local name="$1" logfile status
  shift
  logfile="$LOGDIR/${name}.log"
  echo "[$(date -u +%FT%TZ)] $name" | tee -a "$LOGDIR/pipeline.log"
  set +e
  "$@" 2>&1 | tee "$logfile"
  status=${PIPESTATUS[0]}
  set -e
  if [[ "$status" -ne 0 ]]; then
    echo "[$(date -u +%FT%TZ)] $name FAILED status=$status log=$logfile" | tee -a "$LOGDIR/pipeline.log" >&2
    return "$status"
  fi
  echo "[$(date -u +%FT%TZ)] $name ok" >> "$LOGDIR/pipeline.log"
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
    exit "$status"
  fi
}

first_summary() {
  local dir="$1"
  [[ -d "$dir" ]] || return 0
  find "$dir" -maxdepth 1 -type f -name '*SUM.SAS' | sort | head -n 1
}

init_products_ready() {
  local sum
  sum="$(first_summary "$REAL_INITDIR")"
  [[ -s "$REAL_INITDIR/ccf.cif" && -n "$sum" && -s "$sum" ]]
}

mirror_odf_into_init() {
  is_yes "$LINK_ODF_CONSTITUENTS" || return 0
  [[ -d "$ODFDIR" ]] || return 0

  local f base skip pat
  IFS='|' read -r -a SKIP_PATTERNS <<< "$SKIP_ODF_LINK_PATTERNS"
  for f in "$ODFDIR"/*; do
    [[ -f "$f" ]] || continue
    base="$(basename "$f")"
    skip=0
    for pat in "${SKIP_PATTERNS[@]}"; do
      [[ -n "$pat" && "$base" == $pat ]] && { skip=1; break; }
    done
    [[ "$skip" == "1" ]] && continue
    [[ -e "$INITDIR/$base" ]] || ln -sf "$f" "$INITDIR/$base"
  done
}

write_sas_setup() {
  local sum_real="$1" ccfpath verbosity
  ccfpath="${SAS_CCFPATH:-$SAS_CCFPATH_CONFIG}"
  verbosity="${SAS_VERBOSITY:-$SAS_VERBOSITY_CONFIG}"
  {
    echo "# Generated by reduction_v6/pipeline.sh"
    printf 'export SAS_CCF=%q\n' "$REAL_INITDIR/ccf.cif"
    printf 'export SAS_ODF=%q\n' "$sum_real"
    if [[ -n "$ccfpath" ]]; then
      printf 'export SAS_CCFPATH=%q\n' "$ccfpath"
    fi
    if [[ -n "$verbosity" ]]; then
      printf 'export SAS_VERBOSITY=%q\n' "$verbosity"
    fi
  } | write_if_changed "$REAL_WORKDIR/sas_setup.env"
}

adopt_existing_init() {
  local sum
  sum="$(first_summary "$REAL_INITDIR")"
  [[ -s "$REAL_INITDIR/ccf.cif" && -n "$sum" && -s "$sum" ]] || return 1
  mirror_odf_into_init
  "$PYTHON" "$SCRIPT_DIR/tools.py" rewrite-path "$sum" "$REAL_INITDIR/"
  write_sas_setup "$sum"
}

print_env() {
  need_file "$REAL_WORKDIR/sas_setup.env"
  cat "$REAL_WORKDIR/sas_setup.env"
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

  local sum
  sum="$(first_summary "$INITDIR")"
  [[ -n "$sum" ]] && export SAS_ODF="$sum"
  export SAS_CCF="$INITDIR/ccf.cif"
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

stage_init() {
  if [[ "$FORCE" == "1" ]]; then
    rm -f "$INITDIR/ccf.cif" "$INITDIR"/*SUM.SAS "$REAL_WORKDIR/sas_setup.env"
  elif adopt_existing_init; then
    echo "Skipping init; found existing ccf.cif and SUM.SAS"
    [[ "$PRINT_ENV" == "1" ]] && print_env
    return 0
  fi

  [[ -d "$ODFDIR" ]] || { echo "ODF directory does not exist: $ODFDIR" >&2; exit 1; }

  set_sas_init_env
  need_cmd cifbuild
  need_cmd odfingest

  pushd "$INITDIR" >/dev/null
  runlog init_cifbuild cifbuild
  export SAS_CCF="$INITDIR/ccf.cif"
  runlog init_odfingest odfingest odfdir="$ODFDIR" outdir="$INITDIR"
  popd >/dev/null

  local sum real_sum
  sum="$(first_summary "$INITDIR")"
  [[ -n "$sum" ]] || { echo "odfingest did not create a *SUM.SAS file in $INITDIR" >&2; exit 1; }
  real_sum="$REAL_INITDIR/$(basename "$sum")"

  mirror_odf_into_init
  "$PYTHON" "$SCRIPT_DIR/tools.py" rewrite-path "$sum" "$REAL_INITDIR/"
  write_sas_setup "$real_sum"
  echo "Wrote $REAL_WORKDIR/sas_setup.env"
  [[ "$PRINT_ENV" == "1" ]] && print_env
}

write_manifest() {
  local outfile="$1" root="$2" tmp
  shift 2
  tmp="${outfile}.tmp.$$"
  mkdir -p "$(dirname "$outfile")"

  shopt -s globstar nullglob
  for pattern in "$@"; do
    for f in "$root"/**/$pattern; do
      [[ -f "$f" ]] && readlink -f "$f"
    done
  done | awk '!seen[$0]++' > "$tmp"

  if [[ -s "$tmp" || ! -e "$outfile" ]]; then
    if [[ -e "$outfile" ]] && cmp -s "$tmp" "$outfile"; then
      rm -f "$tmp"
    else
      mv -f "$tmp" "$outfile"
    fi
  else
    rm -f "$tmp"
  fi
}

write_raw_manifests() {
  local inst
  for inst in $DETECTORS; do
    case "$inst" in
      PN) write_manifest "$REPRODIR/manifest/PN_raw.txt" "$REPRODIR" '*EPN*ImagingEvts.ds' '*EPN*ImagingEvts*.FIT*' '*PIEVLI*.FIT*' '*EPN*EVLI*.FIT*' ;;
      M1) write_manifest "$REPRODIR/manifest/M1_raw.txt" "$REPRODIR" '*EMOS1*ImagingEvts.ds' '*EMOS1*ImagingEvts*.FIT*' '*M1EVLI*.FIT*' '*EMOS1*EVLI*.FIT*' ;;
      M2) write_manifest "$REPRODIR/manifest/M2_raw.txt" "$REPRODIR" '*EMOS2*ImagingEvts.ds' '*EMOS2*ImagingEvts*.FIT*' '*M2EVLI*.FIT*' '*EMOS2*EVLI*.FIT*' ;;
    esac
  done
}

manifest_paths_exist() {
  local manifest="$1" path have=0
  [[ -s "$manifest" ]] || return 1
  while IFS= read -r path; do
    [[ -n "$path" ]] || continue
    [[ -s "$path" ]] || return 1
    have=1
  done < "$manifest"
  [[ "$have" == "1" ]]
}

repro_products_ready() {
  local inst
  for inst in $DETECTORS; do
    manifest_paths_exist "$REAL_REPRODIR/manifest/${inst}_raw.txt" || return 1
  done
  [[ -s "$REAL_REPRODIR/atthk.dat" && -s "$REAL_WORKDIR/attitude.env" ]]
}

drop_stale_manifests() {
  local inst manifest
  for inst in $DETECTORS; do
    manifest="$REAL_REPRODIR/manifest/${inst}_raw.txt"
    [[ -e "$manifest" ]] || continue
    manifest_paths_exist "$manifest" || rm -f "$manifest"
  done
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

write_attitude_env() {
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
  } | write_if_changed "$REAL_WORKDIR/attitude.env"
}

stage_repro() {
  if ! init_products_ready; then
    stage_init
  else
    adopt_existing_init
  fi

  mkdir -p "$REPRODIR/manifest"
  if [[ "$FORCE" == "1" ]]; then
    rm -f "$REPRODIR/manifest/"*_raw.txt "$REPRODIR/atthk.dat" "$REAL_WORKDIR/attitude.env"
  else
    drop_stale_manifests
  fi
  write_raw_manifests

  if [[ "$FORCE" != "1" ]] && repro_products_ready; then
    echo "Skipping repro; found existing raw manifests, atthk.dat, and attitude.env"
    return 0
  fi

  set_sas_repro_env
  need_cmd atthkgen

  pushd "$REPRODIR" >/dev/null

  local want_pn=0 want_mos=0 inst pn_ready=0 mos_ready=1
  for inst in $DETECTORS; do
    [[ "$inst" == "PN" ]] && want_pn=1
    if [[ "$inst" == "M1" || "$inst" == "M2" ]]; then
      want_mos=1
      manifest_paths_exist "$REAL_REPRODIR/manifest/${inst}_raw.txt" || mos_ready=0
    fi
  done
  manifest_paths_exist "$REAL_REPRODIR/manifest/PN_raw.txt" && pn_ready=1

  if [[ "$want_pn" == "1" && ( "$FORCE" == "1" || "$pn_ready" == "0" ) ]]; then
    need_cmd epproc
    runlog repro_epproc epproc
  fi
  if [[ "$want_mos" == "1" && ( "$FORCE" == "1" || "$mos_ready" == "0" ) ]]; then
    need_cmd emproc
    runlog repro_emproc emproc
  fi

  write_raw_manifests
  for inst in $DETECTORS; do
    manifest_paths_exist "$REAL_REPRODIR/manifest/${inst}_raw.txt" || {
      echo "No selected $inst imaging event lists found after repro" >&2
      exit 1
    }
  done

  if [[ "$FORCE" == "1" || ! -s "$REPRODIR/atthk.dat" ]]; then
    runlog repro_atthkgen atthkgen atthkset="$REPRODIR/atthk.dat"
  fi
  popd >/dev/null

  write_attitude_env
  echo "Wrote $REAL_WORKDIR/attitude.env"
}

clean_gti_enabled() {
  is_yes "${CLEAN_GTI_ENABLED:-yes}"
}

clean_band_labels() {
  "$PYTHON" "$SCRIPT_DIR/tools.py" clean-band-labels "$CONFIG"
}

clean_expr() {
  "$PYTHON" "$SCRIPT_DIR/tools.py" clean-expr "$CONFIG" "$1" "$2"
}

maps_clean_labels() {
  "$PYTHON" "$SCRIPT_DIR/tools.py" maps-clean-labels "$CONFIG" "$1"
}

fits_table_hdu() {
  "$PYTHON" "$SCRIPT_DIR/tools.py" fits-table-hdu "$@"
}

fits_table_has_column() {
  local out
  out="$("$PYTHON" "$SCRIPT_DIR/tools.py" fits-table-has-column "$@")" || return $?
  [[ "$out" == "yes" ]]
}

flare_lc_expr() {
  local inst="$1" qual
  if [[ "$inst" == "PN" ]]; then
    qual="(FLAG==0)"
  else
    qual="#XMMEA_EM"
  fi
  echo "(PATTERN<=$CLEAN_GTI_PATTERN_MAX)&&(PI in [$CLEAN_GTI_PI_MIN:$CLEAN_GTI_PI_MAX])&&$qual"
}

build_flare_gti() {
  local inst="$1" evt="$2" base="$3"
  local lc gti expr lc_real gti_real
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

  lc_real="$(readlink -f "$lc" 2>/dev/null || echo "$lc")"
  gti_real="$(readlink -f "$gti" 2>/dev/null || echo "$gti")"
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$inst" "$base" "$CLEAN_GTI_RATE_CUT" "$CLEAN_GTI_TIMEBIN" "$CLEAN_GTI_PI_MIN" "$CLEAN_GTI_PI_MAX" "$lc_real" "$gti_real" >> "$GTI_SUMMARY"
  CLEAN_GTI_RESULT="$gti"
}

clean_manifest_valid() {
  local manifest="$1" path
  [[ -f "$manifest" ]] || return 1
  while IFS= read -r path; do
    [[ -n "$path" ]] || continue
    [[ -s "$path" ]] || return 1
  done < "$manifest"
}

clean_manifest_for_detector() {
  local inst="$1" candidate
  for candidate in "$REAL_CLEANDIR/manifest/${inst}_clean_files.txt" "$REAL_CLEANDIR/${inst}_clean_files.txt"; do
    if manifest_paths_exist "$candidate"; then
      echo "$candidate"
      return 0
    fi
  done
  return 1
}

clean_stage_complete() {
  local inst label manifest legacy
  if [[ -s "$REAL_CLEANDIR/.complete" ]]; then
    for inst in $DETECTORS; do
      while IFS= read -r label; do
        [[ -n "$label" ]] || continue
        manifest="$REAL_CLEANDIR/manifest/${inst}_${label}_clean_files.txt"
        clean_manifest_valid "$manifest" || return 1
      done < <(clean_band_labels)
    done
    return 0
  fi

  for inst in $DETECTORS; do
    legacy="$(clean_manifest_for_detector "$inst" || true)"
    [[ -n "$legacy" ]] || return 1
  done
}

write_clean_band_manifest() {
  local inst="$1" label="$2" out="$CLEANDIR/manifest/${inst}_${label}_clean_files.txt"
  mkdir -p "$CLEANDIR/manifest"
  {
    if [[ -d "$CLEANDIR/events/$inst/$label" ]]; then
      find "$CLEANDIR/events/$inst/$label" -maxdepth 1 -type f -name '*_clean.fits' | sort | while IFS= read -r path; do
        readlink -f "$path"
      done
    fi
  } | write_if_changed "$out"
}

write_clean_manifest() {
  local inst="$1" out="$CLEANDIR/manifest/${inst}_clean_files.txt"
  mkdir -p "$CLEANDIR/manifest"
  {
    if [[ -d "$CLEANDIR/events/$inst" ]]; then
      find "$CLEANDIR/events/$inst" -type f -name '*_clean.fits' | sort | while IFS= read -r path; do
        readlink -f "$path"
      done
    fi
  } | write_if_changed "$out"
}

clean_one_event() {
  local inst="$1" evt="$2" label="$3" outfile="$4" gti="$5"
  local expr rows

  if [[ "$FORCE" != "1" && -s "$outfile" ]]; then
    return 0
  fi

  expr="$(clean_expr "$inst" "$label")"
  if [[ -n "$gti" ]]; then
    expr="($expr)&&gti($gti,TIME)"
  fi

  runlog "clean_${inst}_${label}_$(basename "${outfile%.fits}")" \
    evselect table="${evt}:EVENTS" \
      withfilteredset=yes filteredset="$outfile" \
      destruct=yes keepfilteroutput=yes updateexposure=yes writedss=yes \
      expression="$expr"

  rows="$("$PYTHON" "$SCRIPT_DIR/tools.py" fits-rows "$outfile")"
  if [[ "$rows" -le 0 ]]; then
    echo "No events left after cleaning: $outfile" >&2
    rm -f "$outfile"
    return 0
  fi

  attcalc eventset="$outfile" attitudelabel=ahf refpointlabel=pnt withatthkset=yes atthkset="$ATTHKGEN_FILE"
}

stage_clean() {
  local inst label rawlist evt base gti outfile
  local labels=()

  if ! repro_products_ready; then
    stage_repro
  fi

  if [[ "$FORCE" != "1" ]] && clean_stage_complete; then
    echo "Skipping clean; found existing clean products in $REAL_WORKDIR/clean"
    return 0
  fi

  while IFS= read -r label; do
    [[ -n "$label" ]] && labels+=("$label")
  done < <(clean_band_labels)
  [[ "${#labels[@]}" -gt 0 ]] || { echo "No clean bands configured" >&2; exit 1; }

  set_sas_clean_env
  need_cmd evselect
  need_cmd attcalc
  clean_gti_enabled && need_cmd tabgtigen

  if [[ "$FORCE" == "1" ]]; then
    shopt -s nullglob
    rm -f "$CLEANDIR/.complete" "$CLEANDIR/clean_band_filters.tsv" "$CLEANDIR/manifest/"*_clean_files.txt
    rm -f "$CLEANDIR/events/"*/*_clean.fits "$CLEANDIR/events/"*/*/*_clean.fits
    rm -f "$GTIDIR/"*/*_flare_gti.fits "$LIGHTCURVEDIR/"*/*_flare_lc.fits
    shopt -u nullglob
  fi

  mkdir -p "$CLEANDIR/manifest"
  "$PYTHON" "$SCRIPT_DIR/tools.py" clean-band-table "$CONFIG" > "$CLEANDIR/clean_band_filters.tsv"
  if clean_gti_enabled; then
    mkdir -p "$GTIDIR" "$LIGHTCURVEDIR"
    printf 'inst\tevent\trate_cut_ct_s\ttimebin_s\tpi_min\tpi_max\tlightcurve\tgti\n' > "$GTI_SUMMARY"
  fi

  for inst in $DETECTORS; do
    rawlist="$REAL_REPRODIR/manifest/${inst}_raw.txt"
    manifest_paths_exist "$rawlist" || { echo "Missing or stale raw manifest: $rawlist" >&2; exit 1; }

    while IFS= read -r evt; do
      [[ -n "$evt" ]] || continue
      base="$(basename "$evt")"
      base="${base%.*}"
      gti=""
      if clean_gti_enabled; then
        CLEAN_GTI_RESULT=""
        build_flare_gti "$inst" "$evt" "$base"
        gti="$CLEAN_GTI_RESULT"
      fi

      for label in "${labels[@]}"; do
        mkdir -p "$CLEANDIR/events/$inst/$label"
        outfile="$CLEANDIR/events/$inst/$label/${base}_${label}_clean.fits"
        clean_one_event "$inst" "$evt" "$label" "$outfile" "$gti"
      done
    done < "$rawlist"

    for label in "${labels[@]}"; do
      write_clean_band_manifest "$inst" "$label"
    done
    write_clean_manifest "$inst"
  done

  date -u +%FT%TZ > "$CLEANDIR/.complete"
  echo "Wrote $REAL_CLEANDIR/clean_band_filters.tsv"
  echo "Wrote $REAL_CLEANDIR/manifest"
}

maps_stage_complete() {
  local manifest="$REAL_MAPSDIR/maps_manifest.tsv"
  local sources="$REAL_MAPSDIR/source_manifest.tsv"
  local inst label base event counts exp_vig exp_unvig extra
  local src_inst src_base src_bands source_list
  local have=0
  [[ -s "$REAL_MAPSDIR/.complete" && -s "$manifest" && -s "$sources" ]] || return 1
  while IFS=$'\t' read -r inst label base event counts exp_vig exp_unvig extra; do
    [[ "$inst" == "inst" ]] && continue
    [[ -n "$inst" ]] || continue
    have=1
    for path in "$event" "$counts" "$exp_vig" "$exp_unvig"; do
      [[ -s "$path" ]] || return 1
    done
  done < "$manifest"
  while IFS=$'\t' read -r src_inst src_base src_bands source_list; do
    [[ "$src_inst" == "inst" ]] && continue
    [[ -n "$src_inst" ]] || continue
    [[ -s "$source_list" ]] || return 1
  done < "$sources"
  [[ "$have" == "1" ]]
}

write_maps_manifest_header() {
  mkdir -p "$MAPSDIR"
  printf 'inst\tband\tbase\tevent\tcounts\texposure_vig\texposure_unvig\n' > "$MAPSDIR/maps_manifest.tsv"
}

append_maps_manifest() {
  local inst="$1" label="$2" base="$3"
  shift 3
  local paths=() path
  for path in "$@"; do
    paths+=("$(readlink -f "$path" 2>/dev/null || echo "$path")")
  done
  printf '%s\t%s\t%s' "$inst" "$label" "$base" >> "$MAPSDIR/maps_manifest.tsv"
  for path in "${paths[@]}"; do
    printf '\t%s' "$path" >> "$MAPSDIR/maps_manifest.tsv"
  done
  printf '\n' >> "$MAPSDIR/maps_manifest.tsv"
}

label_in_list() {
  local needle="$1" item
  shift
  for item in "$@"; do
    [[ "$item" == "$needle" ]] && return 0
  done
  return 1
}

source_list_path() {
  local inst="$1" base="$2"
  echo "$MAPSDIR/sources/$inst/$base/${base}_source_list.fits"
}

source_local_list_path() {
  local inst="$1" base="$2"
  echo "$MAPSDIR/sources/$inst/$base/${base}_source_local.fits"
}

source_bkg_path() {
  local inst="$1" base="$2" label="$3"
  echo "$MAPSDIR/sources/$inst/$base/${base}_${label}_source_background.fits"
}

build_map_event_from_clean() {
  local inst="$1" label="$2" base="$3" pimin="$4" pimax="$5" outfile="$6"
  shift 6
  local partdir part source current next tmp idx=0
  local -a parts=()
  [[ "$#" -gt 0 ]] || { echo "No clean event inputs for $inst $label $base" >&2; return 1; }
  mkdir -p "$(dirname "$outfile")"
  partdir="$(dirname "$outfile")/parts"
  mkdir -p "$partdir"

  for source in "$@"; do
    [[ -s "$source" ]] || { echo "Missing clean event file: $source" >&2; exit 1; }
    part="$partdir/${base}_${label}_part${idx}.fits"
    if [[ "$FORCE" == "1" || ! -s "$part" ]]; then
      runlog "maps_${label}_${inst}_${base}_part${idx}" \
        evselect table="${source}:EVENTS" \
          withfilteredset=yes filteredset="$part" \
          destruct=yes keepfilteroutput=yes updateexposure=yes writedss=yes \
          expression="(PI in [$pimin:$pimax])"
    fi
    if [[ "$("$PYTHON" "$SCRIPT_DIR/tools.py" fits-rows "$part")" -gt 0 ]]; then
      parts+=("$part")
    else
      rm -f "$part"
    fi
    idx=$((idx + 1))
  done
  [[ "${#parts[@]}" -gt 0 ]] || { echo "No clean events overlap maps band $label for $inst $base" >&2; return 1; }

  if [[ "${#parts[@]}" -eq 1 ]]; then
    if [[ "$FORCE" == "1" || ! -s "$outfile" ]]; then
      runlog "maps_${label}_${inst}_${base}_events" \
        evselect table="${parts[0]}:EVENTS" \
          withfilteredset=yes filteredset="$outfile" \
          destruct=yes keepfilteroutput=yes updateexposure=yes writedss=yes \
          expression="(PI in [$pimin:$pimax])"
    fi
  else
    current="${parts[0]}"
    idx=1
    while [[ "$idx" -lt "${#parts[@]}" ]]; do
      next="${parts[$idx]}"
      tmp="$partdir/${base}_${label}_merge${idx}.fits"
      if [[ "$FORCE" == "1" || ! -s "$tmp" ]]; then
        runlog "maps_${label}_${inst}_${base}_merge${idx}" \
          merge set1="$current" set2="$next" outset="$tmp"
      fi
      current="$tmp"
      idx=$((idx + 1))
    done
    if [[ "$FORCE" == "1" || ! -s "$outfile" ]]; then
      runlog "maps_${label}_${inst}_${base}_events" \
        evselect table="${current}:EVENTS" \
          withfilteredset=yes filteredset="$outfile" \
          destruct=yes keepfilteroutput=yes updateexposure=yes writedss=yes \
          expression="(PI in [$pimin:$pimax])"
    fi
  fi

  if [[ "$("$PYTHON" "$SCRIPT_DIR/tools.py" fits-rows "$outfile")" -le 0 ]]; then
    echo "No events left in maps event file: $outfile" >&2
    rm -f "$outfile"
    return 1
  fi
  attcalc eventset="$outfile" attitudelabel=ahf refpointlabel=pnt withatthkset=yes atthkset="$ATTHKGEN_FILE"
}

drop_maps_from() {
  local step="$1"
  [[ -z "$step" ]] && return 0
  case "$step" in
    grid|counts|exposure|mask|sources) ;;
    *) echo "Unknown maps rerun step: $step" >&2; exit 2 ;;
  esac

  shopt -s globstar nullglob
  rm -f "$MAPSDIR/.complete" "$MAPSDIR/maps_manifest.tsv" "$MAPSDIR/maps_work.tsv"

  case "$step" in
    grid)
      rm -f "$MAPSDIR/grid.env" "$MAPSDIR/grid_summary.tsv" "$MAPSDIR/maps_band_table.tsv" "$MAPSDIR/source_manifest.tsv"
      rm -f "$MAPSDIR"/events/*/*/*_map_events.fits "$MAPSDIR"/events/*/*/parts/*.fits
      rm -f "$MAPSDIR"/**/*_counts.fits "$MAPSDIR"/**/*_exp_vig.fits "$MAPSDIR"/**/*_exp_unvig.fits
      rm -f "$MAPSDIR"/sources/*/*/*.fits
      ;;
    counts)
      rm -f "$MAPSDIR/maps_band_table.tsv" "$MAPSDIR/source_manifest.tsv"
      rm -f "$MAPSDIR"/events/*/*/*_map_events.fits "$MAPSDIR"/events/*/*/parts/*.fits
      rm -f "$MAPSDIR"/**/*_counts.fits "$MAPSDIR"/**/*_exp_vig.fits "$MAPSDIR"/**/*_exp_unvig.fits
      rm -f "$MAPSDIR"/sources/*/*/*.fits
      ;;
    exposure)
      rm -f "$MAPSDIR/maps_band_table.tsv" "$MAPSDIR/source_manifest.tsv"
      rm -f "$MAPSDIR"/**/*_exp_vig.fits "$MAPSDIR"/**/*_exp_unvig.fits
      rm -f "$MAPSDIR"/sources/*/*/*.fits
      ;;
    mask)
      rm -f "$MAPSDIR/source_manifest.tsv"
      rm -f "$MAPSDIR"/sources/*/*/*.fits
      ;;
    sources)
      rm -f "$MAPSDIR/source_manifest.tsv"
      rm -f "$MAPSDIR"/sources/*/*/*.fits
      ;;
  esac
  shopt -u globstar nullglob
}

stage_maps() {
  local inst label desc pimin pimax clean_label manifest clean_evt base outdir records input_index
  local map_event counts exp_vig exp_unvig hard_sources hard_local_sources
  local source_bkg source_label source_fitmethod
  local have_products=0
  local -a source_bands=() clean_inputs=()

  if ! clean_stage_complete; then
    stage_clean
  fi

  if [[ "$FORCE" != "1" && -z "$MAPS_RERUN_FROM" ]] && maps_stage_complete; then
    echo "Skipping maps; found existing maps products in $REAL_WORKDIR/maps"
    return 0
  fi

  set_sas_clean_env
  need_cmd evselect
  need_cmd attcalc
  need_cmd merge
  need_cmd eexpmap
  need_cmd eboxdetect
  need_cmd esplinemap
  read -r -a source_bands <<< "$MAP_SOURCE_DETECTION_BANDS"

  if [[ "$FORCE" == "1" ]]; then
    drop_maps_from grid
  else
    drop_maps_from "$MAPS_RERUN_FROM"
  fi

  if [[ "$FORCE" == "1" || ! -s "$REAL_MAPSDIR/grid.env" ]]; then
    "$PYTHON" "$SCRIPT_DIR/tools.py" maps-grid "$MAPSDIR" "$MAP_BIN_PHYS" "$MAP_PAD_FRAC" "$REAL_CLEANDIR" "$DETECTORS"
  fi
  # shellcheck source=/dev/null
  source "$REAL_MAPSDIR/grid.env"
  "$PYTHON" "$SCRIPT_DIR/tools.py" maps-band-table "$CONFIG" > "$MAPSDIR/maps_band_table.tsv"
  records="$MAPSDIR/maps_work.tsv"
  printf 'inst\tband\tbase\tpimin\tpimax\tevent\tcounts\texposure_vig\texposure_unvig\n' > "$records"
  write_maps_manifest_header

  while IFS=$'\t' read -r label desc pimin pimax; do
    [[ "$label" == "label" ]] && continue
    [[ -n "$label" ]] || continue
    for inst in $DETECTORS; do
      input_index="$MAPSDIR/events/$label/$inst/${label}_${inst}_clean_inputs.tsv"
      mkdir -p "$(dirname "$input_index")"
      : > "$input_index"

      while IFS= read -r clean_label; do
        [[ -n "$clean_label" ]] || continue
        manifest="$REAL_CLEANDIR/manifest/${inst}_${clean_label}_clean_files.txt"
        clean_manifest_valid "$manifest" || { echo "Missing or stale clean manifest: $manifest" >&2; exit 1; }
        while IFS= read -r clean_evt; do
          [[ -n "$clean_evt" ]] || continue
          [[ -s "$clean_evt" ]] || { echo "Missing clean event file: $clean_evt" >&2; exit 1; }
          base="$(basename "$clean_evt")"
          base="${base%_${clean_label}_clean.fits}"
          printf '%s\t%s\t%s\n' "$base" "$clean_label" "$clean_evt" >> "$input_index"
        done < "$manifest"
      done < <(maps_clean_labels "$label")

      while IFS= read -r base; do
        [[ -n "$base" ]] || continue
        outdir="$MAPSDIR/$label/$inst/$base"
        mkdir -p "$outdir"

        map_event="$MAPSDIR/events/$label/$inst/${base}_${label}_map_events.fits"
        counts="$outdir/${base}_${label}_counts.fits"
        exp_vig="$outdir/${base}_${label}_exp_vig.fits"
        exp_unvig="$outdir/${base}_${label}_exp_unvig.fits"

        if [[ "$FORCE" == "1" || ! -s "$map_event" ]]; then
          clean_inputs=()
          while IFS=$'\t' read -r _base _clean_label clean_evt; do
            [[ "$_base" == "$base" ]] || continue
            clean_inputs+=("$clean_evt")
          done < "$input_index"
          build_map_event_from_clean "$inst" "$label" "$base" "$pimin" "$pimax" "$map_event" "${clean_inputs[@]}"
        fi

        if [[ "$FORCE" == "1" || ! -s "$counts" ]]; then
          runlog "maps_${label}_${inst}_${base}_counts" \
            evselect table="${map_event}:EVENTS" \
              withimageset=yes imageset="$counts" \
              xcolumn=X ycolumn=Y imagebinning=binSize \
              ximagebinsize="$MAP_GRID_BIN_PHYS" yimagebinsize="$MAP_GRID_BIN_PHYS" \
              withxranges=yes ximagemin="$MAP_GRID_X_MIN" ximagemax="$MAP_GRID_X_MAX" \
              withyranges=yes yimagemin="$MAP_GRID_Y_MIN" yimagemax="$MAP_GRID_Y_MAX" \
              ignorelegallimits=yes \
              writedss=yes
        fi

        if [[ "$FORCE" == "1" || ! -s "$exp_vig" ]]; then
          runlog "maps_${label}_${inst}_${base}_eexpmap_vig" \
            eexpmap imageset="$counts" \
              attitudeset="$ATTHKGEN_FILE" \
              eventset="$map_event" \
              expimageset="$exp_vig" \
              pimin="$pimin" pimax="$pimax" \
              withvignetting=yes \
              attrebin="$MAP_EEXPMAP_ATTREBIN"
        fi
        if [[ "$FORCE" == "1" || ! -s "$exp_unvig" ]]; then
          runlog "maps_${label}_${inst}_${base}_eexpmap_unvig" \
            eexpmap imageset="$counts" \
              attitudeset="$ATTHKGEN_FILE" \
              eventset="$map_event" \
              expimageset="$exp_unvig" \
              pimin="$pimin" pimax="$pimax" \
              withvignetting=no \
              attrebin="$MAP_EEXPMAP_ATTREBIN"
        fi
        printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$inst" "$label" "$base" "$pimin" "$pimax" "$(readlink -f "$map_event")" "$(readlink -f "$counts")" "$(readlink -f "$exp_vig")" "$(readlink -f "$exp_unvig")" >> "$records"
        append_maps_manifest "$inst" "$label" "$base" "$map_event" "$counts" "$exp_vig" "$exp_unvig"
        have_products=1
      done < <(cut -f1 "$input_index" | sort -u)
    done
  done < "$MAPSDIR/maps_band_table.tsv"

  [[ "$have_products" == "1" ]] || { echo "No maps products were generated; check clean manifests in $REAL_CLEANDIR/manifest" >&2; exit 1; }

  printf 'inst\tbase\tbands\tsource_list\n' > "$MAPSDIR/source_manifest.tsv"
  tail -n +2 "$records" | cut -f1,3 | sort -u | while IFS=$'\t' read -r inst base; do
    [[ -n "$inst" && -n "$base" ]] || continue
    local -a images=()
    local -a exps=()
    local -a exps2=()
    local -a bkgs=()
    local -a pmins=()
    local -a pmaxs=()
    local idx
    local band_line rec_inst rec_label rec_base rec_pimin rec_pimax rec_evt rec_counts rec_exp_vig rec_exp_unvig
    for source_label in "${source_bands[@]}"; do
      band_line=""
      while IFS=$'\t' read -r rec_inst rec_label rec_base rec_pimin rec_pimax rec_evt rec_counts rec_exp_vig rec_exp_unvig; do
        [[ "$rec_inst" == "inst" ]] && continue
        [[ "$rec_inst" == "$inst" && "$rec_base" == "$base" && "$rec_label" == "$source_label" ]] || continue
        band_line=1
        images+=("$rec_counts")
        exps+=("$rec_exp_vig")
        exps2+=("$rec_exp_unvig")
        pmins+=("$rec_pimin")
        pmaxs+=("$rec_pimax")
        break
      done < "$records"
      [[ -n "$band_line" ]] || { echo "Missing source-detection band $source_label for $inst $base" >&2; exit 1; }
    done
    [[ "${#images[@]}" -gt 0 ]] || continue
    hard_local_sources="$(source_local_list_path "$inst" "$base")"
    hard_sources="$(source_list_path "$inst" "$base")"
    mkdir -p "$(dirname "$hard_sources")"
    if [[ "$FORCE" == "1" || ! -s "$hard_local_sources" ]]; then
      runlog "maps_source_${inst}_${base}_ebox_local" \
        eboxdetect imagesets="${images[*]}" boxlistset="$hard_local_sources" \
          expimagesets="${exps[*]}" \
          withexpimage=yes usemap=no \
          likemin="$MAP_EBOX_LIKEMIN_LOCAL" boxsize="$MAP_EBOX_BOXSIZE" nruns="$MAP_EBOX_NRUNS" \
          pimin="${pmins[*]}" pimax="${pmaxs[*]}"
    fi
    bkgs=()
    for ((idx=0; idx<${#source_bands[@]}; idx++)); do
      source_label="${source_bands[$idx]}"
      source_bkg="$(source_bkg_path "$inst" "$base" "$source_label")"
      bkgs+=("$source_bkg")
      if [[ "$FORCE" == "1" || ! -s "$source_bkg" ]]; then
        source_fitmethod="${MAP_BKG_FITMETHOD:-model}"
        if [[ "$source_fitmethod" == "model" ]]; then
          runlog "maps_source_${inst}_${base}_${source_label}_esplinemap" \
            esplinemap imageset="${images[$idx]}" boxlistset="$hard_local_sources" bkgimageset="$source_bkg" \
              fitmethod="$source_fitmethod" \
              withexpimage=yes expimageset="${exps[$idx]}" \
              withexpimage2=yes expimageset2="${exps2[$idx]}" \
              withdetmask=no withootset=no \
              idband="$((idx + 1))" pimin="${pmins[$idx]}" pimax="${pmaxs[$idx]}" \
              mlmin="$MAP_BKG_MLMIN" scut="$MAP_BKG_SCUT" \
              nsplinenodes="$MAP_BKG_NSPLINENODES" excesssigma="$MAP_BKG_EXCESSSIGMA" nfitrun="$MAP_BKG_NFITRUN" \
              snrmin="$MAP_BKG_SNRMIN" smoothsigma="$MAP_BKG_SMOOTHSIGMA"
        else
          runlog "maps_source_${inst}_${base}_${source_label}_esplinemap" \
            esplinemap imageset="${images[$idx]}" boxlistset="$hard_local_sources" bkgimageset="$source_bkg" \
              fitmethod="$source_fitmethod" \
              withexpimage=yes expimageset="${exps[$idx]}" \
              withexpimage2=no \
              withdetmask=no withootset=no \
              idband="$((idx + 1))" pimin="${pmins[$idx]}" pimax="${pmaxs[$idx]}" \
              mlmin="$MAP_BKG_MLMIN" scut="$MAP_BKG_SCUT" \
              nsplinenodes="$MAP_BKG_NSPLINENODES" excesssigma="$MAP_BKG_EXCESSSIGMA" nfitrun="$MAP_BKG_NFITRUN" \
              snrmin="$MAP_BKG_SNRMIN" smoothsigma="$MAP_BKG_SMOOTHSIGMA"
        fi
      fi
    done
    if [[ "$FORCE" == "1" || ! -s "$hard_sources" ]]; then
      runlog "maps_source_${inst}_${base}_ebox_map" \
        eboxdetect imagesets="${images[*]}" boxlistset="$hard_sources" \
          bkgimagesets="${bkgs[*]}" usemap=yes \
          expimagesets="${exps[*]}" withexpimage=yes \
          likemin="$MAP_EBOX_LIKEMIN_LOCAL" boxsize="$MAP_EBOX_BOXSIZE" nruns="$MAP_EBOX_NRUNS" \
          pimin="${pmins[*]}" pimax="${pmaxs[*]}"
    fi
    printf '%s\t%s\t%s\t%s\n' "$inst" "$base" "${source_bands[*]}" "$(readlink -f "$hard_sources")" >> "$MAPSDIR/source_manifest.tsv"
  done

  date -u +%FT%TZ > "$MAPSDIR/.complete"
  echo "Wrote $REAL_MAPSDIR/maps_manifest.tsv"
  echo "Wrote $REAL_MAPSDIR/source_manifest.tsv"
  echo "Wrote $REAL_MAPSDIR/grid.env"
}

cheese_stage_complete() {
  local manifest="$REAL_CHEESEDIR/cheese_manifest.tsv"
  local header has_detmask=0
  local inst label base source_list global_region base_mask detmask cheesemask cheesed_counts cheesed_exp_vig
  local have=0
  [[ -s "$REAL_CHEESEDIR/.complete" && -s "$manifest" ]] || return 1
  IFS= read -r header < "$manifest"
  [[ "$header" == *$'\tdetmask\t'* ]] && has_detmask=1
  if [[ "$has_detmask" == "1" ]]; then
    while IFS=$'\t' read -r inst label base source_list global_region base_mask detmask cheesemask cheesed_counts cheesed_exp_vig; do
      [[ "$inst" == "inst" ]] && continue
      [[ -n "$inst" ]] || continue
      have=1
      for path in "$source_list" "$global_region" "$base_mask" "$cheesemask" "$cheesed_counts" "$cheesed_exp_vig"; do
        [[ -s "$path" ]] || return 1
      done
    done < "$manifest"
  else
    while IFS=$'\t' read -r inst label base source_list global_region base_mask cheesemask cheesed_counts cheesed_exp_vig; do
      [[ "$inst" == "inst" ]] && continue
      [[ -n "$inst" ]] || continue
      have=1
      for path in "$source_list" "$global_region" "$base_mask" "$cheesemask" "$cheesed_counts" "$cheesed_exp_vig"; do
        [[ -s "$path" ]] || return 1
      done
    done < "$manifest"
  fi
  [[ "$have" == "1" ]]
}

write_cheese_manifest_header() {
  mkdir -p "$CHEESEDIR"
  printf 'inst\tband\tbase\tsource_list\tglobal_region\tbase_mask\tcheese_mask\tcheesed_counts\tcheesed_exposure_vig\n' > "$CHEESEDIR/cheese_manifest.tsv"
}

append_cheese_manifest() {
  local inst="$1" label="$2" base="$3"
  shift 3
  local paths=() path
  for path in "$@"; do
    paths+=("$(readlink -f "$path" 2>/dev/null || echo "$path")")
  done
  printf '%s\t%s\t%s' "$inst" "$label" "$base" >> "$CHEESEDIR/cheese_manifest.tsv"
  for path in "${paths[@]}"; do
    printf '\t%s' "$path" >> "$CHEESEDIR/cheese_manifest.tsv"
  done
  printf '\n' >> "$CHEESEDIR/cheese_manifest.tsv"
}

cheese_global_region_path() {
  local inst="$1" base="$2"
  echo "$CHEESEDIR/regions/$inst/$base/${base}_hard_global_region.fits"
}

cheese_base_mask_path() {
  local inst="$1" base="$2"
  echo "$CHEESEDIR/regions/$inst/$base/${base}_hard_base_mask.fits"
}

drop_cheese_from() {
  local step="$1"
  [[ -z "$step" ]] && return 0
  case "$step" in
    regions|mask|images) ;;
    *) echo "Unknown cheese rerun step: $step" >&2; exit 2 ;;
  esac

  shopt -s globstar nullglob
  rm -f "$CHEESEDIR/.complete" "$CHEESEDIR/cheese_manifest.tsv"

  case "$step" in
    regions)
      rm -f "$CHEESEDIR"/regions/*/*/*_hard_global_region.fits
      rm -f "$CHEESEDIR"/regions/*/*/*_hard_base_mask.fits
      rm -f "$CHEESEDIR"/regions/*/*/*_region_tmp.fits
      rm -f "$CHEESEDIR"/**/*_cheesemask.fits
      rm -f "$CHEESEDIR"/**/*_cheesed_counts.fits
      rm -f "$CHEESEDIR"/**/*_cheesed_exp_vig.fits
      ;;
    mask)
      rm -f "$CHEESEDIR"/**/*_cheesemask.fits
      rm -f "$CHEESEDIR"/**/*_cheesed_counts.fits
      rm -f "$CHEESEDIR"/**/*_cheesed_exp_vig.fits
      ;;
    images)
      rm -f "$CHEESEDIR"/**/*_cheesed_counts.fits
      rm -f "$CHEESEDIR"/**/*_cheesed_exp_vig.fits
      ;;
  esac
  shopt -u globstar nullglob
}

stage_cheese() {
  local inst base bands source_list source_table source_expr hard_event hard_counts global_region base_mask tmp_region
  local label map_event counts exp_vig exp_unvig outdir cheesemask cheesed_counts cheesed_exp_vig
  local maps_extra

  if ! maps_stage_complete; then
    stage_maps
  fi

  if [[ "$FORCE" != "1" && -z "$CHEESE_RERUN_FROM" ]] && cheese_stage_complete; then
    echo "Skipping cheese; found existing cheese products in $REAL_WORKDIR/cheese"
    return 0
  fi

  set_sas_clean_env
  need_cmd region
  need_cmd regionmask

  if [[ "$FORCE" == "1" ]]; then
    drop_cheese_from regions
  else
    drop_cheese_from "$CHEESE_RERUN_FROM"
  fi

  write_cheese_manifest_header

  while IFS=$'\t' read -r inst base bands source_list; do
    [[ "$inst" == "inst" ]] && continue
    [[ -n "$inst" && -n "$base" && -n "$source_list" ]] || continue
    [[ -s "$source_list" ]] || { echo "Missing hard-band source list: $source_list" >&2; exit 1; }

    source_table="$(fits_table_hdu "$source_list" SRCLIST)"
    source_expr=""
    if [[ -n "$MAP_SOURCE_DETECTION_IDBAND" ]] && fits_table_has_column "$source_list" "$source_table" ID_BAND; then
      source_expr="ID_BAND == $MAP_SOURCE_DETECTION_IDBAND"
    fi

    hard_event="$MAPSDIR/events/hard/$inst/${base}_hard_map_events.fits"
    hard_counts="$MAPSDIR/hard/$inst/$base/${base}_hard_counts.fits"
    [[ -s "$hard_event" ]] || { echo "Missing hard-band map event file: $hard_event" >&2; exit 1; }
    [[ -s "$hard_counts" ]] || { echo "Missing hard-band counts image: $hard_counts" >&2; exit 1; }

    global_region="$(cheese_global_region_path "$inst" "$base")"
    base_mask="$(cheese_base_mask_path "$inst" "$base")"
    tmp_region="$(dirname "$global_region")/${base}_region_tmp.fits"
    mkdir -p "$(dirname "$global_region")"

    if [[ "$FORCE" == "1" || ! -s "$global_region" ]]; then
      local -a region_args=(
        eventset="$hard_event"
        tempset="$tmp_region"
        srclisttab="${source_list}:${source_table}"
        operationstyle=global
        outunit=xy
        bkgregionset="$global_region"
      )
      if [[ -n "$source_expr" ]]; then
        region_args+=(expression="$source_expr")
      fi
      runlog "cheese_${inst}_${base}_region" \
        region "${region_args[@]}"
      rm -f "$tmp_region"
    fi

    if [[ "$FORCE" == "1" || ! -s "$base_mask" ]]; then
      runlog "cheese_${inst}_${base}_regionmask" \
        regionmask region="${global_region}:REGION" autosize=no \
          withmaskset=yes maskset="$base_mask" \
          whichpixdef=image pixdefset="$hard_counts"
    fi
  done < "$MAPSDIR/source_manifest.tsv"

  while IFS=$'\t' read -r inst label base map_event counts exp_vig exp_unvig maps_extra; do
    [[ "$inst" == "inst" ]] && continue
    [[ -n "$inst" && -n "$label" && -n "$base" ]] || continue
    source_list="$(source_list_path "$inst" "$base")"
    global_region="$(cheese_global_region_path "$inst" "$base")"
    base_mask="$(cheese_base_mask_path "$inst" "$base")"

    [[ -s "$source_list" ]] || { echo "Missing hard-band source list: $source_list" >&2; exit 1; }
    [[ -s "$global_region" ]] || { echo "Missing hard-band global region file: $global_region" >&2; exit 1; }
    [[ -s "$base_mask" ]] || { echo "Missing hard-band base mask: $base_mask" >&2; exit 1; }

    outdir="$CHEESEDIR/$label/$inst/$base"
    mkdir -p "$outdir"
    cheesemask="$outdir/${base}_${label}_cheesemask.fits"
    cheesed_counts="$outdir/${base}_${label}_cheesed_counts.fits"
    cheesed_exp_vig="$outdir/${base}_${label}_cheesed_exp_vig.fits"

    if [[ "$FORCE" == "1" || ! -s "$cheesemask" ]]; then
      runlog "cheese_${label}_${inst}_${base}_mask" \
        cp -f "$base_mask" "$cheesemask"
    fi
    if [[ "$FORCE" == "1" || ! -s "$cheesed_counts" ]]; then
      runlog "cheese_${label}_${inst}_${base}_counts" \
        "$PYTHON" "$SCRIPT_DIR/tools.py" apply-image-mask "$counts" "$cheesemask" "$cheesed_counts" image
    fi
    if [[ "$FORCE" == "1" || ! -s "$cheesed_exp_vig" ]]; then
      runlog "cheese_${label}_${inst}_${base}_expvig" \
        "$PYTHON" "$SCRIPT_DIR/tools.py" apply-image-mask "$exp_vig" "$cheesemask" "$cheesed_exp_vig" image
    fi

    append_cheese_manifest "$inst" "$label" "$base" "$source_list" "$global_region" "$base_mask" "$cheesemask" "$cheesed_counts" "$cheesed_exp_vig"
  done < "$MAPSDIR/maps_manifest.tsv"

  date -u +%FT%TZ > "$CHEESEDIR/.complete"
  echo "Wrote $REAL_CHEESEDIR/cheese_manifest.tsv"
}

file_type() {
  local name="$1"
  case "$name" in
    *SUM.SAS) echo "SUM.SAS" ;;
    ccf.cif) echo "ccf.cif" ;;
    *.FIT|*.FTZ|*.fits|*.fits.gz) echo "FITS" ;;
    *.ASC) echo "ASC" ;;
    MANIFEST.*) echo "MANIFEST" ;;
    *.ds) echo "ds" ;;
    *.*) echo "${name##*.}" ;;
    *) echo "other" ;;
  esac
}

count_lines() {
  [[ -s "$1" ]] && wc -l < "$1" || echo 0
}

qc_init() {
  local qcdir="$WORKDIR/qc/init" files="$WORKDIR/qc/init/output_files.txt"
  local counts="$WORKDIR/qc/init/file_type_counts.txt"
  [[ -d "$INITDIR" ]] || { echo "Missing init output directory: $INITDIR" >&2; exit 1; }
  mkdir -p "$qcdir"

  {
    [[ -f "$REAL_WORKDIR/sas_setup.env" ]] && echo "$REAL_WORKDIR/sas_setup.env"
    [[ -d "$LOGDIR" ]] && find "$LOGDIR" -maxdepth 1 -type f -name 'init_*.log' | sort
    find "$INITDIR" -maxdepth 1 \( -type f -o -type l \) | sort
  } > "$files"

  while IFS= read -r path; do
    file_type "$(basename "$path")"
  done < "$files" | sort | uniq -c | awk '{ printf "%s %s\n", $2, $1 }' > "$counts"

  echo "Wrote $files"
  echo "Wrote $counts"
}

qc_repro() {
  local qcdir="$WORKDIR/qc/repro" files="$WORKDIR/qc/repro/output_files.txt"
  local counts="$WORKDIR/qc/repro/manifest_counts.txt" status="$WORKDIR/qc/repro/status.txt"
  local inst
  [[ -d "$REPRODIR" ]] || { echo "Missing repro output directory: $REPRODIR" >&2; exit 1; }
  mkdir -p "$qcdir"

  {
    [[ -f "$REAL_WORKDIR/attitude.env" ]] && echo "$REAL_WORKDIR/attitude.env"
    [[ -d "$LOGDIR" ]] && find "$LOGDIR" -maxdepth 1 -type f -name 'repro_*.log' | sort
    find "$REPRODIR" -maxdepth 2 \( -type f -o -type l \) | sort
  } > "$files"

  {
    for inst in $DETECTORS; do
      echo "$inst $(count_lines "$REAL_REPRODIR/manifest/${inst}_raw.txt")"
    done
  } > "$counts"

  {
    [[ -s "$REAL_REPRODIR/atthk.dat" ]] && echo "atthk.dat ok" || echo "atthk.dat missing"
    [[ -s "$REAL_WORKDIR/attitude.env" ]] && echo "attitude.env ok" || echo "attitude.env missing"
    for inst in $DETECTORS; do
      manifest_paths_exist "$REAL_REPRODIR/manifest/${inst}_raw.txt" && echo "${inst}_raw.txt ok" || echo "${inst}_raw.txt missing_or_stale"
    done
  } > "$status"

  rm -f "$qcdir"/soft*_mosaic.png "$qcdir"/hard*_mosaic.png
  "$PYTHON" "$SCRIPT_DIR/tools.py" event-mosaic "$REAL_REPRODIR/manifest" "$qcdir" "$DETECTORS" "$QC_SOFT_HARD_SPLIT_EV"

  echo "Wrote $files"
  echo "Wrote $counts"
  echo "Wrote $status"
  echo "Wrote $qcdir/soft_mosaic.png"
  echo "Wrote $qcdir/hard_mosaic.png"
}

qc_clean() {
  local qcdir="$WORKDIR/qc/clean" files="$WORKDIR/qc/clean/output_files.txt"
  local counts="$WORKDIR/qc/clean/manifest_counts.txt" status="$WORKDIR/qc/clean/status.txt"
  local inst label manifest legacy
  [[ -d "$CLEANDIR" ]] || { echo "Missing clean output directory: $CLEANDIR" >&2; exit 1; }
  mkdir -p "$qcdir"

  {
    find "$CLEANDIR" -maxdepth 5 \( -type f -o -type l \) | sort
    [[ -d "$LOGDIR" ]] && find "$LOGDIR" -maxdepth 1 -type f \( -name 'clean_*.log' -o -name 'lc_*.log' -o -name 'gti_*.log' \) | sort
  } > "$files"

  {
    for inst in $DETECTORS; do
      while IFS= read -r label; do
        [[ -n "$label" ]] || continue
        manifest="$REAL_CLEANDIR/manifest/${inst}_${label}_clean_files.txt"
        echo "$inst $label $(count_lines "$manifest")"
      done < <(clean_band_labels)
      legacy="$(clean_manifest_for_detector "$inst" || true)"
      [[ -n "$legacy" ]] && echo "$inst aggregate $(count_lines "$legacy")"
    done
  } > "$counts"

  {
    if clean_stage_complete; then
      echo "clean_stage complete"
    else
      echo "clean_stage missing_or_stale"
    fi
    for inst in $DETECTORS; do
      while IFS= read -r label; do
        [[ -n "$label" ]] || continue
        manifest="$REAL_CLEANDIR/manifest/${inst}_${label}_clean_files.txt"
        clean_manifest_valid "$manifest" && echo "${inst}_${label}_clean_files.txt ok" || echo "${inst}_${label}_clean_files.txt missing_or_stale"
      done < <(clean_band_labels)
      legacy="$(clean_manifest_for_detector "$inst" || true)"
      [[ -n "$legacy" ]] && echo "$(basename "$legacy") legacy_or_aggregate_ok"
    done
  } > "$status"

  [[ -s "$CLEANDIR/clean_band_filters.tsv" ]] && cp -f "$CLEANDIR/clean_band_filters.tsv" "$qcdir/clean_band_filters.tsv"
  [[ -s "$CLEANDIR/flare_gti_summary.tsv" ]] && cp -f "$CLEANDIR/flare_gti_summary.tsv" "$qcdir/flare_gti_summary.tsv"

  rm -f "$qcdir"/soft*_mosaic.png "$qcdir"/hard*_mosaic.png
  "$PYTHON" "$SCRIPT_DIR/tools.py" event-mosaic "$REAL_CLEANDIR" "$qcdir" "$DETECTORS" "$QC_SOFT_HARD_SPLIT_EV"

  echo "Wrote $files"
  echo "Wrote $counts"
  echo "Wrote $status"
  [[ -s "$qcdir/clean_band_filters.tsv" ]] && echo "Wrote $qcdir/clean_band_filters.tsv"
  [[ -s "$qcdir/flare_gti_summary.tsv" ]] && echo "Wrote $qcdir/flare_gti_summary.tsv"
  echo "Wrote $qcdir/soft_mosaic.png"
  echo "Wrote $qcdir/hard_mosaic.png"
}

qc_maps() {
  local qcdir="$WORKDIR/qc/maps" files="$WORKDIR/qc/maps/output_files.txt"
  local counts="$WORKDIR/qc/maps/file_type_counts.txt" status="$WORKDIR/qc/maps/status.txt"
  [[ -d "$MAPSDIR" ]] || { echo "Missing maps output directory: $MAPSDIR" >&2; exit 1; }
  mkdir -p "$qcdir"

  {
    find "$MAPSDIR" -maxdepth 5 \( -type f -o -type l \) | sort
    [[ -d "$LOGDIR" ]] && find "$LOGDIR" -maxdepth 1 -type f -name 'maps_*.log' | sort
  } > "$files"

  while IFS= read -r path; do
    file_type "$(basename "$path")"
  done < "$files" | sort | uniq -c | awk '{ printf "%s %s\n", $2, $1 }' > "$counts"

  {
    if maps_stage_complete; then
      echo "maps_stage complete"
    else
      echo "maps_stage missing_or_stale"
    fi
    [[ -s "$REAL_MAPSDIR/grid.env" ]] && echo "grid.env ok" || echo "grid.env missing"
    [[ -s "$REAL_MAPSDIR/grid_summary.tsv" ]] && echo "grid_summary.tsv ok" || echo "grid_summary.tsv missing"
    [[ -s "$REAL_MAPSDIR/maps_band_table.tsv" ]] && echo "maps_band_table.tsv ok" || echo "maps_band_table.tsv missing"
    [[ -s "$REAL_MAPSDIR/source_manifest.tsv" ]] && echo "source_manifest.tsv $(($(count_lines "$REAL_MAPSDIR/source_manifest.tsv") - 1))" || echo "source_manifest.tsv missing"
    [[ -s "$REAL_MAPSDIR/maps_manifest.tsv" ]] && echo "maps_manifest.tsv $(($(count_lines "$REAL_MAPSDIR/maps_manifest.tsv") - 1))" || echo "maps_manifest.tsv missing"
  } > "$status"

  [[ -s "$MAPSDIR/grid_summary.tsv" ]] && cp -f "$MAPSDIR/grid_summary.tsv" "$qcdir/grid_summary.tsv"
  [[ -s "$MAPSDIR/maps_band_table.tsv" ]] && cp -f "$MAPSDIR/maps_band_table.tsv" "$qcdir/maps_band_table.tsv"
  [[ -s "$MAPSDIR/source_manifest.tsv" ]] && cp -f "$MAPSDIR/source_manifest.tsv" "$qcdir/source_manifest.tsv"
  [[ -s "$MAPSDIR/maps_manifest.tsv" ]] && cp -f "$MAPSDIR/maps_manifest.tsv" "$qcdir/maps_manifest.tsv"
  if [[ -s "$REAL_MAPSDIR/maps_manifest.tsv" ]]; then
    "$PYTHON" "$SCRIPT_DIR/tools.py" maps-qc "$REAL_MAPSDIR/maps_manifest.tsv" "$qcdir"
  fi

  echo "Wrote $files"
  echo "Wrote $counts"
  echo "Wrote $status"
  [[ -s "$qcdir/grid_summary.tsv" ]] && echo "Wrote $qcdir/grid_summary.tsv"
  [[ -s "$qcdir/maps_band_table.tsv" ]] && echo "Wrote $qcdir/maps_band_table.tsv"
  [[ -s "$qcdir/source_manifest.tsv" ]] && echo "Wrote $qcdir/source_manifest.tsv"
  [[ -s "$qcdir/maps_manifest.tsv" ]] && echo "Wrote $qcdir/maps_manifest.tsv"
  [[ -s "$qcdir/maps_qc_summary.tsv" ]] && echo "Wrote $qcdir/*_mosaic.png"
}

qc_cheese() {
  local qcdir="$WORKDIR/qc/cheese" files="$WORKDIR/qc/cheese/output_files.txt"
  local counts="$WORKDIR/qc/cheese/file_type_counts.txt" status="$WORKDIR/qc/cheese/status.txt"
  [[ -d "$CHEESEDIR" ]] || { echo "Missing cheese output directory: $CHEESEDIR" >&2; exit 1; }
  mkdir -p "$qcdir"

  {
    find "$CHEESEDIR" -maxdepth 6 \( -type f -o -type l \) | sort
    [[ -d "$LOGDIR" ]] && find "$LOGDIR" -maxdepth 1 -type f -name 'cheese_*.log' | sort
  } > "$files"

  while IFS= read -r path; do
    file_type "$(basename "$path")"
  done < "$files" | sort | uniq -c | awk '{ printf "%s %s\n", $2, $1 }' > "$counts"

  {
    if cheese_stage_complete; then
      echo "cheese_stage complete"
    else
      echo "cheese_stage missing_or_stale"
    fi
    [[ -s "$REAL_CHEESEDIR/cheese_manifest.tsv" ]] && echo "cheese_manifest.tsv $(($(count_lines "$REAL_CHEESEDIR/cheese_manifest.tsv") - 1))" || echo "cheese_manifest.tsv missing"
  } > "$status"

  [[ -s "$CHEESEDIR/cheese_manifest.tsv" ]] && cp -f "$CHEESEDIR/cheese_manifest.tsv" "$qcdir/cheese_manifest.tsv"
  if [[ -s "$REAL_CHEESEDIR/cheese_manifest.tsv" ]]; then
    "$PYTHON" "$SCRIPT_DIR/tools.py" cheese-qc "$REAL_CHEESEDIR/cheese_manifest.tsv" "$qcdir"
  fi

  echo "Wrote $files"
  echo "Wrote $counts"
  echo "Wrote $status"
  [[ -s "$qcdir/cheese_manifest.tsv" ]] && echo "Wrote $qcdir/cheese_manifest.tsv"
  [[ -s "$qcdir/cheese_qc_summary.tsv" ]] && echo "Wrote $qcdir/*_mosaic.png"
}

run_qc_stage() {
  case "$1" in
    init) qc_init ;;
    repro) qc_repro ;;
    clean) qc_clean ;;
    maps) qc_maps ;;
    cheese) qc_cheese ;;
    all|qc) qc_init; qc_repro; qc_clean; qc_maps; qc_cheese ;;
    *) echo "No QC stage for: $1" >&2; exit 2 ;;
  esac
}

maybe_run_qc() {
  [[ "$RUN_QC" == "1" ]] || return 0
  run_qc_stage "$STAGE"
}

case "$STAGE" in
  init) stage_init; maybe_run_qc ;;
  repro) stage_repro; maybe_run_qc ;;
  clean) stage_clean; maybe_run_qc ;;
  maps) stage_maps; maybe_run_qc ;;
  cheese) stage_cheese; maybe_run_qc ;;
  all) stage_init; stage_repro; stage_clean; stage_maps; stage_cheese; maybe_run_qc ;;
  qc-init) qc_init ;;
  qc-repro) qc_repro ;;
  qc-clean) qc_clean ;;
  qc-maps) qc_maps ;;
  qc-cheese) qc_cheese ;;
  qc) run_qc_stage all ;;
  env) print_env ;;
  *) echo "Unknown stage: $STAGE" >&2; exit 2 ;;
esac
