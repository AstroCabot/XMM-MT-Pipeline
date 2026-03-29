#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 2 ]]; then
  echo "Usage: $0 config.env <stage>" >&2
  exit 1
fi

CONFIG_FILE="$1"
STAGE="$2"
source "$CONFIG_FILE"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
HELPER_DIR="$SCRIPT_DIR/scripts"

: "${WORKDIR:?Set WORKDIR in config.env}"
: "${ODFDIR:?Set ODFDIR in config.env}"

: "${TARGET_ID:=}"
: "${TARGET_ID_TYPE:=smallbody}"
: "${TRACK_INPUT:=}"
: "${FORCE:=0}"
: "${ONLY_INSTS:=PN M1 M2}"

: "${CLEAN_PI_MIN:=200}"
: "${CLEAN_PI_MAX:=12000}"
: "${CLEAN_SKIP_ESPFILT:=no}"
: "${ESPFILT_RANGESCALE_PN:=25.0}"
: "${ESPFILT_RANGESCALE_MOS:=20.0}"
: "${ESPFILT_ALLOWSIGMA:=3.0}"

: "${SRCDET_PI_MIN:=1000}"
: "${SRCDET_PI_MAX:=12000}"
: "${DETECT_MIN_ML:=10.0}"
: "${DETECT_EXCLUDE_TRACK_ARCSEC:=30.0}"
: "${DETECT_DS9_RADIUS_ARCSEC:=20.0}"
: "${DETECT_QC_BANDS:=soft:200:1000;hard:1000:12000}"
: "${DETECT_QC_BIN_ARCSEC:=6.0}"
: "${DETECT_QC_PAD_ARCMIN:=2.0}"

: "${TRACK_STEP:=30s}"
: "${MAX_SEG_BLUR_ARCSEC:=1.0}"
: "${SEGMENT_SNAP_S:=0.0}"
: "${MIN_SEG_S:=0.0}"
: "${COMET_REF_RA:=}"
: "${COMET_REF_DEC:=}"
: "${AHF_INPUT:=}"
: "${PREMERGE_IMAGE_SIZE_DEG:=1.5}"
: "${ATTMOVE_GRANULARITY:=5}"
: "${ATTMOVE_MINSTABLE:=30}"
: "${ATTCALC_IMAGE_SIZE_DEG:=0.6}"

: "${SRC_DX_ARCSEC:=0.0}"
: "${SRC_DY_ARCSEC:=0.0}"
: "${SRC_R_ARCSEC:=500.0}"
: "${BKG_MODE:=annulus}"
: "${BKG_DX_ARCSEC:=0.0}"
: "${BKG_DY_ARCSEC:=0.0}"
: "${BKG_R_ARCSEC:=550.0}"
: "${BKG_RIN_ARCSEC:=500.0}"
: "${BKG_ROUT_ARCSEC:=550.0}"
: "${FIELD_SOURCE_MASK_R_ARCSEC:=20.0}"
: "${CONTAM_SAMPLE_STEP_S:=5.0}"
: "${CONTAM_MIN_GOOD_SPAN_S:=0.0}"
: "${SCIENCE_GTI_POLICY:=full}"
: "${SCIENCE_APPLY_TRAIL_MASK:=yes}"

: "${IMAGE_BANDS:=soft:200:1000;broad:300:2000;hard:1000:12000}"
: "${IMAGE_RADIUS_ARCSEC:=900.0}"
: "${IMAGE_BIN_ARCSEC:=2.0}"
: "${IMAGE_MASK_RADIUS_ARCSEC:=$FIELD_SOURCE_MASK_R_ARCSEC}"
: "${SCIENCE_EEXPMAP_ATTREBIN:=4.0}"
: "${EEXPMAP_ATTREBIN:=$SCIENCE_EEXPMAP_ATTREBIN}"

: "${LC_PI_MIN:=300}"
: "${LC_PI_MAX:=2000}"
: "${LC_BIN_S:=500.0}"
: "${LC_DETXBINS:=5}"
: "${LC_DETYBINS:=5}"
: "${LC_APPLY_ABSOLUTE_CORRECTIONS:=no}"
: "${LC_ALLOW_NOABS_FALLBACK:=yes}"

: "${SPEC_PI_MIN:=300}"
: "${SPEC_PI_MAX:=5000}"
: "${SPEC_DETXBINS:=5}"
: "${SPEC_DETYBINS:=5}"
: "${SPEC_MIN_SAFE_DETBINS:=15}"
: "${SPEC_USE_ODF_ATT:=yes}"
: "${SPECGROUP_MINCOUNTS:=25}"
: "${SPECGROUP_OVERSAMPLE:=3.0}"

: "${SHORTLINK_NAME:=}"
: "${SHORTLINK_KEEP:=no}"

mkdir -p "$WORKDIR"

REAL_WORKDIR="$WORKDIR"
if [[ -n "$SHORTLINK_NAME" ]]; then
  _SHORTLINK="/tmp/_xmm_${SHORTLINK_NAME}"
else
  _SHORTLINK="/tmp/_xmmpipe_2_$$"
fi
if [[ ${#WORKDIR} -gt 60 ]]; then
  rm -f "$_SHORTLINK"
  ln -sf "$WORKDIR" "$_SHORTLINK"
  WORKDIR="$_SHORTLINK"
  _keep="${SHORTLINK_KEEP,,}"
  if [[ "$_keep" != "1" && "$_keep" != "y" && "$_keep" != "yes" && "$_keep" != "true" ]]; then
    trap 'rm -f "$_SHORTLINK"' EXIT
  fi
fi

LOGDIR="$WORKDIR/logs"
mkdir -p "$LOGDIR"

runlog() {
  local name="$1"
  shift
  echo "[$(date -u +%FT%TZ)] $name" | tee -a "$LOGDIR/pipeline.log"
  "$@" 2>&1 | tee "$LOGDIR/${name}.log"
}

need_cmd() {
  command -v "$1" >/dev/null 2>&1 || { echo "Missing required command: $1" >&2; exit 1; }
}

need_file() {
  [[ -f "$1" ]] || { echo "Missing required file: $1" >&2; exit 1; }
}

source_if_exists() {
  local f="$1"
  [[ -f "$f" ]] && source "$f"
}

bool_yes() {
  local v="${1,,}"
  [[ "$v" == "1" || "$v" == "y" || "$v" == "yes" || "$v" == "true" ]]
}

skip_file() {
  local f="$1"
  [[ "$FORCE" != "1" && -s "$f" ]]
}

SKIPLIST="$WORKDIR/clean/skipped_exposures.txt"

is_skipped() {
  local key="$1"
  [[ "$FORCE" != "1" && -f "$SKIPLIST" ]] && grep -qxF "$key" "$SKIPLIST"
}

mark_skipped() {
  local key="$1" reason="$2"
  mkdir -p "$(dirname "$SKIPLIST")"
  if ! grep -qxF "$key" "$SKIPLIST" 2>/dev/null; then
    echo "$key" >> "$SKIPLIST"
  fi
  echo "SKIP [$reason]: $key"
}

fits_rows() {
  python3 - "$1" <<'PY'
import os, sys
from astropy.io import fits
path = sys.argv[1]
if not os.path.exists(path):
    print(-1)
    raise SystemExit(0)
with fits.open(path) as hdul:
    if 'EVENTS' in hdul:
        data = hdul['EVENTS'].data
        print(len(data) if data is not None else 0)
    else:
        for hdu in hdul[1:]:
            if isinstance(hdu, fits.BinTableHDU) and hdu.data is not None:
                print(len(hdu.data))
                break
        else:
            print(-1)
PY
}

fits_filter() {
  python3 - "$1" <<'PY'
import sys
from astropy.io import fits
with fits.open(sys.argv[1]) as hdul:
    for name in ('EVENTS', 0):
        try:
            val = hdul[name].header.get('FILTER', '')
            if val:
                print(val)
                raise SystemExit(0)
        except (KeyError, IndexError):
            pass
    print('')
PY
}

srclist_rows() {
  python3 - "$1" <<'PY'
import os, sys
from astropy.io import fits
path = sys.argv[1]
if not os.path.exists(path):
    print(-1)
    raise SystemExit(0)
with fits.open(path) as hdul:
    for hdu in hdul[1:]:
        if isinstance(hdu, fits.BinTableHDU):
            print(len(hdu.data) if hdu.data is not None else 0)
            break
    else:
        print(-1)
PY
}

setup_sas_env_for_init() {
  export SAS_ODF="$ODFDIR"
  export SAS_CCF="$WORKDIR/init/ccf.cif"
}

setup_sas_env() {
  need_file "$WORKDIR/sas_setup.env"
  source "$WORKDIR/sas_setup.env"
  export SAS_CCF="$WORKDIR/init/ccf.cif"
  local sumsas
  sumsas="$(find "$WORKDIR/init" -maxdepth 1 -type f -name '*SUM.SAS' | sort | head -n 1 || true)"
  [[ -n "$sumsas" ]] && export SAS_ODF="$sumsas"
  if [[ -n "$sumsas" && -f "$sumsas" ]]; then
    sed -i "s|^PATH .*|PATH $REAL_WORKDIR/init/|" "$sumsas"
  fi
}

load_attitude() {
  need_file "$WORKDIR/attitude.env"
  source "$WORKDIR/attitude.env"
  ATTHKGEN_FILE="$WORKDIR/repro/atthk.dat"
  ATTHK_FILE="$WORKDIR/repro/atthk.dat"
}

load_ref_attitude() {
  need_file "$WORKDIR/track/ref_attitude.env"
  source "$WORKDIR/track/ref_attitude.env"
}

first_available() {
  local f
  for f in "$@"; do
    if [[ -n "$f" && -f "$f" ]]; then
      echo "$f"
      return 0
    fi
  done
  echo ""
}

find_and_write_list() {
  local outfile="$1"
  local root="$2"
  shift 2
  python3 - "$outfile" "$root" "$@" <<'PY'
import sys
from pathlib import Path

out = Path(sys.argv[1])
root = Path(sys.argv[2])
patterns = sys.argv[3:]
seen = set()
rows = []
for patt in patterns:
    for p in sorted(root.rglob(patt)):
        if p.is_file():
            s = str(p.resolve())
            if s not in seen:
                seen.add(s)
                rows.append(s)
out.parent.mkdir(parents=True, exist_ok=True)
out.write_text("\n".join(rows) + ("\n" if rows else ""), encoding="utf-8")
print(len(rows))
PY
}

choose_original_ahf() {
  python3 - "$ODFDIR" <<'PY'
import sys
from pathlib import Path

root = Path(sys.argv[1])
patterns = [
    '*SC*ATS.FIT*', '*SC*ATS.FT*', '*ATS.FIT*', '*ATS.FT*',
    '*ATTTSR*.FIT*', '*ATTTSR*.FT*', '*AHF*.FIT*', '*AHF*.FT*',
]
cands = []
seen = set()
for patt in patterns:
    for p in sorted(root.rglob(patt)):
        if p.is_file():
            s = str(p.resolve())
            if s not in seen:
                seen.add(s)
                cands.append(s)
if not cands:
    raise SystemExit(1)
print(cands[0])
PY
}

write_list_from_clean_dir() {
  local inst="$1"
  local outfile="$2"
  python3 - "$outfile" "$WORKDIR/clean/$inst" <<'PY'
import sys
from pathlib import Path
out = Path(sys.argv[1])
root = Path(sys.argv[2])
rows = [str(p.resolve()) for p in sorted(root.glob('*.clean.fits')) if p.is_file()]
out.parent.mkdir(parents=True, exist_ok=True)
out.write_text("\n".join(rows) + ("\n" if rows else ""), encoding='utf-8')
print(len(rows))
PY
}

merge_pairwise_generic() {
  local outfile="$1"
  local imagesize="$2"
  shift 2
  local -a files=("$@")
  [[ ${#files[@]} -ge 1 ]] || return 1
  if [[ ${#files[@]} -eq 1 ]]; then
    cp -f "${files[0]}" "$outfile"
    return 0
  fi
  local tmp1="${files[0]}"
  local tmp2
  local i=1
  local -a tmp_created=()
  while [[ $i -lt ${#files[@]} ]]; do
    tmp2="${files[$i]}"
    local outtmp="$outfile.tmp${i}.fits"
    merge set1="$tmp1" set2="$tmp2" outset="$outtmp" imagesize="$imagesize"
    tmp_created+=("$outtmp")
    tmp1="$outtmp"
    i=$((i+1))
  done
  mv -f "$tmp1" "$outfile"
  local f
  for f in "${tmp_created[@]}"; do
    [[ "$f" == "$outfile" ]] || rm -f "$f"
  done
}

merge_pairwise_cometref() {
  local outfile="$1"
  shift
  local -a files=("$@")
  [[ ${#files[@]} -ge 1 ]] || return 1
  if [[ ${#files[@]} -eq 1 ]]; then
    cp -f "${files[0]}" "$outfile"
    return 0
  fi
  local tmp1="${files[0]}"
  local tmp2
  local i=1
  local -a tmp_created=()
  while [[ $i -lt ${#files[@]} ]]; do
    tmp2="${files[$i]}"
    local outtmp="$outfile.tmp${i}.fits"
    merge set1="$tmp1" set2="$tmp2" outset="$outtmp" \
      withradec=Y ra="$COMET_REF_RA" dec="$COMET_REF_DEC" \
      imagesize="$ATTCALC_IMAGE_SIZE_DEG" columns="INSTID"
    tmp_created+=("$outtmp")
    tmp1="$outtmp"
    i=$((i+1))
  done
  mv -f "$tmp1" "$outfile"
  local f
  for f in "${tmp_created[@]}"; do
    [[ "$f" == "$outfile" ]] || rm -f "$f"
  done
}

selected_insts() {
  echo "$ONLY_INSTS"
}

clean_expr_for_inst() {
  local inst="$1"
  if [[ "$inst" == "PN" ]]; then
    echo "(FLAG==0)&&(PATTERN<=4)&&(PI in [$CLEAN_PI_MIN:$CLEAN_PI_MAX])"
  else
    echo "#XMMEA_EM&&(PATTERN<=12)&&(PI in [$CLEAN_PI_MIN:$CLEAN_PI_MAX])"
  fi
}

merged_clean_evt() {
  local inst="$1"
  echo "$WORKDIR/clean/${inst}_clean_merged.fits"
}

flare_gti_for_inst() {
  local inst="$1"
  echo "$WORKDIR/clean/${inst}_flare_gti.fits"
}

merge_flare_gtis() {
  local out="$1"
  shift
  python3 - "$out" "$@" <<'PYGTI'
import sys
from astropy.io import fits
from astropy.table import Table, vstack
import numpy as np

out_path = sys.argv[1]
gti_files = sys.argv[2:]
tables = []
for f in gti_files:
    with fits.open(f) as hdul:
        for ext in hdul[1:]:
            if hasattr(ext, 'columns') and 'START' in [c.name for c in ext.columns]:
                tables.append(Table(ext.data))
                break
if not tables:
    print("WARNING: no GTI rows found in input files", file=sys.stderr)
    sys.exit(1)
merged = vstack(tables)
max_valid = 1e12  # ~31,700 years in seconds; XMM times are ~8e8
keep = (merged['START'] > 0) & (merged['STOP'] < max_valid)
nrej = len(merged) - int(np.sum(keep))
if nrej > 0:
    print(f"Rejected {nrej} GTI rows with absurd bounds", file=sys.stderr)
merged = merged[keep]
if len(merged) == 0:
    print("WARNING: no valid GTI rows remain after sanitization", file=sys.stderr)
    sys.exit(1)
merged.sort('START')
merged.write(out_path, format='fits', overwrite=True)
print(f"Merged {len(merged)} GTI intervals from {len(gti_files)} files into {out_path}")
PYGTI
}

inst_comet_evt() {
  local inst="$1"
  echo "$WORKDIR/comet/${inst}_comet.fits"
}

inst_analysis_evt() {
  local inst="$1"
  local masked="$WORKDIR/contam/${inst}_science_base.fits"
  if [[ -f "$masked" ]]; then
    echo "$masked"
  else
    inst_comet_evt "$inst"
  fi
}

image_bands_lines() {
  local entry
  IFS=';' read -ra __bands <<< "$IMAGE_BANDS"
  for entry in "${__bands[@]}"; do
    [[ -n "$entry" ]] || continue
    IFS=':' read -r label pimin pimax <<< "$entry"
    if [[ -z "$label" || -z "$pimin" || -z "$pimax" ]]; then
      echo "Bad IMAGE_BANDS entry: $entry" >&2
      exit 1
    fi
    echo "$label $pimin $pimax"
  done
}

detect_qc_bands_lines() {
  local entry
  IFS=';' read -ra __bands <<< "$DETECT_QC_BANDS"
  for entry in "${__bands[@]}"; do
    [[ -n "$entry" ]] || continue
    IFS=':' read -r label pimin pimax <<< "$entry"
    if [[ -z "$label" || -z "$pimin" || -z "$pimax" ]]; then
      echo "Bad DETECT_QC_BANDS entry: $entry" >&2
      exit 1
    fi
    echo "$label $pimin $pimax"
  done
}

srcdet_bands_lines() {
  local entry
  if [[ -n "$SRCDET_BANDS" ]]; then
    IFS=';' read -ra __bands <<< "$SRCDET_BANDS"
    for entry in "${__bands[@]}"; do
      [[ -n "$entry" ]] || continue
      local pimin pimax
      IFS=':' read -r pimin pimax <<< "$entry"
      if [[ -z "$pimin" || -z "$pimax" ]]; then
        echo "Bad SRCDET_BANDS entry: $entry" >&2
        exit 1
      fi
      echo "$pimin $pimax"
    done
  else
    echo "$SRCDET_PI_MIN $SRCDET_PI_MAX"
  fi
}

build_region_env() {
  local eventfile="$1"
  local prefix="$2"
  local shape="$3"
  shift 3
  local outf="$WORKDIR/tmp_${prefix}_region.env"
  python3 "$HELPER_DIR/xmm_event_geometry.py" region --event "$eventfile" --shape "$shape" "$@" > "$outf"
  source "$outf"
  eval "${prefix}_EXPR=\$REGION_EXPR"
  eval "${prefix}_X=\$CENTER_X"
  eval "${prefix}_Y=\$CENTER_Y"
  rm -f "$outf"
}

make_image_grid_if_needed() {
  local refevt
  refevt="$(first_available "$(inst_comet_evt PN)" "$(inst_comet_evt M1)" "$(inst_comet_evt M2)")"
  [[ -n "$refevt" ]] || { echo "No comet-frame event list found for imaging" >&2; exit 1; }
  mkdir -p "$WORKDIR/images"
  if [[ "$FORCE" == "1" || ! -s "$WORKDIR/images/grid.json" ]]; then
    python3 "$HELPER_DIR/xmm_event_geometry.py" image-grid \
      --event "$refevt" \
      --radius-arcsec "$IMAGE_RADIUS_ARCSEC" \
      --bin-arcsec "$IMAGE_BIN_ARCSEC" \
      --center-mode median \
      --out-json "$WORKDIR/images/grid.json" > "$WORKDIR/images/grid.env"
  fi
  source "$WORKDIR/images/grid.env"
}

stage_init() {
  need_cmd cifbuild
  need_cmd odfingest

  if [[ "$FORCE" != "1" && -f "$WORKDIR/sas_setup.env" ]]; then
    echo "Skipping init; found $WORKDIR/sas_setup.env"
    return 0
  fi

  mkdir -p "$WORKDIR/init"
  pushd "$WORKDIR/init" >/dev/null
  setup_sas_env_for_init
  runlog init_cifbuild cifbuild
  runlog init_odfingest odfingest odfdir="$ODFDIR" outdir="$WORKDIR/init"
  local sumsas
  sumsas="$(find "$WORKDIR/init" -maxdepth 1 -type f -name '*SUM.SAS' | sort | head -n 1 || true)"
  [[ -n "$sumsas" ]] || { echo "Could not find *SUM.SAS after odfingest" >&2; exit 1; }
  cat > "$WORKDIR/sas_setup.env" <<EOF
export SAS_CCF="${REAL_WORKDIR}/init/ccf.cif"
export SAS_ODF="${sumsas/$WORKDIR/$REAL_WORKDIR}"
EOF

  local odf_file
  for odf_file in "$ODFDIR"/*; do
    [[ -f "$odf_file" ]] || continue
    case "$(basename "$odf_file")" in
      *.gz.gpg|*.tar.gz|*.tar) continue ;;
    esac
    local odf_base
    odf_base="$(basename "$odf_file")"
    [[ -e "$WORKDIR/init/$odf_base" ]] || ln -sf "$odf_file" "$WORKDIR/init/$odf_base"
  done
  sed -i "s|^PATH .*|PATH $REAL_WORKDIR/init/|" "$sumsas"

  popd >/dev/null
}

stage_repro() {
  need_cmd epproc
  need_cmd emproc
  need_cmd atthkgen
  setup_sas_env

  mkdir -p "$WORKDIR/repro" "$WORKDIR/repro/manifest"
  pushd "$WORKDIR/repro" >/dev/null

  if [[ "$FORCE" == "1" || ! -s "$WORKDIR/repro/manifest/PN_raw.txt" || ! -s "$WORKDIR/repro/manifest/M1_raw.txt" || ! -s "$WORKDIR/repro/manifest/M2_raw.txt" ]]; then
    runlog repro_epproc epproc
    runlog repro_emproc emproc
  fi

  find_and_write_list "$WORKDIR/repro/manifest/PN_raw.txt" "$WORKDIR/repro" '*EPN*ImagingEvts.ds' '*EPN*ImagingEvts*.FIT*' '*PIEVLI*.FIT*' '*EPN*EVLI*.FIT*' >/dev/null
  find_and_write_list "$WORKDIR/repro/manifest/M1_raw.txt" "$WORKDIR/repro" '*EMOS1*ImagingEvts.ds' '*EMOS1*ImagingEvts*.FIT*' '*M1EVLI*.FIT*' '*EMOS1*EVLI*.FIT*' >/dev/null
  find_and_write_list "$WORKDIR/repro/manifest/M2_raw.txt" "$WORKDIR/repro" '*EMOS2*ImagingEvts.ds' '*EMOS2*ImagingEvts*.FIT*' '*M2EVLI*.FIT*' '*EMOS2*EVLI*.FIT*' >/dev/null

  if [[ ! -s "$WORKDIR/repro/manifest/PN_raw.txt" && ! -s "$WORKDIR/repro/manifest/M1_raw.txt" && ! -s "$WORKDIR/repro/manifest/M2_raw.txt" ]]; then
    echo "No EPIC imaging event lists found after repro" >&2
    exit 1
  fi

  local atthk="$WORKDIR/repro/atthk.dat"
  if [[ "$FORCE" == "1" || ! -s "$atthk" ]]; then
    runlog repro_atthkgen atthkgen atthkset="$atthk"
  fi

  local ahf="$AHF_INPUT"
  if [[ -z "$ahf" ]]; then
    if ! ahf="$(choose_original_ahf)"; then
      echo "Could not locate an original spacecraft attitude history file (*SC*ATS.FIT*) in $ODFDIR" >&2
      exit 1
    fi
  fi
  [[ -f "$ahf" ]] || { echo "AHF/ATS file not found: $ahf" >&2; exit 1; }

  cat > "$WORKDIR/attitude.env" <<EOF
ATTHKGEN_FILE="${atthk/$WORKDIR/$REAL_WORKDIR}"
ATTHK_FILE="${atthk/$WORKDIR/$REAL_WORKDIR}"
AHF_FILE="$ahf"
EOF
  popd >/dev/null
}

clean_one_event() {
  local inst="$1"
  local evt="$2"
  local outfile="$3"

  local skipkey="${inst}:$(basename "$evt")"
  if is_skipped "$skipkey"; then
    return 0
  fi

  local filt
  filt="$(fits_filter "$evt")"
  if [[ "${filt,,}" == "closed" ]]; then
    mark_skipped "$skipkey" "Closed filter"
    return 0
  fi

  local base
  base="$(basename "$outfile")"
  base="${base%.fits}"
  local outdir="$WORKDIR/clean/${inst}/${base}.espfilt"
  mkdir -p "$WORKDIR/clean/${inst}" "$outdir"

  if skip_file "$outfile" && [[ "$(fits_rows "$outfile")" -gt 0 ]]; then
    return 0
  fi

  pushd "$outdir" >/dev/null
  rm -f input.fits
  ln -sf "$evt" input.fits

  local espfilt_ok=1
  if bool_yes "$CLEAN_SKIP_ESPFILT"; then
    espfilt_ok=0
    echo "NOTE: espfilt skipped for ${inst}/${base} (CLEAN_SKIP_ESPFILT=yes)"
  elif [[ "$inst" == "PN" ]]; then
    runlog "clean_espfilt_${inst}_${base}" espfilt eventfile=input.fits method=histogram withoot=no rangescale="$ESPFILT_RANGESCALE_PN" allowsigma="$ESPFILT_ALLOWSIGMA" || espfilt_ok=0
  else
    runlog "clean_espfilt_${inst}_${base}" espfilt eventfile=input.fits method=histogram withoot=no rangescale="$ESPFILT_RANGESCALE_MOS" allowsigma="$ESPFILT_ALLOWSIGMA" || espfilt_ok=0
  fi

  if [[ "$espfilt_ok" -eq 1 ]]; then
    local espfilt_gti
    espfilt_gti="$(ls -1 *gti*fits *GTI*fits 2>/dev/null | head -n 1 || true)"
    if [[ -n "$espfilt_gti" ]]; then
      cp -f "$espfilt_gti" "$WORKDIR/clean/${inst}/${base}.flare_gti.fits"
      local _gti_valid
      _gti_valid="$(python3 -c "
import sys
from astropy.io import fits
with fits.open(sys.argv[1]) as h:
    for e in h[1:]:
        if hasattr(e,'columns') and 'START' in [c.name for c in e.columns]:
            ok = sum(1 for s,t in zip(e.data['START'],e.data['STOP'])
                     if float(s) > 0 and float(t) < 1e12)
            print(ok)
            break
" "$WORKDIR/clean/${inst}/${base}.flare_gti.fits" 2>/dev/null || echo "0")"
      if [[ "${_gti_valid:-0}" -eq 0 ]]; then
        echo "WARNING: espfilt produced dummy GTI for ${inst}/${base} (all rows absurd); removing"
        rm -f "$WORKDIR/clean/${inst}/${base}.flare_gti.fits"
      fi
    fi
  else
    echo "WARNING: espfilt failed for ${inst}/${base}; no flare GTI for this exposure"
  fi
  popd >/dev/null

  local expr
  expr="$(clean_expr_for_inst "$inst")"

  evselect table="${evt}:EVENTS" \
    withfilteredset=yes filteredset="$outfile" \
    destruct=yes keepfilteroutput=yes updateexposure=yes writedss=yes \
    expression="$expr"

  if [[ "$(fits_rows "$outfile")" -le 0 ]]; then
    mark_skipped "$skipkey" "empty after cleaning"
    rm -f "$outfile"
    return 0
  fi

  attcalc eventset="$outfile" \
    attitudelabel=ahf \
    refpointlabel=pnt \
    withatthkset=yes atthkset="$ATTHKGEN_FILE"
}

stage_clean() {
  need_cmd espfilt
  need_cmd evselect
  need_cmd attcalc
  need_cmd merge
  setup_sas_env
  load_attitude

  mkdir -p "$WORKDIR/clean"

  local inst rawlist evt base outfile merged listfile
  for inst in $(selected_insts); do
    rawlist="$WORKDIR/repro/manifest/${inst}_raw.txt"
    [[ -f "$rawlist" ]] || continue
    mkdir -p "$WORKDIR/clean/$inst"

    while IFS= read -r evt; do
      [[ -n "$evt" ]] || continue
      base="$(basename "$evt")"
      base="${base%.*}"
      outfile="$WORKDIR/clean/$inst/${base}.clean.fits"
      clean_one_event "$inst" "$evt" "$outfile"
    done < "$rawlist"

    listfile="$WORKDIR/clean/${inst}_clean_files.txt"
    write_list_from_clean_dir "$inst" "$listfile" >/dev/null

    mapfile -t clean_files < "$listfile"
    [[ ${#clean_files[@]} -gt 0 ]] || continue

    merged="$(merged_clean_evt "$inst")"
    if [[ "$FORCE" == "1" || ! -s "$merged" ]]; then
      merge_pairwise_generic "$merged" "$PREMERGE_IMAGE_SIZE_DEG" "${clean_files[@]}"
      attcalc eventset="$merged" \
        attitudelabel=ahf \
        refpointlabel=pnt \
        withatthkset=yes atthkset="$ATTHKGEN_FILE"
    fi

    local merged_gti="$(flare_gti_for_inst "$inst")"
    if [[ "$FORCE" == "1" || ! -s "$merged_gti" ]]; then
      local -a gti_list=()
      local gf
      for gf in "$WORKDIR/clean/$inst"/*.flare_gti.fits; do
        [[ -f "$gf" ]] && gti_list+=("$gf")
      done
      if [[ ${#gti_list[@]} -gt 0 ]]; then
        merge_flare_gtis "$merged_gti" "${gti_list[@]}"
      fi
    fi
  done
}

stage_track() {
  need_cmd python3
  local refevt
  refevt="$(first_available "$(merged_clean_evt PN)" "$(merged_clean_evt M1)" "$(merged_clean_evt M2)")"
  [[ -n "$refevt" ]] || { echo "No merged cleaned event list found; run clean first" >&2; exit 1; }

  mkdir -p "$WORKDIR/track"
  if [[ "$FORCE" != "1" && -s "$WORKDIR/track/comet_track.fits" && -s "$WORKDIR/track/ref_attitude.env" ]]; then
    echo "Skipping track; existing track products found"
    return 0
  fi

  local cmd=(python3 "$HELPER_DIR/make_comet_track_and_segments.py"
    --reference-event "$refevt"
    --track-out "$WORKDIR/track/comet_track.fits"
    --segments-out "$WORKDIR/track/motion_segments.csv"
    --max-blur-arcsec "$MAX_SEG_BLUR_ARCSEC"
    --snap-s "$SEGMENT_SNAP_S"
    --min-seg-s "$MIN_SEG_S")
  if [[ -n "$TRACK_INPUT" ]]; then
    cmd+=(--track-input "$TRACK_INPUT")
  else
    [[ -n "$TARGET_ID" ]] || { echo "Set TARGET_ID or TRACK_INPUT in config.env" >&2; exit 1; }
    cmd+=(--target-id "$TARGET_ID" --id-type "$TARGET_ID_TYPE" --step "$TRACK_STEP")
  fi
  runlog track_build "${cmd[@]}"

  python3 - "$WORKDIR/track/comet_track.fits" "$WORKDIR/track/ref_attitude.env" "$COMET_REF_RA" "$COMET_REF_DEC" <<'PY'
import sys
from astropy.io import fits
import numpy as np
track_fits, out_env, ref_ra_in, ref_dec_in = sys.argv[1:5]
with fits.open(track_fits) as hdul:
    tab = None
    for hdu in hdul[1:]:
        if isinstance(hdu, fits.BinTableHDU) and hdu.data is not None and len(hdu.data) > 0:
            tab = hdu.data
            break
    if tab is None:
        raise RuntimeError("No non-empty track table")
    ra = np.asarray(tab['RA'], dtype=float)
    dec = np.asarray(tab['DEC'], dtype=float)
mid = len(ra) // 2
ref_ra = float(ref_ra_in) if ref_ra_in.strip() else float(ra[mid])
ref_dec = float(ref_dec_in) if ref_dec_in.strip() else float(dec[mid])
with open(out_env, 'w', encoding='utf-8') as fh:
    fh.write(f'COMET_REF_RA="{ref_ra:.10f}"\n')
    fh.write(f'COMET_REF_DEC="{ref_dec:.10f}"\n')
print(f"Reference comet frame center: RA={ref_ra:.10f} DEC={ref_dec:.10f}")
PY
}

stage_detect() {
  need_cmd emosaic_prep
  need_cmd edetect_stack
  need_cmd atthkgen
  need_cmd python3
  setup_sas_env
  load_attitude
  need_file "$WORKDIR/track/comet_track.fits"

  mkdir -p "$WORKDIR/detect" "$WORKDIR/detect_mosaic"

  local pn_evt="$(merged_clean_evt PN)"
  local m1_evt="$(merged_clean_evt M1)"
  local m2_evt="$(merged_clean_evt M2)"
  [[ -f "$pn_evt" || -f "$m1_evt" || -f "$m2_evt" ]] || { echo "No merged cleaned event lists available for detect stage" >&2; exit 1; }

  local stack_atthk="$WORKDIR/detect/stack_atthk.dat"
  if [[ "$FORCE" == "1" || ! -s "$stack_atthk" ]]; then
    ( export SAS_CCF="$SAS_CCF"; export SAS_ODF="$SAS_ODF"; export SAS_ATTITUDE=AHF; atthkgen atthkset="$stack_atthk" )
  fi

  local mosaic_done="$WORKDIR/detect_mosaic/.done"
  if [[ "$FORCE" == "1" || ! -f "$mosaic_done" ]]; then
    local tmpmos="/tmp/_xmmpipe_mosaic_$$"
    rm -rf "$tmpmos"
    mkdir -p "$tmpmos"

    local short_atthk="$tmpmos/atthk.dat"
    cp -f "$stack_atthk" "$short_atthk"
    local mosaic_args=(atthkfile="$short_atthk" pseudoexpid=10)
    local _det_inst _det_evt _det_gti _det_key
    for _det_inst in PN:pn:pnevtfile M1:m1:mos1evtfile M2:m2:mos2evtfile; do
      IFS=: read -r _det_key _det_short _det_param <<< "$_det_inst"
      _det_evt="$(merged_clean_evt "$_det_key")"
      [[ -f "$_det_evt" ]] || continue
      _det_gti="$(flare_gti_for_inst "$_det_key")"
      if [[ -f "$_det_gti" ]]; then
        evselect table="${_det_evt}:EVENTS" \
          withfilteredset=yes filteredset="$tmpmos/${_det_short}.fits" \
          destruct=yes keepfilteroutput=yes updateexposure=yes writedss=yes \
          expression="GTI(${_det_gti},TIME)"
      else
        cp -f "$_det_evt" "$tmpmos/${_det_short}.fits"
      fi
      mosaic_args+=(${_det_param}="$tmpmos/${_det_short}.fits")
    done

    pushd "$tmpmos" >/dev/null
    emosaic_prep "${mosaic_args[@]}"
    popd >/dev/null

    rm -rf "$WORKDIR/detect_mosaic"
    mkdir -p "$WORKDIR/detect_mosaic"
    mv "$tmpmos"/prep_mosaic_* "$WORKDIR/detect_mosaic/" 2>/dev/null || true
    mv "$tmpmos"/gti_positions.ds "$WORKDIR/detect_mosaic/" 2>/dev/null || true
    rm -rf "$tmpmos"
    local _link _target _fixed
    for _link in "$WORKDIR/detect_mosaic"/prep_mosaic_*/*SUM.SAS \
                 "$WORKDIR/detect_mosaic"/prep_mosaic_*/ccf.cif; do
      [[ -L "$_link" ]] || continue
      _target="$(readlink "$_link")"
      _fixed="${_target/$WORKDIR/$REAL_WORKDIR}"
      [[ "$_fixed" != "$_target" ]] && ln -sfn "$_fixed" "$_link"
    done
    touch "$mosaic_done"
  fi

  local tmpdet="/tmp/_xmmpipe_detect_$$"
  rm -rf "$tmpdet"
  mkdir -p "$tmpdet/mosaic" "$tmpdet/detect"

  python3 - "$WORKDIR/detect_mosaic" "$tmpdet" "$WORKDIR/detect/stack_eventsets.txt" "$WORKDIR/detect/pseudoexposure_manifest.tsv" <<'PY'
import csv
import sys, shutil
from pathlib import Path
from astropy.io import fits

root = Path(sys.argv[1])
tmpdet = Path(sys.argv[2])
out_list = Path(sys.argv[3])
out_tsv = Path(sys.argv[4])

cands = []
for p in root.rglob('*'):
    if not p.is_file():
        continue
    b = p.name.lower()
    if b.endswith('.ds') and '_p' in b and 'gti' not in b and 'tmp' not in b:
        cands.append(p)

def read_header_metadata(hdul, hdr, nevts, ngti, p):
    obs = str(hdr.get('OBS_ID', ''))
    exp = str(hdr.get('EXP_ID', hdr.get('EXPIDSTR', '')))
    inst = str(hdr.get('INSTRUME', ''))
    submode = str(hdr.get('SUBMODE', ''))
    datamode = str(hdr.get('DATAMODE', ''))
    filt = str(hdr.get('FILTER', ''))
    tstart = hdr.get('TSTART', '')
    tstop = hdr.get('TSTOP', '')
    ontime = hdr.get('ONTIME', '')
    livetime = hdr.get('LIVETIME', '')
    return {
        'filename': p.name, 'obs_id': obs, 'exp_id': exp, 'instrument': inst,
        'submode': submode, 'datamode': datamode, 'filter': filt,
        'tstart': str(tstart), 'tstop': str(tstop),
        'ontime': str(ontime), 'livetime': str(livetime),
        'events': nevts, 'gti_rows': ngti, 'path': str(p),
    }

rows = []
skipped_rows = []
all_rows = []
for p in cands:
    try:
        with fits.open(p) as hdul:
            hdu = hdul['EVENTS'] if 'EVENTS' in hdul else hdul[1]
            hdr = hdu.header
            nevts = len(hdu.data) if hdu.data is not None else 0
            gti_ext = None
            for i in range(1, 20):
                if hdr.get(f'DSTYP{i}') == 'TIME':
                    ref = hdr.get(f'DSREF{i}', '')
                    if ref.startswith(':'):
                        gti_ext = ref[1:]
                    break
            ngti = 0
            if gti_ext and gti_ext in hdul:
                gd = hdul[gti_ext].data
                ngti = len(gd) if gd is not None else 0
            meta = read_header_metadata(hdul, hdr, nevts, ngti, p)
            if nevts == 0 or ngti == 0:
                meta['status'] = 'skipped'
                skipped_rows.append(meta)
                all_rows.append(meta)
                print(f'SKIP {p.name}: events={nevts} gti_rows={ngti}', file=sys.stderr)
                continue
            obs = meta['obs_id']
            exp = meta['exp_id']
            inst = meta['instrument']
            short = tmpdet / 'mosaic' / p.name
            shutil.copy2(str(p), str(short))
            rows.append((obs, exp, inst, nevts, ngti, str(short)))
            meta['status'] = 'accepted'
            meta['short_path'] = str(short)
            all_rows.append(meta)
    except Exception:
        continue
rows.sort()
out_list.parent.mkdir(parents=True, exist_ok=True)
out_list.write_text('\n'.join(r[5] for r in rows) + ('\n' if rows else ''), encoding='utf-8')
out_tsv.parent.mkdir(parents=True, exist_ok=True)
with out_tsv.open('w', newline='', encoding='utf-8') as fh:
    writer = csv.writer(fh, delimiter='\t')
    writer.writerow(['obs_id', 'exp_id', 'instrument', 'events', 'gti_rows', 'short_path'])
    for row in rows:
        writer.writerow(row)

all_tsv = out_tsv.parent / 'all_pseudoexposures.tsv'
cols = ['filename', 'status', 'obs_id', 'exp_id', 'instrument', 'submode', 'datamode',
        'filter', 'tstart', 'tstop', 'ontime', 'livetime', 'events', 'gti_rows', 'path']
all_rows.sort(key=lambda r: r['filename'])
with all_tsv.open('w', newline='', encoding='utf-8') as fh:
    writer = csv.writer(fh, delimiter='\t')
    writer.writerow(cols)
    for r in all_rows:
        writer.writerow([r.get(c, '') for c in cols])

print(f'{len(rows)} pseudo-exposures accepted for detection', file=sys.stderr)
if skipped_rows:
    print(f'{len(skipped_rows)} pseudo-exposures skipped (empty GTI/events)', file=sys.stderr)
PY

  [[ -s "$WORKDIR/detect/stack_eventsets.txt" ]] || { echo "No pseudo-exposure event lists were created under $WORKDIR/detect_mosaic" >&2; rm -rf "$tmpdet"; exit 1; }
  cp -f "$WORKDIR/detect/stack_eventsets.txt" "$tmpdet/detect/stack_eventsets.txt"
  cp -f "$stack_atthk" "$tmpdet/detect/atthk.dat"
  cp -f "$SAS_ODF" "$tmpdet/detect/SUM.SAS"

  local -a srcdet_pimin=() srcdet_pimax=()
  local _srcline _spimin _spimax
  while read -r _srcline; do
    [[ -n "$_srcline" ]] || continue
    read -r _spimin _spimax <<< "$_srcline"
    srcdet_pimin+=("$_spimin")
    srcdet_pimax+=("$_spimax")
  done < <(srcdet_bands_lines)
  [[ ${#srcdet_pimin[@]} -gt 0 ]] || { echo "No source-detection energy bands configured" >&2; rm -rf "$tmpdet"; exit 1; }
  local srcdet_pimin_arg="${srcdet_pimin[*]}"
  local srcdet_pimax_arg="${srcdet_pimax[*]}"

  pushd "$tmpdet/detect" >/dev/null
  if [[ "$FORCE" == "1" || ! -s "$WORKDIR/detect/stack_srclist.fits" ]]; then
    edetect_stack \
      eventsets="@stack_eventsets.txt" \
      attitudesets="$tmpdet/detect/atthk.dat" \
      summarysets="$tmpdet/detect/SUM.SAS" \
      pimin="$srcdet_pimin_arg" \
      pimax="$srcdet_pimax_arg" \
      mlmin="$DETECT_MIN_ML" \
      srclistset="srclist.fits" \
      informational=all \
      compress=false \
      prefix="stack_" \
      emos_mosaicedset="mosaic" \
      > "$LOGDIR/detect_edetect_stack.log" 2>&1

    cp -f stack_srclist.fits "$WORKDIR/detect/" 2>/dev/null || true
    cp -f stack_*.fits "$WORKDIR/detect/" 2>/dev/null || true
    cp -f mosaic*.fits "$WORKDIR/detect/" 2>/dev/null || true
  fi
  popd >/dev/null
  rm -rf "$tmpdet"

  local -a preferred_detect_bands=("00500 02000" "01000 02000")
  preferred_detect_bands+=("$(printf '%05d %05d' "${srcdet_pimin[0]}" "${srcdet_pimax[0]}")")
  local pair low high cand
  for pair in "${preferred_detect_bands[@]}"; do
    read -r low high <<< "$pair"
    for cand in       "$WORKDIR/detect/stack_mosaic_EPIC_${low}_${high}.fits"       "$WORKDIR/detect/mosaic_EPIC_${low}_${high}.fits"; do
      if [[ -f "$cand" ]]; then
        cp -f "$cand" "$WORKDIR/detect/EPIC_srcdet_mosaic.fits"
        break 2
      fi
    done
  done

  python3 "$HELPER_DIR/curate_emllist_for_region.py"     --input "$WORKDIR/detect/stack_srclist.fits"     --output "$WORKDIR/detect/field_sources_all.fits"     --ds9 "$WORKDIR/detect/field_sources_all.reg"     --min-det-ml "$DETECT_MIN_ML"     --ds9-radius-arcsec "$DETECT_DS9_RADIUS_ARCSEC"

  python3 "$HELPER_DIR/curate_emllist_for_region.py"     --input "$WORKDIR/detect/stack_srclist.fits"     --output "$WORKDIR/detect/field_sources_curated.fits"     --ds9 "$WORKDIR/detect/field_sources_curated.reg"     --track "$WORKDIR/track/comet_track.fits"     --exclude-track-radius-arcsec "$DETECT_EXCLUDE_TRACK_ARCSEC"     --min-det-ml "$DETECT_MIN_ML"     --ds9-radius-arcsec "$DETECT_DS9_RADIUS_ARCSEC"

  local -a qc_events=()
  [[ -f "$pn_evt" ]] && qc_events+=("$pn_evt")
  [[ -f "$m1_evt" ]] && qc_events+=("$m1_evt")
  [[ -f "$m2_evt" ]] && qc_events+=("$m2_evt")
  local qc_line qc_label qc_pimin qc_pimax
  while read -r qc_line; do
    [[ -n "$qc_line" ]] || continue
    read -r qc_label qc_pimin qc_pimax <<< "$qc_line"
    python3 "$HELPER_DIR/make_still_sky_mosaic.py"       --events "${qc_events[@]}"       --pimin "$qc_pimin"       --pimax "$qc_pimax"       --track "$WORKDIR/track/comet_track.fits"       --srclist "$WORKDIR/detect/field_sources_all.fits"       --bin-arcsec "$DETECT_QC_BIN_ARCSEC"       --pad-arcmin "$DETECT_QC_PAD_ARCMIN"       --label "$qc_label"       --out-fits "$WORKDIR/detect/EPIC_${qc_label}_mosaic.fits"       --out-json "$WORKDIR/detect/EPIC_${qc_label}_mosaic.json"
  done < <(detect_qc_bands_lines)

  if [[ -f "$WORKDIR/detect/EPIC_hard_mosaic.fits" ]]; then
    :
  elif [[ -f "$WORKDIR/detect/EPIC_srcdet_mosaic.fits" ]]; then
    cp -f "$WORKDIR/detect/EPIC_srcdet_mosaic.fits" "$WORKDIR/detect/EPIC_hard_mosaic.fits"
  fi
}

stage_comet() {
  need_cmd attmove
  need_cmd attcalc
  need_cmd odfingest
  need_cmd python3
  setup_sas_env
  load_attitude
  load_ref_attitude

  mkdir -p "$WORKDIR/comet"
  local moved_root="$WORKDIR/comet/moved_odf"
  local temp_odf_dir="$moved_root/odf"
  local ingest_dir="$moved_root/ingest"
  local moved_env="$WORKDIR/comet/moved_sas_setup.env"
  local moved_atthk="$WORKDIR/comet/moved_atthk.dat"

  local ahf_input="${AHF_INPUT:-$AHF_FILE}"
  [[ -f "$ahf_input" ]] || { echo "AHF/ATS file not found: $ahf_input" >&2; exit 1; }
  need_file "$ATTHKGEN_FILE"

  local relpath
  relpath="$(python3 - "$ODFDIR" "$ahf_input" <<'PY'
import sys
from pathlib import Path
odf = Path(sys.argv[1]).resolve()
ahf = Path(sys.argv[2]).resolve()
try:
    print(ahf.relative_to(odf))
except Exception:
    matches = list(odf.rglob(ahf.name))
    if not matches:
        raise SystemExit(1)
    print(matches[0].resolve().relative_to(odf))
PY
)" || { echo "Could not determine the relative path of the original ATS/AHF inside ODFDIR" >&2; exit 1; }

  # --- sub-step 1: attmove (create moved ODF symlink tree + shifted AHF) ---
  if [[ "$FORCE" == "1" || ! -s "$temp_odf_dir/$relpath" ]]; then
    rm -rf "$temp_odf_dir"
    mkdir -p "$temp_odf_dir"
    cp -as "${ODFDIR%/}/." "$temp_odf_dir"/
    find "$temp_odf_dir" \( -name '*SUM.SAS' -o -name 'ccf.cif' \) -exec rm -f {} +
    mkdir -p "$(dirname "$temp_odf_dir/$relpath")"
    rm -f "$temp_odf_dir/$relpath"

    attmove \
      input="$ahf_input" \
      output="$temp_odf_dir/$relpath" \
      track="$WORKDIR/track/comet_track.fits" \
      withrefatt=yes \
      refra="$COMET_REF_RA" \
      refdec="$COMET_REF_DEC" \
      granularity="$ATTMOVE_GRANULARITY" \
      minstable="$ATTMOVE_MINSTABLE" \
      creatediagnostics=yes \
      diagfile="$WORKDIR/comet/attmove_diag.fits"
  else
    echo "  [skip] attmove output already exists"
  fi

  # --- sub-step 2: odfingest (re-ingest moved ODF) ---
  local moved_sum
  moved_sum="$(find "$ingest_dir" -maxdepth 1 -type f -name '*SUM.SAS' 2>/dev/null | sort | head -n 1 || true)"
  if [[ "$FORCE" == "1" || -z "$moved_sum" || ! -s "$moved_sum" ]]; then
    rm -rf "$ingest_dir"
    mkdir -p "$ingest_dir"
    odfingest odfdir="$temp_odf_dir" outdir="$ingest_dir"
    moved_sum="$(find "$ingest_dir" -maxdepth 1 -type f -name '*SUM.SAS' | sort | head -n 1 || true)"
    [[ -n "$moved_sum" ]] || { echo "odfingest did not create a moved *SUM.SAS under $ingest_dir" >&2; exit 1; }
    sed -i "s|^PATH .*|PATH ${temp_odf_dir/$WORKDIR/$REAL_WORKDIR}/|" "$moved_sum"
  else
    echo "  [skip] odfingest SUM.SAS already exists"
  fi

  # --- sub-step 3: build_moved_atthk.py ---
  if [[ "$FORCE" == "1" || ! -s "$moved_atthk" ]]; then
    python3 "$HELPER_DIR/build_moved_atthk.py" \
      --input-atthk "$ATTHKGEN_FILE" \
      --track "$WORKDIR/track/comet_track.fits" \
      --output "$moved_atthk" \
      --ref-ra "$COMET_REF_RA" \
      --ref-dec "$COMET_REF_DEC" \
      --report-json "$WORKDIR/comet/moved_atthk.json"
  else
    echo "  [skip] moved_atthk already exists"
  fi

  # --- sub-step 4: write moved_sas_setup.env ---
  if [[ "$FORCE" == "1" || ! -s "$moved_env" ]]; then
    local real_moved_env="${moved_env/$WORKDIR/$REAL_WORKDIR}"
    local real_ccf="${SAS_CCF/$WORKDIR/$REAL_WORKDIR}"
    local real_sum="${moved_sum/$WORKDIR/$REAL_WORKDIR}"
    local real_odf="${temp_odf_dir/$WORKDIR/$REAL_WORKDIR}"
    local real_atthk="${moved_atthk/$WORKDIR/$REAL_WORKDIR}"
    cat > "$moved_env" <<EOF
export SAS_CCF="$real_ccf"
export SAS_ODF="$real_sum"
export SAS_ATTITUDE="AHF"
export MOVED_ODF_DIR="$real_odf"
export MOVED_AHF_FILE="$real_odf/$relpath"
export MOVED_ATTHK_FILE="$real_atthk"
EOF
  else
    echo "  [skip] moved_sas_setup.env already exists"
  fi

  source "$moved_env"
  export SAS_CCF="$WORKDIR/init/ccf.cif"
  local cur_moved_sum
  cur_moved_sum="$(find "$WORKDIR/comet/moved_odf/ingest" -maxdepth 1 -type f -name '*SUM.SAS' | sort | head -n 1 || true)"
  if [[ -n "$cur_moved_sum" && -f "$cur_moved_sum" ]]; then
    export SAS_ODF="$cur_moved_sum"
    sed -i "s|^PATH .*|PATH $REAL_WORKDIR/comet/moved_odf/odf/|" "$cur_moved_sum"
  fi
  export MOVED_ATTHK_FILE="$WORKDIR/comet/moved_atthk.dat"
  cp -f "$MOVED_ATTHK_FILE" "$WORKDIR/comet/comet_atthk.fits"

  local inst infile outfile
  for inst in $(selected_insts); do
    infile="$(merged_clean_evt "$inst")"
    [[ -f "$infile" ]] || continue
    outfile="$(inst_comet_evt "$inst")"
    if skip_file "$outfile" && [[ "$(fits_rows "$outfile")" -gt 0 ]]; then
      continue
    fi
    cp -f "$infile" "$outfile"
    attcalc \
      eventset="$outfile" \
      attitudelabel=ahf \
      refpointlabel=user \
      nominalra="$COMET_REF_RA" \
      nominaldec="$COMET_REF_DEC" \
      setpnttouser=yes \
      imagesize="$ATTCALC_IMAGE_SIZE_DEG" \
      withatthkset=yes \
      atthkset="$MOVED_ATTHK_FILE"
  done
}

stage_contam() {
  need_cmd python3
  load_ref_attitude
  need_file "$WORKDIR/track/comet_track.fits"
  local contam_srclist="$WORKDIR/detect/field_sources_all.fits"
  [[ -f "$contam_srclist" ]] || contam_srclist="$WORKDIR/detect/field_sources_curated.fits"
  need_file "$contam_srclist"
  mkdir -p "$WORKDIR/contam"
  make_image_grid_if_needed

  local inst needs_mask=0
  local have_all=1
  for inst in $(selected_insts); do
    [[ -f "$(inst_comet_evt "$inst")" ]] || continue
    [[ -f "$WORKDIR/contam/${inst}_science_base.fits" ]] || have_all=0
  done
  if [[ "$FORCE" != "1" && -s "$WORKDIR/contam/science_gti.fits" && -s "$WORKDIR/contam/contamination_summary.json" && $have_all -eq 1 ]]; then
    echo "Skipping contam; found contamination products and science-base event lists"
    return 0
  fi

  local bkg_args=()
  if [[ "${BKG_MODE,,}" == "circle" ]]; then
    bkg_args=(--bkg-shape circle --bkg-dx-arcsec "$BKG_DX_ARCSEC" --bkg-dy-arcsec "$BKG_DY_ARCSEC" --bkg-r-arcsec "$BKG_R_ARCSEC")
  else
    bkg_args=(--bkg-shape annulus --bkg-dx-arcsec "$BKG_DX_ARCSEC" --bkg-dy-arcsec "$BKG_DY_ARCSEC" --bkg-rin-arcsec "$BKG_RIN_ARCSEC" --bkg-rout-arcsec "$BKG_ROUT_ARCSEC")
  fi

  python3 "$HELPER_DIR/build_contamination_products.py"     --track "$WORKDIR/track/comet_track.fits"     --srclist "$contam_srclist"     --output-dir "$WORKDIR/contam"     --science-policy "$SCIENCE_GTI_POLICY"     --sample-step-s "$CONTAM_SAMPLE_STEP_S"     --mask-radius-arcsec "$FIELD_SOURCE_MASK_R_ARCSEC"     --min-good-span-s "$CONTAM_MIN_GOOD_SPAN_S"     --src-dx-arcsec "$SRC_DX_ARCSEC"     --src-dy-arcsec "$SRC_DY_ARCSEC"     --src-r-arcsec "$SRC_R_ARCSEC"     --summary-json "$WORKDIR/contam/contamination_summary.json"     "${bkg_args[@]}"

  if bool_yes "$SCIENCE_APPLY_TRAIL_MASK"; then
    python3 "$HELPER_DIR/trail_mask_tools.py" make       --grid-json "$WORKDIR/images/grid.json"       --track "$WORKDIR/track/comet_track.fits"       --srclist "$contam_srclist"       --mask-radius-arcsec "$IMAGE_MASK_RADIUS_ARCSEC"       --out-mask "$WORKDIR/contam/EPIC_trail_mask.fits"       --out-json "$WORKDIR/contam/EPIC_trail_mask.json"
    needs_mask=1
  else
    needs_mask=0
  fi

  for inst in $(selected_insts); do
    local in_evt="$(inst_comet_evt "$inst")"
    [[ -f "$in_evt" ]] || continue
    local out_evt="$WORKDIR/contam/${inst}_science_base.fits"
    if [[ $needs_mask -eq 1 ]]; then
      python3 "$HELPER_DIR/trail_mask_tools.py" apply         --event "$in_evt"         --mask "$WORKDIR/contam/EPIC_trail_mask.fits"         --grid-json "$WORKDIR/images/grid.json"         --output "$out_evt"         --report-json "$WORKDIR/contam/${inst}_science_base.json"
    else
      cp -f "$in_evt" "$out_evt"
    fi
  done

  cat > "$WORKDIR/contam/science_selection.env" <<EOF
SCIENCE_GTI="$REAL_WORKDIR/contam/science_gti.fits"
SCIENCE_GTI_POLICY="$SCIENCE_GTI_POLICY"
SCIENCE_TRAIL_MASK="$REAL_WORKDIR/contam/EPIC_trail_mask.fits"
SCIENCE_APPLY_TRAIL_MASK="$SCIENCE_APPLY_TRAIL_MASK"
SCIENCE_SOURCE_LIST="${contam_srclist/$WORKDIR/$REAL_WORKDIR}"
EOF
}

stage_image() {
  need_cmd evselect
  need_cmd eexpmap
  need_cmd python3
  setup_sas_env
  source_if_exists "$WORKDIR/comet/moved_sas_setup.env"
  load_ref_attitude
  mkdir -p "$WORKDIR/images"
  make_image_grid_if_needed

  local moved_ahf="$WORKDIR/comet/moved_atthk.dat"
  need_file "$moved_ahf"

  local _img_inst _img_src _img_gti _img_filt
  for _img_inst in $(selected_insts); do
    _img_src="$(inst_comet_evt "$_img_inst")"
    [[ -f "$_img_src" ]] || continue
    _img_gti="$(flare_gti_for_inst "$_img_inst")"
    _img_filt="$WORKDIR/images/${_img_inst}_filt_events.fits"
    if [[ -f "$_img_gti" && ( "$FORCE" == "1" || ! -s "$_img_filt" ) ]]; then
      evselect table="${_img_src}:EVENTS" \
        withfilteredset=yes filteredset="$_img_filt" \
        destruct=yes keepfilteroutput=yes updateexposure=yes writedss=yes \
        expression="GTI(${_img_gti},TIME)"
    fi
  done

  local line label pimin pimax inst evt img exp
  while read -r line; do
    [[ -n "$line" ]] || continue
    read -r label pimin pimax <<< "$line"
    local banddir="$WORKDIR/images/$label"
    mkdir -p "$banddir"
    local -a counts=() expos=()
    for inst in $(selected_insts); do
      local _filt="$WORKDIR/images/${inst}_filt_events.fits"
      if [[ -f "$_filt" ]]; then
        evt="$_filt"
      else
        evt="$(inst_comet_evt "$inst")"
      fi
      [[ -f "$evt" ]] || continue
      img="$banddir/${inst}_counts.fits"
      exp="$banddir/${inst}_exp.fits"
      if [[ "$FORCE" == "1" || ! -s "$img" ]]; then
        evselect table="${evt}:EVENTS" \
          withimageset=yes imageset="$img" \
          xcolumn=X ycolumn=Y imagebinning=binSize \
          ximagebinsize="$BIN_PHYS" yimagebinsize="$BIN_PHYS" \
          withxranges=yes ximagemin="$X_MIN_PHYS" ximagemax="$X_MAX_PHYS" \
          withyranges=yes yimagemin="$Y_MIN_PHYS" yimagemax="$Y_MAX_PHYS" \
          writedss=yes \
          expression="PI in [$pimin:$pimax]"
      fi
      if [[ "$FORCE" == "1" || ! -s "$exp" ]]; then
        eexpmap imageset="$img" \
          attitudeset="$moved_ahf" \
          eventset="$evt" \
          expimageset="$exp" \
          pimin="$pimin" pimax="$pimax" \
          attrebin="$SCIENCE_EEXPMAP_ATTREBIN"
      fi
      counts+=("$img")
      expos+=("$exp")
    done

    [[ ${#counts[@]} -gt 0 ]] || continue

    local combine_args=(
      --counts "${counts[@]}"
      --exposure "${expos[@]}"
      --out-counts "$banddir/EPIC_counts.fits"
      --out-exposure "$banddir/EPIC_exp.fits"
      --out-rate "$banddir/EPIC_rate.fits"
      --out-json "$banddir/EPIC_image_summary.json"
    )
    local mask_srclist="$WORKDIR/detect/field_sources_all.fits"
    [[ -f "$mask_srclist" ]] || mask_srclist="$WORKDIR/detect/field_sources_curated.fits"
    if [[ -f "$mask_srclist" && "$(srclist_rows "$mask_srclist")" -gt 0 ]]; then
      combine_args+=(
        --grid-json "$WORKDIR/images/grid.json"
        --track "$WORKDIR/track/comet_track.fits"
        --srclist "$mask_srclist"
        --mask-radius-arcsec "$IMAGE_MASK_RADIUS_ARCSEC"
        --out-mask "$banddir/EPIC_trail_mask.fits"
        --out-rate-masked "$banddir/EPIC_rate_masked.fits"
      )
    fi
    python3 "$HELPER_DIR/combine_epic_images.py" "${combine_args[@]}"
    if bool_yes "$IMAGE_WRITE_NET_PRODUCTS"; then
      python3 "$HELPER_DIR/image_postproc.py" science \
        --banddir "$banddir" \
        --sigma-px "$IMAGE_BKG_SMOOTH_SIGMA_PX" \
        --out-json "$banddir/EPIC_clean_summary.json"
    fi
  done < <(image_bands_lines)
}

stage_lcurve() {
  need_cmd evselect
  need_cmd epiclccorr
  need_cmd elcbuild
  need_cmd python3
  setup_sas_env
  need_file "$WORKDIR/contam/science_gti.fits"

  mkdir -p "$WORKDIR/lcurve"
  local -a abs_lcs=() rel_lcs=()
  local inst evt lc_evt src_raw bkg_raw corr_abs corr_rel mode_txt chosen mode status
  local summary_tsv="$WORKDIR/lcurve/lc_summary.tsv"
  printf 'inst\tmode\tfile\n' > "$summary_tsv"

  for inst in $(selected_insts); do
    evt="$(inst_analysis_evt "$inst")"
    [[ -f "$evt" ]] || continue

    lc_evt="$WORKDIR/lcurve/${inst}_lc_events.fits"
    if [[ "$FORCE" == "1" || ! -s "$lc_evt" ]]; then
      evselect table="${evt}:EVENTS" \
        withfilteredset=yes filteredset="$lc_evt" \
        destruct=yes keepfilteroutput=yes updateexposure=yes writedss=yes \
        expression="GTI($WORKDIR/contam/science_gti.fits,TIME)&&(PI in [$LC_PI_MIN:$LC_PI_MAX])"
    fi
    if [[ "$(fits_rows "$lc_evt")" -le 0 ]]; then
      echo "No events left for $inst light curve after science GTI and PI cut; skipping" >&2
      continue
    fi

    build_region_env "$lc_evt" SRC circle --dx "$SRC_DX_ARCSEC" --dy "$SRC_DY_ARCSEC" --r "$SRC_R_ARCSEC"
    if [[ "${BKG_MODE,,}" == "circle" ]]; then
      build_region_env "$lc_evt" BKG circle --dx "$BKG_DX_ARCSEC" --dy "$BKG_DY_ARCSEC" --r "$BKG_R_ARCSEC"
    else
      build_region_env "$lc_evt" BKG annulus --dx "$BKG_DX_ARCSEC" --dy "$BKG_DY_ARCSEC" --rin "$BKG_RIN_ARCSEC" --rout "$BKG_ROUT_ARCSEC"
    fi

    src_raw="$WORKDIR/lcurve/${inst}_src_raw.fits"
    bkg_raw="$WORKDIR/lcurve/${inst}_bkg_raw.fits"
    corr_abs="$WORKDIR/lcurve/${inst}_corr_abs.fits"
    corr_rel="$WORKDIR/lcurve/${inst}_corr_relonly.fits"
    mode_txt="$WORKDIR/lcurve/${inst}_corr_mode.txt"

    if [[ "$FORCE" == "1" || ! -s "$src_raw" ]]; then
      evselect table="${lc_evt}:EVENTS" \
        withrateset=yes rateset="$src_raw" \
        maketimecolumn=yes makeratecolumn=yes \
        timebinsize="$LC_BIN_S" writedss=yes \
        expression="$SRC_EXPR"
    fi
    if [[ "$FORCE" == "1" || ! -s "$bkg_raw" ]]; then
      evselect table="${lc_evt}:EVENTS" \
        withrateset=yes rateset="$bkg_raw" \
        maketimecolumn=yes makeratecolumn=yes \
        timebinsize="$LC_BIN_S" writedss=yes \
        expression="$BKG_EXPR"
    fi

    chosen=""
    mode=""
    rm -f "$WORKDIR/lcurve/${inst}_corr.fits"
    if [[ "$FORCE" == "1" || ( ! -s "$corr_abs" && ! -s "$corr_rel" ) ]]; then
      rm -f "$corr_abs" "$corr_rel" "$mode_txt"
      local log1="$WORKDIR/lcurve/${inst}_epiclccorr.log"
      set +e
      epiclccorr srctslist="$src_raw" \
        eventlist="$lc_evt" \
        withbkgset=yes bkgtslist="$bkg_raw" \
        outset="$corr_abs" \
        applyabsolutecorrections="$LC_APPLY_ABSOLUTE_CORRECTIONS" \
        withsourcepos=yes sourcecoords=pos sourcex="$SRC_X" sourcey="$SRC_Y" \
        detxbins="$LC_DETXBINS" detybins="$LC_DETYBINS" \
        > "$log1" 2>&1
      status=$?
      set -e

      if [[ $status -eq 0 && -s "$corr_abs" ]]; then
        if [[ "${LC_APPLY_ABSOLUTE_CORRECTIONS,,}" == "yes" ]]; then
          mode="abs"
          chosen="$corr_abs"
        else
          mode="relonly"
          mv -f "$corr_abs" "$corr_rel"
          chosen="$corr_rel"
        fi
      elif [[ "${LC_APPLY_ABSOLUTE_CORRECTIONS,,}" == "yes" ]] && bool_yes "$LC_ALLOW_NOABS_FALLBACK" && grep -Eq 'ZeroAreaCorrection|ZeroArfgenArea|wrongArfgenResults' "$log1"; then
        echo "epiclccorr absolute correction failed for $inst; retrying with applyabsolutecorrections=no" >&2
        rm -f "$corr_abs"
        local log2="$WORKDIR/lcurve/${inst}_epiclccorr_noabs.log"
        set +e
        epiclccorr srctslist="$src_raw" \
          eventlist="$lc_evt" \
          withbkgset=yes bkgtslist="$bkg_raw" \
          outset="$corr_rel" \
          applyabsolutecorrections=no \
          withsourcepos=yes sourcecoords=pos sourcex="$SRC_X" sourcey="$SRC_Y" \
          detxbins="$LC_DETXBINS" detybins="$LC_DETYBINS" \
          > "$log2" 2>&1
        status=$?
        set -e
        if [[ $status -eq 0 && -s "$corr_rel" ]]; then
          mode="relonly"
          chosen="$corr_rel"
        else
          echo "WARNING: epiclccorr failed for $inst even with applyabsolutecorrections=no; using raw light curve. See $log1 and $log2" >&2
          mode="raw"
          chosen="$src_raw"
        fi
      else
        echo "WARNING: epiclccorr failed for $inst; using raw light curve. See $log1" >&2
        mode="raw"
        chosen="$src_raw"
      fi
      printf 'mode=%s\nfile=%s\n' "$mode" "$chosen" > "$mode_txt"
    else
      if [[ -s "$corr_abs" ]]; then
        mode="abs"
        chosen="$corr_abs"
      fi
      if [[ -s "$corr_rel" ]]; then
        if [[ -z "$chosen" || "${LC_APPLY_ABSOLUTE_CORRECTIONS,,}" != "yes" ]]; then
          mode="relonly"
          chosen="$corr_rel"
        fi
      fi
      if [[ -z "$chosen" && -s "$src_raw" ]]; then
        mode="raw"
        chosen="$src_raw"
      fi
    fi

    [[ -n "$chosen" && -s "$chosen" ]] || continue
    ln -sfn "$(basename "$chosen")" "$WORKDIR/lcurve/${inst}_corr.fits"
    printf '%s\t%s\t%s\n' "$inst" "$mode" "$chosen" >> "$summary_tsv"
    if [[ "$mode" == "abs" ]]; then
      abs_lcs+=("$chosen")
    elif [[ "$mode" == "relonly" ]]; then
      rel_lcs+=("$chosen")
    fi
  done

  local preferred_combined=""
  local preferred_mode=""
  if [[ ${#abs_lcs[@]} -ge 1 ]]; then
    if [[ ${#abs_lcs[@]} -ge 2 ]]; then
      local sets_csv_abs
      sets_csv_abs="$(IFS=,; echo "${abs_lcs[*]}")"
      if [[ "$FORCE" == "1" || ! -s "$WORKDIR/lcurve/EPIC_bundle_lc_abs.fits" ]]; then
        elcbuild sets="$sets_csv_abs" outset="$WORKDIR/lcurve/EPIC_bundle_lc_abs.fits" || echo "WARNING: elcbuild (abs) failed; skipping bundle" >&2
      fi
    fi
    if python3 "$HELPER_DIR/combine_epic_lightcurves.py" \
      --inputs "${abs_lcs[@]}" \
      --out-fits "$WORKDIR/lcurve/EPIC_total_corr_abs_lc.fits" \
      --out-csv "$WORKDIR/lcurve/EPIC_total_corr_abs_lc.csv"; then
      preferred_combined="$WORKDIR/lcurve/EPIC_total_corr_abs_lc.fits"
      preferred_mode="abs"
    else
      echo "WARNING: combine_epic_lightcurves.py (abs) failed; per-instrument LCs still available" >&2
    fi
  fi

  if [[ ${#rel_lcs[@]} -ge 1 ]]; then
    if [[ ${#rel_lcs[@]} -ge 2 ]]; then
      local sets_csv_rel
      sets_csv_rel="$(IFS=,; echo "${rel_lcs[*]}")"
      if [[ "$FORCE" == "1" || ! -s "$WORKDIR/lcurve/EPIC_bundle_lc_relonly.fits" ]]; then
        elcbuild sets="$sets_csv_rel" outset="$WORKDIR/lcurve/EPIC_bundle_lc_relonly.fits" || echo "WARNING: elcbuild (rel) failed; skipping bundle" >&2
      fi
    fi
    if python3 "$HELPER_DIR/combine_epic_lightcurves.py" \
      --inputs "${rel_lcs[@]}" \
      --out-fits "$WORKDIR/lcurve/EPIC_total_corr_relonly_lc.fits" \
      --out-csv "$WORKDIR/lcurve/EPIC_total_corr_relonly_lc.csv"; then
      if [[ -z "$preferred_combined" ]]; then
        preferred_combined="$WORKDIR/lcurve/EPIC_total_corr_relonly_lc.fits"
        preferred_mode="relonly"
      fi
    else
      echo "WARNING: combine_epic_lightcurves.py (rel) failed; per-instrument LCs still available" >&2
    fi
  fi

  if [[ -n "$preferred_combined" && -s "$preferred_combined" ]]; then
    ln -sfn "$(basename "$preferred_combined")" "$WORKDIR/lcurve/EPIC_total_corr_lc.fits"
    if [[ "$preferred_combined" == *.fits ]]; then
      local csv_pref="${preferred_combined%.fits}.csv"
      if [[ -f "$csv_pref" ]]; then
        ln -sfn "$(basename "$csv_pref")" "$WORKDIR/lcurve/EPIC_total_corr_lc.csv"
      fi
    fi
    printf 'mode=%s\nfile=%s\n' "$preferred_mode" "$preferred_combined" > "$WORKDIR/lcurve/EPIC_total_corr_mode.txt"
  else
    echo "WARNING: No combined light curve produced; per-instrument raw LCs are in $WORKDIR/lcurve/" >&2
  fi
}

stage_spectrum() {
  need_cmd evselect
  need_cmd especget
  need_cmd epicspeccombine
  need_cmd specgroup
  need_cmd python3
  setup_sas_env
  source_if_exists "$WORKDIR/comet/moved_sas_setup.env"
  need_file "$WORKDIR/contam/science_gti.fits"

  mkdir -p "$WORKDIR/spectra"
  local -a srcspecs=() bkgspecs=() arfs=() rmfs=()
  local detx_use="$SPEC_DETXBINS"
  local dety_use="$SPEC_DETYBINS"
  (( detx_use < SPEC_MIN_SAFE_DETBINS )) && detx_use="$SPEC_MIN_SAFE_DETBINS"
  (( dety_use < SPEC_MIN_SAFE_DETBINS )) && dety_use="$SPEC_MIN_SAFE_DETBINS"

  local inst evt spec_evt srcspec bkgspec arf rmf
  for inst in $(selected_insts); do
    evt="$(inst_analysis_evt "$inst")"
    [[ -f "$evt" ]] || continue

    spec_evt="$WORKDIR/spectra/${inst}_spec_events.fits"
    if [[ "$FORCE" == "1" || ! -s "$spec_evt" ]]; then
      local spec_gti_expr="GTI($WORKDIR/contam/science_gti.fits,TIME)"
      local spec_flare_gti="$(flare_gti_for_inst "$inst")"
      if [[ -f "$spec_flare_gti" ]]; then
        spec_gti_expr="${spec_gti_expr}&&GTI(${spec_flare_gti},TIME)"
      fi
      evselect table="${evt}:EVENTS" \
        withfilteredset=yes filteredset="$spec_evt" \
        destruct=yes keepfilteroutput=yes updateexposure=yes writedss=yes \
        expression="${spec_gti_expr}&&(PI in [$SPEC_PI_MIN:$SPEC_PI_MAX])"
    fi
    if [[ "$(fits_rows "$spec_evt")" -le 0 ]]; then
      echo "No events left for $inst spectrum after science GTI and PI cut; skipping" >&2
      continue
    fi

    build_region_env "$spec_evt" SRC circle --dx "$SRC_DX_ARCSEC" --dy "$SRC_DY_ARCSEC" --r "$SRC_R_ARCSEC"
    if [[ "${BKG_MODE,,}" == "circle" ]]; then
      build_region_env "$spec_evt" BKG circle --dx "$BKG_DX_ARCSEC" --dy "$BKG_DY_ARCSEC" --r "$BKG_R_ARCSEC"
    else
      build_region_env "$spec_evt" BKG annulus --dx "$BKG_DX_ARCSEC" --dy "$BKG_DY_ARCSEC" --rin "$BKG_RIN_ARCSEC" --rout "$BKG_ROUT_ARCSEC"
    fi

    srcspec="$WORKDIR/spectra/${inst}_src_spec.fits"
    bkgspec="$WORKDIR/spectra/${inst}_bkg_spec.fits"
    arf="$WORKDIR/spectra/${inst}_src.arf"
    rmf="$WORKDIR/spectra/${inst}_src.rmf"

    if [[ "$FORCE" == "1" || ! -s "$srcspec" || ! -s "$bkgspec" || ! -s "$arf" || ! -s "$rmf" ]]; then
      especget \
        table="${spec_evt}:EVENTS" \
        withfilestem=no \
        srcspecset="$srcspec" \
        bckspecset="$bkgspec" \
        witharfset=yes srcarfset="$arf" \
        withrmfset=no srcrmfset="$rmf" \
        srcexp="$SRC_EXPR" \
        backexp="$BKG_EXPR" \
        extendedsource=yes \
        useodfatt="$SPEC_USE_ODF_ATT" \
        withsourcepos=yes sourcecoords=pos sourcex="$SRC_X" sourcey="$SRC_Y" \
        detxbins="$detx_use" detybins="$dety_use"
    fi

    srcspecs+=("$srcspec")
    bkgspecs+=("$bkgspec")
    arfs+=("$arf")
    rmfs+=("$rmf")
  done

  if [[ ${#srcspecs[@]} -ge 2 ]]; then
    if epicspeccombine \
      pha="$(printf '%s ' "${srcspecs[@]}")" \
      bkg="$(printf '%s ' "${bkgspecs[@]}")" \
      arf="$(printf '%s ' "${arfs[@]}")" \
      rmf="$(printf '%s ' "${rmfs[@]}")" \
      filepha="$WORKDIR/spectra/EPIC_src_combined.fits" \
      filebkg="$WORKDIR/spectra/EPIC_bkg_combined.fits" \
      filersp="$WORKDIR/spectra/EPIC_rsp_combined.fits" \
      allowHEdiff=yes; then
      specgroup \
        spectrumset="$WORKDIR/spectra/EPIC_src_combined.fits" \
        groupedset="$WORKDIR/spectra/EPIC_src_combined_grp.fits" \
        backgndset="$WORKDIR/spectra/EPIC_bkg_combined.fits" \
        rmfset="$WORKDIR/spectra/EPIC_rsp_combined.fits" \
        mincounts="$SPECGROUP_MINCOUNTS" \
        oversample="$SPECGROUP_OVERSAMPLE" \
        addfilenames=yes
    else
      echo "WARNING: epicspeccombine failed (likely mismatched RMF energy grids); per-instrument spectra still available for simultaneous fitting" >&2
    fi
  elif [[ ${#srcspecs[@]} -eq 1 ]]; then
    cp -f "${srcspecs[0]}" "$WORKDIR/spectra/EPIC_src_combined.fits"
    cp -f "${bkgspecs[0]}" "$WORKDIR/spectra/EPIC_bkg_combined.fits"
    cp -f "${rmfs[0]}" "$WORKDIR/spectra/EPIC_rsp_combined.fits"
    specgroup \
      spectrumset="$WORKDIR/spectra/EPIC_src_combined.fits" \
      groupedset="$WORKDIR/spectra/EPIC_src_combined_grp.fits" \
      backgndset="$WORKDIR/spectra/EPIC_bkg_combined.fits" \
      rmfset="$WORKDIR/spectra/EPIC_rsp_combined.fits" \
      arfset="${arfs[0]}" \
      mincounts="$SPECGROUP_MINCOUNTS" \
      oversample="$SPECGROUP_OVERSAMPLE" \
      addfilenames=yes
  else
    echo "No spectra were produced." >&2
    exit 1
  fi
}

prepare_instid_event() {
  local inst="$1"
  local infile="$2"
  local outfile="$3"
  python3 "$HELPER_DIR/add_instid_column.py" --input "$infile" --output "$outfile" --inst "$inst"
}

stage_merge() {
  need_cmd merge
  need_cmd python3
  need_cmd evselect
  setup_sas_env
  load_ref_attitude
  mkdir -p "$WORKDIR/final"

  local inst evt base_evt scievt tagged
  local -a full_files=() sci_files=()
  for inst in $(selected_insts); do
    evt="$(inst_comet_evt "$inst")"
    [[ -f "$evt" ]] || continue

    tagged="$WORKDIR/final/${inst}_full_instid.fits"
    if [[ "$FORCE" == "1" || ! -s "$tagged" ]]; then
      prepare_instid_event "$inst" "$evt" "$tagged"
    fi
    full_files+=("$tagged")

    base_evt="$(inst_analysis_evt "$inst")"
    [[ -f "$base_evt" ]] || continue
    scievt="$WORKDIR/final/${inst}_science.fits"
    if [[ "$FORCE" == "1" || ! -s "$scievt" ]]; then
      evselect table="${base_evt}:EVENTS" \
        withfilteredset=yes filteredset="$scievt" \
        destruct=yes keepfilteroutput=yes updateexposure=yes writedss=yes \
        expression="GTI($WORKDIR/contam/science_gti.fits,TIME)"
    fi
    if [[ "$(fits_rows "$scievt")" -gt 0 ]]; then
      tagged="$WORKDIR/final/${inst}_science_instid.fits"
      if [[ "$FORCE" == "1" || ! -s "$tagged" ]]; then
        prepare_instid_event "$inst" "$scievt" "$tagged"
      fi
      sci_files+=("$tagged")
    fi
  done

  if [[ ${#full_files[@]} -ge 1 ]]; then
    merge_pairwise_cometref "$WORKDIR/final/EPIC_comet_merged_full.fits" "${full_files[@]}"
  fi
  if [[ ${#sci_files[@]} -ge 1 ]]; then
    merge_pairwise_cometref "$WORKDIR/final/EPIC_comet_merged_sciencegti.fits" "${sci_files[@]}"
  fi
}

stage_qc() {
  need_cmd python3
  "$SCRIPT_DIR/xmm_comet_run_qc.sh" "$CONFIG_FILE"
}

run_stage() {
  local _stage="$1"
  echo ""
  echo "================================================================"
  echo "  STAGE: ${_stage}"
  echo "  $(date -u +%FT%TZ)"
  echo "================================================================"
  case "$_stage" in
    init) stage_init ;;
    repro) stage_repro ;;
    clean) stage_clean ;;
    track) stage_track ;;
    detect) stage_detect ;;
    comet) stage_comet ;;
    contam) stage_contam ;;
    image) stage_image ;;
    lcurve) stage_lcurve ;;
    spectrum) stage_spectrum ;;
    merge) stage_merge ;;
    qc) stage_qc ;;
    science)
      for _s in clean track detect comet contam image lcurve spectrum merge; do
        run_stage "$_s"
      done
      return
      ;;
    all)
      for _s in init repro clean track detect comet contam image lcurve spectrum merge qc; do
        run_stage "$_s"
      done
      return
      ;;
    *)
      echo "Unknown stage: $_stage" >&2
      exit 1
      ;;
  esac
  echo "--- ${_stage} done ($(date -u +%FT%TZ)) ---"
}

run_stage "$STAGE"
