#!/usr/bin/env bash
set -euo pipefail

###############################################################################
# Build per-instrument still-sky and comet-frame QC mosaic products.
#
# Produces:
#   $WORKDIR/qc/mosaics/stillsky/{band}/{inst}_{counts,exp}.fits
#   $WORKDIR/qc/mosaics/comet/{band}/{inst}_{counts,exp,rate,rate_masked}.fits
#
# Uses only SAS commands (evselect, eexpmap) and the combine helper.
#
# Usage:
#   bash xmm_qc_build_mosaics.sh config.env
#
# Requires: pipeline stages init, repro, clean, track, detect, comet, contam,
#           and image to have completed first.
###############################################################################

if [[ $# -lt 1 ]]; then
  echo "Usage: $0 config.env" >&2
  exit 1
fi

CONFIG_FILE="$1"
source "$CONFIG_FILE"
: "${WORKDIR:?Set WORKDIR in config.env}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
FORCE="${FORCE:-0}"
ONLY_INSTS="${ONLY_INSTS:-PN M1 M2}"

source_if_exists() { [[ -f "$1" ]] && source "$1"; }
need_cmd() { command -v "$1" >/dev/null 2>&1 || { echo "Missing: $1" >&2; exit 1; }; }
skip_file() { [[ "$FORCE" != "1" && -s "$1" ]]; }

need_cmd evselect
need_cmd eexpmap
need_cmd python3

# Load SAS env
source_if_exists "$WORKDIR/sas_setup.env"
[[ -n "${SAS_CCF:-}" && -n "${SAS_ODF:-}" ]] || { echo "Missing SAS setup" >&2; exit 1; }
export SAS_CCF="$WORKDIR/init/ccf.cif"
SUMSAS="$(find "$WORKDIR/init" -maxdepth 1 -type f -name '*SUM.SAS' | sort | head -n 1 || true)"
[[ -n "$SUMSAS" ]] && export SAS_ODF="$SUMSAS"

# Load attitude files
source_if_exists "$WORKDIR/attitude.env"
ORIG_ATTHK="$WORKDIR/repro/atthk.dat"
MOVED_ATTHK="$WORKDIR/comet/moved_atthk.dat"
[[ -f "$ORIG_ATTHK" ]] || { echo "Missing original attitude: $ORIG_ATTHK" >&2; exit 1; }

# Load image grid
source_if_exists "$WORKDIR/images/grid.env"
[[ -n "${BIN_PHYS:-}" ]] || { echo "Missing grid.env; run image stage first" >&2; exit 1; }

# Resolve bands from DETECT_QC_BANDS (soft+hard by default for still-sky)
parse_bands() {
  local spec="$1"
  local entry
  IFS=';' read -ra __bands <<< "$spec"
  for entry in "${__bands[@]}"; do
    [[ -n "$entry" ]] || continue
    IFS=':' read -r label pimin pimax <<< "$entry"
    [[ -n "$label" && -n "$pimin" && -n "$pimax" ]] || continue
    echo "$label $pimin $pimax"
  done
}

MOSAICS_ROOT="$WORKDIR/qc/mosaics"
STILLSKY_ROOT="$MOSAICS_ROOT/stillsky"
COMET_ROOT="$MOSAICS_ROOT/comet"

###############################################################################
# Part 1: Still-sky per-instrument count images and exposure maps
###############################################################################
echo "=== Building still-sky per-instrument mosaics ==="

while read -r label pimin pimax; do
  [[ -n "$label" ]] || continue
  banddir="$STILLSKY_ROOT/$label"
  mkdir -p "$banddir"

  counts_list=()
  expos_list=()

  for inst in $ONLY_INSTS; do
    # Use merged cleaned event list (still-sky frame)
    evt="$WORKDIR/clean/${inst}_clean_merged.fits"
    [[ -f "$evt" ]] || continue

    img="$banddir/${inst}_counts.fits"
    exp="$banddir/${inst}_exp.fits"

    # Count image in still-sky frame
    if ! skip_file "$img"; then
      evselect table="${evt}:EVENTS" \
        withimageset=yes imageset="$img" \
        xcolumn=X ycolumn=Y imagebinning=binSize \
        ximagebinsize="$BIN_PHYS" yimagebinsize="$BIN_PHYS" \
        withxranges=yes ximagemin="$X_MIN_PHYS" ximagemax="$X_MAX_PHYS" \
        withyranges=yes yimagemin="$Y_MIN_PHYS" yimagemax="$Y_MAX_PHYS" \
        writedss=yes \
        expression="PI in [$pimin:$pimax]"
    fi

    # Exposure map with ORIGINAL (non-moved) attitude → still-sky frame
    if ! skip_file "$exp"; then
      eexpmap imageset="$img" \
        attitudeset="$ORIG_ATTHK" \
        eventset="$evt" \
        expimageset="$exp" \
        pimin="$pimin" pimax="$pimax" \
        attrebin="${EEXPMAP_ATTREBIN:-0.020626481}"
    fi

    counts_list+=("$img")
    expos_list+=("$exp")
  done

  # Combine into EPIC totals
  if [[ ${#counts_list[@]} -gt 0 ]]; then
    if ! skip_file "$banddir/EPIC_counts.fits"; then
      python3 "$SCRIPT_DIR/combine_epic_images.py" \
        --counts "${counts_list[@]}" \
        --exposure "${expos_list[@]}" \
        --out-counts "$banddir/EPIC_counts.fits" \
        --out-exposure "$banddir/EPIC_exp.fits" \
        --out-rate "$banddir/EPIC_rate.fits" \
        --out-json "$banddir/EPIC_summary.json"
    fi
  fi

  echo "  still-sky $label: ${#counts_list[@]} instruments"
done < <(parse_bands "${DETECT_QC_BANDS:-soft:200:1000;hard:1000:12000}")

###############################################################################
# Part 2: Comet-frame per-instrument images (counts, exposure, rate)
#          These mostly already exist in WORKDIR/images — we produce rate maps
#          and masked rate maps per instrument.
###############################################################################
echo "=== Building comet-frame per-instrument rate images ==="

MASK_SRCLIST="$WORKDIR/detect/field_sources_all.fits"
[[ -f "$MASK_SRCLIST" ]] || MASK_SRCLIST="$WORKDIR/detect/field_sources_curated.fits"

while read -r label pimin pimax; do
  [[ -n "$label" ]] || continue
  src_banddir="$WORKDIR/images/$label"
  dstdir="$COMET_ROOT/$label"
  mkdir -p "$dstdir"

  counts_list=()
  expos_list=()

  for inst in $ONLY_INSTS; do
    cimg="$src_banddir/${inst}_counts.fits"
    eimg="$src_banddir/${inst}_exp.fits"
    [[ -f "$cimg" && -f "$eimg" ]] || continue

    # Symlink the existing per-inst counts and exposure into our QC tree
    ln -sfn "$(readlink -f "$cimg")" "$dstdir/${inst}_counts.fits"
    ln -sfn "$(readlink -f "$eimg")" "$dstdir/${inst}_exp.fits"

    counts_list+=("$dstdir/${inst}_counts.fits")
    expos_list+=("$dstdir/${inst}_exp.fits")
  done

  # Build per-instrument rate + masked rate images
  if [[ ${#counts_list[@]} -gt 0 ]]; then
    combine_args=(
      --counts "${counts_list[@]}"
      --exposure "${expos_list[@]}"
      --out-counts "$dstdir/EPIC_counts.fits"
      --out-exposure "$dstdir/EPIC_exp.fits"
      --out-rate "$dstdir/EPIC_rate.fits"
      --out-json "$dstdir/EPIC_summary.json"
    )
    if [[ -f "$MASK_SRCLIST" && -f "$WORKDIR/images/grid.json" && -f "$WORKDIR/track/comet_track.fits" ]]; then
      combine_args+=(
        --grid-json "$WORKDIR/images/grid.json"
        --track "$WORKDIR/track/comet_track.fits"
        --srclist "$MASK_SRCLIST"
        --mask-radius-arcsec "${IMAGE_MASK_RADIUS_ARCSEC:-${FIELD_SOURCE_MASK_R_ARCSEC:-20.0}}"
        --out-mask "$dstdir/EPIC_trail_mask.fits"
        --out-rate-masked "$dstdir/EPIC_rate_masked.fits"
      )
    fi
    if ! skip_file "$dstdir/EPIC_rate.fits"; then
      python3 "$SCRIPT_DIR/combine_epic_images.py" "${combine_args[@]}"
    fi
  fi

  echo "  comet $label: ${#counts_list[@]} instruments"
done < <(parse_bands "${IMAGE_BANDS:-soft:200:1000;broad:300:2000;hard:1000:12000}")

###############################################################################
# Part 3: Let the Python helper do source masking, background estimation,
#          and per-instrument rate computation.
###############################################################################
echo "=== Building masked products and background maps ==="

python3 "$SCRIPT_DIR/build_qc_mosaics.py" \
  --config "$CONFIG_FILE" \
  --stillsky-root "$STILLSKY_ROOT" \
  --comet-root "$COMET_ROOT" \
  --detect-srclist "$MASK_SRCLIST" \
  --grid-json "$WORKDIR/images/grid.json" \
  --track "$WORKDIR/track/comet_track.fits" \
  --mask-radius-arcsec "${IMAGE_MASK_RADIUS_ARCSEC:-${FIELD_SOURCE_MASK_R_ARCSEC:-20.0}}"

echo "=== QC mosaics complete ==="
echo "  Still-sky: $STILLSKY_ROOT"
echo "  Comet:     $COMET_ROOT"
