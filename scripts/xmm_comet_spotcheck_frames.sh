#!/usr/bin/env bash
set -euo pipefail

###############################################################################
# Build a small set of SAS-corrected comet-frame spot-check frames.
#
# Purpose:
#   generate a few representative short-interval counts / exposure / rate
#   images so you can visually check that the comet-frame attitude transform,
#   exposure correction, and stationary-source trail masking behave as expected.
#
# Usage:
#   ./xmm_comet_spotcheck_frames.sh config.env [band_label]
#
# Resume behaviour:
#   existing outputs are reused unless FORCE=1 is set.
#
# Optional environment overrides:
#   FORCE=1
#   ONLY_INSTS="PN M1"
#   SPOTCHECK_SEGMENTS="1,17,33"   # explicit segids from motion_segments.csv
###############################################################################

if [[ $# -lt 1 ]]; then
  echo "Usage: $0 config.env [band_label]"
  exit 1
fi

CONFIG_FILE="$1"
BAND_LABEL="${2:-}"
source "$CONFIG_FILE"

: "${WORKDIR:?Set WORKDIR in config.env}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
FORCE="${FORCE:-0}"
ONLY_INSTS="${ONLY_INSTS:-PN M1 M2}"

source_if_exists() {
  local f="$1"
  [[ -f "$f" ]] && source "$f"
}

need_cmd() {
  command -v "$1" >/dev/null 2>&1 || { echo "Missing required command: $1" >&2; exit 1; }
}

skip_file() {
  local f="$1"
  [[ "$FORCE" != "1" && -s "$f" ]]
}

first_existing() {
  local f
  for f in "$@"; do
    if [[ -f "$f" ]]; then
      echo "$f"
      return 0
    fi
  done
  echo ""
}

# Need SAS env and already-built pipeline products.
source_if_exists "$WORKDIR/sas_setup.env"
[[ -n "${SAS_CCF:-}" && -n "${SAS_ODF:-}" ]] || { echo "Missing SAS setup; run pipeline init/repro first" >&2; exit 1; }

MOVED_ATTHK="$WORKDIR/comet/moved_atthk.dat"
[[ -f "$MOVED_ATTHK" ]] || { echo "Missing $MOVED_ATTHK; run pipeline comet stage first" >&2; exit 1; }
[[ -f "$WORKDIR/track/motion_segments.csv" ]] || { echo "Missing $WORKDIR/track/motion_segments.csv; run pipeline track stage first" >&2; exit 1; }

need_cmd evselect
need_cmd eexpmap
need_cmd python3

# Resolve image band.
python3 - "$CONFIG_FILE" "$BAND_LABEL" <<'PY' > "$WORKDIR/qc_band.env"
import os, shlex, subprocess, sys
cfg, band = sys.argv[1:3]
cmd = ["bash", "-lc", f"set -a; source {shlex.quote(cfg)} >/dev/null 2>&1; env -0"]
proc = subprocess.run(cmd, capture_output=True, check=True)
env = {}
for chunk in proc.stdout.split(b"\x00"):
    if not chunk or b"=" not in chunk:
        continue
    k, v = chunk.split(b"=", 1)
    env[k.decode()] = v.decode()
image_bands = env.get("IMAGE_BANDS", "broad:300:2000")
bands = []
for entry in image_bands.split(";"):
    entry = entry.strip()
    if not entry:
        continue
    label, pimin, pimax = entry.split(":")
    bands.append((label.strip(), int(pimin), int(pimax)))
if not bands:
    raise RuntimeError("No IMAGE_BANDS configured")
if band:
    matches = [b for b in bands if b[0] == band]
    if not matches:
        raise RuntimeError(f"Requested band label {band!r} not found in IMAGE_BANDS")
    label, pimin, pimax = matches[0]
else:
    label, pimin, pimax = bands[0]
print(f'SPOT_BAND_LABEL="{label}"')
print(f'SPOT_PIMIN="{pimin}"')
print(f'SPOT_PIMAX="{pimax}"')
PY
source "$WORKDIR/qc_band.env"

# Ensure image grid exists (reuse science pipeline grid if already present).
mkdir -p "$WORKDIR/qc/spotchecks/$SPOT_BAND_LABEL"
if [[ ! -f "$WORKDIR/images/grid.env" || ! -f "$WORKDIR/images/grid.json" ]]; then
  REF_EVT="$(first_existing "$WORKDIR/comet/PN_comet.fits" "$WORKDIR/comet/M1_comet.fits" "$WORKDIR/comet/M2_comet.fits")"
  [[ -n "$REF_EVT" ]] || { echo "No comet-frame event list found" >&2; exit 1; }
  python3 "$SCRIPT_DIR/xmm_event_geometry.py" image-grid \
    --event "$REF_EVT" \
    --radius-arcsec "${IMAGE_RADIUS_ARCSEC:-300.0}" \
    --bin-arcsec "${IMAGE_BIN_ARCSEC:-2.0}" \
    --out-json "$WORKDIR/images/grid.json" > "$WORKDIR/images/grid.env"
fi
source "$WORKDIR/images/grid.env"

# Choose representative segments: first, middle, last by default.
python3 - "$WORKDIR/track/motion_segments.csv" "${SPOTCHECK_SEGMENTS:-}" <<'PY' > "$WORKDIR/qc_spot_segments.tsv"
import csv, sys
csv_path, manual = sys.argv[1:3]
rows = []
with open(csv_path, 'r', newline='', encoding='utf-8') as fh:
    for row in csv.DictReader(fh):
        rows.append(row)
if not rows:
    raise RuntimeError("motion_segments.csv is empty")
if manual.strip():
    wanted = {int(x.strip()) for x in manual.split(',') if x.strip()}
    keep = [r for r in rows if int(float(r['segid'])) in wanted]
else:
    idxs = sorted(set([0, len(rows)//2, len(rows)-1]))
    keep = [rows[i] for i in idxs]
for r in keep:
    print(f"{int(float(r['segid']))}\t{float(r['tstart']):.6f}\t{float(r['tstop']):.6f}\t{r.get('utc_start','')}\t{r.get('utc_stop','')}")
PY

while IFS=$'\t' read -r segid tstart tstop utc_start utc_stop; do
  [[ -n "$segid" ]] || continue
  segdir="$WORKDIR/qc/spotchecks/$SPOT_BAND_LABEL/seg${segid}"
  mkdir -p "$segdir"

  meta="$segdir/segment_info.txt"
  if [[ ! -f "$meta" || "$FORCE" == "1" ]]; then
    cat > "$meta" <<EOF
segid=${segid}
tstart=${tstart}
tstop=${tstop}
utc_start=${utc_start}
utc_stop=${utc_stop}
band_label=${SPOT_BAND_LABEL}
pimin=${SPOT_PIMIN}
pimax=${SPOT_PIMAX}
EOF
  fi

  counts_list=()
  expos_list=()

  for inst in $ONLY_INSTS; do
    evt="$WORKDIR/comet/${inst}_comet.fits"
    [[ -f "$evt" ]] || continue
    img="$segdir/${inst}_counts.fits"
    exp="$segdir/${inst}_exp.fits"

    if ! skip_file "$img"; then
      evselect table="${evt}:EVENTS" \
        withimageset=yes imageset="$img" \
        xcolumn=X ycolumn=Y imagebinning=binSize \
        ximagebinsize="$BIN_PHYS" yimagebinsize="$BIN_PHYS" \
        withxranges=yes ximagemin="$X_MIN_PHYS" ximagemax="$X_MAX_PHYS" \
        withyranges=yes yimagemin="$Y_MIN_PHYS" yimagemax="$Y_MAX_PHYS" \
        writedss=yes \
        expression="(TIME in [${tstart}:${tstop}])&&(PI in [${SPOT_PIMIN}:${SPOT_PIMAX}])"
    fi

    if ! skip_file "$exp"; then
      eexpmap imageset="$img" \
        attitudeset="$MOVED_ATTHK" \
        eventset="$evt" \
        expimageset="$exp" \
        pimin="$SPOT_PIMIN" pimax="$SPOT_PIMAX" \
        attrebin="${EEXPMAP_ATTREBIN:-0.020626481}"
    fi

    counts_list+=("$img")
    expos_list+=("$exp")
  done

  if [[ ${#counts_list[@]} -eq 0 ]]; then
    echo "No comet event files available for seg${segid}; skipping" >&2
    continue
  fi

  combine_args=(
    --counts "${counts_list[@]}"
    --exposure "${expos_list[@]}"
    --out-counts "$segdir/EPIC_counts.fits"
    --out-exposure "$segdir/EPIC_exp.fits"
    --out-rate "$segdir/EPIC_rate.fits"
  )
  MASK_SRCLIST="$WORKDIR/detect/field_sources_all.fits"
  [[ -f "$MASK_SRCLIST" ]] || MASK_SRCLIST="$WORKDIR/detect/field_sources_curated.fits"
  if [[ -f "$MASK_SRCLIST" ]]; then
    combine_args+=(
      --grid-json "$WORKDIR/images/grid.json"
      --track "$WORKDIR/track/comet_track.fits"
      --track-tmin-sec "$tstart"
      --track-tmax-sec "$tstop"
      --srclist "$MASK_SRCLIST"
      --mask-radius-arcsec "${IMAGE_MASK_RADIUS_ARCSEC:-${FIELD_SOURCE_MASK_R_ARCSEC:-20.0}}"
      --out-mask "$segdir/EPIC_trail_mask.fits"
      --out-rate-masked "$segdir/EPIC_rate_masked.fits"
      --out-json "$segdir/EPIC_image_summary.json"
    )
  fi
  if ! skip_file "$segdir/EPIC_rate.fits"; then
    python3 "$SCRIPT_DIR/combine_epic_images.py" "${combine_args[@]}"
  fi

done < "$WORKDIR/qc_spot_segments.tsv"

cat <<EOF
Spot-check frames written under:
  $WORKDIR/qc/spotchecks/$SPOT_BAND_LABEL
Representative segments came from:
  $WORKDIR/track/motion_segments.csv
EOF
