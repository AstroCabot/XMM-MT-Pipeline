#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 1 ]]; then
  echo "Usage: $0 config.env" >&2
  exit 1
fi

CONFIG_FILE="$1"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$CONFIG_FILE"

: "${WORKDIR:?Set WORKDIR in config.env}"
mkdir -p "$WORKDIR/qc"

# Resolve a default spot-check band from IMAGE_BANDS unless overridden.
QC_SPOTCHECK_BAND="${QC_SPOTCHECK_BAND:-}"
if [[ -z "$QC_SPOTCHECK_BAND" ]]; then
  IFS=';' read -ra __bands <<< "${IMAGE_BANDS:-broad:300:2000}"
  for __entry in "${__bands[@]}"; do
    [[ -n "$__entry" ]] || continue
    IFS=':' read -r QC_SPOTCHECK_BAND _ _ <<< "$__entry"
    [[ -n "$QC_SPOTCHECK_BAND" ]] && break
  done
fi

QC_COMPARE_HORIZONS="${QC_COMPARE_HORIZONS:-no}"
HORIZONS_FLAG=()
case "${QC_COMPARE_HORIZONS,,}" in
  1|y|yes|true) HORIZONS_FLAG=(--horizons) ;;
  *) HORIZONS_FLAG=() ;;
esac

python3 "$SCRIPT_DIR/scripts/xmm_qc_manifest.py" manifest \
  --config "$CONFIG_FILE" \
  --outdir "$WORKDIR/qc"

# ---------------------------------------------------------------------------
# Quick comet-frame diagnostic (counts-only, no eexpmap — runs in seconds).
# Only requires comet/ stage outputs to exist.
# ---------------------------------------------------------------------------
if ls "$WORKDIR/comet/"*_comet.fits 1>/dev/null 2>&1; then
  # Parse PI range from first IMAGE_BANDS entry (e.g. "soft:200:1000")
  IFS=':' read -r _qd_name _qd_pimin _qd_pimax <<< "${IMAGE_BANDS%%\;*}"
  _qd_pimin="${_qd_pimin:-200}"
  _qd_pimax="${_qd_pimax:-2000}"

  : "${QC_QUICKDIAG_NBINS:=10}"
  : "${QC_QUICKDIAG_BIN_ARCSEC:=8.0}"

  python3 "$SCRIPT_DIR/scripts/xmm_comet_quick_diag.py" "$WORKDIR" \
    --pimin "$_qd_pimin" --pimax "$_qd_pimax" \
    --nbins "$QC_QUICKDIAG_NBINS" \
    --bin-arcsec "$QC_QUICKDIAG_BIN_ARCSEC" \
    --radius-arcsec "${IMAGE_RADIUS_ARCSEC:-900}" \
    --src-r-arcsec "${SRC_R_ARCSEC:-500}" || true
  echo "  Quick diagnostic: $WORKDIR/qc/quick_diag/"
else
  echo "  Skipping quick diagnostic (no comet-frame event lists yet)"
fi

python3 "$SCRIPT_DIR/scripts/xmm_comet_quick_checks.py" \
  --config "$CONFIG_FILE" \
  --outdir "$WORKDIR/qc" \
  "${HORIZONS_FLAG[@]}" \
  clean exposures track detect contam image lcurve spectrum

# Region-support reports help diagnose geometry-driven correction failures.
for inst in EPIC PN M1 M2; do
  python3 "$SCRIPT_DIR/scripts/xmm_comet_quick_checks.py" \
    --config "$CONFIG_FILE" \
    --outdir "$WORKDIR/qc" \
    --inst "$inst" \
    region-support || true
done

# Build still-sky and comet-frame mosaic products (evselect + eexpmap + Python).
"$SCRIPT_DIR/scripts/xmm_qc_build_mosaics.sh" "$CONFIG_FILE" || true

# QC visualisation of the mosaic products.
python3 "$SCRIPT_DIR/scripts/xmm_comet_quick_checks.py" \
  --config "$CONFIG_FILE" \
  --outdir "$WORKDIR/qc" \
  mosaics || true

# Comet-frame spot checks are the most important visual check of the moving-frame logic.
if [[ -f "$WORKDIR/track/motion_segments.csv" && -f "$WORKDIR/comet/moved_atthk.dat" ]]; then
  FORCE="${FORCE:-0}" "$SCRIPT_DIR/scripts/xmm_comet_spotcheck_frames.sh" "$CONFIG_FILE" "$QC_SPOTCHECK_BAND" || true
  python3 "$SCRIPT_DIR/scripts/xmm_comet_quick_checks.py" \
    --config "$CONFIG_FILE" \
    --outdir "$WORKDIR/qc" \
    --band "$QC_SPOTCHECK_BAND" \
    spotcheck || true
fi

python3 "$SCRIPT_DIR/scripts/xmm_qc_manifest.py" index \
  --config "$CONFIG_FILE" \
  --outdir "$WORKDIR/qc"

echo "QC outputs written under $WORKDIR/qc"
echo "Start here:"
echo "  $WORKDIR/qc/STAGE_QC.md"
echo "Key images:"
echo "  $WORKDIR/qc/00_manifest.png"
echo "  $WORKDIR/qc/00_clean.png"
echo "  $WORKDIR/qc/00b_exposures.png"
echo "  $WORKDIR/qc/01_track.png"
echo "  $WORKDIR/qc/02_detect.png"
echo "  $WORKDIR/qc/03_contam.png"
echo "  $WORKDIR/qc/04_image.png"
echo "  $WORKDIR/qc/05_lcurve.png"
echo "  $WORKDIR/qc/06_spectrum.png"
echo "  $WORKDIR/qc/07_spotchecks_${QC_SPOTCHECK_BAND}.png"
echo "  $WORKDIR/qc/08_mosaics.png"
echo "  $WORKDIR/qc/quick_diag/mosaic.png"
echo "  $WORKDIR/qc/quick_diag/comparison.png"
