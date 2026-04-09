#!/bin/bash
set -e
BASEDIR="/mnt/c/Users/cabot/OneDrive/Documents/Research/Cambridge/X3I/reduction_v4"
SCRIPTS="$BASEDIR/scripts"
cd "$BASEDIR/output"
export PYTHONPATH="$SCRIPTS:$PYTHONPATH"
python "$SCRIPTS/xmm_comet_quick_checks.py" --config "$BASEDIR/my_comet.env" "$@"
