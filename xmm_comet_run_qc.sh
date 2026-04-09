#!/usr/bin/env bash
set -euo pipefail
[[ $# -ge 1 ]] || { echo "Usage: $0 config.env" >&2; exit 1; }
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
python3 "$SCRIPT_DIR/qc/runner.py" "$1"
