#!/usr/bin/env bash
# Run both segment builds. Study paths and parameters are read from JSON.
set -euo pipefail
ROOT="$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)"
cd "$ROOT"
CONFIG="${1:-config/analysis.json}"
THREADS="${2:-4}"
python scripts/run_nextstrain.py --config "$CONFIG" --segment HA --threads "$THREADS"
python scripts/run_nextstrain.py --config "$CONFIG" --segment NA --threads "$THREADS"
printf 'HA and NA builds complete. iTOL upload and figure styling are performed separately.\n'
