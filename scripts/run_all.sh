#!/usr/bin/env bash
set -euo pipefail

mkdir -p figures results

python scripts/create_figure6_clade_structure.py \
  --data-dir data \
  --outdir figures \
  --results-dir results \
  --basename Figure6_clade_structure \
  --formats png pdf svg \
  --dpi 600

python scripts/summarize_ha_na_clade_separation.py \
  --data-dir data \
  --outdir results

printf "Analysis complete. Outputs written to figures/ and results/.\n"
