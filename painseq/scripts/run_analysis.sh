#!/usr/bin/env bash
# Run one (or all) of the painseq analysis scripts with the correct module env.
#   bash scripts/run_analysis.sh 03            # run analysis/03_*.R
#   bash scripts/run_analysis.sh all           # run 01..04 in order
#
# The .extralib shim supplies libssl.so.10 / libcrypto.so.10, which R's curl
# package needs on this cluster. Without it, anything that touches curl (e.g.
# plotly widget bundling) fails with "libssl.so.10: cannot open shared object file".
set -euo pipefail
PROJ=/hpf/projects/msalter/sghazis/painseq
APP=/hpf/projects/msalter/sghazis/rnaseq_TUY35595

module load R/4.5.2 libarchive/3.8.1 2>/dev/null || true
export R_LIBS_USER="$APP/.Rlib"
export LD_LIBRARY_PATH="$APP/.extralib:${LD_LIBRARY_PATH:-}"
cd "$PROJ"

run_one() {
  local s
  s=$(ls analysis/"$1"_*.R 2>/dev/null | head -1)
  [ -n "$s" ] || { echo "no analysis script matching '$1'" >&2; exit 1; }
  echo "=== $s ==="
  Rscript "$s" 2>&1 | tee "logs/$(basename "${s%.R}").log"
}

if [ "${1:-all}" = "all" ]; then
  for n in 01 02 03 04; do run_one "$n"; done
else
  run_one "$1"
fi
