#!/usr/bin/env bash
# Wait for the REAL count matrices, then run the full analysis on them and
# rebuild the STAGING site (local only — never auto-deploys to prod).
# Logs to logs/real_analysis.log.
set -uo pipefail
PROJ="/hpf/projects/msalter/sghazis/rnaseq_TUY35595"
cd "$PROJ"
export R_LIBS_USER="$PROJ/.Rlib"
unset TUY_DEV_DATA                       # use the REAL results/star_salmon
module load R/4.5.1 libarchive/3.8.1 2>/dev/null || true

G="results/star_salmon/salmon.merged.gene_counts.tsv"
T="results/star_salmon/salmon.merged.transcript_counts.tsv"

echo "[$(date +%T)] waiting for count matrices…"
while [ ! -s "$G" ] || [ ! -s "$T" ]; do sleep 20; done
# let the merge finish writing (size stable for two checks)
prev=0; while :; do cur=$(stat -c%s "$G" 2>/dev/null)$(stat -c%s "$T" 2>/dev/null)
  [ "$cur" = "$prev" ] && break; prev=$cur; sleep 15; done
echo "[$(date +%T)] matrices present: $(wc -l < "$G") gene rows, $(wc -l < "$T") tx rows"

for s in 20_deseq2.R 30_edger.R 10_pca.R 40_assemble.R 50_gsea.R; do
  echo "[$(date +%T)] $s"; Rscript "analysis/$s" || { echo "$s FAILED"; exit 1; }
done

echo "[$(date +%T)] rebuilding staging site with REAL data…"
bash deploy/publish.sh staging || { echo "build FAILED"; exit 1; }
echo "[$(date +%T)] REAL ANALYSIS + STAGING BUILD DONE"