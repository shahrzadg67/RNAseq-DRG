#!/usr/bin/env bash
# Wait for the package install to finish, then smoke-test the analysis chain on
# synthetic data (DESeq2 -> edgeR -> PCA -> assemble). Verifies the code end-to-end
# before the real count matrices land. Logs to logs/smoke_test.log.
set -uo pipefail
PROJ="/hpf/projects/msalter/sghazis/rnaseq_TUY35595"
cd "$PROJ"
export R_LIBS_USER="$PROJ/.Rlib"
export TUY_DEV_DATA="$PROJ/dev_data"
module load R/4.5.1 2>/dev/null || true

echo "[$(date +%T)] waiting for install_packages.R to finish…"
while pgrep -f install_packages.R >/dev/null; do sleep 30; done
echo "[$(date +%T)] install finished. Checking key packages…"
Rscript -e 'for (p in c("DESeq2","edgeR","clusterProfiler","plotly","DT","shinylive","org.Mm.eg.db")) cat(sprintf("%-16s %s\n", p, requireNamespace(p, quietly=TRUE)))'

echo "[$(date +%T)] generating synthetic data…"
Rscript analysis/make_synthetic_data.R || { echo "synthetic-data FAILED"; exit 1; }

for s in 20_deseq2.R 30_edger.R 10_pca.R 40_assemble.R; do
  echo "[$(date +%T)] running $s …"
  Rscript "analysis/$s" || { echo "$s FAILED"; exit 1; }
done

echo "[$(date +%T)] artifacts produced:"
ls -lh app/data/*.rds 2>/dev/null
echo "[$(date +%T)] SMOKE TEST OK"
