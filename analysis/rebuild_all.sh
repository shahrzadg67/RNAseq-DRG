#!/usr/bin/env bash
# Full rebuild once the ~/.R/Makevars compiler/linker fix is in place:
# (re)install the whole R stack, then smoke-test the analysis chain on synthetic
# data. Safe to re-run — install steps skip already-installed packages.
# Logs to logs/rebuild_all.log.
set -uo pipefail
PROJ="/hpf/projects/msalter/sghazis/rnaseq_TUY35595"
cd "$PROJ"
export R_LIBS_USER="$PROJ/.Rlib"
export TUY_DEV_DATA="$PROJ/dev_data"
module load R/4.5.1 2>/dev/null || true

echo "[$(date +%T)] CXX check: $(R CMD config CXX)"
echo "[$(date +%T)] installing full stack (CRAN + Bioc + shinylive)…"
Rscript analysis/install_packages.R || echo "install_packages returned non-zero (will verify below)"

echo "[$(date +%T)] verifying required packages…"
Rscript -e '.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))
need <- c("DESeq2","edgeR","limma","tximport","clusterProfiler","org.Mm.eg.db",
          "enrichplot","pathview","EnhancedVolcano","fgsea","AnnotationDbi",
          "ggplot2","dplyr","tidyr","plotly","DT","shiny","bslib","shinylive",
          "stringi","fs","digest")
ok <- vapply(need, requireNamespace, logical(1), quietly=TRUE)
cat(sprintf("%-16s %s\n", need, ok))
miss <- need[!ok]; if (length(miss)) { cat("MISSING:", paste(miss, collapse=", "), "\n"); quit(status=2) }
cat("ALL REQUIRED PACKAGES OK\n")' || { echo "PACKAGES STILL MISSING — stopping before smoke test"; exit 2; }

echo "[$(date +%T)] generating synthetic data + running analysis…"
Rscript analysis/make_synthetic_data.R || { echo "synthetic-data FAILED"; exit 1; }
for s in 20_deseq2.R 30_edger.R 10_pca.R 40_assemble.R; do
  echo "[$(date +%T)] $s"; Rscript "analysis/$s" || { echo "$s FAILED"; exit 1; }
done
echo "[$(date +%T)] artifacts:"; ls -lh app/data/*.rds 2>/dev/null
echo "[$(date +%T)] REBUILD_ALL OK"
