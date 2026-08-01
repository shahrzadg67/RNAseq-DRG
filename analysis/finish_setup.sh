#!/usr/bin/env bash
# Wait for the main install to finish, SERIALLY repair any packages that failed
# to compile in the parallel pass, then smoke-test the analysis chain on
# synthetic data. Logs to logs/finish_setup.log.
set -uo pipefail
PROJ="/hpf/projects/msalter/sghazis/rnaseq_TUY35595"
cd "$PROJ"
export R_LIBS_USER="$PROJ/.Rlib"
export TUY_DEV_DATA="$PROJ/dev_data"
module load R/4.5.1 2>/dev/null || true

echo "[$(date +%T)] waiting for install_packages.R to finish…"
while pgrep -f install_packages.R >/dev/null; do sleep 30; done
echo "[$(date +%T)] main install finished."

echo "[$(date +%T)] SERIAL repair pass (Ncpus=1) for any missing packages…"
Rscript - <<'RS'
.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))
options(Ncpus = 1, repos = c(CRAN = "https://cloud.r-project.org"))
have <- function(p) requireNamespace(p, quietly = TRUE)
cran <- c("stringi","digest","fastmap","sourcetools","farver","isoband",
          "dplyr","tidyr","tibble","ggplot2","plotly","DT","shiny","bslib",
          "cowplot","scales","yaml","jsonlite","matrixStats","htmlwidgets",
          "purrr","stringr","readr","httpuv")
miss <- cran[!vapply(cran, have, logical(1))]
cat("Missing CRAN:", if(length(miss)) paste(miss, collapse=", ") else "none", "\n")
for (p in miss) { cat(">> installing", p, "\n"); try(install.packages(p)) }
if (!have("BiocManager")) install.packages("BiocManager")
bioc <- c("DESeq2","edgeR","limma","tximport","clusterProfiler","org.Mm.eg.db",
          "enrichplot","pathview","EnhancedVolcano","apeglm","fgsea","AnnotationDbi")
bmiss <- bioc[!vapply(bioc, have, logical(1))]
cat("Missing Bioc:", if(length(bmiss)) paste(bmiss, collapse=", ") else "none", "\n")
if (length(bmiss)) BiocManager::install(bmiss, update=FALSE, ask=FALSE, Ncpus=1)
if (!have("shinylive")) install.packages("shinylive")
# final report
allp <- c(cran, bioc, "shinylive")
st <- data.frame(pkg=allp, ok=vapply(allp, have, logical(1)))
print(st[!st$ok,,drop=FALSE]); if(all(st$ok)) cat("ALL PACKAGES OK\n")
RS

echo "[$(date +%T)] generating synthetic data…"
Rscript analysis/make_synthetic_data.R || { echo "synthetic-data FAILED"; exit 1; }

for s in 20_deseq2.R 30_edger.R 10_pca.R 40_assemble.R; do
  echo "[$(date +%T)] running $s …"
  Rscript "analysis/$s" || { echo "$s FAILED"; exit 1; }
done

echo "[$(date +%T)] artifacts:"; ls -lh app/data/*.rds 2>/dev/null
echo "[$(date +%T)] FINISH_SETUP OK"
