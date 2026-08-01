#!/usr/bin/env Rscript
# Orchestrate the full downstream analysis. Run after the pipeline produces
# results/star_salmon/salmon.merged.*_counts.tsv (or with TUY_DEV_DATA set).
#   module load R/4.5.1
#   export R_LIBS_USER=/hpf/projects/msalter/sghazis/rnaseq_TUY35595/.Rlib
#   Rscript analysis/run_all.R
A <- "/hpf/projects/msalter/sghazis/rnaseq_TUY35595/analysis"
steps <- c("20_deseq2.R", "30_edger.R", "10_pca.R", "40_assemble.R", "50_gsea.R")
for (s in steps) {
  message("\n############ ", s, " ############")
  t0 <- proc.time()[3]
  source(file.path(A, s), local = new.env())
  message(sprintf("   (%.1fs)", proc.time()[3] - t0))
}
message("\nALL DONE. Artifacts in app/data/ and app/www/.")
