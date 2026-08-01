#!/usr/bin/env Rscript
# Orchestrate the BATCH-CORRECTED downstream analysis (NP/Sham only, replicate
# modelled/removed as batch). Runs alongside — and never overwrites — the
# original analysis. Run after the pipeline has produced
# results/star_salmon/salmon.merged.*_counts.tsv (or with TUY_DEV_DATA set).
#   module load R/4.5.2 libarchive/3.8.1
#   export R_LIBS_USER=/hpf/projects/msalter/sghazis/rnaseq_TUY35595/.Rlib
#   Rscript analysis/run_all_bc.R
A <- "/hpf/projects/msalter/sghazis/rnaseq_TUY35595/analysis"
steps <- c("21_deseq2_bc.R", "31_edger_bc.R", "11_pca_bc.R", "41_assemble_bc.R", "51_gsea_bc.R")
for (s in steps) {
  message("\n############ ", s, " ############")
  t0 <- proc.time()[3]
  source(file.path(A, s), local = new.env())
  message(sprintf("   (%.1fs)", proc.time()[3] - t0))
}
message("\nALL DONE (batch-corrected). Artifacts: app/data/*_bc.rds and app/www/{gsea_bc,pathview_bc}/.")
