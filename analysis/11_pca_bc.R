#!/usr/bin/env Rscript
# ============================================================================
# 11_pca_bc.R  —  BATCH-CORRECTED PCA (gene + transcript level), NP/Sham only.
# SNI is dropped entirely. The replicate batch effect is removed from the VST
# matrix with limma::removeBatchEffect BEFORE top-variable-feature selection
# and prcomp. The `design = model.matrix(~ group)` argument tells
# removeBatchEffect which biological structure to PRESERVE (group encodes
# condition x sex) while it strips the replicate-attributable variance.
# Output keys are plain "gene"/"transcript" (no SNI suffix); each variant
# matches the existing pca.rds per-variant schema. Writes app/data/pca_bc.rds
# — the original pca.rds (4 with/without-SNI variants) is untouched.
#
# NOTE: removeBatchEffect output is for VISUALISATION ONLY. Differential
# expression models the batch internally via the design term (21/31_*_bc.R);
# corrected values are never fed into DE.
# ============================================================================
source("/hpf/projects/msalter/sghazis/rnaseq_TUY35595/analysis/_common.R")
suppressPackageStartupMessages({ library(DESeq2); library(matrixStats); library(limma) })

meta_all <- load_metadata()
NTOP  <- 1000      # top variable features for PCA (à la DESeq2::plotPCA)
NPC   <- 10        # how many PCs to export

vst_matrix <- function(level, samples) {
  cd  <- load_counts(level)
  mat <- cd$mat[, intersect(samples, colnames(cd$mat)), drop = FALSE]
  meta <- droplevels(meta_all[colnames(mat), ])
  dds <- DESeqDataSetFromMatrix(mat, meta, design = ~ 1)
  dds <- dds[rowSums(counts(dds) >= 10) >= 2, ]
  v <- tryCatch(vst(dds, blind = TRUE),
                error = function(e) varianceStabilizingTransformation(dds, blind = TRUE))
  list(vst = assay(v), meta = meta)
}

run_pca_bc <- function(level) {
  samples <- rownames(meta_all)[meta_all$condition != "SNI"]   # NP + Sham only
  vm   <- vst_matrix(level, samples)
  x    <- vm$vst
  meta <- vm$meta
  meta$replicate <- factor(meta$replicate)                     # CRITICAL: batch as factor
  # remove replicate batch effect, preserving group (= condition x sex) structure
  d0 <- model.matrix(~ group, data = meta)
  x  <- limma::removeBatchEffect(x, batch = meta$replicate, design = d0)

  rv <- rowVars(x)
  x  <- x[order(rv, decreasing = TRUE)[seq_len(min(NTOP, nrow(x)))], ]
  pca <- prcomp(t(x), center = TRUE, scale. = FALSE)
  k <- min(NPC, ncol(pca$x))
  ve <- (pca$sdev^2 / sum(pca$sdev^2))[seq_len(k)]
  coords <- data.frame(sample_id = rownames(pca$x),
                       meta[rownames(pca$x), c("condition","sex","group","replicate")],
                       pca$x[, seq_len(k), drop = FALSE], row.names = NULL,
                       check.names = FALSE)
  list(coords = coords,
       var_explained = data.frame(PC = paste0("PC", seq_len(k)),
                                   pct = round(100 * ve, 2),
                                   cumpct = round(100 * cumsum(ve), 2)),
       level = level, include_sni = FALSE, n_features = nrow(x))
}

pca <- list()
for (lv in LEVELS) {
  message("PCA (batch-corrected) :: ", lv)
  pca[[lv]] <- run_pca_bc(lv)
}
save_artifact(pca, "pca_bc")
message("PCA (batch-corrected) done -> app/data/pca_bc.rds  (keys: ",
        paste(names(pca), collapse = ", "), ")")
