#!/usr/bin/env Rscript
# ============================================================================
# 10_pca.R  —  PCA + scree/elbow for gene & transcript level, computed BOTH
# with SNI and without SNI. Exports tidy PC coordinates (PC1..PCk per sample)
# and variance-explained so the app can offer a free PC-pair selector and an
# elbow plot without recomputing anything.
# ============================================================================
source("/hpf/projects/msalter/sghazis/rnaseq_TUY35595/analysis/_common.R")
suppressPackageStartupMessages({ library(DESeq2); library(matrixStats) })

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

run_pca <- function(level, include_sni) {
  samples <- rownames(meta_all)
  if (!include_sni) samples <- samples[meta_all$condition != "SNI"]
  vm <- vst_matrix(level, samples)
  x  <- vm$vst
  rv <- rowVars(x)
  x  <- x[order(rv, decreasing = TRUE)[seq_len(min(NTOP, nrow(x)))], ]
  pca <- prcomp(t(x), center = TRUE, scale. = FALSE)
  k <- min(NPC, ncol(pca$x))
  ve <- (pca$sdev^2 / sum(pca$sdev^2))[seq_len(k)]
  coords <- data.frame(sample_id = rownames(pca$x),
                       vm$meta[rownames(pca$x), c("condition","sex","group","replicate")],
                       pca$x[, seq_len(k), drop = FALSE], row.names = NULL,
                       check.names = FALSE)
  list(coords = coords,
       var_explained = data.frame(PC = paste0("PC", seq_len(k)),
                                   pct = round(100 * ve, 2),
                                   cumpct = round(100 * cumsum(ve), 2)),
       level = level, include_sni = include_sni, n_features = nrow(x))
}

pca <- list()
for (lv in LEVELS) for (sni in c(TRUE, FALSE)) {
  key <- sprintf("%s_%s", lv, if (sni) "withSNI" else "noSNI")
  message("PCA :: ", key)
  pca[[key]] <- run_pca(lv, sni)
}
save_artifact(pca, "pca")
message("PCA done -> app/data/pca.rds  (keys: ", paste(names(pca), collapse=", "), ")")
