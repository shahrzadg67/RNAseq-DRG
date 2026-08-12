#!/usr/bin/env Rscript
# ============================================================================
# 14_pca_uncorr_noxy.R  —  UNCORRECTED PCA with the SEX CHROMOSOMES REMOVED.
# Same recipe as 10_pca.R (VST -> top-1000 variable -> prcomp, NO batch
# correction), but chrX/chrY features are dropped from the count matrix first
# (list from analysis/cache/xy_features.rds). Computes gene & transcript ×
# withSNI & noSNI, and APPENDS them to app/data/pca.rds under keys
# "<level>_<sni>_noXY" so the uncorrected PCA tab can offer a "remove chrX/Y"
# toggle. The original 4 variants in pca.rds are kept unchanged.
# ============================================================================
source("/hpf/projects/msalter/sghazis/rnaseq_TUY35595/analysis/_common.R")
suppressPackageStartupMessages({ library(DESeq2); library(matrixStats) })

meta_all <- load_metadata()
xy   <- read_cache("xy_features")        # list(gene=<ids>, transcript=<ids>)
NTOP <- 1000; NPC <- 10

vst_matrix_noxy <- function(level, samples) {
  cd  <- load_counts(level)
  mat <- cd$mat[!rownames(cd$mat) %in% xy[[level]], , drop = FALSE]   # drop chrX/Y
  mat <- mat[, intersect(samples, colnames(mat)), drop = FALSE]
  meta <- droplevels(meta_all[colnames(mat), ])
  dds <- DESeqDataSetFromMatrix(mat, meta, design = ~ 1)
  dds <- dds[rowSums(counts(dds) >= 10) >= 2, ]
  v <- tryCatch(vst(dds, blind = TRUE),
                error = function(e) varianceStabilizingTransformation(dds, blind = TRUE))
  list(vst = assay(v), meta = meta)
}

run_pca <- function(level, include_sni) {          # UNCORRECTED (no removeBatchEffect)
  samples <- rownames(meta_all)
  if (!include_sni) samples <- samples[meta_all$condition != "SNI"]
  vm <- vst_matrix_noxy(level, samples); x <- vm$vst
  x  <- x[order(rowVars(x), decreasing = TRUE)[seq_len(min(NTOP, nrow(x)))], ]
  pca <- prcomp(t(x), center = TRUE, scale. = FALSE); k <- min(NPC, ncol(pca$x))
  ve <- (pca$sdev^2 / sum(pca$sdev^2))[seq_len(k)]
  coords <- data.frame(sample_id = rownames(pca$x),
                       vm$meta[rownames(pca$x), c("condition","sex","group","replicate")],
                       pca$x[, seq_len(k), drop = FALSE], row.names = NULL, check.names = FALSE)
  list(coords = coords,
       var_explained = data.frame(PC = paste0("PC", seq_len(k)),
                                   pct = round(100 * ve, 2), cumpct = round(100 * cumsum(ve), 2)),
       level = level, include_sni = include_sni, n_features = nrow(x))
}

pca <- readRDS(file.path(APP_DATA, "pca.rds"))     # existing 4 uncorrected variants
for (lv in LEVELS) for (sni in c(TRUE, FALSE)) {
  key <- sprintf("%s_%s_noXY", lv, if (sni) "withSNI" else "noSNI")
  message("uncorrected noXY PCA :: ", key)
  pca[[key]] <- run_pca(lv, sni)
}
saveRDS(pca, file.path(APP_DATA, "pca.rds"))
message("pca.rds now has keys: ", paste(names(pca), collapse = ", "))
