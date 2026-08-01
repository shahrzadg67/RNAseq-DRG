#!/usr/bin/env Rscript
# ============================================================================
# 40_assemble.R  —  combine DESeq2 + edgeR into app-ready artifacts:
#   de_long.rds      : one row per feature x contrast x level x engine (stats)
#   expr.rds         : per-level VST matrix for ALL samples (incl SNI) + meta,
#                      for the gene/transcript expression boxplot lookup
#   search_index.rds : feature_id <-> symbol/gene_id/level for typeahead
# Cutoffs are NOT baked in — full per-feature stats are stored so the app can
# re-threshold live.
# ============================================================================
source("/hpf/projects/msalter/sghazis/rnaseq_TUY35595/analysis/_common.R")
suppressPackageStartupMessages({ library(DESeq2); library(matrixStats) })

de   <- read_cache("deseq2")
ed   <- read_cache("edger")
meta_all <- load_metadata()

## ---- 1. long DE table (both engines) — SLIM, app-only columns --------------
# Keep only what the volcano / table / concordance need; round to shrink the
# bundle (full stats remain in the cached intermediates + per-contrast CSVs).
collect <- function(src) {
  do.call(rbind, unlist(lapply(LEVELS, function(lv) src[[lv]]$results), recursive = FALSE))
}
de_long <- rbind(collect(de), collect(ed))
de_long <- data.frame(
  feature_id = de_long$feature_id,
  symbol     = de_long$symbol,
  log2FC     = round(de_long$log2FC, 3),
  padj       = signif(de_long$padj, 4),
  contrast   = factor(de_long$contrast),
  level      = factor(de_long$level),
  engine     = factor(de_long$engine),
  stringsAsFactors = FALSE)
rownames(de_long) <- NULL
saveRDS(de_long, file.path(APP_DATA, "de_long.rds"), compress = "xz")
message("de_long rows: ", nrow(de_long),
        "  size(MB): ", round(file.size(file.path(APP_DATA,"de_long.rds"))/1e6, 1))

## ---- 2. all-sample expression (incl SNI) for boxplots ----------------------
expr <- list()
for (lv in LEVELS) {
  cd  <- load_counts(lv)
  mat <- cd$mat
  meta <- droplevels(meta_all[colnames(mat), ])
  dds <- DESeqDataSetFromMatrix(mat, meta, design = ~ 1)
  dds <- dds[rowSums(counts(dds) >= 10) >= 2, ]
  v <- tryCatch(vst(dds, blind = TRUE),
                error = function(e) varianceStabilizingTransformation(dds, blind = TRUE))
  expr[[lv]] <- list(vst = assay(v), meta = meta)
}
saveRDS(expr, file.path(APP_DATA, "expr.rds"))

## ---- 3. search index (symbol / id typeahead) -------------------------------
idx <- do.call(rbind, lapply(LEVELS, function(lv) {
  r <- de[[lv]]$results[[1]]                     # any contrast has the full feature set
  data.frame(level = lv, feature_id = r$feature_id,
             gene_id = r$gene_id, symbol = r$symbol,
             label = ifelse(is.na(r$symbol), r$feature_id,
                            paste0(r$symbol, " (", r$feature_id, ")")),
             stringsAsFactors = FALSE)
}))
idx <- idx[!duplicated(paste(idx$level, idx$feature_id)), ]
saveRDS(idx, file.path(APP_DATA, "search_index.rds"))
message("search_index rows: ", nrow(idx))
message("assemble done -> de_long.rds, expr.rds, search_index.rds")
