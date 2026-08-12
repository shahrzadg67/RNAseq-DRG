#!/usr/bin/env Rscript
# ============================================================================
# 41_assemble_bc.R  —  assemble BATCH-CORRECTED app artifacts from the *_bc
# caches. Parallel to 40_assemble.R but every output gets a _bc suffix so the
# originals (de_long.rds, expr.rds, search_index.rds) are never touched.
#   de_long_bc.rds      : one row per feature x contrast x level x engine
#   expr_bc.rds         : per-level VST for the 12 NP/Sham samples, with the
#                         replicate batch effect REMOVED, for the boxplot lookup
#   search_index_bc.rds : feature_id <-> symbol/gene_id/level typeahead
# ============================================================================
source("/hpf/projects/msalter/sghazis/rnaseq_TUY35595/analysis/_common.R")
suppressPackageStartupMessages({ library(limma) })

de <- read_cache("deseq2_bc")
ed <- read_cache("edger_bc")

## ---- 1. long DE table (both engines) — SLIM, app-only columns --------------
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
saveRDS(de_long, file.path(APP_DATA, "de_long_bc.rds"), compress = "xz")
message("de_long_bc rows: ", nrow(de_long),
        "  size(MB): ", round(file.size(file.path(APP_DATA,"de_long_bc.rds"))/1e6, 1))

## ---- 2. batch-corrected expression (NP/Sham only) for boxplots -------------
# Reuse the 12-sample blind VST from the bc DESeq2 cache and strip the
# replicate effect (preserving group), so boxplots show corrected expression.
expr <- list()
for (lv in LEVELS) {
  v    <- de[[lv]]$vst
  meta <- de[[lv]]$meta
  meta$replicate <- factor(meta$replicate)
  d0 <- model.matrix(~ group, data = meta)
  v  <- limma::removeBatchEffect(v, batch = meta$replicate, design = d0)
  expr[[lv]] <- list(vst = v, meta = meta)
}
saveRDS(expr, file.path(APP_DATA, "expr_bc.rds"))
message("expr_bc done (", paste(LEVELS, collapse = ", "), ")")

## ---- 3. search index (symbol / id typeahead) ------------------------------
idx <- do.call(rbind, lapply(LEVELS, function(lv) {
  r <- de[[lv]]$results[[1]]                     # any contrast has the full feature set
  data.frame(level = lv, feature_id = r$feature_id,
             gene_id = r$gene_id, symbol = r$symbol,
             label = ifelse(is.na(r$symbol), r$feature_id,
                            paste0(r$symbol, " (", r$feature_id, ")")),
             stringsAsFactors = FALSE)
}))
idx <- idx[!duplicated(paste(idx$level, idx$feature_id)), ]
saveRDS(idx, file.path(APP_DATA, "search_index_bc.rds"))
message("search_index_bc rows: ", nrow(idx))
message("assemble (batch-corrected) done -> de_long_bc.rds, expr_bc.rds, search_index_bc.rds")
