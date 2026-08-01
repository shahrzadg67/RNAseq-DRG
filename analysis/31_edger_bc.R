#!/usr/bin/env Rscript
# ============================================================================
# 31_edger_bc.R  —  BATCH-CORRECTED edgeR/limma-voom differential expression.
# Identical to 30_edger.R except the design adds `replicate` as an additive
# batch covariate:   model.matrix(~ 0 + group + replicate)
# Only the 4 group columns are renamed to the group levels; the 2 replicate
# columns keep their names. The contrast matrix is built zero-padded over the
# FULL set of design columns (the old `contrast_vec(cn)[colnames(design)]`
# would return NA for the replicate columns and corrupt the contrasts).
# Writes a SEPARATE cache (edger_bc) — the original edger.rds is untouched.
# ============================================================================
source("/hpf/projects/msalter/sghazis/rnaseq_TUY35595/analysis/_common.R")
suppressPackageStartupMessages({ library(edgeR); library(limma) })

meta_all <- load_metadata()

run_level <- function(level) {
  message("\n=== edgeR/voom (batch-corrected) :: ", level, " level ===")
  cd  <- load_counts(level)
  mat <- cd$mat
  de_samp <- intersect(rownames(meta_all)[meta_all$group %in% DE_GROUPS], colnames(mat))
  meta <- droplevels(meta_all[de_samp, ])
  meta$group     <- factor(meta$group, levels = DE_GROUPS)
  meta$replicate <- factor(meta$replicate)              # CRITICAL: batch as factor, not integer
  mat <- mat[, de_samp, drop = FALSE]

  # annotation parallel to DESeq2
  if (level == "gene") {
    ann <- annotate_genes(rownames(mat)); ann$feature_id <- ann$gene_id
  } else {
    t2g <- load_tx2gene()
    gid <- if (!is.null(t2g)) t2g$gene_id[match(rownames(mat), t2g$tx_id)]
           else cd$annot$gene_id[match(rownames(mat), cd$annot$feature_id)]
    g_ann <- annotate_genes(gid)
    ann <- data.frame(feature_id = rownames(mat), gene_id = gid,
                      ensembl = g_ann$ensembl, symbol = g_ann$symbol,
                      entrez = g_ann$entrez, stringsAsFactors = FALSE)
  }

  y <- DGEList(counts = mat, group = meta$group, genes = ann)
  keep <- filterByExpr(y, group = meta$group)
  y <- y[keep, , keep.lib.sizes = FALSE]
  y <- calcNormFactors(y, method = "TMM")
  message(sprintf("  features kept: %d", nrow(y)))

  design <- model.matrix(~ 0 + group + replicate, data = meta)
  # rename ONLY the group columns (groupNP_F.. -> NP_F..); leave replicate2/3 as-is
  gcols <- paste0("group", levels(meta$group))
  colnames(design)[match(gcols, colnames(design))] <- levels(meta$group)
  stopifnot(qr(design)$rank == ncol(design))           # full-rank guard
  v    <- voom(y, design)
  vfit <- lmFit(v, design)

  # Contrast matrix over the FULL design columns: group columns carry the
  # contrast vector, replicate (batch) columns stay 0.
  cmat <- sapply(names(CONTRASTS), function(cn) {
    col <- setNames(numeric(ncol(design)), colnames(design))
    cv  <- contrast_vec(cn)                              # named over DE_GROUPS
    col[names(cv)] <- cv
    col
  })
  rownames(cmat) <- colnames(design)
  stopifnot(identical(rownames(cmat), colnames(design)))
  # replicate rows must be exactly 0 (batch is adjusted for, not contrasted)
  rep_rows <- grep("^replicate", rownames(cmat))
  stopifnot(all(cmat[rep_rows, ] == 0))
  efit <- eBayes(contrasts.fit(vfit, contrasts = cmat))

  res_list <- lapply(names(CONTRASTS), function(cn) {
    tt <- topTable(efit, coef = cn, number = Inf, sort.by = "P")
    data.frame(
      feature_id = tt$feature_id, gene_id = tt$gene_id, symbol = tt$symbol,
      entrez = tt$entrez, baseMean = 2^tt$AveExpr,
      log2FC = tt$logFC, lfcSE = NA_real_, stat = tt$t,
      pvalue = tt$P.Value, padj = tt$adj.P.Val,
      contrast = cn, level = level, engine = "edgeR",
      stringsAsFactors = FALSE)
  })
  names(res_list) <- names(CONTRASTS)

  list(results = res_list,
       logcpm = cpm(y, log = TRUE),          # for QC / expression views
       norm.factors = y$samples$norm.factors,
       lib.size = y$samples$lib.size,
       meta = meta, annot = ann)
}

ed <- lapply(LEVELS, run_level); names(ed) <- LEVELS
save_cache(ed, "edger_bc")    # heavy intermediate — not bundled into the app

for (lv in LEVELS) for (cn in names(CONTRASTS)) {
  write_table_artifact(ed[[lv]]$results[[cn]], sprintf("DE_edgeR_bc_%s_%s", lv, cn))
  cat(sprintf("  %-10s %-18s sig(padj<.05): %d\n", lv, cn,
              sum(ed[[lv]]$results[[cn]]$padj < 0.05, na.rm = TRUE)))
}
message("edgeR (batch-corrected) done -> analysis/cache/edger_bc.rds")
