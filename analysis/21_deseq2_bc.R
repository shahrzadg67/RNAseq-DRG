#!/usr/bin/env Rscript
# ============================================================================
# 21_deseq2_bc.R  —  BATCH-CORRECTED DESeq2 differential expression.
# Identical to 20_deseq2.R except the design adds `replicate` as an additive
# batch covariate:   ~ 0 + group + replicate
# SNI is already excluded (DE uses DE_GROUPS = NP/Sham only); after dropping
# SNI the 12 samples are balanced across replicate 1/2/3, so the design is
# full rank (4 group + 2 replicate = 6 coefs, 12 samples). The contrast
# vectors are still pure group-mean contrasts — the replicate coefficients
# (named `replicate*`, NOT `group*`) stay 0 in to_full(), so the comparison
# is now marginal over (adjusted for) the replicate batch effect.
# Writes a SEPARATE cache (deseq2_bc) — the original deseq2.rds is untouched.
# ============================================================================
source("/hpf/projects/msalter/sghazis/rnaseq_TUY35595/analysis/_common.R")
suppressPackageStartupMessages({ library(DESeq2) })

meta_all <- load_metadata()

run_level <- function(level) {
  message("\n=== DESeq2 (batch-corrected) :: ", level, " level ===")
  cd <- load_counts(level)
  mat <- cd$mat

  # keep only NP + Sham samples that are present in the matrix
  de_samp <- rownames(meta_all)[meta_all$group %in% DE_GROUPS]
  de_samp <- intersect(de_samp, colnames(mat))
  meta <- droplevels(meta_all[de_samp, ])
  meta$group     <- factor(meta$group, levels = DE_GROUPS)
  meta$replicate <- factor(meta$replicate)              # CRITICAL: batch as factor, not integer
  mat <- mat[, de_samp, drop = FALSE]

  # design full-rank guard (fails loudly if metadata ever becomes unbalanced)
  stopifnot(qr(model.matrix(~ 0 + group + replicate, data = meta))$rank ==
              nlevels(meta$group) + nlevels(meta$replicate) - 1)

  dds <- DESeqDataSetFromMatrix(countData = mat, colData = meta,
                                design = ~ 0 + group + replicate)
  keep <- rowSums(counts(dds) >= 10) >= 2
  dds <- dds[keep, ]
  message(sprintf("  features kept: %d / %d", sum(keep), length(keep)))
  dds <- DESeq(dds)

  rn <- resultsNames(dds)                       # groupNP_F.. + replicate_2_vs_1, replicate_3_vs_1
  stopifnot(all(paste0("group", DE_GROUPS) %in% rn))   # contrast machinery relies on this
  # build a contrast numeric vector aligned to resultsNames; replicate coefs stay 0
  to_full <- function(vec) {
    full <- setNames(numeric(length(rn)), rn)
    full[paste0("group", names(vec))] <- vec
    full
  }

  # annotation: gene-level uses gene_id; transcript-level maps tx -> parent gene
  if (level == "gene") {
    ann <- annotate_genes(rownames(dds))
    ann$feature_id <- ann$gene_id
  } else {
    t2g <- load_tx2gene()
    gid <- if (!is.null(t2g)) t2g$gene_id[match(rownames(dds), t2g$tx_id)]
           else cd$annot$gene_id[match(rownames(dds), cd$annot$feature_id)]
    g_ann <- annotate_genes(gid)
    ann <- data.frame(feature_id = rownames(dds), gene_id = gid,
                      ensembl = g_ann$ensembl, symbol = g_ann$symbol,
                      entrez = g_ann$entrez, stringsAsFactors = FALSE)
  }

  res_list <- lapply(names(CONTRASTS), function(cn) {
    r <- as.data.frame(results(dds, contrast = to_full(contrast_vec(cn))))
    r$feature_id <- rownames(r)
    r <- merge(r, ann, by = "feature_id", all.x = TRUE, sort = FALSE)
    data.frame(
      feature_id = r$feature_id, gene_id = r$gene_id, symbol = r$symbol,
      entrez = r$entrez, baseMean = r$baseMean,
      log2FC = r$log2FoldChange, lfcSE = r$lfcSE, stat = r$stat,
      pvalue = r$pvalue, padj = r$padj,
      contrast = cn, level = level, engine = "DESeq2",
      stringsAsFactors = FALSE)[order(r$pvalue), ]
  })
  names(res_list) <- names(CONTRASTS)

  # VST for PCA / expression boxplots (computed on the DE sample set)
  vsd <- tryCatch(vst(dds, blind = TRUE),
                  error = function(e) varianceStabilizingTransformation(dds, blind = TRUE))
  list(results = res_list,
       vst = assay(vsd),
       norm = counts(dds, normalized = TRUE),
       annot = ann, meta = meta)
}

de <- lapply(LEVELS, run_level); names(de) <- LEVELS
save_cache(de, "deseq2_bc")   # heavy intermediate — not bundled into the app

# flat CSVs per level x contrast for download / inspection
for (lv in LEVELS) for (cn in names(CONTRASTS)) {
  write_table_artifact(de[[lv]]$results[[cn]], sprintf("DE_DESeq2_bc_%s_%s", lv, cn))
  cat(sprintf("  %-10s %-18s sig(padj<.05): %d\n", lv, cn,
              sum(de[[lv]]$results[[cn]]$padj < 0.05, na.rm = TRUE)))
}
message("DESeq2 (batch-corrected) done -> analysis/cache/deseq2_bc.rds")
