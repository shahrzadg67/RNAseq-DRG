#!/usr/bin/env Rscript
# ============================================================================
# 08_deseq2_analysis.R  —  Differential expression for rnaseq_TUY35595
# ----------------------------------------------------------------------------
# Loads the nf-core/rnaseq gene-level Salmon counts + metadata.csv, builds a
# DESeq2 dataset with Sham as the reference condition, and runs the key
# contrasts (NP vs Sham, SNI vs Sham), controlling for sex.
#
# Design: 3 conditions (NP / SNI / Sham) x 2 sexes (M / F).
#   NP   : n=3 per sex   (solid)
#   Sham : n=3 per sex   (solid)  <- reference / control
#   SNI  : n=1 per sex   (NO within-group replication -- see WARNING below)
#
# Run AFTER the pipeline finishes, from the project root:
#   cd /hpf/projects/msalter/sghazis/rnaseq_TUY35595
#   Rscript scripts/08_deseq2_analysis.R
# (Needs R with DESeq2; e.g. `module load R` then install if absent.)
# ============================================================================

suppressPackageStartupMessages({
  library(DESeq2)
})

## ---- paths -----------------------------------------------------------------
proj      <- "/hpf/projects/msalter/sghazis/rnaseq_TUY35595"
counts_f  <- file.path(proj, "results/star_salmon/salmon.merged.gene_counts.tsv")
meta_f    <- file.path(proj, "metadata.csv")
outdir    <- file.path(proj, "deseq2")
dir.create(outdir, showWarnings = FALSE)

stopifnot(file.exists(counts_f), file.exists(meta_f))

## ---- 1. load gene counts ---------------------------------------------------
# nf-core star_salmon matrix: cols = gene_id, gene_name, then one col per sample.
raw <- read.delim(counts_f, header = TRUE, check.names = FALSE)
gene_id   <- raw$gene_id
gene_name <- raw$gene_name
mat <- as.matrix(raw[, !(colnames(raw) %in% c("gene_id", "gene_name"))])
rownames(mat) <- gene_id
# Salmon gives estimated (fractional) counts -> DESeq2 needs integers.
mat <- round(mat)
mode(mat) <- "integer"

## ---- 2. load + align metadata ----------------------------------------------
meta <- read.csv(meta_f, stringsAsFactors = FALSE)
rownames(meta) <- meta$sample_id
# reorder metadata to match the count-matrix columns exactly
meta <- meta[colnames(mat), , drop = FALSE]
stopifnot(identical(rownames(meta), colnames(mat)))

# factors; Sham is the reference level so log2FC reads as "<condition> vs Sham"
meta$condition <- relevel(factor(meta$condition), ref = "Sham")
meta$sex       <- factor(meta$sex)

## ---- 3. build DESeq2 dataset -----------------------------------------------
# Design controls for sex, then tests condition. (Drop + sex if sexes are to be
# analysed separately, or use ~ sex * condition to test sex:condition interaction.)
dds <- DESeqDataSetFromMatrix(countData = mat,
                              colData   = meta,
                              design    = ~ sex + condition)

# light pre-filter: keep genes with >= 10 reads across >= the smallest group (2)
keep <- rowSums(counts(dds) >= 10) >= 2
dds  <- dds[keep, ]
message(sprintf("Genes retained after filter: %d / %d", sum(keep), length(keep)))

dds <- DESeq(dds)
saveRDS(dds, file.path(outdir, "dds.rds"))

## ---- 4. contrasts vs Sham --------------------------------------------------
name_map <- data.frame(gene_id = gene_id, gene_name = gene_name)

write_res <- function(dds, contrast, tag) {
  res <- results(dds, contrast = contrast)
  res <- res[order(res$padj), ]
  df  <- as.data.frame(res)
  df$gene_id   <- rownames(df)
  df$gene_name <- name_map$gene_name[match(df$gene_id, name_map$gene_id)]
  df <- df[, c("gene_id", "gene_name",
               "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")]
  out <- file.path(outdir, sprintf("DE_%s.csv", tag))
  write.csv(df, out, row.names = FALSE)
  message(sprintf("%-14s  sig(padj<0.05): %d   -> %s",
                  tag, sum(df$padj < 0.05, na.rm = TRUE), out))
}

write_res(dds, c("condition", "NP",  "Sham"), "NP_vs_Sham")   # n=3 vs n=3 (solid)
write_res(dds, c("condition", "SNI", "Sham"), "SNI_vs_Sham")  # n=1/sex (see WARNING)

# WARNING: SNI has a single replicate per sex. DESeq2 will still report a result
# by borrowing the shared dispersion estimate, but SNI contrasts have no real
# within-group replication -> treat p-values/padj as exploratory only.

## ---- 5. QC: VST, PCA, sample distances -------------------------------------
vsd <- vst(dds, blind = TRUE)
saveRDS(vsd, file.path(outdir, "vsd.rds"))

pca <- plotPCA(vsd, intgroup = c("condition", "sex"), returnData = TRUE)
write.csv(pca, file.path(outdir, "pca_coords.csv"), row.names = FALSE)

pdf(file.path(outdir, "QC_plots.pdf"), width = 7, height = 6)
print(plotPCA(vsd, intgroup = c("condition", "sex")))
plotDispEsts(dds)
dev.off()

# normalized counts for downstream / sanity checks
write.csv(counts(dds, normalized = TRUE),
          file.path(outdir, "normalized_counts.csv"))

message("\nDone. Outputs in: ", outdir)
