#!/usr/bin/env Rscript
# ============================================================================
# run_all_noxy.R  —  BATCH-CORRECTED analysis with the SEX CHROMOSOMES REMOVED.
# Same design as the *_bc pipeline (NP/Sham only, replicate modelled as an
# additive batch covariate; replicate removed via removeBatchEffect for PCA/
# boxplots), but chrX and chrY features are dropped from the count matrix first
# (list from analysis/cache/xy_features.rds, built by 60_annotate.R).
# All outputs get a _noxy suffix — nothing existing is touched.
#   caches : deseq2_noxy, edger_noxy
#   app/data: pca_noxy.rds, de_long_noxy.rds, expr_noxy.rds,
#             search_index_noxy.rds, gsea_noxy.rds
#   app/www : gsea_noxy/<level>/<contrast>/<cat>/*.png, pathview_noxy/<contrast>/*
# Run via sbatch (scripts/run_noxy.sbatch) so it never competes with the app.
# ============================================================================
source("/hpf/projects/msalter/sghazis/rnaseq_TUY35595/analysis/_common.R")
suppressPackageStartupMessages({
  library(DESeq2); library(edgeR); library(limma); library(matrixStats)
})
options(bitmapType = "cairo")

meta_all <- load_metadata()
xy <- read_cache("xy_features")                 # list(gene=<ids>, transcript=<ids>)

# count matrix with chrX/chrY features removed
counts_noxy <- function(level) {
  cd  <- load_counts(level)
  keep <- !rownames(cd$mat) %in% xy[[level]]
  message(sprintf("  %s: dropped %d sex-chr features (%d -> %d)",
                  level, sum(!keep), nrow(cd$mat), sum(keep)))
  list(mat = cd$mat[keep, , drop = FALSE], annot = cd$annot)
}

# annotation parallel to the other pipelines
annot_for <- function(level, ids, cd) {
  if (level == "gene") { a <- annotate_genes(ids); a$feature_id <- a$gene_id; return(a) }
  t2g <- load_tx2gene()
  gid <- if (!is.null(t2g)) t2g$gene_id[match(ids, t2g$tx_id)]
         else cd$annot$gene_id[match(ids, cd$annot$feature_id)]
  g <- annotate_genes(gid)
  data.frame(feature_id = ids, gene_id = gid, ensembl = g$ensembl,
             symbol = g$symbol, entrez = g$entrez, stringsAsFactors = FALSE)
}

# ---------------------------------------------------------------- DESeq2 (bc) --
run_deseq2 <- function(level) {
  message("\n=== DESeq2 (noXY) :: ", level, " ===")
  cd <- counts_noxy(level); mat <- cd$mat
  de_samp <- intersect(rownames(meta_all)[meta_all$group %in% DE_GROUPS], colnames(mat))
  meta <- droplevels(meta_all[de_samp, ])
  meta$group <- factor(meta$group, levels = DE_GROUPS)
  meta$replicate <- factor(meta$replicate)
  mat <- mat[, de_samp, drop = FALSE]
  stopifnot(qr(model.matrix(~ 0 + group + replicate, meta))$rank ==
              nlevels(meta$group) + nlevels(meta$replicate) - 1)
  dds <- DESeqDataSetFromMatrix(mat, meta, ~ 0 + group + replicate)
  dds <- dds[rowSums(counts(dds) >= 10) >= 2, ]
  dds <- DESeq(dds)
  rn <- resultsNames(dds); stopifnot(all(paste0("group", DE_GROUPS) %in% rn))
  to_full <- function(vec){ full <- setNames(numeric(length(rn)), rn); full[paste0("group", names(vec))] <- vec; full }
  ann <- annot_for(level, rownames(dds), cd)
  res_list <- lapply(names(CONTRASTS), function(cn){
    r <- as.data.frame(results(dds, contrast = to_full(contrast_vec(cn)))); r$feature_id <- rownames(r)
    r <- merge(r, ann, by = "feature_id", all.x = TRUE, sort = FALSE)
    data.frame(feature_id=r$feature_id, gene_id=r$gene_id, symbol=r$symbol, entrez=r$entrez,
               baseMean=r$baseMean, log2FC=r$log2FoldChange, lfcSE=r$lfcSE, stat=r$stat,
               pvalue=r$pvalue, padj=r$padj, contrast=cn, level=level, engine="DESeq2",
               stringsAsFactors=FALSE)[order(r$pvalue), ]
  }); names(res_list) <- names(CONTRASTS)
  vsd <- tryCatch(vst(dds, blind=TRUE), error=function(e) varianceStabilizingTransformation(dds, blind=TRUE))
  list(results=res_list, vst=assay(vsd), norm=counts(dds, normalized=TRUE), annot=ann, meta=meta)
}

# ---------------------------------------------------------------- edgeR (bc) --
run_edger <- function(level) {
  message("\n=== edgeR (noXY) :: ", level, " ===")
  cd <- counts_noxy(level); mat <- cd$mat
  de_samp <- intersect(rownames(meta_all)[meta_all$group %in% DE_GROUPS], colnames(mat))
  meta <- droplevels(meta_all[de_samp, ])
  meta$group <- factor(meta$group, levels = DE_GROUPS)
  meta$replicate <- factor(meta$replicate)
  mat <- mat[, de_samp, drop = FALSE]
  ann <- annot_for(level, rownames(mat), cd)
  y <- DGEList(counts=mat, group=meta$group, genes=ann)
  y <- y[filterByExpr(y, group=meta$group), , keep.lib.sizes=FALSE]
  y <- calcNormFactors(y, method="TMM")
  design <- model.matrix(~ 0 + group + replicate, meta)
  gcols <- paste0("group", levels(meta$group))
  colnames(design)[match(gcols, colnames(design))] <- levels(meta$group)
  stopifnot(qr(design)$rank == ncol(design))
  v <- voom(y, design); vfit <- lmFit(v, design)
  cmat <- sapply(names(CONTRASTS), function(cn){ col <- setNames(numeric(ncol(design)), colnames(design))
    cv <- contrast_vec(cn); col[names(cv)] <- cv; col })
  rownames(cmat) <- colnames(design)
  stopifnot(all(cmat[grep("^replicate", rownames(cmat)), ] == 0))
  efit <- eBayes(contrasts.fit(vfit, cmat))
  res_list <- lapply(names(CONTRASTS), function(cn){
    tt <- topTable(efit, coef=cn, number=Inf, sort.by="P")
    data.frame(feature_id=tt$feature_id, gene_id=tt$gene_id, symbol=tt$symbol, entrez=tt$entrez,
               baseMean=2^tt$AveExpr, log2FC=tt$logFC, lfcSE=NA_real_, stat=tt$t,
               pvalue=tt$P.Value, padj=tt$adj.P.Val, contrast=cn, level=level, engine="edgeR",
               stringsAsFactors=FALSE)
  }); names(res_list) <- names(CONTRASTS)
  list(results=res_list, meta=meta, annot=ann)
}

# ------------------------------------------------------------------- PCA (bc) --
run_pca <- function(level) {
  cd <- counts_noxy(level); samples <- rownames(meta_all)[meta_all$condition != "SNI"]
  mat <- cd$mat[, intersect(samples, colnames(cd$mat)), drop=FALSE]
  meta <- droplevels(meta_all[colnames(mat), ]); meta$replicate <- factor(meta$replicate)
  dds <- DESeqDataSetFromMatrix(mat, meta, ~ 1); dds <- dds[rowSums(counts(dds) >= 10) >= 2, ]
  v <- tryCatch(vst(dds, blind=TRUE), error=function(e) varianceStabilizingTransformation(dds, blind=TRUE))
  x <- limma::removeBatchEffect(assay(v), batch=meta$replicate, design=model.matrix(~ group, meta))
  x <- x[order(rowVars(x), decreasing=TRUE)[seq_len(min(1000, nrow(x)))], ]
  pca <- prcomp(t(x), center=TRUE, scale.=FALSE); k <- min(10, ncol(pca$x))
  ve <- (pca$sdev^2/sum(pca$sdev^2))[seq_len(k)]
  list(coords=data.frame(sample_id=rownames(pca$x),
         meta[rownames(pca$x), c("condition","sex","group","replicate")],
         pca$x[,seq_len(k),drop=FALSE], row.names=NULL, check.names=FALSE),
       var_explained=data.frame(PC=paste0("PC",seq_len(k)), pct=round(100*ve,2), cumpct=round(100*cumsum(ve),2)),
       level=level, include_sni=FALSE, n_features=nrow(x))
}

# ============================ run DE + PCA + assemble ============================
de <- lapply(LEVELS, run_deseq2); names(de) <- LEVELS; save_cache(de, "deseq2_noxy")
ed <- lapply(LEVELS, run_edger);  names(ed) <- LEVELS; save_cache(ed, "edger_noxy")
for (lv in LEVELS) for (cn in names(CONTRASTS)) {
  write_table_artifact(de[[lv]]$results[[cn]], sprintf("DE_DESeq2_noxy_%s_%s", lv, cn))
  write_table_artifact(ed[[lv]]$results[[cn]], sprintf("DE_edgeR_noxy_%s_%s", lv, cn))
}
pca <- setNames(lapply(LEVELS, run_pca), LEVELS); save_artifact(pca, "pca_noxy")

collect <- function(src) do.call(rbind, unlist(lapply(LEVELS, function(lv) src[[lv]]$results), recursive=FALSE))
dl <- rbind(collect(de), collect(ed))
de_long <- data.frame(feature_id=dl$feature_id, symbol=dl$symbol, log2FC=round(dl$log2FC,3),
                      padj=signif(dl$padj,4), contrast=factor(dl$contrast), level=factor(dl$level),
                      engine=factor(dl$engine), stringsAsFactors=FALSE)
saveRDS(de_long, file.path(APP_DATA,"de_long_noxy.rds"), compress="xz")
expr <- setNames(lapply(LEVELS, function(lv){ meta <- de[[lv]]$meta; meta$replicate <- factor(meta$replicate)
  list(vst=limma::removeBatchEffect(de[[lv]]$vst, batch=meta$replicate, design=model.matrix(~group,meta)), meta=meta) }), LEVELS)
saveRDS(expr, file.path(APP_DATA,"expr_noxy.rds"))
idx <- do.call(rbind, lapply(LEVELS, function(lv){ r <- de[[lv]]$results[[1]]
  data.frame(level=lv, feature_id=r$feature_id, gene_id=r$gene_id, symbol=r$symbol,
             label=ifelse(is.na(r$symbol), r$feature_id, paste0(r$symbol," (",r$feature_id,")")), stringsAsFactors=FALSE) }))
idx <- idx[!duplicated(paste(idx$level, idx$feature_id)), ]
saveRDS(idx, file.path(APP_DATA,"search_index_noxy.rds"))
message("DE + PCA + assemble (noXY) done")

# ================================== GSEA ========================================
source(file.path(PROJ,"analysis","gsea_core_noxy.R"))
message("\nALL DONE (noXY).")
