#!/usr/bin/env Rscript
# ============================================================================
# 60_annotate.R  —  build a gene/transcript annotation table with SYMBOL, full
# DESCRIPTION (gene name), and CHROMOSOME, for hover tooltips + tables in the
# app, and to drive the sex-chromosome (chrX/chrY) removal downstream.
#   -> app/data/gene_annot.rds : data.frame(level, feature_id, gene_id, symbol,
#                                            description, chr)
#   -> analysis/cache/xy_features.rds : list(gene=<ids>, transcript=<ids>) on chrX/Y
# ============================================================================
source("/hpf/projects/msalter/sghazis/rnaseq_TUY35595/analysis/_common.R")
suppressPackageStartupMessages({ library(org.Mm.eg.db); library(AnnotationDbi) })

map1 <- function(ens, col) {
  ok <- !is.na(ens) & ens %in% keys(org.Mm.eg.db, "ENSEMBL")
  out <- rep(NA_character_, length(ens))
  if (any(ok)) out[ok] <- suppressMessages(
    mapIds(org.Mm.eg.db, keys = ens[ok], column = col, keytype = "ENSEMBL",
           multiVals = "first"))
  out
}

annot_level <- function(level) {
  cd  <- load_counts(level)
  fid <- rownames(cd$mat)
  if (level == "gene") {
    gid <- fid
  } else {
    t2g <- load_tx2gene()
    gid <- if (!is.null(t2g)) t2g$gene_id[match(fid, t2g$tx_id)]
           else cd$annot$gene_id[match(fid, cd$annot$feature_id)]
  }
  ens <- sub("\\..*$", "", gid)
  data.frame(level = level, feature_id = fid, gene_id = gid,
             symbol      = map1(ens, "SYMBOL"),
             description = map1(ens, "GENENAME"),
             chr         = map1(ens, "CHR"),
             stringsAsFactors = FALSE)
}

ann <- do.call(rbind, lapply(LEVELS, annot_level))
saveRDS(ann, file.path(APP_DATA, "gene_annot.rds"))
message("gene_annot rows: ", nrow(ann),
        "  (chr known: ", sum(!is.na(ann$chr)), ")")

# chrX / chrY feature lists per level (for the sex-chromosome-removed pipeline)
xy <- lapply(LEVELS, function(lv) {
  a <- ann[ann$level == lv, ]
  a$feature_id[a$chr %in% c("X", "Y")]
})
names(xy) <- LEVELS
save_cache(xy, "xy_features")
for (lv in LEVELS) {
  a <- ann[ann$level == lv, ]
  message(sprintf("  %-10s total=%d  chrX=%d  chrY=%d  -> removing %d",
                  lv, nrow(a), sum(a$chr=="X", na.rm=TRUE), sum(a$chr=="Y", na.rm=TRUE),
                  length(xy[[lv]])))
}
message("annotate done -> app/data/gene_annot.rds, analysis/cache/xy_features.rds")
