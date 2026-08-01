#!/usr/bin/env Rscript
# ============================================================================
# 92_all_png.R  —  render EVERY PCA and volcano as a static PNG for BOTH the
# original (uncorrected, with SNI) and batch-corrected analyses, into one
# obvious folder: plots_png/ . PNGs open natively in VS Code (image viewer) —
# no server, no proxy, immune to the OnDemand session hopping nodes.
#   plots_png/original/{PCA,volcano}/*.png
#   plots_png/batch_corrected/{PCA,volcano}/*.png
# ============================================================================
source("/hpf/projects/msalter/sghazis/rnaseq_TUY35595/analysis/_common.R")
suppressPackageStartupMessages({ library(ggplot2) })
options(bitmapType = "cairo")
have_repel <- requireNamespace("ggrepel", quietly = TRUE)

ROOT <- file.path(PROJ, "plots_png")
mk <- function(...) { d <- file.path(ROOT, ...); dir.create(d, recursive = TRUE, showWarnings = FALSE); d }

pca_png <- function(co, ve, cb, title, file) {
  p <- ggplot(co, aes(PC1, PC2, colour = .data[[cb]], label = sample_id)) +
    geom_point(size = 4)
  p <- if (have_repel) p + ggrepel::geom_text_repel(size = 3, show.legend = FALSE, max.overlaps = 30)
       else p + geom_text(size = 2.6, vjust = -0.8, show.legend = FALSE)
  p <- p + labs(title = title, x = sprintf("PC1 (%.1f%%)", ve$pct[1]),
                y = sprintf("PC2 (%.1f%%)", ve$pct[2])) + theme_bw(base_size = 13)
  ggsave(file, p, width = 7.5, height = 6, dpi = 150)
}

volcano_png <- function(d, title, file) {
  d <- d[!is.na(d$padj), ]
  d$sig <- d$padj < 0.05 & abs(d$log2FC) >= 1
  ns <- sum(d$sig)
  p <- ggplot(d, aes(log2FC, -log10(padj), colour = sig)) +
    geom_point(size = 0.6, alpha = 0.5) +
    scale_colour_manual(values = c("FALSE"="#b0b8bf","TRUE"="#e74c3c"), guide = "none") +
    geom_vline(xintercept = c(-1,1), linetype = "dotted", colour = "#888") +
    geom_hline(yintercept = -log10(0.05), linetype = "dotted", colour = "#888") +
    labs(title = title, subtitle = sprintf("%d significant (padj<0.05, |log2FC|>=1)", ns),
         x = "log2 fold-change", y = "-log10 adj.p") + theme_bw(base_size = 12)
  ggsave(file, p, width = 7, height = 5.5, dpi = 150)
}

render_set <- function(pca_rds, de_rds, sub) {
  pca <- readRDS(pca_rds); de <- readRDS(de_rds)
  dP <- mk(sub, "PCA"); dV <- mk(sub, "volcano")
  for (key in names(pca)) {
    x <- pca[[key]]
    for (cb in c("group","condition","sex"))
      pca_png(x$coords, x$var_explained, cb,
              sprintf("PCA (%s) — %s, by %s", sub, key, cb),
              file.path(dP, sprintf("PCA_%s_by_%s.png", key, cb)))
  }
  for (en in c("DESeq2","edgeR")) for (lv in c("gene","transcript")) for (cn in names(CONTRASTS)) {
    d <- de[de$engine==en & de$level==lv & de$contrast==cn, ]
    if (!nrow(d)) next
    volcano_png(d, sprintf("%s — %s — %s (%s)", en, lv, cn, sub),
                file.path(dV, sprintf("volcano_%s_%s_%s.png", en, lv, cn)))
  }
  message(sub, ": PCA=", length(list.files(dP)), " volcano=", length(list.files(dV)))
}

render_set(file.path(APP_DATA,"pca.rds"),    file.path(APP_DATA,"de_long.rds"),    "original")
render_set(file.path(APP_DATA,"pca_bc.rds"), file.path(APP_DATA,"de_long_bc.rds"), "batch_corrected")
message("\nAll PNGs -> ", ROOT)
