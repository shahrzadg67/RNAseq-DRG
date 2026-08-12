#!/usr/bin/env Rscript
# ============================================================================
# 90_render_bc_figures.R  —  render static PNG figures for the BATCH-CORRECTED
# analysis (the app shows these interactively; this writes downloadable PNGs).
# Also assembles a single deliverables folder with DE tables + all PNGs.
#   batch_corrected_results/
#     DE_tables/        (the 20 DESeq2/edgeR CSVs, copied)
#     PCA_plots/        (scatter coloured by group/condition/sex + scree)
#     volcano_plots/    (DESeq2 & edgeR volcano per level x contrast)
#     GSEA_plots/       (running-score PNGs, copied from app/www/gsea_bc)
# ============================================================================
source("/hpf/projects/msalter/sghazis/rnaseq_TUY35595/analysis/_common.R")
suppressPackageStartupMessages({ library(ggplot2) })
options(bitmapType = "cairo")

OUT  <- file.path(PROJ, "batch_corrected_results")
dirs <- c("DE_tables", "PCA_plots", "volcano_plots", "GSEA_plots")
for (d in dirs) dir.create(file.path(OUT, d), recursive = TRUE, showWarnings = FALSE)

## ---- 1. copy DE tables -----------------------------------------------------
csvs <- list.files(file.path(CACHE, "tables"), pattern = "_bc_.*\\.csv$", full.names = TRUE)
file.copy(csvs, file.path(OUT, "DE_tables"), overwrite = TRUE)
message("DE tables copied: ", length(csvs))

## ---- 2. PCA plots ----------------------------------------------------------
pca <- readRDS(file.path(APP_DATA, "pca_bc.rds"))
pal <- function(v) setNames(scales::hue_pal()(length(unique(v))), sort(unique(as.character(v))))
for (lv in names(pca)) {
  d  <- pca[[lv]]; co <- d$coords; ve <- d$var_explained
  vx <- ve$pct[1]; vy <- ve$pct[2]
  for (cb in c("group","condition","sex")) {
    p <- ggplot(co, aes(PC1, PC2, colour = .data[[cb]], label = sample_id)) +
      geom_point(size = 4) +
      ggrepel::geom_text_repel(size = 3, show.legend = FALSE, max.overlaps = 20) +
      labs(title = sprintf("PCA (batch-corrected) — %s level, coloured by %s", lv, cb),
           subtitle = "SNI excluded · replicate batch effect removed (limma::removeBatchEffect)",
           x = sprintf("PC1 (%.1f%%)", vx), y = sprintf("PC2 (%.1f%%)", vy)) +
      theme_bw(base_size = 13)
    ggsave(file.path(OUT, "PCA_plots", sprintf("PCA_bc_%s_PC1-PC2_by_%s.png", lv, cb)),
           p, width = 7.5, height = 6, dpi = 150)
  }
  # scree / elbow
  ve$PC <- factor(ve$PC, levels = ve$PC)
  ps <- ggplot(ve, aes(PC, pct)) +
    geom_col(fill = "#18bc9c") +
    geom_line(aes(y = cumpct, group = 1), colour = "#e74c3c") +
    geom_point(aes(y = cumpct), colour = "#e74c3c") +
    labs(title = sprintf("Variance explained per PC (batch-corrected) — %s", lv),
         x = "Principal component", y = "% variance (bars) / cumulative % (line)") +
    theme_bw(base_size = 13)
  ggsave(file.path(OUT, "PCA_plots", sprintf("PCA_bc_%s_scree.png", lv)),
         ps, width = 7.5, height = 5, dpi = 150)
}
message("PCA plots written")

## ---- 3. volcano plots (DESeq2 & edgeR, per level x contrast) ---------------
de_long <- readRDS(file.path(APP_DATA, "de_long_bc.rds"))
PADJ <- 0.05; LFC <- 1
for (en in c("DESeq2","edgeR")) for (lv in c("gene","transcript")) for (cn in names(CONTRASTS)) {
  d <- de_long[de_long$engine==en & de_long$level==lv & de_long$contrast==cn, ]
  if (!nrow(d)) next
  d$sig <- !is.na(d$padj) & d$padj < PADJ & abs(d$log2FC) >= LFC
  ns <- sum(d$sig)
  p <- ggplot(d, aes(log2FC, -log10(padj), colour = sig)) +
    geom_point(size = 0.7, alpha = 0.5) +
    scale_colour_manual(values = c("FALSE"="#b0b8bf","TRUE"="#e74c3c"), guide = "none") +
    geom_vline(xintercept = c(-LFC, LFC), linetype = "dotted", colour = "#888") +
    geom_hline(yintercept = -log10(PADJ), linetype = "dotted", colour = "#888") +
    labs(title = sprintf("%s — %s level — %s (batch-corrected)", en, lv, cn),
         subtitle = sprintf("%d significant at padj<%.2f & |log2FC|>=%g", ns, PADJ, LFC),
         x = "log2 fold-change", y = "-log10 adj.p") +
    theme_bw(base_size = 12)
  ggsave(file.path(OUT, "volcano_plots", sprintf("volcano_bc_%s_%s_%s.png", en, lv, cn)),
         p, width = 7, height = 5.5, dpi = 150)
}
message("volcano plots written")

## ---- 4. copy GSEA running-score PNGs (preserve level/contrast/category tree) -
gsrc <- file.path(PROJ, "app", "www", "gsea_bc")
if (dir.exists(gsrc)) {
  pngs <- list.files(gsrc, pattern = "\\.png$", recursive = TRUE)
  for (rel in pngs) {
    dst <- file.path(OUT, "GSEA_plots", rel)
    dir.create(dirname(dst), recursive = TRUE, showWarnings = FALSE)
    file.copy(file.path(gsrc, rel), dst, overwrite = TRUE)
  }
  message("GSEA PNGs copied: ", length(pngs))
}
pvsrc <- file.path(PROJ, "app", "www", "pathview_bc")
if (dir.exists(pvsrc)) {
  pv <- list.files(pvsrc, pattern = "\\.png$", recursive = TRUE)
  if (length(pv)) {
    dir.create(file.path(OUT, "GSEA_plots", "pathview_KEGG"), showWarnings = FALSE)
    file.copy(file.path(pvsrc, pv), file.path(OUT, "GSEA_plots", "pathview_KEGG"), overwrite = TRUE)
    message("pathview KEGG PNGs copied: ", length(pv))
  }
}
message("\nDeliverables assembled -> ", OUT)
