#!/usr/bin/env Rscript
# ============================================================================
# 93_pub_figs.R  —  pre-render publication-quality PCA + volcano figures as
# PNG (300 dpi) + SVG + PDF, using the lab conventions (Male=blue #021893,
# Female=red #941200; Sham=black circle, MINP=orange #E67E22 triangle; NP shown
# as MINP), for BOTH the uncorrected and corrected analyses. These static files
# back the "Download PNG/SVG/PDF" buttons and work in any deployment (shinylive
# or server). Output: app/www/figs/ .
# ============================================================================
source("/hpf/projects/msalter/sghazis/rnaseq_TUY35595/analysis/_common.R")
suppressPackageStartupMessages({ library(ggplot2); library(ragg) })

FIG <- file.path(PROJ, "app", "www", "figs"); dir.create(FIG, recursive = TRUE, showWarnings = FALSE)

# publication theme: rectangular border, no gridlines, larger fonts
TP <- theme_bw(base_size = 14) +
  theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.9),
        panel.grid = element_blank(), axis.text = element_text(colour = "black"),
        plot.title = element_text(face = "bold"))

# --- lab conventions (mirror app/R/helpers.R) -------------------------------
SEX_COLORS  <- c(F = "#941200", M = "#021893")
COND_COLORS <- c(Sham = "#000000", NP = "#E67E22", SNI = "#7F8C8D")
GG_SHAPES   <- c(Sham = 16, MINP = 17, SNI = 15)          # circle / triangle / square
GROUP_SEX   <- c(NP_F="F", NP_M="M", Sham_F="F", Sham_M="M", SNI_F="F", SNI_M="M")
disp <- function(x) gsub("NP", "MINP", as.character(x))
col_map <- function(vals, by) {
  vals <- as.character(vals)
  if (by == "sex")        { lv <- intersect(c("F","M"), vals); m <- setNames(unname(SEX_COLORS[lv]), lv) }
  else if (by=="condition"){ lv <- intersect(names(COND_COLORS), vals); m <- setNames(unname(COND_COLORS[lv]), lv) }
  else                     { lv <- intersect(names(GROUP_SEX), vals); m <- setNames(unname(SEX_COLORS[GROUP_SEX[lv]]), lv) }  # colour groups by sex
  setNames(unname(m), disp(names(m)))
}

save3 <- function(gg, base, w = 7, h = 5.5) {
  ragg::agg_png(paste0(base, ".png"), width = w, height = h, units = "in", res = 300); print(gg); dev.off()
  grDevices::svg(paste0(base, ".svg"), width = w, height = h); print(gg); dev.off()
  grDevices::cairo_pdf(paste0(base, ".pdf"), width = w, height = h); print(gg); dev.off()
}

pca_gg <- function(co, ve, by, title) {
  cols <- col_map(co[[by]], by)
  co$.col <- factor(disp(co[[by]]), levels = names(cols))
  co$.shp <- factor(disp(co$condition), levels = intersect(names(GG_SHAPES), disp(co$condition)))
  ggplot(co, aes(PC1, PC2, colour = .col, shape = .shp)) +
    geom_point(size = 3.4) +
    scale_colour_manual(values = cols, name = tools::toTitleCase(by)) +
    scale_shape_manual(values = GG_SHAPES[levels(co$.shp)], name = "Condition") +
    labs(title = title, x = sprintf("PC1 (%.1f%%)", ve$pct[1]), y = sprintf("PC2 (%.1f%%)", ve$pct[2])) +
    TP
}

volcano_gg <- function(d, title) {
  d <- d[!is.na(d$padj), ]
  d$sig <- d$padj < 0.05 & abs(d$log2FC) >= 1
  if (sum(!d$sig) > 8000) d <- d[sort(c(which(d$sig), sample(which(!d$sig), 8000))), ]
  ns <- sum(d$sig)
  ggplot(d, aes(log2FC, -log10(padj), colour = sig)) +
    geom_point(size = 0.7, alpha = 0.5) +
    scale_colour_manual(values = c("FALSE"="#B0B8BF","TRUE"="#E67E22"), guide = "none") +
    geom_vline(xintercept = c(-1,1), linetype = "dotted", colour = "#888") +
    geom_hline(yintercept = -log10(0.05), linetype = "dotted", colour = "#888") +
    labs(title = title, subtitle = sprintf("%d significant (padj<0.05, |log2FC|>=1)", ns),
         x = "log2 fold-change", y = "-log10 adj.p") +
    TP
}

## ---- PCA figures -----------------------------------------------------------
render_pca <- function(rds, tag, keyfun) {
  pca <- readRDS(file.path(APP_DATA, rds))
  for (k in names(pca)) for (by in c("sex","condition","group")) {
    d <- pca[[k]]
    save3(pca_gg(d$coords, d$var_explained, by,
                 sprintf("PCA (%s) — %s, by %s", tag, disp(k), by)),
          file.path(FIG, keyfun(k, by)), w = 7.5, h = 6)
  }
}
render_pca("pca.rds",     "uncorrected", function(k, by) sprintf("PCA_uncorrected_%s_by_%s", k, by))     # k = gene_withSNI etc.
render_pca("pca_noxy.rds","corrected",   function(k, by) sprintf("PCA_corrected_%s_by_%s", k, by))       # k = gene/transcript
message("PCA figures done")

## ---- Volcano figures -------------------------------------------------------
render_volc <- function(rds, tag) {
  de <- readRDS(file.path(APP_DATA, rds))
  for (en in c("DESeq2","edgeR")) for (lv in c("gene","transcript")) for (cn in names(CONTRASTS)) {
    d <- de[de$engine==en & de$level==lv & de$contrast==cn, ]; if (!nrow(d)) next
    save3(volcano_gg(d, sprintf("%s — %s — %s (%s)", en, lv, disp(cn), tag)),
          file.path(FIG, sprintf("volcano_%s_%s_%s_%s", tag, en, lv, cn)), w = 7, h = 5.5)
  }
}
render_volc("de_long.rds",      "uncorrected")
render_volc("de_long_noxy.rds", "corrected")
message("Volcano figures done -> ", FIG, "  (", length(list.files(FIG)), " files)")
