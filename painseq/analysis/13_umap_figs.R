# =====================================================================
# 13_umap_figs.R — draw the scanpy embeddings in the app's house style.
#
# Reads data/umap_<name>.csv written by 11_qc_umap.py. Every column needed for
# plotting is already in that CSV, so no join back to cell_meta is required —
# which matters, because GEO's cell names are not unique (see 10_export).
#
# Out: figs/umap_<name>_{celltype,injury,time,panelscore}.png/pdf
#      figs/qc_<name>.png/pdf
# =====================================================================
suppressMessages({ library(ggplot2) })
source("/hpf/projects/msalter/sghazis/painseq/analysis/00_common.R")

# categorical palette — colourblind-safe, distinct at 20 levels
PAL20 <- c("#4E79A7","#F28E2B","#E15759","#76B7B2","#59A14F","#EDC948","#B07AA1","#FF9DA7",
           "#9C755F","#BAB0AC","#1F77B4","#D62728","#2CA02C","#9467BD","#8C564B","#E377C2",
           "#7F7F7F","#BCBD22","#17BECF","#AEC7E8")

fmt_time <- function(h) ifelse(h == 0, "0", ifelse(h < 168, paste0(h, "h"), paste0(h / 24, "d")))

umap_cat <- function(d, col, title, pal = NULL, label_centroids = FALSE, legend_cols = 1) {
  d$.g <- factor(d[[col]])
  cols <- if (!is.null(pal) && all(levels(d$.g) %in% names(pal))) pal[levels(d$.g)]
          else setNames(rep(PAL20, length.out = nlevels(d$.g)), levels(d$.g))
  g <- ggplot(d[sample(nrow(d)), ], aes(UMAP1, UMAP2, colour = .g)) +
    geom_point(size = 0.12, alpha = 0.55, shape = 16) +
    scale_colour_manual(values = cols, name = NULL,
                        guide = guide_legend(override.aes = list(size = 3, alpha = 1), ncol = legend_cols)) +
    labs(title = title, x = "UMAP 1", y = "UMAP 2") +
    theme_pub(base_size = 13) +
    theme(axis.text = element_blank(), axis.ticks = element_blank(),
          legend.text = element_text(size = 8.5), legend.key.height = unit(11, "pt"))
  if (label_centroids) {
    cen <- aggregate(cbind(UMAP1, UMAP2) ~ .g, data = d, FUN = median)
    g <- g + geom_text(data = cen, aes(UMAP1, UMAP2, label = .g), inherit.aes = FALSE,
                       size = 2.9, fontface = "bold", colour = "grey10")
  }
  g
}

umap_num <- function(d, col, title, lab) {
  v <- d[[col]]; lim <- stats::quantile(v, c(0.01, 0.99), na.rm = TRUE)
  d$.v <- pmax(pmin(v, lim[2]), lim[1])
  ggplot(d[order(d$.v), ], aes(UMAP1, UMAP2, colour = .v)) +
    geom_point(size = 0.12, alpha = 0.6, shape = 16) +
    scale_colour_gradient2(low = DIVERGE_HEX[1], mid = DIVERGE_HEX[2], high = DIVERGE_HEX[3],
                           midpoint = 0, name = lab) +
    labs(title = title, x = "UMAP 1", y = "UMAP 2") +
    theme_pub(base_size = 13) +
    theme(axis.text = element_blank(), axis.ticks = element_blank(),
          legend.key.height = unit(24, "pt"))
}

sv <- function(g, stem, w = 8.5, h = 7) {
  ggsave(file.path(FIGS, paste0(stem, ".png")), g, width = w, height = h, dpi = 300, bg = "white")
  ggsave(file.path(FIGS, paste0(stem, ".pdf")), g, width = w, height = h, device = cairo_pdf)
  msg("  wrote figs/", stem, ".{png,pdf}")
}

run <- function(name, label) {
  f <- file.path(DATA, paste0("umap_", name, ".csv"))
  if (!file.exists(f)) { msg("SKIP ", name, " — ", f, " not found"); return(invisible(NULL)) }
  d <- read.csv(f, stringsAsFactors = FALSE)
  msg("=== ", name, ": ", nrow(d), " cells ===")
  integ <- unique(d$integration)[1]
  # Un-integrated is the deliberate choice here, not a shortfall: `sample` is
  # nested within injury x time, so integrating on it would delete the condition
  # effect. See the rationale in 11_qc_umap.py.
  note <- if (identical(integ, "harmony")) "Harmony-integrated across samples"
          else "un-integrated by design"
  msg("  integration: ", integ)

  d$time_lab <- factor(fmt_time(d$time), levels = fmt_time(sort(unique(d$time))))
  d$inj      <- factor(d$inj, levels = intersect(INJ_ORDER, unique(d$inj)))

  sv(umap_cat(d, "celltype", sprintf("%s — cell type  (%s)", label, note),
              label_centroids = TRUE, legend_cols = 1), paste0("umap_", name, "_celltype"), 9.5, 7.5)
  sv(umap_cat(d, "inj", sprintf("%s — injury model  (%s)", label, note)),
     paste0("umap_", name, "_injury"))
  sv(umap_cat(d, "time_lab", sprintf("%s — time after injury  (%s)", label, note)),
     paste0("umap_", name, "_time"))
  if ("geno" %in% names(d) && length(unique(d$geno)) > 1)
    sv(umap_cat(d, "geno", sprintf("%s — genotype  (%s)", label, note)),
       paste0("umap_", name, "_genotype"))
  sv(umap_num(d, "panel_score",
              sprintf("%s\n47-gene conserved injury panel score", label), "panel\nscore"),
     paste0("umap_", name, "_panelscore"), 9.5, 7.5)

  # QC: per-sample distributions
  q <- data.frame(sample = d$sample, inj = d$inj,
                  UMI = d$total_counts, genes = d$n_genes_by_counts, pct_mt = d$pct_counts_mt)
  ql <- rbind(
    data.frame(sample = q$sample, inj = q$inj, metric = "UMI per nucleus",   value = q$UMI),
    data.frame(sample = q$sample, inj = q$inj, metric = "Genes per nucleus", value = q$genes),
    data.frame(sample = q$sample, inj = q$inj, metric = "% mitochondrial",   value = q$pct_mt))
  gq <- ggplot(ql, aes(x = reorder(sample, as.integer(inj)), y = value, fill = inj)) +
    geom_violin(scale = "width", linewidth = 0.15, colour = "grey30") +
    facet_wrap(~ metric, ncol = 1, scales = "free_y") +
    scale_fill_manual(values = setNames(PAL20[seq_len(nlevels(ql$inj))], levels(ql$inj)), name = NULL) +
    labs(title = sprintf("%s — QC by sample", label),
         subtitle = paste(strwrap(paste(
           "Single-NUCLEUS data: median %mito is well under 1, so a mitochondrial filter is not",
           "informative here. The authors also pre-filtered, so our QC removes almost nothing —",
           "it is confirmatory, not corrective."), width = 130), collapse = "\n"),
         x = NULL, y = NULL) +
    theme_pub(base_size = 11) +
    theme(axis.text.x = element_text(angle = 90, vjust = .5, hjust = 1, size = 5.5),
          plot.subtitle = element_text(colour = "grey30", size = 8.5),
          strip.background = element_rect(fill = "grey93", colour = "black"),
          strip.text = element_text(face = "bold"))
  sv(gq, paste0("qc_", name), 13, 9)

  qs <- file.path(DATA, paste0("qc_summary_", name, ".csv"))
  if (file.exists(qs)) print(read.csv(qs))
}

run("c57",  "GSE154659 C57 injury atlas")
run("atf3", "GSE154659 Atf3-WT / Atf3-KO")
msg("done")
