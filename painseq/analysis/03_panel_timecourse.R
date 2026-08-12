# =====================================================================
# 03_panel_timecourse.R — the 47-gene conserved injury panel across
# injury model x time, from whole-DRG pseudobulk.
#
# Panel A (primary)   : log2FC vs Naive          -> figs/panelA_log2fc_timecourse.{png,pdf,html}
# Panel B (companion) : row z-score of mean logCPM -> figs/panelB_zscore_timecourse.{png,pdf,html}
# Sensitivity         : Panel A on male C57 cells only
#
# Model: limma on pseudobulk logCPM, ~0 + group + sex + geno, each of the 24
# injury x time groups contrasted against Naive_0. Moderated variance lets the
# n=1 groups (CFA-48h, ScNT-6h, ScNT-24h) still yield an estimate; those columns
# are marked on the figure with a dagger.
# =====================================================================
suppressMessages({ library(Matrix); library(edgeR); library(limma); library(ggplot2) })
source("/hpf/projects/msalter/sghazis/painseq/analysis/00_common.R")

CLIP_L2FC <- 3      # log2FCs in this panel reach ~5; clip so mid-range structure stays visible
CLIP_Z    <- 2      # matches the app's existing 47-marker heatmap

# hours -> compact readable label
fmt_time <- function(h) {
  ifelse(h == 0, "0",
  ifelse(h < 168, paste0(h, "h"),
         paste0(h / 24, "d")))
}

# =====================================================================
# Fit: pseudobulk -> logCPM + per-group log2FC vs Naive
# =====================================================================
fit_timecourse <- function(pb, label, use_covariates = TRUE) {
  msg("--- fitting: ", label, " (", ncol(pb$counts), " samples) ---")
  cts <- pb$counts; md <- pb$meta
  md$group <- factor(md$group, levels = unique(md$group[order(
    match(md$inj, INJ_ORDER), md$time)]))

  y <- DGEList(counts = cts, group = md$group)
  keep <- filterByExpr(y, group = md$group)
  msg("  filterByExpr keeps ", sum(keep), "/", length(keep), " genes")
  y <- y[keep, , keep.lib.sizes = FALSE]
  y <- calcNormFactors(y, method = "TMM")
  lcpm <- cpm(y, log = TRUE, prior.count = 1)

  # design: group means, plus sex/genotype covariates when they actually vary
  vary <- function(v) length(unique(v)) > 1
  cov_terms <- if (use_covariates) c(if (vary(md$sex)) "sex", if (vary(md$geno)) "geno") else character(0)
  form <- as.formula(paste("~0 + group", if (length(cov_terms)) paste("+", paste(cov_terms, collapse = " + ")) else ""))
  design <- model.matrix(form, data = md)
  colnames(design) <- make.names(colnames(design))
  # drop aliased columns so lmFit does not fail on a rank-deficient design
  qrd <- qr(design)
  if (qrd$rank < ncol(design)) {
    keepcol <- qrd$pivot[seq_len(qrd$rank)]
    msg("  design rank-deficient; dropping: ", paste(setdiff(colnames(design), colnames(design)[keepcol]), collapse = ", "))
    design <- design[, keepcol, drop = FALSE]
  }
  msg("  design: ", ncol(design), " columns [", paste(cov_terms, collapse = "+"), "]")

  fit <- lmFit(lcpm, design)
  grp_levels <- levels(md$group)
  ref <- paste0("group", make.names("Naive_0"))
  stopifnot(ref %in% colnames(design))
  targets <- setdiff(paste0("group", make.names(grp_levels)), ref)
  targets <- intersect(targets, colnames(design))

  cm <- makeContrasts(contrasts = paste0(targets, " - ", ref), levels = design)
  colnames(cm) <- sub("^group", "", targets)
  fit2 <- eBayes(contrasts.fit(fit, cm))

  l2fc <- sapply(colnames(cm), function(cn) topTable(fit2, coef = cn, number = Inf, sort.by = "none")$logFC)
  padj <- sapply(colnames(cm), function(cn) topTable(fit2, coef = cn, number = Inf, sort.by = "none")$adj.P.Val)
  rownames(l2fc) <- rownames(padj) <- rownames(lcpm)

  # Naive is the reference: log2FC is 0 by construction. Keep it as a visual anchor.
  naive_col <- matrix(0, nrow(l2fc), 1, dimnames = list(rownames(l2fc), make.names("Naive_0")))
  l2fc <- cbind(naive_col, l2fc)
  padj <- cbind(matrix(NA_real_, nrow(padj), 1, dimnames = list(rownames(padj), make.names("Naive_0"))), padj)

  # group-mean logCPM, for Panel B
  gm <- sapply(grp_levels, function(g) rowMeans(lcpm[, md$group == g, drop = FALSE]))
  colnames(gm) <- make.names(grp_levels)

  # column metadata, in fixed biological order
  cmeta <- unique(md[, c("group", "inj", "time")])
  cmeta$col <- make.names(as.character(cmeta$group))
  cmeta <- cmeta[order(match(cmeta$inj, INJ_ORDER), cmeta$time), ]
  ns <- table(md$group)
  cmeta$n_samples <- as.integer(ns[as.character(cmeta$group)])
  cmeta$label <- fmt_time(cmeta$time)
  rownames(cmeta) <- NULL

  list(l2fc = l2fc[, cmeta$col, drop = FALSE],
       padj = padj[, cmeta$col, drop = FALSE],
       gmean = gm[, cmeta$col, drop = FALSE],
       cmeta = cmeta, kept = rownames(lcpm), label = label)
}

# =====================================================================
# Assemble the panel gene matrix, keeping genes that were filtered out or are
# absent as explicit NA rows rather than dropping them silently.
# =====================================================================
panel_matrix <- function(f, pg, what = c("l2fc", "z")) {
  what <- match.arg(what)
  rows <- pg[order(pg$role != "panel", -pg$mean_sni), ]     # panel by effect size, refs last
  M <- matrix(NA_real_, nrow(rows), ncol(f$l2fc), dimnames = list(rows$gene, colnames(f$l2fc)))
  src <- if (what == "l2fc") f$l2fc else f$gmean
  hit <- rows$matrix_sym %in% rownames(src)
  M[hit, ] <- src[rows$matrix_sym[hit], , drop = FALSE]
  if (what == "z") {
    M <- t(scale(t(M)))                                     # row z-score across the columns
    M[!is.finite(M)] <- NA_real_
  }
  P <- matrix(NA_real_, nrow(rows), ncol(M), dimnames = dimnames(M))
  P[hit, ] <- f$padj[rows$matrix_sym[hit], , drop = FALSE]
  rows$in_model <- hit
  list(M = M, padj = P, rows = rows)
}

# =====================================================================
# Static ggplot heatmap — facet per injury model gives the block separators
# =====================================================================
heat_gg <- function(pm, cmeta, title, subtitle, fill_lab, clip) {
  M <- pm$M; rows <- pm$rows
  d <- expand.grid(gene = factor(rownames(M), levels = rev(rownames(M))),
                   col  = factor(colnames(M), levels = colnames(M)), stringsAsFactors = FALSE)
  d$value <- as.vector(M)
  d$padj  <- as.vector(pm$padj)
  j <- match(d$col, cmeta$col)
  d$inj   <- factor(as.character(cmeta$inj)[j], levels = INJ_ORDER)
  # x must key on the unique column id, NOT the time label. Time labels repeat
  # across injury models (12h occurs in ScNT and SpNT but not Crush), so a factor
  # built from labels gets one shared global level order and scrambles the times
  # inside each facet. Key on `col` and relabel at the scale.
  d$col   <- factor(d$col, levels = cmeta$col)
  d$vclip <- pmax(pmin(d$value, clip), -clip)
  d$sig   <- !is.na(d$padj) & d$padj < 0.05
  xlabs   <- setNames(cmeta$label, cmeta$col)

  # gene labels: mark low-detection genes and the two reference rows
  lab <- setNames(rows$gene, rows$gene)
  lab[rows$low_detection %in% TRUE] <- paste0(rows$gene[rows$low_detection %in% TRUE], " *")
  lab[rows$role == "reference"] <- paste0(rows$gene[rows$role == "reference"], " ‡")

  ggplot(d, aes(x = col, y = gene, fill = vclip)) +
    geom_tile(colour = "grey92", linewidth = 0.15) +
    geom_point(data = subset(d, sig), aes(x = col, y = gene), inherit.aes = FALSE,
               size = 0.45, colour = "grey15", alpha = 0.85) +
    facet_grid(~ inj, scales = "free_x", space = "free_x", switch = "x",
               labeller = labeller(inj = INJ_SHORT)) +
    scale_x_discrete(labels = function(v) xlabs[v]) +
    scale_fill_gradient2(low = DIVERGE_HEX[1], mid = DIVERGE_HEX[2], high = DIVERGE_HEX[3],
                         midpoint = 0, limits = c(-clip, clip), na.value = "grey88",
                         name = fill_lab) +
    scale_y_discrete(labels = function(g) lab[g]) +
    labs(title = title, subtitle = paste(strwrap(subtitle, width = 118), collapse = "\n"),
         x = NULL, y = NULL) +
    theme_pub(base_size = 11) +
    theme(axis.text.x = element_text(angle = 90, vjust = .5, hjust = 1, size = 8),
          axis.text.y = element_text(size = 7.2),
          plot.subtitle = element_text(colour = "grey30", size = 8.5, lineheight = 1.15),
          strip.background = element_rect(fill = "grey93", colour = "black"),
          strip.text = element_text(face = "bold", size = 8.5),
          panel.spacing.x = unit(2.4, "pt"),
          legend.position = "right", legend.key.height = unit(26, "pt"))
}

# =====================================================================
# Interactive plotly heatmap — matches the app's existing 47-marker heatmap
# =====================================================================
heat_plotly <- function(pm, cmeta, title, zlab, clip) {
  M <- pm$M; rows <- pm$rows
  g <- rownames(M); n <- ncol(M)
  hov <- matrix("", nrow(M), ncol(M))
  for (i in seq_len(nrow(M))) for (jj in seq_len(ncol(M)))
    hov[i, jj] <- sprintf("<b>%s</b><br>%s %s<br>%s: %s<br>adj.p: %s<br>n samples: %d",
      g[i], cmeta$inj[jj], cmeta$label[jj], zlab,
      ifelse(is.na(M[i, jj]), "n/a", sprintf("%+.2f", M[i, jj])),
      ifelse(is.na(pm$padj[i, jj]), "-", format.pval(pm$padj[i, jj], digits = 2)),
      cmeta$n_samples[jj])

  bnd <- cumsum(rle(as.character(cmeta$inj))$lengths)
  bnd <- bnd[-length(bnd)]
  mid <- tapply(seq_len(n), as.character(cmeta$inj), mean)
  mid <- mid[INJ_ORDER[INJ_ORDER %in% cmeta$inj]]

  p <- plotly::plot_ly(x = seq_len(n), y = g, z = pmax(pmin(M, clip), -clip),
        type = "heatmap", colorscale = CPCA_DIVERGE, zmid = 0, zmin = -clip, zmax = clip,
        text = hov, hoverinfo = "text",
        colorbar = list(title = list(text = gsub(" ", "<br>", zlab))))
  shp <- lapply(bnd, function(b) list(type = "line", x0 = b + .5, x1 = b + .5,
        y0 = -.5, y1 = length(g) - .5, line = list(color = "black", width = 2.5)))
  plotly::layout(p, title = title, shapes = shp,
      xaxis = list(tickmode = "array", tickvals = unname(mid), ticktext = names(mid), tickangle = -35),
      yaxis = list(autorange = "reversed", categoryorder = "array", categoryarray = rev(g),
                   tickfont = list(size = 9))) |> ply_pub()
}

save_fig <- function(gg, ply, stem, w = 11, h = 9) {
  ggsave(file.path(FIGS, paste0(stem, ".png")), gg, width = w, height = h, dpi = 300, bg = "white")
  ggsave(file.path(FIGS, paste0(stem, ".pdf")), gg, width = w, height = h, device = cairo_pdf)
  msg("  wrote figs/", stem, ".{png,pdf}")
  # The interactive copy is a bonus. plotly::partial_bundle() pulls in curl, which
  # needs libssl.so.10 (see scripts/run_analysis.sh for the .extralib shim) — so
  # skip the bundling step and never let an HTML failure kill the run.
  if (!is.null(ply) && requireNamespace("htmlwidgets", quietly = TRUE)) {
    ok <- tryCatch({
      htmlwidgets::saveWidget(ply, file.path(FIGS, paste0(stem, ".html")), selfcontained = TRUE)
      TRUE
    }, error = function(e) { msg("  (html skipped: ", conditionMessage(e), ")"); FALSE })
    if (ok) msg("  wrote figs/", stem, ".html")
  }
}

# =====================================================================
# Run
# =====================================================================
pb <- readRDS(file.path(DATA, "pb_sample.rds"))
pg <- readRDS(file.path(DATA, "panel_genes.rds"))

fit_all <- fit_timecourse(pb, "all cells", use_covariates = TRUE)

npanel <- sum(pg$role == "panel")
cap_n1 <- {
  n1 <- fit_all$cmeta[fit_all$cmeta$n_samples == 1, ]
  if (nrow(n1)) paste0("  n=1: ", paste(paste(n1$inj, n1$label), collapse = ", ")) else ""
}
sub_common <- sprintf(
  paste("GSE154659 mouse DRG snRNA-seq, whole-DRG pseudobulk, all %s cells, %d samples.",
        "Dot = adj.p<0.05. * = detected in <0.5%% of cells (treat with caution).",
        "grey = gene absent from the matrix or below the expression filter.",
        "‡ = RAG archetype shown for reference, not one of the 47.",
        "Naive is the reference column, so it is 0 by construction.%s"),
  format(sum(pb$meta$n_cells), big.mark = ","), nrow(pb$meta), cap_n1)

# ---- Panel A: log2FC vs Naive -------------------------------------------
pmA <- panel_matrix(fit_all, pg, "l2fc")
ggA <- heat_gg(pmA, fit_all$cmeta,
  "Conserved injury panel across injury model and time",
  sub_common, "log2FC\nvs Naive", CLIP_L2FC)
plA <- heat_plotly(pmA, fit_all$cmeta,
  "Conserved injury panel - log2FC vs Naive", "log2FC vs Naive", CLIP_L2FC)
save_fig(ggA, plA, "panelA_log2fc_timecourse")

# ---- Panel B: row z-score of mean logCPM ---------------------------------
pmB <- panel_matrix(fit_all, pg, "z")
ggB <- heat_gg(pmB, fit_all$cmeta,
  "Conserved injury panel - expression pattern (z-scored)",
  sub_common, "z-score\n(logCPM)", CLIP_Z)
plB <- heat_plotly(pmB, fit_all$cmeta,
  "Conserved injury panel - z-scored logCPM", "z-score", CLIP_Z)
save_fig(ggB, plB, "panelB_zscore_timecourse")

# ---- Sensitivity: male C57 only ------------------------------------------
keep_s <- pb$meta$sex == "male" & pb$meta$geno == "C57"
pb_m <- list(counts = pb$counts[, keep_s, drop = FALSE], meta = pb$meta[keep_s, ])
fit_m <- fit_timecourse(pb_m, "male C57 only", use_covariates = FALSE)
pmA_m <- panel_matrix(fit_m, pg, "l2fc")
ggA_m <- heat_gg(pmA_m, fit_m$cmeta,
  "Conserved injury panel - log2FC vs Naive (male C57 only, sensitivity)",
  sprintf("Sensitivity check: %s male C57 cells, no sex/genotype covariates",
          format(sum(pb_m$meta$n_cells), big.mark = ",")), "log2FC\nvs Naive", CLIP_L2FC)
save_fig(ggA_m, NULL, "panelA_log2fc_timecourse_maleC57")

# concordance between the two versions (verification step 5)
shared_c <- intersect(colnames(pmA$M), colnames(pmA_m$M))
r_sens <- cor(as.vector(pmA$M[, shared_c]), as.vector(pmA_m$M[, shared_c]), use = "complete.obs")
msg("SENSITIVITY: all-cells vs male-C57-only log2FC, Pearson r = ", round(r_sens, 4),
    " over ", length(shared_c), " shared columns")

saveRDS(list(all = fit_all, male = fit_m, panelA = pmA, panelB = pmB,
             sensitivity_r = r_sens), file.path(DATA, "panel_timecourse.rds"))
msg("wrote data/panel_timecourse.rds")
