# =====================================================================
# 04_panel_celltype.R — which cell type expresses the conserved injury panel?
#
# Fits each cell type INDEPENDENTLY on its own (sample x celltype) pseudobulk,
# contrasting every injury x time group against that cell type's own Naive.
# Fitting per cell type (rather than one model with a celltype term) means each
# cell type gets its own library sizes, its own expression filter and its own
# dispersion — a satellite glia unit and an NF2 unit are not forced to share a
# mean-variance relationship.
#
# Fig C : 49 genes x 20 cell types, log2FC of axotomy peak (72h + 7d, pooled
#         over Crush/ScNT/SpNT) vs Naive, computed WITHIN each cell type
# Fig D : the same, restricted to the 9 neuronal subtypes
# Fig E : small multiples — one mini heatmap per cell type, genes x time, SpNT
# Table : per gene, its peak cell type / timepoint / log2FC  -> data/panel_celltype_peak.csv
# =====================================================================
suppressMessages({ library(Matrix); library(edgeR); library(limma); library(ggplot2) })
source("/hpf/projects/msalter/sghazis/painseq/analysis/00_common.R")

CLIP_L2FC   <- 3
PEAK_TIMES  <- c(72L, 168L)          # the axotomy peak window
FACET_MODEL <- "SpNT"                # closest analogue to SNI, for Fig E

fmt_time <- function(h) ifelse(h == 0, "0", ifelse(h < 168, paste0(h, "h"), paste0(h / 24, "d")))

pbsc <- readRDS(file.path(DATA, "pb_sample_celltype.rds"))
pg   <- readRDS(file.path(DATA, "panel_genes.rds"))
rows <- pg[order(pg$role != "panel", -pg$mean_sni), ]

# =====================================================================
# Per-cell-type fit: log2FC of every injury x time group vs that type's Naive
# =====================================================================
fit_one_celltype <- function(ct) {
  sel <- pbsc$meta$celltype == ct & pbsc$meta$n_cells >= MIN_CELLS
  md  <- pbsc$meta[sel, ]; cts <- pbsc$counts[, sel, drop = FALSE]
  if (sum(md$inj == "Naive") < 2) return(NULL)          # need a baseline with replication

  md$group <- factor(md$group)
  ok <- names(which(table(md$group) >= 1))
  md$group <- factor(as.character(md$group), levels = ok)

  y <- DGEList(counts = cts, group = md$group)
  keep <- filterByExpr(y, group = md$group)
  if (sum(keep) < 500) return(NULL)
  y <- calcNormFactors(y[keep, , keep.lib.sizes = FALSE], method = "TMM")
  lcpm <- cpm(y, log = TRUE, prior.count = 1)

  design <- model.matrix(~0 + group, data = md)
  colnames(design) <- make.names(colnames(design))
  qrd <- qr(design)
  if (qrd$rank < ncol(design)) design <- design[, qrd$pivot[seq_len(qrd$rank)], drop = FALSE]
  ref <- paste0("group", make.names("Naive_0"))
  if (!ref %in% colnames(design)) return(NULL)
  targets <- setdiff(colnames(design), ref)
  if (!length(targets)) return(NULL)

  fit <- lmFit(lcpm, design)
  cm  <- makeContrasts(contrasts = paste0(targets, " - ", ref), levels = design)
  colnames(cm) <- sub("^group", "", targets)
  fit2 <- eBayes(contrasts.fit(fit, cm))
  l2fc <- sapply(colnames(cm), function(cn) topTable(fit2, coef = cn, number = Inf, sort.by = "none")$logFC)
  padj <- sapply(colnames(cm), function(cn) topTable(fit2, coef = cn, number = Inf, sort.by = "none")$adj.P.Val)
  rownames(l2fc) <- rownames(padj) <- rownames(lcpm)

  cmeta <- unique(md[, c("group", "inj", "time")])
  cmeta$col <- make.names(as.character(cmeta$group))
  cmeta <- cmeta[cmeta$col %in% colnames(l2fc), ]
  cmeta <- cmeta[order(match(cmeta$inj, INJ_ORDER), cmeta$time), ]
  list(l2fc = l2fc, padj = padj, cmeta = cmeta, n_units = nrow(md),
       n_cells = sum(md$n_cells))
}

cts_all <- names(CELLTYPE_FULL)
msg("fitting ", length(cts_all), " cell types independently")
fits <- setNames(vector("list", length(cts_all)), cts_all)
for (ct in cts_all) {
  f <- tryCatch(fit_one_celltype(ct), error = function(e) { msg("  ", ct, ": FAILED - ", conditionMessage(e)); NULL })
  fits[[ct]] <- f
  msg("  ", formatC(ct, width = -18), if (is.null(f)) "skipped (too few cells / no Naive baseline)"
      else sprintf("%d units, %s cells, %d contrasts", f$n_units, format(f$n_cells, big.mark = ","), ncol(f$l2fc)))
}
fits <- Filter(Negate(is.null), fits)
saveRDS(fits, file.path(DATA, "panel_celltype_fits.rds"))

# =====================================================================
# Fig C / D — genes x cell types, axotomy peak vs Naive
# =====================================================================
peak_cols <- function(cmeta) cmeta$col[cmeta$inj %in% AXOTOMY & cmeta$time %in% PEAK_TIMES]

M_ct <- matrix(NA_real_, nrow(rows), length(fits), dimnames = list(rows$gene, names(fits)))
P_ct <- M_ct
for (ct in names(fits)) {
  f  <- fits[[ct]]; pc <- peak_cols(f$cmeta)
  if (!length(pc)) next
  hit <- rows$matrix_sym %in% rownames(f$l2fc)
  M_ct[hit, ct] <- rowMeans(f$l2fc[rows$matrix_sym[hit], pc, drop = FALSE], na.rm = TRUE)
  # summarise significance across the peak columns by the smallest adj.p
  P_ct[hit, ct] <- apply(f$padj[rows$matrix_sym[hit], pc, drop = FALSE], 1, min, na.rm = TRUE)
}
P_ct[!is.finite(P_ct)] <- NA_real_

heat_ct <- function(M, P, title, subtitle, w, h, stem, cluster_rows = FALSE) {
  ord <- if (cluster_rows) cluster_order(M, TRUE, FALSE)$row else seq_len(nrow(M))
  M <- M[ord, , drop = FALSE]; P <- P[ord, , drop = FALSE]
  rr <- rows[ord, ]
  d <- expand.grid(gene = factor(rownames(M), levels = rev(rownames(M))),
                   ct = factor(colnames(M), levels = colnames(M)), stringsAsFactors = FALSE)
  d$value <- as.vector(M); d$padj <- as.vector(P)
  d$vclip <- pmax(pmin(d$value, CLIP_L2FC), -CLIP_L2FC)
  d$sig   <- !is.na(d$padj) & d$padj < 0.05
  d$neuronal <- d$ct %in% NEURONAL_TYPES

  lab <- setNames(rr$gene, rr$gene)
  lab[rr$low_detection %in% TRUE] <- paste0(rr$gene[rr$low_detection %in% TRUE], " *")
  lab[rr$role == "reference"]     <- paste0(rr$gene[rr$role == "reference"], " ‡")

  gg <- ggplot(d, aes(x = ct, y = gene, fill = vclip)) +
    geom_tile(colour = "grey92", linewidth = 0.15) +
    geom_point(data = subset(d, sig), size = 0.5, colour = "grey15", alpha = .85, show.legend = FALSE) +
    scale_fill_gradient2(low = DIVERGE_HEX[1], mid = DIVERGE_HEX[2], high = DIVERGE_HEX[3],
                         midpoint = 0, limits = c(-CLIP_L2FC, CLIP_L2FC), na.value = "grey88",
                         name = "log2FC\nvs Naive") +
    scale_y_discrete(labels = function(g) lab[g]) +
    labs(title = title, subtitle = paste(strwrap(subtitle, width = 110), collapse = "\n"),
         x = NULL, y = NULL) +
    theme_pub(base_size = 11) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8.5),
          axis.text.y = element_text(size = 7.2),
          plot.subtitle = element_text(colour = "grey30", size = 8.5, lineheight = 1.15),
          legend.key.height = unit(26, "pt"))
  ggsave(file.path(FIGS, paste0(stem, ".png")), gg, width = w, height = h, dpi = 300, bg = "white")
  ggsave(file.path(FIGS, paste0(stem, ".pdf")), gg, width = w, height = h, device = cairo_pdf)
  msg("  wrote figs/", stem, ".{png,pdf}")
}

sub_ct <- sprintf(paste(
  "Axotomy peak (%s, pooled over Crush/ScNT/SpNT) vs Naive, fitted WITHIN each cell type.",
  "Dot = adj.p<0.05. * = detected in <0.5%% of cells. ‡ = RAG archetype, not one of the 47.",
  "grey = gene below that cell type's expression filter, or the cell type lacks that contrast."),
  paste(fmt_time(PEAK_TIMES), collapse = " + "))

heat_ct(M_ct, P_ct, "Conserved injury panel by cell type", sub_ct, 10.5, 9, "figC_panel_by_celltype")

neur <- intersect(NEURONAL_TYPES, colnames(M_ct))
heat_ct(M_ct[, neur, drop = FALSE], P_ct[, neur, drop = FALSE],
        "Conserved injury panel by sensory-neuron subtype",
        paste("Neuronal subtypes only.", sub_ct), 7.5, 9, "figD_panel_by_neuron_subtype")

# =====================================================================
# Fig E — small multiples: genes x time, one panel per cell type, SpNT
# =====================================================================
long <- do.call(rbind, lapply(names(fits), function(ct) {
  f  <- fits[[ct]]
  cm <- f$cmeta[f$cmeta$inj == FACET_MODEL, ]
  if (!nrow(cm)) return(NULL)
  hit <- rows$matrix_sym %in% rownames(f$l2fc)
  v <- matrix(NA_real_, nrow(rows), nrow(cm), dimnames = list(rows$gene, cm$col))
  p <- v
  v[hit, ] <- f$l2fc[rows$matrix_sym[hit], cm$col, drop = FALSE]
  p[hit, ] <- f$padj[rows$matrix_sym[hit], cm$col, drop = FALSE]
  data.frame(gene = rep(rownames(v), ncol(v)),
             col  = rep(cm$col, each = nrow(v)),
             time = rep(cm$time, each = nrow(v)),
             celltype = ct, value = as.vector(v), padj = as.vector(p),
             stringsAsFactors = FALSE)
}))
long$gene     <- factor(long$gene, levels = rev(rows$gene))
long$celltype <- factor(long$celltype, levels = intersect(names(CELLTYPE_FULL), unique(long$celltype)))
long$tlab     <- factor(fmt_time(long$time), levels = fmt_time(sort(unique(long$time))))
long$vclip    <- pmax(pmin(long$value, CLIP_L2FC), -CLIP_L2FC)

ggE <- ggplot(long, aes(x = tlab, y = gene, fill = vclip)) +
  geom_tile(colour = "grey93", linewidth = 0.1) +
  facet_wrap(~ celltype, nrow = 3) +
  scale_fill_gradient2(low = DIVERGE_HEX[1], mid = DIVERGE_HEX[2], high = DIVERGE_HEX[3],
                       midpoint = 0, limits = c(-CLIP_L2FC, CLIP_L2FC), na.value = "grey88",
                       name = "log2FC\nvs Naive") +
  labs(title = sprintf("Conserved injury panel — %s time-course within each cell type", FACET_MODEL),
       subtitle = paste(strwrap(paste(
         "Each panel is one cell type fitted independently against its own Naive baseline.",
         "grey = below that cell type's expression filter or no unit with >= 10 cells."), width = 150), collapse = "\n"),
       x = NULL, y = NULL) +
  theme_pub(base_size = 10) +
  theme(axis.text.x = element_text(angle = 90, vjust = .5, hjust = 1, size = 6.5),
        axis.text.y = element_text(size = 5.4),
        strip.background = element_rect(fill = "grey93", colour = "black"),
        strip.text = element_text(face = "bold", size = 8),
        plot.subtitle = element_text(colour = "grey30", size = 8),
        panel.spacing = unit(3, "pt"))
ggsave(file.path(FIGS, "figE_panel_celltype_timecourse.png"), ggE, width = 20, height = 13, dpi = 300, bg = "white")
ggsave(file.path(FIGS, "figE_panel_celltype_timecourse.pdf"), ggE, width = 20, height = 13, device = cairo_pdf)
msg("  wrote figs/figE_panel_celltype_timecourse.{png,pdf}")

# =====================================================================
# Peak table — where and when is each gene maximal?
# =====================================================================
peak <- do.call(rbind, lapply(seq_len(nrow(rows)), function(i) {
  g <- rows$gene[i]; sym <- rows$matrix_sym[i]
  best <- NULL
  for (ct in names(fits)) {
    f <- fits[[ct]]
    if (!sym %in% rownames(f$l2fc)) next
    v <- f$l2fc[sym, f$cmeta$col]; p <- f$padj[sym, f$cmeta$col]
    k <- which.max(v)
    if (!length(k) || !is.finite(v[k])) next
    if (is.null(best) || v[k] > best$log2FC)
      best <- data.frame(Gene = g, role = rows$role[i],
                         peak_celltype = ct, peak_celltype_full = unname(CELLTYPE_FULL[ct]),
                         peak_injury = as.character(f$cmeta$inj[k]), peak_time_h = f$cmeta$time[k],
                         peak_time = fmt_time(f$cmeta$time[k]),
                         log2FC = round(unname(v[k]), 3), adj_p = signif(unname(p[k]), 3),
                         pct_cells_detected = rows$pct_cells[i],
                         stringsAsFactors = FALSE)
  }
  if (is.null(best)) data.frame(Gene = g, role = rows$role[i], peak_celltype = NA, peak_celltype_full = NA,
      peak_injury = NA, peak_time_h = NA, peak_time = NA, log2FC = NA, adj_p = NA,
      pct_cells_detected = rows$pct_cells[i], stringsAsFactors = FALSE) else best
}))
peak <- peak[order(-peak$log2FC), ]
write.csv(peak, file.path(DATA, "panel_celltype_peak.csv"), row.names = FALSE)
msg("wrote data/panel_celltype_peak.csv")
print(utils::head(peak[peak$role == "panel", c("Gene","peak_celltype","peak_injury","peak_time","log2FC","adj_p")], 20))
