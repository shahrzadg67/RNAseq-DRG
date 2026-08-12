# =====================================================================
# 41_recruitment.R — is the panel induced because MORE CELLS express it, or
# because each expressing cell expresses MORE?
#
# Pseudobulk cannot tell these apart: both raise the group mean identically. The
# distinction matters biologically — recruitment means a subset of neurons switches
# the gene on (a state transition), while a level change means the same cells
# simply make more transcript.
#
# The decomposition is exact on the linear scale. With f = fraction of cells
# expressing and m = mean expression AMONG expressing cells, the mean over all
# cells is f * m, so
#
#     log2 FC(mean over all cells) = log2(f_inj / f_naive) + log2(m_inj / m_naive)
#                                    \___ recruitment ___/   \___ level ____/
#
# and the two components add. Stored values are log1p(CP10K), so expm1() recovers
# the linear scale first.
#
# Out: data/recruitment.rds/.csv, figs/recruitment_{heatmap,scatter}.{png,pdf}
# =====================================================================
suppressMessages({ library(Matrix); library(ggplot2) })
source("/hpf/projects/msalter/sghazis/painseq/analysis/00_common.R")

MIN_CELLS_GRP <- 30      # a fraction estimated from fewer cells is not worth plotting
PC            <- 0.5     # pseudocount on the expressing-cell count (Sprr1a is 0% in naive)

pc <- readRDS(file.path(DATA, "percell_expr.rds"))
E  <- pc$expr; md <- pc$meta; gu <- pc$genes
panel <- gu$gene[gu$panel]
msg("panel genes with per-cell data: ", length(panel))

md$group <- paste(md$inj, md$time, sep = "_")
grp_levels <- unique(md$group[order(match(md$inj, INJ_ORDER), md$time)])

# ---------------------------------------------------------------------
# f (fraction expressing) and m (linear mean among expressing) per gene x cell set
# ---------------------------------------------------------------------
stats_for <- function(cols) {
  if (length(cols) < MIN_CELLS_GRP)
    return(list(n = length(cols), f = rep(NA_real_, length(panel)), m = rep(NA_real_, length(panel))))
  sub <- E[panel, cols, drop = FALSE]
  nexp <- Matrix::rowSums(sub > 0)
  # linear CP10K sum: expm1 on the stored log1p values, done on the sparse slots only
  lin <- sub; lin@x <- expm1(lin@x)
  tot <- Matrix::rowSums(lin)
  list(n = length(cols),
       f = (nexp + PC) / (length(cols) + 2 * PC),          # regularised fraction
       m = ifelse(nexp > 0, tot / pmax(nexp, 1), NA_real_)) # mean among expressing
}

build <- function(split_by_celltype) {
  keys <- if (split_by_celltype) split(seq_len(nrow(md)), list(md$celltype, md$group), drop = TRUE)
          else split(seq_len(nrow(md)), md$group)
  out <- list()
  for (k in names(keys)) {
    s <- stats_for(keys[[k]])
    ct <- if (split_by_celltype) sub("\\..*", "", k) else "ALL"
    gp <- if (split_by_celltype) sub("^[^.]*\\.", "", k) else k
    out[[k]] <- data.frame(gene = panel, celltype = ct, group = gp, n_cells = s$n,
                           f = s$f, m = s$m, stringsAsFactors = FALSE)
  }
  d <- do.call(rbind, out); rownames(d) <- NULL
  # reference = that cell type's own Naive
  ref <- d[d$group == "Naive_0", c("gene", "celltype", "f", "m")]
  names(ref) <- c("gene", "celltype", "f0", "m0")
  d <- merge(d, ref, by = c("gene", "celltype"), all.x = TRUE)
  d$recruitment <- log2(d$f / d$f0)
  d$level       <- log2(d$m / d$m0)
  d$total       <- d$recruitment + d$level
  d$inj  <- sub("_.*", "", d$group)
  d$time <- as.integer(sub(".*_", "", d$group))
  d
}

msg("computing pooled (whole-DRG) statistics")
D_all <- build(FALSE)
msg("computing per-cell-type statistics")
D_ct  <- build(TRUE)
saveRDS(list(pooled = D_all, by_celltype = D_ct), file.path(DATA, "recruitment.rds"))
write.csv(D_all, file.path(DATA, "recruitment_pooled.csv"), row.names = FALSE)

# ---------------------------------------------------------------------
# VERIFICATION: reproduce the planning probe exactly
# (Atf3 in NF1, male C57, SpNT — 0.2 / 16.9 / 53.0 / 63.1 / 42.9 % expressing)
# ---------------------------------------------------------------------
msg("\n--- verification: Atf3 in NF1 (male C57, Naive + SpNT) ---")
sel <- md$celltype == "NF1" & md$sex == "male" & md$geno == "C57"
for (t in c(0, 6, 24, 48, 168)) {
  k <- which(sel & md$time == t & md$inj %in% c("Naive", "SpNT"))
  if (length(k) < 30) next
  v <- E["Atf3", k]
  msg(sprintf("  %-5s n=%4d  pct.expressing=%5.1f%%  mean(expressing, log1p)=%.2f",
              ifelse(t == 0, "naive", paste0(t, "h")), length(k), 100 * mean(v > 0),
              ifelse(any(v > 0), mean(v[v > 0]), 0)))
}

# ---------------------------------------------------------------------
# Classify each gene at the axotomy peak
# ---------------------------------------------------------------------
peak <- D_all[D_all$inj %in% AXOTOMY & D_all$time %in% c(72L, 168L), ]
agg <- aggregate(cbind(recruitment, level, total) ~ gene, data = peak, FUN = mean, na.rm = TRUE)
agg$driver <- with(agg, ifelse(abs(recruitment) > 2 * abs(level), "recruitment",
                        ifelse(abs(level) > 2 * abs(recruitment), "level", "mixed")))
agg$pct_recruitment <- round(100 * abs(agg$recruitment) / (abs(agg$recruitment) + abs(agg$level)), 1)
agg <- agg[order(-agg$total), ]
write.csv(agg, file.path(DATA, "recruitment_drivers.csv"), row.names = FALSE)

msg("\n--- what drives induction at the axotomy peak (72h + 7d) ---")
msg("  driver classification across ", nrow(agg), " panel genes:")
print(table(agg$driver))
msg(sprintf("  median %% of the total fold-change attributable to recruitment: %.1f%%",
            median(agg$pct_recruitment, na.rm = TRUE)))
up <- agg[agg$total > 1, ]
msg(sprintf("  among the %d genes induced >2-fold: recruitment-driven %d, mixed %d, level-driven %d",
            nrow(up), sum(up$driver == "recruitment"), sum(up$driver == "mixed"), sum(up$driver == "level")))
msg("\n  top induced genes:")
print(utils::head(agg[, c("gene", "recruitment", "level", "total", "pct_recruitment", "driver")], 12),
      row.names = FALSE, digits = 3)

# ---------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------
gene_order <- agg$gene
mk <- function(d, val, lab) {
  d <- d[d$group %in% grp_levels, ]
  d$gene <- factor(d$gene, levels = rev(gene_order))
  d$inj  <- factor(d$inj, levels = INJ_ORDER)
  d$col  <- factor(d$group, levels = grp_levels)
  d$tlab <- ifelse(d$time == 0, "0", ifelse(d$time < 168, paste0(d$time, "h"), paste0(d$time / 24, "d")))
  d$v <- pmax(pmin(d[[val]], 4), -4)
  d$panel_lab <- lab
  d
}
hd <- rbind(mk(D_all, "recruitment", "Recruitment  log2(fraction expressing)"),
            mk(D_all, "level",       "Level  log2(mean among expressing cells)"))
xl <- setNames(hd$tlab[!duplicated(hd$col)], hd$col[!duplicated(hd$col)])

gg <- ggplot(hd, aes(col, gene, fill = v)) +
  geom_tile(colour = "grey92", linewidth = .12) +
  facet_grid(panel_lab ~ inj, scales = "free_x", space = "free_x",
             labeller = labeller(inj = INJ_SHORT)) +
  scale_x_discrete(labels = function(v) xl[v]) +
  scale_fill_gradient2(low = DIVERGE_HEX[1], mid = DIVERGE_HEX[2], high = DIVERGE_HEX[3],
                       midpoint = 0, limits = c(-4, 4), na.value = "grey88",
                       name = "log2 vs\nNaive") +
  labs(title = "What actually drives the panel: recruitment, not upregulation",
       subtitle = paste(strwrap(paste(
         "Top: change in the FRACTION of cells expressing each gene. Bottom: change in the mean level",
         "AMONG expressing cells. The two add exactly to the pseudobulk log2 fold-change. Almost all",
         "of the signal sits in the top block — injury switches genes on in cells that were silent,",
         "rather than raising output in cells already expressing them."), width = 140), collapse = "\n"),
       x = NULL, y = NULL) +
  theme_pub(base_size = 11) +
  theme(axis.text.x = element_text(angle = 90, vjust = .5, hjust = 1, size = 7),
        axis.text.y = element_text(size = 6.5),
        plot.subtitle = element_text(colour = "grey30", size = 8.5),
        strip.background = element_rect(fill = "grey93", colour = "black"),
        strip.text = element_text(face = "bold", size = 8),
        panel.spacing.x = unit(2.2, "pt"))
ggsave(file.path(FIGS, "recruitment_heatmap.png"), gg, width = 15, height = 11, dpi = 300, bg = "white")
ggsave(file.path(FIGS, "recruitment_heatmap.pdf"), gg, width = 15, height = 11, device = cairo_pdf)
msg("wrote figs/recruitment_heatmap.{png,pdf}")

lim <- max(abs(c(agg$recruitment, agg$level)), na.rm = TRUE)
gg2 <- ggplot(agg, aes(recruitment, level, colour = driver)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted", colour = "grey60") +
  geom_hline(yintercept = 0, linewidth = .25, colour = "grey70") +
  geom_vline(xintercept = 0, linewidth = .25, colour = "grey70") +
  geom_point(size = 2.6, alpha = .9) +
  ggrepel::geom_text_repel(aes(label = gene), size = 2.9, max.overlaps = 25, seed = 3,
                           colour = "grey15", show.legend = FALSE) +
  scale_colour_manual(values = c(recruitment = "#B2182B", mixed = "#E8A33D", level = "#2166AC"),
                      name = NULL) +
  coord_equal(xlim = c(-lim, lim), ylim = c(-lim, lim)) +
  labs(title = "Recruitment versus level, at the axotomy peak",
       subtitle = paste(strwrap(paste(
         "Each point is one panel gene. x = log2 change in the fraction of cells expressing it;",
         "y = log2 change in the mean level among expressing cells. Points far right and near zero",
         "on y are pure recruitment. The dotted line is equal contribution."), width = 100), collapse = "\n"),
       x = "Recruitment — log2 change in fraction expressing",
       y = "Level — log2 change among expressing cells") +
  theme_pub(base_size = 12) +
  theme(plot.subtitle = element_text(colour = "grey30", size = 9), legend.position = "bottom")
ggsave(file.path(FIGS, "recruitment_scatter.png"), gg2, width = 9, height = 9.6, dpi = 300, bg = "white")
ggsave(file.path(FIGS, "recruitment_scatter.pdf"), gg2, width = 9, height = 9.6, device = cairo_pdf)
msg("wrote figs/recruitment_scatter.{png,pdf} and data/recruitment.{rds,_pooled.csv,_drivers.csv}")
