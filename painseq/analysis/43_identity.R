# =====================================================================
# 43_identity.R — does the conserved injury panel mark loss of subtype identity?
#
# Renthal et al.'s central published claim is that injured sensory neurons LOSE
# their subtype identity and converge on a common injured state. If the 47-gene
# panel marks that transition, then within injured neurons the panel score should
# rise as subtype identity falls — a negative per-cell correlation.
#
# Identity score: each cell is scored on ITS OWN subtype's markers, where markers
# were derived in 40_extract_percell.R from NAIVE cells only by detection-rate
# difference. Deriving them from naive cells matters: markers picked using injured
# cells would already encode the answer.
#
# The score is standardised within subtype against that subtype's naive
# distribution, so "identity = 0" means a typical naive cell of that subtype and
# negative values mean the cell has drifted away from it.
#
# Out: data/identity.rds/.csv, figs/identity_{scatter,trajectory}.{png,pdf}
# =====================================================================
suppressMessages({ library(Matrix); library(ggplot2) })
source("/hpf/projects/msalter/sghazis/painseq/analysis/00_common.R")

pc   <- readRDS(file.path(DATA, "percell_expr.rds"))
mark <- readRDS(file.path(DATA, "subtype_markers.rds"))
E    <- pc$expr; md <- pc$meta
u    <- read.csv(file.path(DATA, "umap_c57.csv"), stringsAsFactors = FALSE)
stopifnot(nrow(u) == nrow(md))
md$panel_score <- u$panel_score
md$inj <- factor(as.character(md$inj), levels = INJ_ORDER)

# ---- per-cell identity score --------------------------------------------
# mean expression of the cell's own-subtype markers, minus the mean over all
# markers of all other subtypes (so a cell drifting toward a generic state, or
# toward another subtype, both score low)
all_marks <- unique(unlist(mark))
md$identity <- NA_real_
for (ct in names(mark)) {
  k <- which(md$celltype == ct)
  if (!length(k)) next
  own <- intersect(mark[[ct]], rownames(E))
  oth <- setdiff(intersect(all_marks, rownames(E)), own)
  if (length(own) < 3) next
  md$identity[k] <- Matrix::colMeans(E[own, k, drop = FALSE]) -
                    Matrix::colMeans(E[oth, k, drop = FALSE])
}
# standardise within subtype against that subtype's NAIVE cells
for (ct in unique(md$celltype)) {
  k <- which(md$celltype == ct); n <- which(md$celltype == ct & md$inj == "Naive")
  if (length(n) < 30 || all(is.na(md$identity[k]))) next
  mu <- mean(md$identity[n], na.rm = TRUE); sdv <- sd(md$identity[n], na.rm = TRUE)
  if (!is.finite(sdv) || sdv == 0) next
  md$identity[k] <- (md$identity[k] - mu) / sdv
}
msg("identity score computed for ", sum(!is.na(md$identity)), " of ", nrow(md), " cells")

# keep the row position in cell_meta order so downstream artifacts can align by
# position rather than by cell name (GEO names are not unique)
md$row <- seq_len(nrow(md))
neur <- md[md$celltype %in% NEURONAL_TYPES & !is.na(md$identity), ]
msg("neurons with both scores: ", nrow(neur))

# ---- correlation, per injury x time --------------------------------------
cc <- do.call(rbind, lapply(split(seq_len(nrow(neur)), list(neur$inj, neur$time), drop = TRUE),
  function(k) {
    if (length(k) < 100) return(NULL)
    data.frame(inj = as.character(neur$inj[k[1]]), time = neur$time[k[1]], n = length(k),
               r = cor(neur$identity[k], neur$panel_score[k]),
               mean_identity = mean(neur$identity[k]),
               mean_panel = mean(neur$panel_score[k]), stringsAsFactors = FALSE)
  }))
rownames(cc) <- NULL
cc <- cc[order(match(cc$inj, INJ_ORDER), cc$time), ]
write.csv(cc, file.path(DATA, "identity_by_group.csv"), row.names = FALSE)

msg("\n--- identity vs panel score, per injury x time (neurons) ---")
print(cc, row.names = FALSE, digits = 3)

ax <- neur[neur$inj %in% AXOTOMY & neur$time >= 24, ]
nv <- neur[neur$inj == "Naive", ]
r_ax <- cor(ax$identity, ax$panel_score); r_nv <- cor(nv$identity, nv$panel_score)
msg(sprintf("\n  pooled: naive r = %+.3f (n=%d) | axotomy >=24h r = %+.3f (n=%d)",
            r_nv, nrow(nv), r_ax, nrow(ax)))
msg(sprintf("  mean identity score: naive %+.2f -> axotomy >=24h %+.2f (%.1f SD of the naive distribution)",
            mean(nv$identity), mean(ax$identity), mean(ax$identity) - mean(nv$identity)))
if (r_ax < -0.2) {
  msg("  => CONFIRMED: within injured neurons, the higher the panel score the further the")
  msg("     cell has drifted from its own subtype identity. The panel marks the loss of")
  msg("     subtype identity that Renthal et al. describe.")
} else if (r_ax > 0.2) {
  msg("  => The correlation is POSITIVE, which contradicts the expectation. Panel-high")
  msg("     cells retain subtype identity; the panel is not a de-differentiation marker.")
} else {
  msg("  => NULL RESULT: essentially no per-cell relationship between panel score and")
  msg("     subtype identity. The panel tracks injury without tracking identity loss, so")
  msg("     the two are separable processes. Worth reporting as-is.")
}

saveRDS(list(by_group = cc, r_naive = r_nv, r_axotomy = r_ax,
             cells = neur[, c("row", "cell", "celltype", "inj", "time", "identity", "panel_score")]),
        file.path(DATA, "identity.rds"), compress = "xz")

# ---- figures --------------------------------------------------------------
pl <- neur[neur$inj %in% c("Naive", "SpNT", "ScNT", "Crush"), ]
pl$tlab <- factor(ifelse(pl$time == 0, "naive",
                  ifelse(pl$time < 168, paste0(pl$time, "h"), paste0(pl$time / 24, "d"))),
                  levels = c("naive", "6h", "12h", "24h", "36h", "48h", "72h",
                             "7d", "14d", "28d", "60d", "90d"))
keep_t <- c("naive", "24h", "72h", "7d", "28d")
pl <- pl[pl$tlab %in% keep_t, ]; pl$tlab <- droplevels(pl$tlab)

g1 <- ggplot(pl, aes(identity, panel_score)) +
  geom_hex(bins = 55) +
  geom_smooth(method = "lm", se = FALSE, colour = "#B2182B", linewidth = .7, formula = y ~ x) +
  facet_grid(inj ~ tlab) +
  scale_fill_gradient(low = "#EEF2F6", high = "#1F3B54", trans = "log10", name = "cells") +
  labs(title = "Does the panel mark loss of subtype identity?",
       subtitle = paste(strwrap(paste(
         "Each hexagon is a bin of sensory neurons. x = how strongly the cell still expresses its OWN",
         "subtype's naive markers, in SD of that subtype's naive distribution. y = 47-gene panel score.",
         "A downward slope means panel-high cells have drifted from their subtype identity."), width = 128), collapse = "\n"),
       x = "Subtype identity score (SD from that subtype's naive mean)", y = "Panel score") +
  theme_pub(base_size = 11) +
  theme(plot.subtitle = element_text(colour = "grey30", size = 8.5),
        strip.background = element_rect(fill = "grey93", colour = "black"),
        strip.text = element_text(face = "bold", size = 8.5))
ggsave(file.path(FIGS, "identity_scatter.png"), g1, width = 14, height = 9, dpi = 300, bg = "white")
ggsave(file.path(FIGS, "identity_scatter.pdf"), g1, width = 14, height = 9, device = cairo_pdf)
msg("wrote figs/identity_scatter.{png,pdf}")

cc$inj <- factor(cc$inj, levels = INJ_ORDER)
g2 <- ggplot(cc[cc$inj != "Naive", ], aes(pmax(time, 3), mean_identity, colour = inj)) +
  geom_hline(yintercept = 0, linetype = "dotted", colour = "grey45") +
  geom_line(linewidth = .8) + geom_point(size = 2) +
  scale_x_log10(breaks = c(6, 24, 72, 168, 672, 2160), labels = c("6h","24h","72h","7d","28d","90d")) +
  scale_colour_manual(values = c(Crush = "#4E79A7", ScNT = "#59A14F", SpNT = "#B2182B",
                                 CFA = "#F28E2B", Paclitaxel = "#9C755F"), name = NULL) +
  labs(title = "Subtype identity over time after injury",
       subtitle = "Mean identity score of sensory neurons. Zero = a typical naive cell of the same subtype.",
       x = "Time after injury (log scale)", y = "Mean subtype identity (SD from naive)") +
  theme_pub(base_size = 12) +
  theme(plot.subtitle = element_text(colour = "grey30", size = 9), legend.position = "bottom")
ggsave(file.path(FIGS, "identity_trajectory.png"), g2, width = 9, height = 6.5, dpi = 300, bg = "white")
ggsave(file.path(FIGS, "identity_trajectory.pdf"), g2, width = 9, height = 6.5, device = cairo_pdf)
msg("wrote figs/identity_trajectory.{png,pdf} and data/identity.{rds,_by_group.csv}")
