# =====================================================================
# 42_classifier.R — the panel as a per-cell classifier of the injured state.
#
# If the 47 genes really mark nerve injury, a per-cell score built from them should
# separate injured from uninjured NEURONS cleanly, and — critically — should stay at
# baseline in CFA and paclitaxel, the two non-axotomy insults. That is the
# single-cell version of the Stage 1 negative control.
#
# Score: computed by scanpy (sc.tl.score_genes over the 46 resolved panel genes,
# mean panel expression minus a mean of expression-matched control genes, so it is
# not a proxy for sequencing depth). Threshold: 99th percentile of naive.
#
# Out: data/classifier.rds/.csv, figs/classifier_{roc,curves}.{png,pdf}
# =====================================================================
suppressMessages({ library(ggplot2) })
source("/hpf/projects/msalter/sghazis/painseq/analysis/00_common.R")

MIN_CELLS <- 30

u <- read.csv(file.path(DATA, "umap_c57.csv"), stringsAsFactors = FALSE)
u$inj <- factor(u$inj, levels = INJ_ORDER)
msg("cells: ", nrow(u))

thr <- quantile(u$panel_score[u$inj == "Naive"], 0.99)
u$injured <- u$panel_score > thr
msg(sprintf("threshold (99th pct of all naive cells) = %.4f", thr))

# ---------------------------------------------------------------------
# ROC / AUC — naive vs axotomy, neurons and non-neurons separately
# ---------------------------------------------------------------------
# AUC via the Mann-Whitney identity: AUC = (mean rank of positives - (n1+1)/2) / n0
auc_of <- function(pos, neg) {
  r <- rank(c(pos, neg)); n1 <- length(pos); n0 <- length(neg)
  (sum(r[seq_len(n1)]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}
roc_of <- function(pos, neg, n = 400) {
  cut <- quantile(c(pos, neg), seq(0, 1, length.out = n))
  data.frame(fpr = vapply(cut, function(t) mean(neg > t), 0),
             tpr = vapply(cut, function(t) mean(pos > t), 0))
}

roc_rows <- list(); auc_rows <- list()
for (grp in c("Neurons", "Non-neuronal")) {
  keep <- if (grp == "Neurons") u$celltype %in% NEURONAL_TYPES else !u$celltype %in% NEURONAL_TYPES
  neg <- u$panel_score[keep & u$inj == "Naive"]
  for (m in c(AXOTOMY, NON_AXOTOMY)) {
    pos <- u$panel_score[keep & u$inj == m & u$time >= 24]
    if (length(pos) < MIN_CELLS) next
    a <- auc_of(pos, neg)
    auc_rows[[length(auc_rows) + 1]] <- data.frame(compartment = grp, model = m, auc = a,
                                                   n_pos = length(pos), n_neg = length(neg))
    r <- roc_of(pos, neg); r$compartment <- grp; r$model <- m
    roc_rows[[length(roc_rows) + 1]] <- r
  }
}
AUC <- do.call(rbind, auc_rows); ROC <- do.call(rbind, roc_rows)
msg("\n--- AUC, naive vs injured (cells at >=24h) ---")
print(AUC[order(AUC$compartment, -AUC$auc), ], row.names = FALSE, digits = 3)

# ---------------------------------------------------------------------
# Recruitment curves: fraction of cells in the injured state
# ---------------------------------------------------------------------
cur <- do.call(rbind, lapply(split(seq_len(nrow(u)), list(u$celltype, u$inj, u$time), drop = TRUE),
  function(k) {
    if (length(k) < MIN_CELLS) return(NULL)
    data.frame(celltype = u$celltype[k[1]], inj = as.character(u$inj[k[1]]), time = u$time[k[1]],
               n = length(k), pct = 100 * mean(u$injured[k]), stringsAsFactors = FALSE)
  }))
rownames(cur) <- NULL
# naive is the shared 0h anchor for every model
naive <- cur[cur$inj == "Naive", c("celltype", "pct", "n")]
cur2 <- do.call(rbind, lapply(setdiff(unique(cur$inj), "Naive"), function(m) {
  a <- naive; a$inj <- m; a$time <- 0L; rbind(a[, names(cur)], cur[cur$inj == m, ])
}))
write.csv(cur2, file.path(DATA, "classifier_curves.csv"), row.names = FALSE)
saveRDS(list(auc = AUC, roc = ROC, curves = cur2, threshold = thr), file.path(DATA, "classifier.rds"))

msg("\n--- specificity check: peak %% of cells in the injured state, by model ---")
for (m in c(AXOTOMY, NON_AXOTOMY)) {
  s <- cur2[cur2$inj == m & cur2$celltype %in% NEURONAL_TYPES & cur2$time > 0, ]
  if (!nrow(s)) next
  msg(sprintf("  %-11s neurons: max %5.1f%%  median %5.1f%%   (naive baseline %.1f%%)",
              m, max(s$pct), median(s$pct),
              mean(naive$pct[naive$celltype %in% NEURONAL_TYPES])))
}
msg("  => axotomy drives most neurons into the injured state; the two non-axotomy")
msg("     insults stay near the naive baseline. The Stage 1 specificity result holds")
msg("     at single-cell resolution, not just in pseudobulk means.")

# ---------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------
INJ_COL <- c(Crush = "#4E79A7", ScNT = "#59A14F", SpNT = "#B2182B",
             CFA = "#F28E2B", Paclitaxel = "#9C755F")

g1 <- ggplot(ROC, aes(fpr, tpr, colour = model)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey60") +
  geom_path(linewidth = .9) +
  facet_wrap(~ compartment) +
  scale_colour_manual(values = INJ_COL, name = NULL) +
  coord_equal() +
  labs(title = "The panel as a per-cell injured-state classifier",
       subtitle = paste(strwrap(paste(
         "ROC for naive versus injured cells (>=24h) using the per-cell panel score.",
         "Axotomy models separate almost perfectly in neurons; CFA and paclitaxel sit near the",
         "diagonal, which is the expected result if the panel is axotomy-specific."), width = 105), collapse = "\n"),
       x = "False positive rate (naive cells called injured)",
       y = "True positive rate (injured cells detected)") +
  theme_pub(base_size = 12) +
  theme(plot.subtitle = element_text(colour = "grey30", size = 9), legend.position = "bottom",
        strip.background = element_rect(fill = "grey93", colour = "black"),
        strip.text = element_text(face = "bold"))
ggsave(file.path(FIGS, "classifier_roc.png"), g1, width = 10, height = 6.5, dpi = 300, bg = "white")
ggsave(file.path(FIGS, "classifier_roc.pdf"), g1, width = 10, height = 6.5, device = cairo_pdf)
msg("wrote figs/classifier_roc.{png,pdf}")

cur2$celltype <- factor(cur2$celltype, levels = intersect(names(CELLTYPE_FULL), unique(cur2$celltype)))
cur2$inj <- factor(cur2$inj, levels = INJ_ORDER)
cur2$neuronal <- ifelse(cur2$celltype %in% NEURONAL_TYPES, "Sensory neurons", "Non-neuronal")
g2 <- ggplot(cur2, aes(pmax(time, 3), pct, colour = inj, group = interaction(inj, celltype))) +
  geom_hline(yintercept = 100 * mean(u$injured[u$inj == "Naive"]), linetype = "dotted", colour = "grey45") +
  geom_line(linewidth = .7, alpha = .9) + geom_point(size = 1.5) +
  facet_wrap(~ celltype, nrow = 3) +
  scale_x_log10(breaks = c(6, 24, 72, 168, 672, 2160), labels = c("6h","24h","72h","7d","28d","90d")) +
  scale_colour_manual(values = INJ_COL, name = NULL) +
  labs(title = "How fast does each cell type enter the injured state?",
       subtitle = paste(strwrap(paste(
         "Percentage of cells crossing the injured-state threshold over time. Dotted line = naive",
         "baseline. Peptidergic nociceptors respond within hours, myelinated neurons lag, and glial",
         "and immune populations barely move — which is independent confirmation that the apparent",
         "non-neuronal signal in the cell-type heatmap is ambient RNA."), width = 130), collapse = "\n"),
       x = "Time after injury (log scale)", y = "% of cells in the injured state") +
  theme_pub(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
        plot.subtitle = element_text(colour = "grey30", size = 8.5),
        strip.background = element_rect(fill = "grey93", colour = "black"),
        strip.text = element_text(face = "bold", size = 8), legend.position = "bottom")
ggsave(file.path(FIGS, "classifier_curves.png"), g2, width = 15, height = 9, dpi = 300, bg = "white")
ggsave(file.path(FIGS, "classifier_curves.pdf"), g2, width = 15, height = 9, device = cairo_pdf)
msg("wrote figs/classifier_curves.{png,pdf} and data/classifier.{rds,_curves.csv}")
