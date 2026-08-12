# =====================================================================
# 45_atf3_percell.R — the Atf3-KO experiment, done at the resolution it deserves.
#
# WHY THIS EXISTS. 32_atf3_dependency.R analysed this arm as whole-neuron
# pseudobulk and concluded "underpowered": with n=3 samples per genotype x state,
# the genotype x injury interaction reached significance for only 3 of 25 induced
# genes. That is a real limit of the SAMPLE-level test — but it is the wrong test
# for this dataset.
#
# 41_recruitment.R showed that injury induction is overwhelmingly RECRUITMENT:
# genes switch on in cells that were silent, rather than rising in cells already
# expressing them. Recruitment is a per-CELL quantity. With 17,665 neurons, the
# fraction of cells expressing a gene is estimated to within a fraction of a
# percent, even though the sample-level test has n=3. So the right question is not
# "is the interaction significant" but:
#
#     Does losing Atf3 stop neurons from ENTERING the injured state?
#
# Three analyses, all per-cell:
#   1. Injured-state entry — fraction of neurons crossing the panel-score threshold,
#      WT vs KO, per subtype. A proportion test over thousands of cells.
#   2. Per-gene recruitment — does the KO recruit fewer cells for each panel gene?
#   3. Per-subtype breakdown — the previous script pooled all 9 neuronal subtypes.
#
# ALSO records the allele caveat quantitatively: Atf3 mRNA is not reduced in the KO.
#
# Out: data/atf3_percell.rds/.csv, figs/atf3_{entry,recruitment}.{png,pdf}
# =====================================================================
suppressMessages({ library(Matrix); library(ggplot2) })
source("/hpf/projects/msalter/sghazis/painseq/analysis/00_common.R")

MIN_CELLS <- 40

pc  <- readRDS(file.path(DATA, "percell_expr_atf3.rds"))
E   <- pc$expr; md <- pc$meta; gu <- pc$genes
panel <- gu$gene[gu$panel]
md$geno  <- factor(md$geno, levels = c("Atf3-WT", "Atf3-KO"))
md$state <- factor(ifelse(md$inj == "Naive", "Naive", paste0("Crush ", md$time, "h")),
                   levels = c("Naive", "Crush 36h", "Crush 168h"))
msg("cells: ", nrow(md), " | ", length(panel), " panel genes | subtypes: ",
    paste(sort(unique(md$celltype)), collapse = ", "))
print(table(md$geno, md$state))

# =====================================================================
# 0. The allele caveat, stated in numbers
# =====================================================================
msg("\n--- Atf3 transcript itself ---")
for (g in c("Atf3-WT", "Atf3-KO")) for (s in levels(md$state)) {
  k <- which(md$geno == g & md$state == s)
  if (!length(k)) next
  v <- E["Atf3", k]
  msg(sprintf("  %-8s %-11s n=%5d  pct.expressing=%5.1f%%  mean(expressing)=%.2f",
              g, s, length(k), 100 * mean(v > 0), ifelse(any(v > 0), mean(v[v > 0]), 0)))
}
msg("  => Atf3 mRNA is present, and induced, in the KO. The allele still yields a")
msg("     transcript that 3'-end snRNA-seq counts (an internal-exon deletion leaves the")
msg("     3' end intact), so transcript level cannot verify the knockout. Everything")
msg("     below is conditional on the published genotype labels being correct.")

# =====================================================================
# 1. Injured-state entry — per-cell panel score
# ---------------------------------------------------------------------
# Score each cell the same way scanpy did for the C57 atlas: mean expression of the
# panel genes minus the mean of a matched random background set, so the score is not
# just a proxy for sequencing depth.
# =====================================================================
set.seed(5)
bg_pool <- setdiff(rownames(E), panel)
bg <- sample(bg_pool, min(length(bg_pool), 5 * length(panel)))
score <- Matrix::colMeans(E[panel, , drop = FALSE]) - Matrix::colMeans(E[bg, , drop = FALSE])
md$score <- as.numeric(score)

# threshold from WT naive, exactly as for the C57 atlas
thr <- quantile(md$score[md$geno == "Atf3-WT" & md$state == "Naive"], 0.99)
md$injured <- md$score > thr
msg(sprintf("\n--- injured-state entry (threshold = 99th pct of WT naive = %.3f) ---", thr))

entry <- do.call(rbind, lapply(split(seq_len(nrow(md)), list(md$geno, md$state, md$celltype), drop = TRUE),
  function(k) {
    if (length(k) < MIN_CELLS) return(NULL)
    data.frame(geno = md$geno[k[1]], state = md$state[k[1]], celltype = md$celltype[k[1]],
               n = length(k), n_injured = sum(md$injured[k]),
               pct = 100 * mean(md$injured[k]), stringsAsFactors = FALSE)
  }))
rownames(entry) <- NULL

# pooled over subtypes, per genotype x state, with a proportion test
pooled <- do.call(rbind, lapply(split(seq_len(nrow(md)), list(md$geno, md$state), drop = TRUE),
  function(k) data.frame(geno = md$geno[k[1]], state = md$state[k[1]], n = length(k),
                         n_injured = sum(md$injured[k]), pct = 100 * mean(md$injured[k]))))
rownames(pooled) <- NULL
print(pooled, row.names = FALSE, digits = 3)

for (s in c("Crush 36h", "Crush 168h")) {
  w <- pooled[pooled$state == s, ]
  if (nrow(w) != 2) next
  tt <- prop.test(w$n_injured, w$n)
  msg(sprintf("  %s: WT %.1f%% vs KO %.1f%% injured  (diff %+.1f pp, p = %.3g, n = %d/%d cells)",
              s, w$pct[w$geno == "Atf3-WT"], w$pct[w$geno == "Atf3-KO"],
              w$pct[w$geno == "Atf3-KO"] - w$pct[w$geno == "Atf3-WT"],
              tt$p.value, w$n[w$geno == "Atf3-WT"], w$n[w$geno == "Atf3-KO"]))
}
msg("  NOTE: a proportion test over thousands of cells treats cells as independent,")
msg("  which overstates significance when they come from 1-3 mice. Read the effect")
msg("  SIZE (percentage-point difference) as the result; the p-value is optimistic.")

# =====================================================================
# 2. Per-gene recruitment, WT vs KO
# =====================================================================
frac_of <- function(k, genes) (Matrix::rowSums(E[genes, k, drop = FALSE] > 0) + 0.5) / (length(k) + 1)
rec <- do.call(rbind, lapply(c("Crush 36h", "Crush 168h"), function(s) {
  kw <- which(md$geno == "Atf3-WT" & md$state == s); kk <- which(md$geno == "Atf3-KO" & md$state == s)
  nw <- which(md$geno == "Atf3-WT" & md$state == "Naive"); nk <- which(md$geno == "Atf3-KO" & md$state == "Naive")
  if (!length(kw) || !length(kk)) return(NULL)
  data.frame(gene = panel, state = s,
             wt_rec = log2(frac_of(kw, panel) / frac_of(nw, panel)),
             ko_rec = log2(frac_of(kk, panel) / frac_of(nk, panel)),
             wt_pct = 100 * frac_of(kw, panel), ko_pct = 100 * frac_of(kk, panel),
             stringsAsFactors = FALSE)
}))
rec$blocked <- rec$wt_rec - rec$ko_rec          # recruitment lost in the KO
rec <- rec[order(-rec$blocked), ]
write.csv(rec, file.path(DATA, "atf3_percell_recruitment.csv"), row.names = FALSE)

msg("\n--- per-gene recruitment, WT vs KO ---")
for (s in unique(rec$state)) {
  r <- rec[rec$state == s, ]
  ind <- r$wt_rec > 1
  msg(sprintf("  %s: of %d genes recruited in WT, median recruitment retained in KO = %.0f%%",
              s, sum(ind), 100 * median(pmax(r$ko_rec[ind], 0) / r$wt_rec[ind], na.rm = TRUE)))
  msg(sprintf("    most Atf3-blocked: %s",
              paste(utils::head(r$gene[ind][order(-r$blocked[ind])], 8), collapse = ", ")))
}

saveRDS(list(entry = entry, pooled = pooled, recruitment = rec, threshold = thr,
             score = md[, c("cell", "geno", "state", "celltype", "score", "injured")]),
        file.path(DATA, "atf3_percell.rds"), compress = "xz")

# =====================================================================
# Figures
# =====================================================================
entry$celltype <- factor(entry$celltype, levels = intersect(names(CELLTYPE_FULL), unique(entry$celltype)))
gg <- ggplot(entry, aes(state, pct, fill = geno)) +
  geom_col(position = position_dodge(preserve = "single"), width = .72) +
  geom_text(aes(label = sprintf("%.0f", pct)), position = position_dodge(width = .72),
            vjust = -0.3, size = 2.5, colour = "grey20") +
  facet_wrap(~ celltype, nrow = 2) +
  scale_fill_manual(values = c(`Atf3-WT` = "#4E79A7", `Atf3-KO` = "#E15759"), name = NULL) +
  labs(title = "Does losing Atf3 stop neurons entering the injured state?",
       subtitle = paste(strwrap(paste(
         "Percentage of neurons crossing the injured-state threshold (99th percentile of WT naive),",
         "by subtype. Estimated per cell from 17,665 neurons, so it does not depend on the n=3",
         "sample-level interaction test that the pseudobulk analysis found underpowered."), width = 125), collapse = "\n"),
       x = NULL, y = "% of neurons in the injured state") +
  theme_pub(base_size = 11) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 8),
        plot.subtitle = element_text(colour = "grey30", size = 8.5),
        strip.background = element_rect(fill = "grey93", colour = "black"),
        strip.text = element_text(face = "bold", size = 8.5),
        legend.position = "bottom")
ggsave(file.path(FIGS, "atf3_entry.png"), gg, width = 12, height = 7.5, dpi = 300, bg = "white")
ggsave(file.path(FIGS, "atf3_entry.pdf"), gg, width = 12, height = 7.5, device = cairo_pdf)
msg("wrote figs/atf3_entry.{png,pdf}")

r2 <- rec[rec$state == "Crush 168h", ]
lim <- range(c(r2$wt_rec, r2$ko_rec), na.rm = TRUE)
gg2 <- ggplot(r2, aes(wt_rec, ko_rec)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey45") +
  geom_hline(yintercept = 0, linewidth = .2, colour = "grey75") +
  geom_vline(xintercept = 0, linewidth = .2, colour = "grey75") +
  geom_point(aes(colour = blocked > 1), size = 2.6, alpha = .9) +
  ggrepel::geom_text_repel(data = subset(r2, wt_rec > 1 | blocked > 1),
                           aes(label = gene), size = 2.9, max.overlaps = 24, seed = 2, colour = "grey15") +
  scale_colour_manual(values = c(`FALSE` = "#4E79A7", `TRUE` = "#B2182B"),
                      labels = c("recruited normally", "recruitment blocked in KO"), name = NULL) +
  coord_equal(xlim = lim, ylim = lim) +
  labs(title = "Recruitment of panel genes at 7 days: Atf3-WT vs Atf3-KO",
       subtitle = paste(strwrap(paste(
         "log2 change in the FRACTION of neurons expressing each gene, crush versus naive, within",
         "each genotype. Points below the dashed line recruit fewer cells without Atf3."), width = 100), collapse = "\n"),
       x = "Recruitment in Atf3-WT (log2 fraction expressing)",
       y = "Recruitment in Atf3-KO") +
  theme_pub(base_size = 12) +
  theme(plot.subtitle = element_text(colour = "grey30", size = 9), legend.position = "bottom")
ggsave(file.path(FIGS, "atf3_recruitment.png"), gg2, width = 9, height = 9.6, dpi = 300, bg = "white")
ggsave(file.path(FIGS, "atf3_recruitment.pdf"), gg2, width = 9, height = 9.6, device = cairo_pdf)
msg("wrote figs/atf3_recruitment.{png,pdf} and data/atf3_percell.{rds,_recruitment.csv}")
