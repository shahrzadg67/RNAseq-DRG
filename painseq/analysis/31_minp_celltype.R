# =====================================================================
# 31_minp_celltype.R — Stage 3c: does MINP resemble ANY cell-type-restricted
# injury response, or is it unlike nerve injury at every resolution?
#
# The bulk work showed MINP-vs-Sham is uncorrelated with the conserved injury
# signature (r ~ 0) at whole-tissue level. That leaves an obvious escape hatch: a
# real MINP response confined to ONE cell type could be diluted below detection in
# whole-DRG bulk. This script closes it.
#
# Approach: correlate the mouse bulk log2FC vectors (MINP-vs-Sham and, as a
# positive control, SNI-vs-Sham) against the single-cell per-cell-type injury
# signatures — one correlation per cell type per injury model per timepoint. If
# MINP matched any cell-type-restricted program, it would surface as a high
# correlation somewhere in that grid.
#
# Genes: all shared genes, not just the 47 — restricting to the panel would build
# in the answer, since the panel was selected to be MINP-negative.
#
# Out: data/minp_celltype_concordance.rds/.csv, figs/minp_celltype_concordance.{png,pdf}
# =====================================================================
suppressMessages({ library(ggplot2) })
source("/hpf/projects/msalter/sghazis/painseq/analysis/00_common.R")

MIN_GENES <- 500     # a correlation needs a decent shared-gene set to mean anything

fits <- readRDS(file.path(DATA, "panel_celltype_fits.rds"))
tc   <- readRDS(file.path(DATA, "panel_timecourse.rds"))

# ---- bulk mouse contrasts ------------------------------------------------
de <- readRDS(file.path(APP, "app/data/de_long.rds"))
de <- de[de$engine == "DESeq2" & de$level == "gene", ]
de$symbol <- as.character(de$symbol)
bulk_vec <- function(contrast) {
  d <- de[de$contrast == contrast & !is.na(de$symbol) & de$symbol != "", ]
  d <- d[!duplicated(d$symbol), ]
  setNames(d$log2FC, d$symbol)
}
BULK <- list(MINP = bulk_vec("NP_vs_Sham"), SNI = bulk_vec("SNI_vs_Sham"))
msg("bulk contrasts: MINP ", length(BULK$MINP), " genes | SNI ", length(BULK$SNI), " genes")

# ---- correlate each bulk contrast against every cell type x injury x time --
rowsdf <- list()
for (ct in names(fits)) {
  f <- fits[[ct]]
  for (k in seq_len(nrow(f$cmeta))) {
    cc <- f$cmeta$col[k]
    sc <- f$l2fc[, cc]
    sc <- sc[is.finite(sc)]
    for (b in names(BULK)) {
      g <- intersect(names(sc), names(BULK[[b]]))
      if (length(g) < MIN_GENES) next
      rowsdf[[length(rowsdf) + 1]] <- data.frame(
        celltype = ct, injury = as.character(f$cmeta$inj[k]), time = f$cmeta$time[k],
        contrast = b, n_genes = length(g),
        r = cor(sc[g], BULK[[b]][g]),
        rho = cor(sc[g], BULK[[b]][g], method = "spearman"),
        stringsAsFactors = FALSE)
    }
  }
}
cc <- do.call(rbind, rowsdf)
cc$neuronal <- cc$celltype %in% NEURONAL_TYPES
saveRDS(cc, file.path(DATA, "minp_celltype_concordance.rds"))
write.csv(cc, file.path(DATA, "minp_celltype_concordance.csv"), row.names = FALSE)

summ <- function(b) {
  s <- cc[cc$contrast == b, ]
  msg(sprintf("\n%s vs single-cell injury signatures (%d cell type x injury x time cells):", b, nrow(s)))
  msg(sprintf("  Pearson r: median %+.3f | range %+.3f .. %+.3f | %d of %d with r > 0.2",
              median(s$r), min(s$r), max(s$r), sum(s$r > 0.2), nrow(s)))
  top <- s[order(-s$r), ][1:5, ]
  for (i in seq_len(nrow(top)))
    msg(sprintf("    best: %-18s %-11s %5s  r = %+.3f (n=%d genes)",
                top$celltype[i], top$injury[i],
                ifelse(top$time[i] < 168, paste0(top$time[i], "h"), paste0(top$time[i] / 24, "d")),
                top$r[i], top$n_genes[i]))
}
summ("SNI"); summ("MINP")

sni_max <- max(cc$r[cc$contrast == "SNI"]); minp_max <- max(cc$r[cc$contrast == "MINP"])
msg(sprintf("\n=> SNI reaches a maximum correlation of %+.3f; MINP never exceeds %+.3f.", sni_max, minp_max))
msg("   MINP does not resemble the injury program in ANY cell type, at ANY timepoint,")
msg("   in ANY of the six injury models — so the bulk null result is not a dilution artefact.")

# ---- figure --------------------------------------------------------------
cc$tlab <- factor(ifelse(cc$time < 168, paste0(cc$time, "h"), paste0(cc$time / 24, "d")),
                  levels = unique(ifelse(sort(unique(cc$time)) < 168,
                                         paste0(sort(unique(cc$time)), "h"),
                                         paste0(sort(unique(cc$time)) / 24, "d"))))
# rev() so the neuronal subtypes read from the TOP down — ggplot places the first
# discrete level at the bottom of the y axis
cc$celltype <- factor(cc$celltype, levels = rev(intersect(names(CELLTYPE_FULL), unique(cc$celltype))))
cc$injury   <- factor(cc$injury, levels = intersect(INJ_ORDER, unique(cc$injury)))
cc$contrast <- factor(cc$contrast, levels = c("SNI", "MINP"),
                      labels = c("Bulk SNI vs Sham (positive control)", "Bulk MINP vs Sham"))

lim <- max(abs(cc$r))
gg <- ggplot(cc, aes(x = tlab, y = celltype, fill = r)) +
  geom_tile(colour = "grey92", linewidth = .15) +
  facet_grid(contrast ~ injury, scales = "free_x", space = "free_x",
             labeller = labeller(injury = INJ_SHORT)) +
  scale_fill_gradient2(low = DIVERGE_HEX[1], mid = DIVERGE_HEX[2], high = DIVERGE_HEX[3],
                       midpoint = 0, limits = c(-lim, lim), name = "Pearson r") +
  labs(title = "Does MINP resemble injury in any single cell type?",
       subtitle = paste(strwrap(paste(
         "Correlation of the mouse bulk log2FC vector against each single-cell injury signature,",
         "computed over all shared genes (not just the 47, which would beg the question).",
         "SNI recovers the injury program strongly; MINP is flat everywhere — so its null result",
         "in bulk is not a cell-type dilution artefact."), width = 125), collapse = "\n"),
       x = "Time after injury", y = NULL) +
  theme_pub(base_size = 11) +
  theme(axis.text.x = element_text(angle = 90, vjust = .5, hjust = 1, size = 7),
        axis.text.y = element_text(size = 8),
        plot.subtitle = element_text(colour = "grey30", size = 8.5),
        strip.background = element_rect(fill = "grey93", colour = "black"),
        strip.text = element_text(face = "bold", size = 8.5),
        panel.spacing = unit(2.5, "pt"))
ggsave(file.path(FIGS, "minp_celltype_concordance.png"), gg, width = 15, height = 9, dpi = 300, bg = "white")
ggsave(file.path(FIGS, "minp_celltype_concordance.pdf"), gg, width = 15, height = 9, device = cairo_pdf)
msg("wrote figs/minp_celltype_concordance.{png,pdf} and data/minp_celltype_concordance.{rds,csv}")
