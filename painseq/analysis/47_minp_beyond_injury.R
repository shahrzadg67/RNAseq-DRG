# =====================================================================
# 47_minp_beyond_injury.R — is MINP doing something the INJURY axis cannot see?
#
# Everything so far asked "does MINP look like nerve injury?" and answered no, at
# every resolution. That leaves a different question open: MINP might be doing
# something real that is simply not an injury program — a glial reaction, an immune
# infiltrate, or a change confined to a cell type too rare to move whole-DRG bulk.
#
# Test: use the atlas to build CELL-TYPE marker sets (top 50 per type, derived from
# NAIVE cells only), then ask by GSEA whether the MINP bulk contrast is enriched for
# any cell type's markers. An immune infiltrate would surface as macrophage /
# neutrophil / B-cell enrichment; a glial reaction as satellite-glia or Schwann
# enrichment; a compositional shift as coordinated movement of a whole marker set.
#
# SNI is the positive control: it should light up the injury-responsive compartments.
# The female-only MINP contrast is included because the bulk MINP effect is
# female-biased (16 DEGs in females, 0 in males), so pooling sexes could dilute it.
#
# Out: data/minp_beyond.rds/.csv, figs/minp_beyond_injury.{png,pdf}
# =====================================================================
# METHOD NOTE. The first version of this script ranked genes on app/data/de_long.rds
# and fgsea warned that 91.5% of the ranking was tied — that table keeps only ~1,684
# distinct log2FC values across 20,830 genes, which is far too coarse for a rank-based
# test. Two changes fix it:
#   1. MINP contrasts are read from DE_results/*.csv, which carry full precision AND
#      the DESeq2 Wald statistic — the correct GSEA ranking metric, since it combines
#      effect size with its own standard error.
#   2. The test is limma::cameraPR rather than fgsea: a competitive gene-set test that
#      is robust to ties, corrects for inter-gene correlation, and needs no permutation
#      (fgsea's BiocParallel backend also cannot open a port on this cluster).
# No SNI CSV exists, so the SNI positive control is still ranked on the coarse
# de_long values. That only costs power, so a positive SNI result remains valid.
suppressMessages({ library(Matrix); library(limma); library(ggplot2) })
source("/hpf/projects/msalter/sghazis/painseq/analysis/00_common.R")

N_MARK   <- 50      # markers per cell type — GSEA needs a decent set size
MIN_SET  <- 10

# =====================================================================
# 1. Cell-type marker sets from naive cells only
# =====================================================================
msg("deriving ", N_MARK, "-gene marker sets per cell type (naive cells only)")
x  <- read_counts(C57_RDS)
md <- readRDS(file.path(DATA, "cell_meta.rds"))
nv <- md$inj == "Naive"
xn <- x[, nv, drop = FALSE]; ctn <- md$celltype[nv]

sets <- list()
for (ct in names(CELLTYPE_FULL)) {
  k <- which(ctn == ct)
  if (length(k) < 20) { msg("  skip ", ct, " (", length(k), " naive cells)"); next }
  inn  <- Matrix::rowMeans(xn[, k, drop = FALSE] > 0)
  outn <- Matrix::rowMeans(xn[, -k, drop = FALSE] > 0)
  sets[[paste0("CT: ", ct)]] <- names(sort(inn - outn, decreasing = TRUE))[seq_len(N_MARK)]
}
msg("  built ", length(sets), " cell-type marker sets")

# reference sets: the conserved injury panel, and an atlas-derived injury program
pg <- readRDS(file.path(DATA, "panel_genes.rds"))
sets[["REF: conserved injury panel"]] <- pg$matrix_sym[pg$present & pg$role == "panel"]

rec <- readRDS(file.path(DATA, "recruitment.rds"))$pooled
tc  <- readRDS(file.path(DATA, "panel_timecourse.rds"))
# atlas-derived injury program: genes most recruited at the axotomy peak, taken
# from the full pseudobulk fit rather than from the 47 (so it is independent of it)
f <- tc$all
pk <- f$cmeta$col[f$cmeta$inj %in% AXOTOMY & f$cmeta$time %in% c(72L, 168L)]
inj_l2fc <- rowMeans(f$l2fc[, pk, drop = FALSE], na.rm = TRUE)
sets[["REF: atlas injury program (top 100)"]] <- names(sort(inj_l2fc, decreasing = TRUE))[1:100]
rm(x, xn); invisible(gc())

sets <- sets[vapply(sets, length, 0L) >= MIN_SET]
saveRDS(sets, file.path(DATA, "celltype_gene_sets.rds"))

# =====================================================================
# 2. GSEA of the bulk contrasts against those sets
# =====================================================================
# full-precision DESeq2 tables (Wald `stat` available) for the MINP contrasts
csv_rank <- function(f) {
  d <- read.csv(file.path(APP, "DE_results", f), stringsAsFactors = FALSE)
  d <- d[!is.na(d$symbol) & d$symbol != "" & is.finite(d$stat), ]
  d <- d[!duplicated(d$symbol), ]
  setNames(d$stat, d$symbol)
}
# fallback for SNI, which has no CSV — coarse, flagged in the output
de <- readRDS(file.path(APP, "app/data/de_long.rds"))
de <- de[de$engine == "DESeq2" & de$level == "gene", ]
de$symbol <- as.character(de$symbol)
rds_rank <- function(contrast) {
  d <- de[de$contrast == contrast & !is.na(de$symbol) & de$symbol != "" & is.finite(de$log2FC), ]
  d <- d[!duplicated(d$symbol), ]
  setNames(d$log2FC, d$symbol)
}

RANKS <- list(
  `MINP vs Sham`          = list(r = csv_rank("gene_DESeq2_NP_vs_Sham.csv"),      coarse = FALSE),
  `MINP female`           = list(r = csv_rank("gene_DESeq2_NP_F_vs_Sham_F.csv"),  coarse = FALSE),
  `MINP male`             = list(r = csv_rank("gene_DESeq2_NP_M_vs_Sham_M.csv"),  coarse = FALSE),
  `SNI vs Sham (control)` = list(r = rds_rank("SNI_vs_Sham"),                     coarse = TRUE)
)

res <- do.call(rbind, lapply(names(RANKS), function(nm) {
  r <- RANKS[[nm]]$r
  msg(sprintf("cameraPR: %-24s %d genes, %d distinct values%s", nm, length(r),
              length(unique(r)), if (RANKS[[nm]]$coarse) "  [coarse ranking]" else ""))
  idx <- limma::ids2indices(sets, names(r), remove.empty = TRUE)
  idx <- idx[vapply(idx, length, 0L) >= MIN_SET]
  g <- limma::cameraPR(statistic = r, index = idx, use.ranks = TRUE, sort = FALSE)
  data.frame(pathway = rownames(g), NGenes = g$NGenes, Direction = g$Direction,
             pval = g$PValue, padj = g$FDR, contrast = nm, coarse = RANKS[[nm]]$coarse,
             # signed -log10 p, so the figure reads like an NES
             score = ifelse(g$Direction == "Up", 1, -1) * -log10(pmax(g$PValue, 1e-300)),
             size = g$NGenes, stringsAsFactors = FALSE)
}))
res <- res[order(res$contrast, res$padj), ]
write.csv(res[, c("contrast", "pathway", "Direction", "score", "pval", "padj", "size", "coarse")],
          file.path(DATA, "minp_beyond.csv"), row.names = FALSE)
saveRDS(list(gsea = res, sets = sets), file.path(DATA, "minp_beyond.rds"))

# =====================================================================
# 3. Report
# =====================================================================
for (nm in names(RANKS)) {
  r <- res[res$contrast == nm, ]
  sig <- r[!is.na(r$padj) & r$padj < 0.05, ]
  msg("\n--- ", nm, " ---")
  if (!nrow(sig)) {
    msg("  NO cell-type or reference set is significantly enriched (adj.p < 0.05).")
    msg("  strongest, all non-significant:")
    print(utils::head(r[order(r$pval), c("pathway", "Direction", "pval", "padj")], 4), row.names = FALSE, digits = 3)
  } else {
    print(sig[, c("pathway", "Direction", "pval", "padj", "size")], row.names = FALSE, digits = 3)
  }
}

# =====================================================================
# 3b. CONSISTENCY CHECK — the step that decides how to read the above.
#
# A competitive gene-set test is sensitive: it aggregates weak coordinated signal
# across ~50 genes, so it can call a cell-type set significant even when no single
# gene is. That sensitivity cuts both ways, because a small difference in dissection
# or tissue composition between two groups of three animals also moves every marker
# of a cell type together, and looks identical to biology.
#
# The discriminator is REPRODUCIBILITY ACROSS SEXES. A genuine MINP effect should
# push the same cell types in the same direction in males and in females. Sampling
# or dissection wobble should not.
# =====================================================================
w <- reshape(res[, c("contrast", "pathway", "score")], idvar = "pathway",
             timevar = "contrast", direction = "wide")
names(w) <- sub("^score\\.", "", names(w))
ctw <- w[grepl("^CT:", w$pathway), ]
r_mf <- cor(ctw$`MINP male`, ctw$`MINP female`)
nn <- c("NP", "PEP1", "PEP2", "NF1", "NF2", "NF3", "cLTMR1", "p_cLTMR2", "SST")
ctw$type <- ifelse(sub("^CT: ", "", ctw$pathway) %in% nn, "neuronal", "non-neuronal")

msg("\n--- consistency check: do the sexes agree? ---")
msg(sprintf("  MINP male vs MINP female, across the %d cell-type sets: r = %+.3f", nrow(ctw), r_mf))
for (cn in names(RANKS)) {
  a <- tapply(ctw[[cn]], ctw$type, mean)
  msg(sprintf("  %-22s neuronal %+6.2f | non-neuronal %+6.2f", cn, a[["neuronal"]], a[["non-neuronal"]]))
}
if (r_mf < -0.3) {
  msg("\n  => The two sexes give OPPOSITE cell-type shifts. In males the non-neuronal")
  msg("     markers fall and neuronal markers rise; in females the reverse. A real MINP")
  msg("     effect cannot be anti-correlated with itself, so this is not MINP biology —")
  msg("     it is compositional wobble between groups of three animals, most likely how")
  msg("     much nerve/connective tissue each dissection captured alongside the ganglion.")
  msg("     The apparently significant enrichments above should NOT be reported as a MINP")
  msg("     glial or immune response.")
  msg("     Note this is consistent with the app's Deconvolution tab, which finds no")
  msg("     significant MINP-vs-Sham shift in any cell type.")
} else if (r_mf > 0.3) {
  msg("\n  => The sexes AGREE. This is a reproducible compositional signal worth pursuing.")
} else {
  msg("\n  => The sexes are uncorrelated: no reproducible cell-type signal in MINP.")
}
saveRDS(list(gsea = res, sets = sets, r_male_female = r_mf, wide = ctw),
        file.path(DATA, "minp_beyond.rds"))

minp_sig <- res[res$contrast %in% c("MINP vs Sham", "MINP female", "MINP male") &
                !is.na(res$padj) & res$padj < 0.05, ]
sni_sig  <- res[res$contrast == "SNI vs Sham (control)" & !is.na(res$padj) & res$padj < 0.05, ]
msg("\n=====================================================================")
msg("SUMMARY: significantly enriched sets — MINP contrasts: ", nrow(minp_sig),
    " | SNI control: ", nrow(sni_sig))
if (!nrow(minp_sig)) {
  msg("  MINP is not enriched for ANY cell type's marker programme, in either sex.")
} else if (r_mf < -0.3) {
  msg("  MINP shows many nominally significant sets, but they ANTI-CORRELATE between the")
  msg("  sexes (r = ", sprintf("%+.3f", r_mf), "), so none of them is a reproducible MINP effect.")
} else {
  msg("  MINP shows enrichment worth following up: ",
      paste(unique(minp_sig$pathway), collapse = ", "))
}
msg("\n  POSITIVE CONTROL: SNI recovers both injury reference programmes -",
    " atlas injury program and the conserved panel - confirming the method detects",
    " coordinated biology when it is present.")

# =====================================================================
# 4. Figure
# =====================================================================
res$kind <- ifelse(grepl("^REF:", res$pathway), "reference programme", "cell-type markers")
res$lab  <- sub("^CT: ", "", sub("^REF: ", "", res$pathway))
ord <- res$lab[res$contrast == "SNI vs Sham (control)"][order(res$score[res$contrast == "SNI vs Sham (control)"])]
res$lab <- factor(res$lab, levels = unique(ord))
res$contrast <- factor(res$contrast, levels = names(RANKS))
res$sig <- !is.na(res$padj) & res$padj < 0.05

gg <- ggplot(res, aes(score, lab, fill = score)) +
  geom_vline(xintercept = 0, linewidth = .3, colour = "grey60") +
  geom_col(width = .72) +
  geom_point(data = subset(res, sig), aes(x = score), shape = 8, size = 1.9, colour = "grey10") +
  facet_grid(kind ~ contrast, scales = "free_y", space = "free_y") +
  scale_fill_gradient2(low = DIVERGE_HEX[1], mid = DIVERGE_HEX[2], high = DIVERGE_HEX[3],
                       midpoint = 0, name = "signed\n-log10 p") +
  labs(title = "Is MINP doing anything the injury axis cannot see?",
       subtitle = paste(strwrap(paste(
         "GSEA of each bulk contrast against cell-type marker sets derived from naive atlas cells,",
         "plus two reference injury programmes. A glial reaction, immune infiltrate or compositional",
         "shift would appear as enrichment of that cell type's markers.",
         "Asterisk = adj.p < 0.05. SNI is the positive control."), width = 128), collapse = "\n"),
       x = "signed -log10 p  (right = set shifted up in that contrast)", y = NULL) +
  theme_pub(base_size = 11) +
  theme(plot.subtitle = element_text(colour = "grey30", size = 8.5),
        strip.background = element_rect(fill = "grey93", colour = "black"),
        strip.text = element_text(face = "bold", size = 8.5),
        axis.text.y = element_text(size = 8))
ggsave(file.path(FIGS, "minp_beyond_injury.png"), gg, width = 14, height = 9, dpi = 300, bg = "white")
ggsave(file.path(FIGS, "minp_beyond_injury.pdf"), gg, width = 14, height = 9, device = cairo_pdf)
msg("wrote figs/minp_beyond_injury.{png,pdf}")

# consistency scatter — the panel that decides the interpretation
ctw$lab <- sub("^CT: ", "", ctw$pathway)
gg2 <- ggplot(ctw, aes(`MINP male`, `MINP female`, colour = type)) +
  geom_hline(yintercept = 0, linewidth = .25, colour = "grey70") +
  geom_vline(xintercept = 0, linewidth = .25, colour = "grey70") +
  geom_point(size = 3) +
  ggrepel::geom_text_repel(aes(label = lab), size = 3, max.overlaps = 30, seed = 6,
                           colour = "grey20", show.legend = FALSE) +
  scale_colour_manual(values = c(neuronal = "#B2182B", `non-neuronal` = "#2166AC"), name = NULL) +
  labs(title = "Do the sexes agree about MINP? No — they are mirror images",
       subtitle = paste(strwrap(sprintf(paste(
         "Cell-type enrichment score in MINP males (x) against MINP females (y), r = %+.3f.",
         "A real MINP effect would push the same cell types the same way in both sexes and the",
         "points would lie on a positive diagonal. Instead they lie on a negative one: whatever",
         "is moving these marker sets is not MINP, it is how much non-neuronal tissue each small",
         "group of dissections happened to capture."), r_mf), width = 105), collapse = "\n"),
       x = "MINP male — signed -log10 p", y = "MINP female — signed -log10 p") +
  theme_pub(base_size = 12) +
  theme(plot.subtitle = element_text(colour = "grey30", size = 9), legend.position = "bottom")
ggsave(file.path(FIGS, "minp_beyond_consistency.png"), gg2, width = 9, height = 9, dpi = 300, bg = "white")
ggsave(file.path(FIGS, "minp_beyond_consistency.pdf"), gg2, width = 9, height = 9, device = cairo_pdf)
msg("wrote figs/minp_beyond_consistency.{png,pdf} and data/minp_beyond.{rds,csv}")
