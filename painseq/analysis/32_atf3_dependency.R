# =====================================================================
# 32_atf3_dependency.R — Stage 3d: is the conserved injury panel Atf3-dependent?
#
# Atf3 is the master regeneration-associated transcription factor and is the
# second-strongest gene in the 47 by mean SNI log2FC. The Atf3-WT / Atf3-KO arm of
# GSE154659 (crush at 36h and 168h, neurons only, male) is a direct loss-of-function
# test: if the panel is downstream of Atf3, its induction should be blunted in the KO.
#
# Design per genotype: Crush(36h, 168h) vs Naive(0h), fitted on neuron pseudobulk.
# The Atf3-dependence of each gene is the interaction — how much smaller the KO
# induction is than the WT induction.
#
# Out: data/atf3_dependency.rds, data/atf3_dependency.csv,
#      figs/atf3_dependency.{png,pdf}
# =====================================================================
suppressMessages({ library(edgeR); library(limma); library(ggplot2) })
source("/hpf/projects/msalter/sghazis/painseq/analysis/00_common.R")

pb <- readRDS(file.path(DATA, "pb_atf3_sample.rds"))
pg <- readRDS(file.path(DATA, "panel_genes.rds"))
rows <- pg[order(pg$role != "panel", -pg$mean_sni), ]

md <- pb$meta
md$geno  <- factor(md$geno, levels = c("Atf3-WT", "Atf3-KO"))
md$state <- factor(ifelse(md$inj == "Naive", "Naive", "Crush"), levels = c("Naive", "Crush"))
md$grp   <- factor(paste(md$geno, md$state, sep = "."),
                   levels = c("Atf3-WT.Naive", "Atf3-WT.Crush", "Atf3-KO.Naive", "Atf3-KO.Crush"))
msg("samples per group:"); print(table(md$grp))

y <- DGEList(counts = pb$counts, group = md$grp)
keep <- filterByExpr(y, group = md$grp)
msg("filterByExpr keeps ", sum(keep), "/", length(keep), " genes")
y <- calcNormFactors(y[keep, , keep.lib.sizes = FALSE], method = "TMM")
lcpm <- cpm(y, log = TRUE, prior.count = 1)

design <- model.matrix(~0 + grp, data = md)
colnames(design) <- make.names(sub("^grp", "", colnames(design)))
fit <- lmFit(lcpm, design)
cm <- makeContrasts(
  WT_crush  = Atf3.WT.Crush - Atf3.WT.Naive,
  KO_crush  = Atf3.KO.Crush - Atf3.KO.Naive,
  # interaction: how much of the WT induction is lost in the KO
  Atf3_dep  = (Atf3.WT.Crush - Atf3.WT.Naive) - (Atf3.KO.Crush - Atf3.KO.Naive),
  levels = design)
fit2 <- eBayes(contrasts.fit(fit, cm))

grab <- function(cn) {
  t <- topTable(fit2, coef = cn, number = Inf, sort.by = "none")
  data.frame(gene = rownames(lcpm), l2fc = t$logFC, padj = t$adj.P.Val, stringsAsFactors = FALSE)
}
W <- grab("WT_crush"); K <- grab("KO_crush"); D <- grab("Atf3_dep")

# ---------------------------------------------------------------------
# SANITY CHECK — and a caveat that changes how everything below is read.
# ---------------------------------------------------------------------
# If the KO were a transcript-null, Atf3 mRNA would be absent or at least reduced.
# It is not: Atf3 is induced by crush MORE strongly in the KO than in the WT. The
# Renthal KO is an allele whose transcript is still produced and still counted by
# 3'-end snRNA-seq (and Atf3 negatively autoregulates, so removing functional
# protein can raise its own mRNA). Consequence: Atf3 mRNA is NOT a validity readout
# for this knockout, and we cannot confirm from these data alone that Atf3 protein
# function is lost. Everything below assumes the published genotype labels are right.
if ("Atf3" %in% W$gene) {
  wa <- W$l2fc[W$gene == "Atf3"]; ka <- K$l2fc[K$gene == "Atf3"]
  msg(sprintf("SANITY Atf3 mRNA: WT crush %+.2f (adj.p=%.2g) | KO crush %+.2f (adj.p=%.2g)",
              wa, W$padj[W$gene == "Atf3"], ka, K$padj[K$gene == "Atf3"]))
  if (ka >= wa)
    msg("  ** Atf3 mRNA is NOT reduced in the KO — it is induced at least as strongly as WT.\n",
        "     The knockout therefore cannot be verified from transcript level; the allele\n",
        "     still produces a counted transcript. Read all 'Atf3-dependence' below as\n",
        "     conditional on the published genotype labels being correct.")
}

i <- match(rows$matrix_sym, W$gene)
out <- data.frame(
  Gene = rows$gene, role = rows$role,
  WT_crush_l2fc = round(W$l2fc[i], 3), WT_crush_padj = signif(W$padj[i], 3),
  KO_crush_l2fc = round(K$l2fc[i], 3), KO_crush_padj = signif(K$padj[i], 3),
  Atf3_dependence = round(D$l2fc[i], 3), dependence_padj = signif(D$padj[i], 3),
  pct_cells_detected = rows$pct_cells, stringsAsFactors = FALSE)
# fraction of the WT induction retained in the KO (only meaningful if WT induces it)
out$pct_retained_in_KO <- ifelse(!is.na(out$WT_crush_l2fc) & out$WT_crush_l2fc > 1,
                                 round(100 * out$KO_crush_l2fc / out$WT_crush_l2fc, 1), NA)
out <- out[order(-out$Atf3_dependence), ]
write.csv(out, file.path(DATA, "atf3_dependency.csv"), row.names = FALSE)
saveRDS(list(table = out, W = W, K = K, D = D, design = md), file.path(DATA, "atf3_dependency.rds"))

ind <- !is.na(out$WT_crush_l2fc) & out$WT_crush_l2fc > 1 & out$WT_crush_padj < 0.05
msg("\npanel genes induced by crush in Atf3-WT (log2FC>1, adj.p<0.05): ", sum(ind), "/", nrow(out))
dep <- ind & !is.na(out$dependence_padj) & out$dependence_padj < 0.05 & out$Atf3_dependence > 0
msg("of those, significantly Atf3-DEPENDENT (WT induction > KO induction): ", sum(dep))
msg("median retention of the WT induction in the KO: ",
    round(median(out$pct_retained_in_KO[ind], na.rm = TRUE), 1), "%")
msg("\nPOWER: the interaction contrast has n=3 pseudobulk samples per genotype x state")
msg("  (Naive rep1-3; Crush = 36h rep1 + 168h rep1,rep2), so it is badly underpowered.")
msg("  Point estimates are informative; the near-absence of significant interaction")
msg("  p-values reflects sample size, NOT evidence that the panel is Atf3-independent.")
msg("  Treat the retention percentages as the readable signal and the p-values as weak.")
msg("\nmost Atf3-dependent panel genes:")
print(utils::head(out[ind, c("Gene", "WT_crush_l2fc", "KO_crush_l2fc", "Atf3_dependence",
                             "dependence_padj", "pct_retained_in_KO")], 15), row.names = FALSE)
msg("\nleast Atf3-dependent (induced normally without Atf3):")
print(utils::tail(out[ind, c("Gene", "WT_crush_l2fc", "KO_crush_l2fc", "Atf3_dependence",
                             "dependence_padj", "pct_retained_in_KO")], 8), row.names = FALSE)

# ---- figure: WT vs KO induction, one point per panel gene ----------------
d <- out[!is.na(out$WT_crush_l2fc), ]
d$sig_dep <- !is.na(d$dependence_padj) & d$dependence_padj < 0.05 & d$Atf3_dependence > 0
lim <- range(c(d$WT_crush_l2fc, d$KO_crush_l2fc), na.rm = TRUE)
gg <- ggplot(d, aes(WT_crush_l2fc, KO_crush_l2fc)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey45") +
  geom_hline(yintercept = 0, linewidth = .2, colour = "grey75") +
  geom_vline(xintercept = 0, linewidth = .2, colour = "grey75") +
  geom_point(aes(colour = sig_dep, shape = role), size = 2.6, alpha = .9) +
  ggrepel::geom_text_repel(
    data = subset(d, WT_crush_l2fc > 1.5 | sig_dep),
    aes(label = Gene), size = 3, max.overlaps = 22, seed = 1, colour = "grey15") +
  scale_colour_manual(values = c(`FALSE` = "#4E79A7", `TRUE` = "#B2182B"),
                      labels = c("not significant", "Atf3-dependent (adj.p<0.05)"), name = NULL) +
  scale_shape_manual(values = c(panel = 16, reference = 17), name = NULL) +
  coord_equal(xlim = lim, ylim = lim) +
  labs(title = "Is the conserved injury panel Atf3-dependent?",
       subtitle = paste(strwrap(paste(
         "Crush vs Naive induction in Atf3-WT (x) against Atf3-KO (y), mouse DRG neurons,",
         "GSE154659. Points on the dashed line induce normally without Atf3; points below it",
         "lose induction in the knockout. Triangles are the two RAG reference genes."), width = 105), collapse = "\n"),
       x = "log2FC crush vs naive — Atf3-WT", y = "log2FC crush vs naive — Atf3-KO") +
  theme_pub(base_size = 12) +
  theme(plot.subtitle = element_text(colour = "grey30", size = 9), legend.position = "bottom")
ggsave(file.path(FIGS, "atf3_dependency.png"), gg, width = 9, height = 9.8, dpi = 300, bg = "white")
ggsave(file.path(FIGS, "atf3_dependency.pdf"), gg, width = 9, height = 9.8, device = cairo_pdf)
msg("wrote figs/atf3_dependency.{png,pdf} and data/atf3_dependency.{rds,csv}")
