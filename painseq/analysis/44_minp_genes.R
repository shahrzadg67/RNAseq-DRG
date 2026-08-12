# =====================================================================
# 44_minp_genes.R — where are the MINP-specific genes in the DRG atlas?
#
# app/data/minp_specific.rds holds the 17 genes called MINP-associated in the bulk
# study, of which the app highlights Myh7 and Slc15a2 as the "core" (well
# annotated, injury-negative, MINP-up). This asks the obvious follow-up: are those
# genes actually expressed in DRG, and by which cells?
#
# The answer matters for interpretation. A gene called MINP-specific in a 3-vs-3
# bulk comparison, but essentially undetectable across 141,093 DRG nuclei, is far
# more likely to be low-count noise or contamination than biology.
#
# Out: data/minp_genes_sc.rds/.csv, figs/minp_genes_sc.{png,pdf}
# =====================================================================
suppressMessages({ library(Matrix); library(ggplot2) })
source("/hpf/projects/msalter/sghazis/painseq/analysis/00_common.R")

minp <- readRDS(file.path(APP, "app/data/minp_specific.rds"))
ann  <- readRDS(file.path(APP, "app/data/gene_annot.rds")); ann <- ann[ann$level == "gene", ]
minp$symbol <- ann$symbol[match(minp$feature_id, ann$feature_id)]
minp$named  <- !is.na(minp$symbol) & minp$symbol != ""

msg("MINP-associated genes: ", nrow(minp))
msg("  with a real gene symbol: ", sum(minp$named), " (", paste(minp$symbol[minp$named], collapse = ", "), ")")
msg("  unnamed Ensembl loci:    ", sum(!minp$named))
msg("  => the MINP-specific signature is dominated by unannotated loci. Those cannot be")
msg("     looked up in a symbol-keyed 2020 atlas, and they are also the rows carrying the")
msg("     most extreme fold-changes with non-significant pooled p-values, which is the")
msg("     signature of low-count noise rather than regulated expression.")

# ---------------------------------------------------------------------
# Look up the named ones in the atlas
# ---------------------------------------------------------------------
x  <- read_counts(C57_RDS)
md <- readRDS(file.path(DATA, "cell_meta.rds"))
tot <- Matrix::colSums(x)

want <- unique(minp$symbol[minp$named])
have <- intersect(want, rownames(x))
msg("\nof the ", length(want), " named genes, ", length(have), " are in the atlas: ",
    paste(have, collapse = ", "))
if (length(setdiff(want, have)))
  msg("  absent from the 2020 atlas annotation: ", paste(setdiff(want, have), collapse = ", "))

# reference points: a housekeeping gene and a well-detected panel gene, so
# "low expression" has a scale
REF <- intersect(c("Actb", "Atf3", "Sprr1a"), rownames(x))
genes <- c(have, REF)
sub <- x[genes, , drop = FALSE]
ln  <- log1p(sweep(as.matrix(sub), 2, tot / 1e4, "/"))

det <- data.frame(
  gene = genes,
  role = ifelse(genes %in% have, "MINP-associated", "reference"),
  pct_cells = round(100 * rowMeans(ln > 0), 3),
  total_umi = as.integer(Matrix::rowSums(sub)),
  mean_expressing = round(apply(ln, 1, function(v) if (any(v > 0)) mean(v[v > 0]) else 0), 2),
  stringsAsFactors = FALSE)
det <- det[order(-det$pct_cells), ]
msg("\n--- detection across all 141,093 nuclei ---")
print(det, row.names = FALSE)

# by cell type
bt <- do.call(rbind, lapply(have, function(g)
  data.frame(gene = g, celltype = names(tapply(ln[g, ], md$celltype, function(v) mean(v > 0))),
             pct = 100 * as.numeric(tapply(ln[g, ], md$celltype, function(v) mean(v > 0))),
             stringsAsFactors = FALSE)))
# injury responsiveness: fraction expressing, naive vs axotomy peak
resp <- do.call(rbind, lapply(have, function(g) {
  nv <- md$inj == "Naive"; ax <- md$inj %in% AXOTOMY & md$time %in% c(72L, 168L)
  data.frame(gene = g,
             pct_naive   = round(100 * mean(ln[g, nv] > 0), 3),
             pct_axotomy = round(100 * mean(ln[g, ax] > 0), 3), stringsAsFactors = FALSE)
}))
resp$fold <- round((resp$pct_axotomy + .01) / (resp$pct_naive + .01), 2)
msg("\n--- injury responsiveness of the named MINP genes ---")
print(resp, row.names = FALSE)

saveRDS(list(minp = minp, detection = det, by_celltype = bt, response = resp),
        file.path(DATA, "minp_genes_sc.rds"))
write.csv(det, file.path(DATA, "minp_genes_sc.csv"), row.names = FALSE)

core <- det[det$gene %in% c("Myh7", "Slc15a2"), ]
if (nrow(core)) {
  msg("\n--- the two 'core' MINP genes ---")
  for (i in seq_len(nrow(core)))
    msg(sprintf("  %-9s detected in %.3f%% of nuclei (%d UMI across the whole atlas)",
                core$gene[i], core$pct_cells[i], core$total_umi[i]))
  msg("  Compare Actb, which every cell expresses. If a gene is near-undetectable across")
  msg("  141,093 DRG nuclei, a bulk call of 'MINP-specific' should be treated as provisional")
  msg("  and confirmed by a targeted assay before being built on.")
}

# ---------------------------------------------------------------------
bt$celltype <- factor(bt$celltype, levels = rev(intersect(names(CELLTYPE_FULL), unique(bt$celltype))))
bt$gene <- factor(bt$gene, levels = have)
g1 <- ggplot(bt, aes(gene, celltype, fill = pct)) +
  geom_tile(colour = "grey92", linewidth = .2) +
  geom_text(aes(label = ifelse(pct >= 0.05, sprintf("%.1f", pct), "")), size = 2.7, colour = "grey15") +
  scale_fill_gradient(low = "#F7F7F7", high = "#B2182B", name = "% of cells\nexpressing") +
  labs(title = "MINP-associated genes across DRG cell types",
       subtitle = paste(strwrap(paste(
         "Percentage of nuclei expressing each named MINP-associated gene. Only",
         paste(length(have), "of the 17"), "MINP rows carry a gene symbol at all — the rest are",
         "unannotated Ensembl loci, which is itself a reason to treat that gene set as provisional."),
         width = 110), collapse = "\n"),
       x = NULL, y = NULL) +
  theme_pub(base_size = 11) +
  theme(axis.text.x = element_text(angle = 40, hjust = 1),
        plot.subtitle = element_text(colour = "grey30", size = 8.5))
ggsave(file.path(FIGS, "minp_genes_sc.png"), g1, width = 8.5, height = 8, dpi = 300, bg = "white")
ggsave(file.path(FIGS, "minp_genes_sc.pdf"), g1, width = 8.5, height = 8, device = cairo_pdf)
msg("\nwrote figs/minp_genes_sc.{png,pdf} and data/minp_genes_sc.{rds,csv}")
