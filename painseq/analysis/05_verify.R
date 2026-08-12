# =====================================================================
# 05_verify.R — the checks that decide whether the Stage 1 result is believable.
#
#  3. POSITIVE control : canonical RAGs must be induced by axotomy at 24h-7d.
#  4. NEGATIVE control : the panel must be near-flat in CFA and Paclitaxel
#                        (non-axotomy insults). This is the single-cell analogue
#                        of the bulk "MINP-negative" property.
#  5. SENSITIVITY      : all-cells vs male-C57-only (computed in 03).
#  6. BULK cross-check : single-cell SpNT-7d log2FC vs the app's mouse
#                        SNI_vs_Sham bulk log2FC, same genes.
#  7. AMBIENT check    : is the panel neuron-specific, or does it appear "induced"
#                        in every cell type (the signature of ambient RNA leakage
#                        in snRNA-seq)?
#
# Writes data/verification_report.txt
# =====================================================================
suppressMessages({ library(edgeR); library(limma) })
source("/hpf/projects/msalter/sghazis/painseq/analysis/00_common.R")

tc    <- readRDS(file.path(DATA, "panel_timecourse.rds"))
fits  <- readRDS(file.path(DATA, "panel_celltype_fits.rds"))
pg    <- readRDS(file.path(DATA, "panel_genes.rds"))
f     <- tc$all
cm    <- f$cmeta
out   <- character(0)
say   <- function(...) { s <- paste0(...); message(s); out <<- c(out, s) }

say("=====================================================================")
say(" VERIFICATION — 47-gene conserved injury panel in GSE154659")
say("=====================================================================")

panel_syms <- pg$matrix_sym[pg$role == "panel" & pg$present]
panel_syms <- intersect(panel_syms, rownames(f$l2fc))
say("\npanel genes carried into the model: ", length(panel_syms), " of 47",
    "  (Crybg1 absent from matrix; others dropped by filterByExpr)")

# ---------------------------------------------------------------------
# 3. Positive control
# ---------------------------------------------------------------------
say("\n--- 3. POSITIVE CONTROL: canonical RAGs, axotomy 24h-7d ---")
rags <- c("Atf3", "Gal", "Npy", "Sox11", "Gap43", "Sprr1a")
pc   <- cm$col[cm$inj %in% AXOTOMY & cm$time %in% c(24L, 72L, 168L)]
for (g in rags) {
  if (!g %in% rownames(f$l2fc)) { say(sprintf("  %-8s NOT IN MODEL", g)); next }
  v <- f$l2fc[g, pc]; p <- f$padj[g, pc]
  say(sprintf("  %-8s mean log2FC %+5.2f  (range %+5.2f..%+5.2f)  sig in %d/%d columns",
              g, mean(v), min(v), max(v), sum(p < 0.05, na.rm = TRUE), length(pc)))
}
rag_ok <- all(vapply(intersect(rags, rownames(f$l2fc)),
                     function(g) mean(f$l2fc[g, pc]) > 1, logical(1)))
say("  => PASS: every canonical RAG has mean log2FC > 1 across axotomy 24h-7d: ", rag_ok)

# ---------------------------------------------------------------------
# 4. Negative control — axotomy vs non-axotomy
# ---------------------------------------------------------------------
say("\n--- 4. NEGATIVE CONTROL: panel response in axotomy vs non-axotomy ---")
ax_cols  <- cm$col[cm$inj %in% AXOTOMY & cm$time %in% c(24L, 72L, 168L)]
non_cols <- cm$col[cm$inj %in% NON_AXOTOMY]
summ <- function(cols, lab) {
  V <- f$l2fc[panel_syms, cols, drop = FALSE]
  P <- f$padj[panel_syms, cols, drop = FALSE]
  say(sprintf("  %-26s n_cols=%2d  mean log2FC %+.2f  mean |log2FC| %.2f  sig-up (p<.05 & FC>1): %4.1f%% of cells",
              lab, length(cols), mean(V, na.rm = TRUE), mean(abs(V), na.rm = TRUE),
              100 * mean(P < 0.05 & V > 1, na.rm = TRUE)))
  invisible(list(V = V, P = P))
}
A <- summ(ax_cols,  "axotomy (24h/72h/7d)")
N <- summ(non_cols, "non-axotomy (CFA, Pacli)")
say(sprintf("  => axotomy induces the panel %.1fx more strongly than non-axotomy (mean log2FC ratio)",
            mean(A$V, na.rm = TRUE) / max(mean(N$V, na.rm = TRUE), 1e-6)))
say("  per non-axotomy column, genes significantly UP (adj.p<0.05 & log2FC>1):")
for (cc in non_cols) {
  v <- f$l2fc[panel_syms, cc]; p <- f$padj[panel_syms, cc]
  up <- panel_syms[which(p < 0.05 & v > 1)]
  say(sprintf("    %-14s %2d / %d  %s", cc, length(up), length(panel_syms),
              if (length(up)) paste(up, collapse = ", ") else "(none)"))
}

# ---------------------------------------------------------------------
# 6. Cross-check against the bulk mouse SNI result
# ---------------------------------------------------------------------
say("\n--- 6. BULK CROSS-CHECK: single-cell SpNT-7d vs bulk mouse SNI_vs_Sham ---")
de <- readRDS(file.path(APP, "app/data/de_long.rds"))
de <- de[de$engine == "DESeq2" & de$level == "gene" & de$contrast == "SNI_vs_Sham", ]
de$symbol <- as.character(de$symbol)
sp_col <- cm$col[cm$inj == "SpNT" & cm$time == 168L]
if (length(sp_col) == 1) {
  sc_v <- f$l2fc[panel_syms, sp_col]
  bulk <- de$log2FC[match(panel_syms, de$symbol)]
  ok   <- is.finite(sc_v) & is.finite(bulk)
  r    <- cor(sc_v[ok], bulk[ok])
  rs   <- cor(sc_v[ok], bulk[ok], method = "spearman")
  say(sprintf("  n=%d panel genes matched | Pearson r = %.3f | Spearman rho = %.3f", sum(ok), r, rs))
  say(sprintf("  same-direction: %d/%d (%.0f%%)", sum(sign(sc_v[ok]) == sign(bulk[ok])), sum(ok),
              100 * mean(sign(sc_v[ok]) == sign(bulk[ok]))))
  say("  => a positive r confirms the gene mapping is sound; these are different")
  say("     studies and platforms, so a moderate r is the expected result.")
}

# ---------------------------------------------------------------------
# 7. Ambient-RNA / cell-type-specificity check
# ---------------------------------------------------------------------
say("\n--- 7. AMBIENT-RNA CHECK: is the panel neuron-specific? ---")
say("  In snRNA-seq, transcripts induced to very high levels in one cell type leak")
say("  into the ambient pool and inflate apparent induction everywhere. If the panel")
say("  is genuinely neuronal, the neuronal log2FC should clearly exceed non-neuronal.")
peak_cols_of <- function(fc) fc$cmeta$col[fc$cmeta$inj %in% AXOTOMY & fc$cmeta$time %in% c(72L, 168L)]
ct_l2fc <- sapply(names(fits), function(ct) {
  fc <- fits[[ct]]; pcx <- peak_cols_of(fc)
  s  <- intersect(panel_syms, rownames(fc$l2fc))
  v  <- rep(NA_real_, length(panel_syms)); names(v) <- panel_syms
  if (length(pcx) && length(s)) v[s] <- rowMeans(fc$l2fc[s, pcx, drop = FALSE], na.rm = TRUE)
  v
})
neur <- intersect(NEURONAL_TYPES, colnames(ct_l2fc))
nonn <- setdiff(colnames(ct_l2fc), neur)
mn <- rowMeans(ct_l2fc[, neur, drop = FALSE], na.rm = TRUE)
mo <- rowMeans(ct_l2fc[, nonn, drop = FALSE], na.rm = TRUE)
say(sprintf("  mean panel log2FC: neuronal %.2f (%d types) vs non-neuronal %.2f (%d types)",
            mean(mn, na.rm = TRUE), length(neur), mean(mo, na.rm = TRUE), length(nonn)))
say(sprintf("  genes with neuronal > non-neuronal: %d/%d",
            sum(mn > mo, na.rm = TRUE), sum(is.finite(mn) & is.finite(mo))))
say("  marquee genes (mean log2FC at axotomy peak, unclipped):")
for (g in intersect(c("Npy","Atf3","Gal","Ecel1","Sprr1a","Sox11","Gap43"), rownames(ct_l2fc))) {
  say(sprintf("    %-8s neuronal %+5.2f | non-neuronal %+5.2f | immune (Macrophage/B/Neutrophil) %+5.2f",
      g, mean(ct_l2fc[g, neur], na.rm = TRUE), mean(ct_l2fc[g, nonn], na.rm = TRUE),
      mean(ct_l2fc[g, intersect(c("Macrophage","B cell","Neutrophil"), colnames(ct_l2fc))], na.rm = TRUE)))
}
say("  CAVEAT: a strongly positive non-neuronal value for a neuron-restricted gene")
say("  (e.g. Npy, Gal) indicates ambient leakage, NOT expression in that cell type.")
say("  Read Fig C column-wise for relative ranking, not as absolute per-type induction.")

# ---------------------------------------------------------------------
# Injury-induced cell types
# ---------------------------------------------------------------------
say("\n--- NOTE: cell types with no Naive baseline ---")
say("  ", paste(setdiff(names(CELLTYPE_FULL), names(fits)), collapse = ", "))
say("  These are injury-induced populations: they do not exist in naive tissue, so")
say("  no 'vs Naive' contrast is definable and they are absent from Fig C/D/E.")
say("  This is exactly why a naive-only deconvolution reference cannot see them.")

say("\n--- 5. SENSITIVITY (from 03) ---")
say(sprintf("  all-cells vs male-C57-only log2FC, Pearson r = %.4f", tc$sensitivity_r))

writeLines(out, file.path(DATA, "verification_report.txt"))
message("\nwrote data/verification_report.txt")
