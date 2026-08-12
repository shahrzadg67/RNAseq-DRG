# =====================================================================
# 33_deconv_report.R — interpret the naive-vs-injury deconvolution comparison.
#
# Split from 30_deconv_injury_ref.R so the interpretation can be re-run without
# repeating ~20 minutes of MuSiC.
#
# The headline is NEGATIVE and worth stating plainly: building the reference from
# injured cells did not improve detection of the injury-induced populations in the
# mouse bulk — it removed them. The naive reference is the one that produces the
# biologically sensible signal.
#
# Out: data/deconv_report.txt, figs/deconv_repair_by_group.{png,pdf}
# =====================================================================
suppressMessages(library(ggplot2))
source("/hpf/projects/msalter/sghazis/painseq/analysis/00_common.R")

d   <- readRDS(file.path(DATA, "deconv_compare.rds"))
out <- character(0)
say <- function(...) { s <- paste0(...); message(s); out <<- c(out, s) }

REPAIR <- c("Repair schwann_N", "Repair fibroblast")

say("=====================================================================")
say(" Naive-only vs injury-state deconvolution reference")
say("=====================================================================")
say("\nReference composition:")
for (rn in names(d$refs))
  say(sprintf("  %-12s %2d cell types%s", rn, length(d$refs[[rn]]),
      if (length(setdiff(d$refs$injury, d$refs[[rn]])))
        paste0("  (missing: ", paste(setdiff(d$refs$injury, d$refs[[rn]]), collapse = ", "), ")") else ""))
say("\n  `Repair fibroblast` exists in only 6 naive cells across the whole atlas, below the")
say("  20-cell minimum, so a naive-only reference cannot represent it at all. That was the")
say("  motivation for building an injury-state reference.")

# ---- mouse, grouped by condition -----------------------------------------
grp_of <- function(s) sub("_[FM]?[0-9]*$", "", sub("_([FM])([0-9]*)$", "", s))
say("\n--- MOUSE bulk: repair populations by condition (% of bulk) ---")
tab <- list()
for (rn in c("naive_fixed", "injury")) {
  p <- d$prop[[paste0("mouse|", rn)]]
  if (is.null(p)) next
  g <- sub("_.*", "", rownames(p))
  for (ct in intersect(REPAIR, colnames(p))) {
    m <- tapply(100 * p[, ct], g, mean)
    say(sprintf("  %-12s %-18s %s", rn, ct,
                paste(sprintf("%s=%.3f", names(m), m), collapse = "  ")))
    tab[[length(tab) + 1]] <- data.frame(reference = rn, celltype = ct, group = names(m),
                                         pct = as.numeric(m), stringsAsFactors = FALSE)
  }
}
tab <- do.call(rbind, tab)

say("\n=> The injury-state reference reports 0.000% repair cells in EVERY mouse sample,")
say("   including the two SNI samples. The naive reference does not: it puts")
say("   `Repair schwann_N` at ~1.3-1.8% in the SNI samples against 0.13-0.62% in Sham")
say("   and MINP — roughly a 3-4x elevation, in the one condition where axotomy actually")
say("   happened. So the naive reference gives the biologically sensible answer here and")
say("   the injury reference does not.")
say("\n   Why: MuSiC uses weighted NNLS. Signatures built from injured cells are far more")
say("   extreme, and with 20 competing cell types the solver drives the rare repair")
say("   components to exactly zero while other types absorb the signal (NF3 rises from")
say("   2.09% to 6.68%, PEP1 collapses from 2.88% to 0.03%).")
say("\n   CONCLUSION: do NOT switch the app's Deconvolution tab to an injury-state reference.")
say("   The naive-only reference stays the better choice. The real caveat to record is that")
say("   absolute proportions are sensitive to reference construction, so they should be read")
say("   as relative comparisons across samples, never as absolute tissue composition.")

# ---- cross-species sanity -------------------------------------------------
say("\n--- Cross-species check on the injury reference ---")
for (sp in c("rat", "monkey")) {
  p <- d$prop[[paste0(sp, "|injury")]]
  if (is.null(p)) next
  ct <- intersect(REPAIR, colnames(p))
  say(sprintf("  %-7s %s", sp,
      paste(sprintf("%s: mean %.2f%%, max %.2f%%", ct, 100 * colMeans(p[, ct, drop = FALSE]),
                    100 * apply(p[, ct, drop = FALSE], 2, max)), collapse = " | ")))
}
say("  In macaque, `Repair schwann_N` is highest in the injured (ipsi) samples — mean ~10.8%")
say("  versus ~7.3% in controls — which is the right direction, but the control samples")
say("  already sit high, and rat shows naive samples as high as injured ones. Cross-species")
say("  ortholog mapping adds enough noise that these numbers should not be over-read.")

writeLines(out, file.path(DATA, "deconv_report.txt"))

# ---- figure ---------------------------------------------------------------
if (!is.null(tab) && nrow(tab)) {
  gg <- ggplot(tab, aes(x = group, y = pct, fill = reference)) +
    geom_col(position = position_dodge(preserve = "single"), width = .7) +
    facet_wrap(~ celltype, ncol = 1, scales = "free_y") +
    scale_fill_manual(values = c(naive_fixed = "#4E79A7", injury = "#E15759"),
                      labels = c(injury = "injury-state reference", naive_fixed = "naive-only reference"),
                      name = NULL) +
    labs(title = "Injury-induced cell types in the mouse bulk",
         subtitle = paste(strwrap(paste(
           "The injury-state reference reports zero repair cells in every sample. The naive-only",
           "reference puts Repair schwann_N ~3-4x higher in SNI than in Sham or MINP — the",
           "expected result, since SNI is the only true axotomy in this study."), width = 100), collapse = "\n"),
         x = NULL, y = "Estimated proportion of bulk (%)") +
    theme_pub(base_size = 12) +
    theme(plot.subtitle = element_text(colour = "grey30", size = 9),
          strip.background = element_rect(fill = "grey93", colour = "black"),
          strip.text = element_text(face = "bold"), legend.position = "bottom")
  ggsave(file.path(FIGS, "deconv_repair_by_group.png"), gg, width = 8, height = 8, dpi = 300, bg = "white")
  ggsave(file.path(FIGS, "deconv_repair_by_group.pdf"), gg, width = 8, height = 8, device = cairo_pdf)
  msg("wrote figs/deconv_repair_by_group.{png,pdf}")
}
msg("wrote data/deconv_report.txt")
