# =====================================================================
# 00_common.R — shared paths, barcode parser, gene panel, house styling.
# Sourced by every other script in analysis/.
#
# Dataset: Renthal et al. 2020, mouse DRG snRNA-seq injury atlas (GEO GSE154659).
#   C57 atlas  : 25,105 genes x 141,093 cells, 66 samples, 6 injury models,
#                12 timepoints, 20 author-assigned cell types.
#   Atf3 WT/KO : 24,597 genes x 17,665 cells, 12 samples, neurons only.
# =====================================================================

PROJ  <- "/hpf/projects/msalter/sghazis/painseq"
APP   <- "/hpf/projects/msalter/sghazis/rnaseq_TUY35595"
DATA  <- file.path(PROJ, "data")
FIGS  <- file.path(PROJ, "figs")

# The C57 counts already exist DECOMPRESSED under the app's refs/ (downloaded for
# the deconvolution). Read that one — it takes a plain readRDS.
# The painseq/*.RDS.gz copies are double-gzipped (an outer gzip around an already
# gzip-compressed RDS), so plain readRDS fails with "unknown input format"; they
# need readRDS(gzcon(gzfile(f, "rb"))). read_counts() handles both.
C57_RDS  <- file.path(APP,  "refs/scref/GSE154659_C57_Raw_counts.RDS")
C57_GZ   <- file.path(PROJ, "GSE154659_C57_Raw_counts.RDS.gz")
ATF3_GZ  <- file.path(PROJ, "GSE154659_Atf3_WT_KO_Raw_counts.RDS.gz")

read_counts <- function(path) {
  stopifnot(file.exists(path))
  if (grepl("\\.gz$", path)) {
    con <- gzcon(gzfile(path, "rb")); on.exit(close(con))
    readRDS(con)
  } else readRDS(path)
}

# =====================================================================
# Barcode parser
# ---------------------------------------------------------------------
# Format: sex_genotype_injury_time_rep_CELLTYPE_barcode
#   e.g.  male_C57_Crush_168_rep1_NF2_bcFJOM
#
# IMPORTANT: the cell type is fields 6 .. (n-1) joined by "_", NOT the
# second-to-last field. Three cell types contain an underscore —
# p_cLTMR2, Schwann_M, Schwann_N, Repair schwann_N — and the
# second-to-last-field shortcut used in rnaseq_TUY35595/analysis/80_deconv.R:60
# silently truncates them to cLTMR2 / M / N / schwann. That mangled nomenclature
# is what ended up in the app's CELLTYPE_FULL and NEURONAL_TYPES
# (rnaseq_TUY35595/app/R/helpers.R:288-310). Do not reproduce it here.
# =====================================================================
parse_barcodes <- function(cn) {
  p <- strsplit(cn, "_", fixed = TRUE)
  n <- lengths(p)
  stopifnot(all(n >= 7))
  data.frame(
    cell     = cn,
    sex      = vapply(p, `[`, "", 1L),
    geno     = vapply(p, `[`, "", 2L),
    inj      = vapply(p, `[`, "", 3L),
    time     = as.integer(vapply(p, `[`, "", 4L)),
    rep      = vapply(p, `[`, "", 5L),
    celltype = mapply(function(v, k) paste(v[6:(k - 1L)], collapse = "_"), p, n, USE.NAMES = FALSE),
    bc       = mapply(function(v, k) v[k], p, n, USE.NAMES = FALSE),
    stringsAsFactors = FALSE
  )
}

# =====================================================================
# Design constants (verified against the matrix, not assumed)
# =====================================================================
INJ_ORDER <- c("Naive", "Crush", "ScNT", "SpNT", "CFA", "Paclitaxel")

# short strip labels — "Paclitaxel" overflows a 2-column facet strip
INJ_SHORT <- c(Naive = "Naive", Crush = "Crush", ScNT = "ScNT",
               SpNT = "SpNT", CFA = "CFA", Paclitaxel = "Pacli")

INJ_FULL <- c(
  Naive      = "Naive (uninjured)",
  Crush      = "Sciatic nerve crush (regeneration-permissive)",
  ScNT       = "Sciatic nerve transection (axotomy)",
  SpNT       = "Spinal nerve transection (axotomy)",
  CFA        = "Complete Freund's adjuvant (inflammation, no axotomy)",
  Paclitaxel = "Paclitaxel (chemotherapy neuropathy, no axotomy)"
)

# The two non-axotomy insults. The panel should be near-SILENT here if it is
# genuinely axotomy-specific — this is the single-cell analogue of MINP-negative.
NON_AXOTOMY <- c("CFA", "Paclitaxel")
AXOTOMY     <- c("Crush", "ScNT", "SpNT")

# Corrected 20-type nomenclature (supersedes the app's mangled CELLTYPE_FULL).
CELLTYPE_FULL <- c(
  NP                  = "Non-peptidergic nociceptor",
  PEP1                = "Peptidergic nociceptor 1",
  PEP2                = "Peptidergic nociceptor 2",
  NF1                 = "Myelinated A-fibre LTMR / neurofilament 1",
  NF2                 = "Myelinated A-fibre LTMR / neurofilament 2",
  NF3                 = "Proprioceptor / neurofilament 3",
  cLTMR1              = "C-fibre low-threshold mechanoreceptor 1",
  p_cLTMR2            = "Putative C-fibre LTMR 2",
  SST                 = "Somatostatin+ pruriceptor",
  Satglia             = "Satellite glia",
  Schwann_M           = "Myelinating Schwann cell",
  Schwann_N           = "Non-myelinating Schwann cell",
  `Repair schwann_N`  = "Repair Schwann cell (injury-induced)",
  Fibroblast          = "Fibroblast",
  `Repair fibroblast` = "Repair fibroblast (injury-induced)",
  Endothelial         = "Endothelial cell",
  Pericyte            = "Pericyte",
  Macrophage          = "Macrophage",
  `B cell`            = "B cell",
  Neutrophil          = "Neutrophil"
)

NEURONAL_TYPES <- c("NP", "PEP1", "PEP2", "NF1", "NF2", "NF3", "cLTMR1", "p_cLTMR2", "SST")

# Cell types that exist only after injury — invisible to a naive-only reference.
INJURY_INDUCED_TYPES <- c("Repair schwann_N", "Repair fibroblast")

# Minimum cells for a (sample x celltype) pseudobulk unit to be trusted.
MIN_CELLS <- 10

# =====================================================================
# The 47-gene conserved injury panel
# ---------------------------------------------------------------------
# Source: rnaseq_TUY35595/app/data/conserved_injury_markers.rds, built by
# app/build_injury_markers.R — genes significantly UP (padj<0.05, log2FC>1) in
# nerve injury in >=2 of {mouse SNI, rat SNI, macaque ipsi} AND flat in mouse
# MINP (padj>0.5, |log2FC|<0.5).
#
# Symbol reconciliation against this 2020-era matrix:
#   Insyn2a -> present as its old symbol Fam196a
#   Crybg1  -> absent entirely (old symbol Aim1; matrix has only Aim1l = Crybg2)
# => 46 of 47 usable.
#
# Sprr1a and Sox11 are the RAG archetypes behind the app's "Sprr1a-negative"
# framing but are NOT in the 47 — Sprr1a only because it has no macaque ortholog
# (app/R/tab_cross_pca.R:156-157), so it cannot survive the three-way inner join.
# Both are well detected here (16.6% / 27.4% of cells), so we carry them as
# explicitly flagged reference rows rather than silently mixing them in.
# =====================================================================
PANEL_ALIAS <- c(Insyn2a = "Fam196a")
REFERENCE_GENES <- c("Sprr1a", "Sox11")

load_panel <- function() {
  mk <- readRDS(file.path(APP, "app/data/conserved_injury_markers.rds"))
  stopifnot(nrow(mk) == 47L, "Gene" %in% names(mk))
  data.frame(
    gene       = mk$Gene,
    matrix_sym = ifelse(mk$Gene %in% names(PANEL_ALIAS), PANEL_ALIAS[mk$Gene], mk$Gene),
    mean_sni   = mk$`Mean SNI log2FC`,
    n_sig      = mk$`n species sig`,
    minp_l2fc  = mk$`MINP log2FC`,
    role       = "panel",
    stringsAsFactors = FALSE
  )
}

# Panel rows + the two flagged reference rows, in plot order
# (panel sorted by mean SNI log2FC descending, references appended at the bottom).
panel_rows <- function(rownames_available) {
  p <- load_panel()
  p <- p[order(-p$mean_sni), ]
  ref <- data.frame(gene = REFERENCE_GENES, matrix_sym = REFERENCE_GENES,
                    mean_sni = NA_real_, n_sig = NA_integer_, minp_l2fc = NA_real_,
                    role = "reference", stringsAsFactors = FALSE)
  out <- rbind(p, ref)
  out$present <- out$matrix_sym %in% rownames_available
  out
}

# =====================================================================
# House styling — mirrors rnaseq_TUY35595/app/R/helpers.R so these figures sit
# next to the app's existing 47-marker heatmap without looking foreign.
# =====================================================================
CPCA_DIVERGE <- list(list(0, "#2166AC"), list(0.5, "#F7F7F7"), list(1, "#B2182B"))
DIVERGE_HEX  <- c("#2166AC", "#F7F7F7", "#B2182B")

theme_pub <- function(base_size = 15) {
  ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      panel.border  = ggplot2::element_rect(colour = "black", fill = NA, linewidth = 0.9),
      panel.grid    = ggplot2::element_blank(),
      axis.ticks    = ggplot2::element_line(colour = "black"),
      axis.text     = ggplot2::element_text(colour = "black"),
      plot.title    = ggplot2::element_text(face = "bold"),
      plot.subtitle = ggplot2::element_text(colour = "grey30"),
      legend.key    = ggplot2::element_blank())
}

ply_pub <- function(p, fontsize = 14) {
  ax <- list(showgrid = FALSE, zeroline = FALSE, showline = TRUE, mirror = TRUE,
             linecolor = "black", linewidth = 1.3, ticks = "outside", tickcolor = "black")
  plotly::layout(p, font = list(size = fontsize), xaxis = ax, yaxis = ax)
}

cluster_order <- function(m, do_row = TRUE, do_col = FALSE) {
  ro <- seq_len(nrow(m)); co <- seq_len(ncol(m))
  safe <- function(x) { x[!is.finite(x)] <- 0; x }
  if (isTRUE(do_row) && nrow(m) > 2)
    ro <- tryCatch(stats::hclust(stats::dist(safe(m)))$order, error = function(e) ro)
  if (isTRUE(do_col) && ncol(m) > 2)
    co <- tryCatch(stats::hclust(stats::dist(t(safe(m))))$order, error = function(e) co)
  list(row = ro, col = co)
}

`%||%` <- function(a, b) if (is.null(a)) b else a

msg <- function(...) { cat(format(Sys.time(), "[%H:%M:%S] "), ..., "\n", sep = ""); flush.console() }
