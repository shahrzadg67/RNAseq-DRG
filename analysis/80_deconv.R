# =====================================================================
# 80_deconv.R — scRNA-seq -> bulk deconvolution (MuSiC)
#
# Estimates the cell-type composition of every bulk RNA-seq sample using a
# published mouse DRG single-cell atlas as the reference, then compares the
# estimated proportions across ALL study contrasts (not only MINP vs Sham).
#
# Reference : Renthal et al. 2020, mouse DRG snRNA-seq (GEO GSE154659, C57).
#             Cell-type label is encoded in each cell barcode:
#             sex_strain_injury_time_rep_CELLTYPE_barcode
#             We use the NAIVE cells only -> clean baseline signatures, not
#             injury-reprogrammed transcriptomes.
# Bulk      : results/star_salmon/salmon.merged.gene_counts.tsv (raw counts;
#             MuSiC needs counts, not VST). Matched to the reference by SYMBOL.
# Method    : MuSiC::music_prop (weighted NNLS).
#
# Output    : app/data/deconv.rds — proportions (long + wide), sample meta,
#             per-cell-type stats for every contrast, and reference provenance.
# =====================================================================
suppressMessages({
  library(Matrix)
  library(SingleCellExperiment)
  library(MuSiC)
})

PROJ    <- "/hpf/projects/msalter/sghazis/rnaseq_TUY35595"
REF_RDS <- file.path(PROJ, "refs/scref/GSE154659_C57_Raw_counts.RDS")
BULK    <- file.path(PROJ, "results/star_salmon/salmon.merged.gene_counts.tsv")
OUT     <- file.path(PROJ, "app/data/deconv.rds")
MIN_CELLS_PER_TYPE <- 20   # drop cell types with too few naive cells for a stable signature

# ---- study design (mirror app/R/helpers.R) -------------------------------
CONTRAST_GROUPS <- list(
  NP_vs_Sham       = c("NP_F","NP_M","Sham_F","Sham_M"),
  NP_M_vs_Sham_M   = c("NP_M","Sham_M"),
  NP_F_vs_Sham_F   = c("NP_F","Sham_F"),
  Sham_M_vs_Sham_F = c("Sham_M","Sham_F"),
  NP_M_vs_NP_F     = c("NP_M","NP_F"),
  SNI_vs_Sham      = c("SNI_F","SNI_M","Sham_F","Sham_M"),
  SNI_vs_NP        = c("SNI_F","SNI_M","NP_F","NP_M"),
  SNI_vs_Sham_M    = c("SNI_F","SNI_M","Sham_M"),
  SNI_vs_Sham_F    = c("SNI_F","SNI_M","Sham_F"),
  SNI_vs_NP_M      = c("SNI_F","SNI_M","NP_M"),
  SNI_vs_NP_F      = c("SNI_F","SNI_M","NP_F")
)

# =====================================================================
# 1. Reference: build a Naive-cell SingleCellExperiment
# =====================================================================
cat(">> reading reference:", REF_RDS, "\n"); flush.console()
ref <- readRDS(REF_RDS)                         # dgCMatrix genes x cells (symbols x barcodes)
stopifnot(inherits(ref, "sparseMatrix"))
cat("   reference:", nrow(ref), "genes x", ncol(ref), "cells\n")

sp       <- strsplit(colnames(ref), "_")
celltype <- vapply(sp, function(z) z[length(z) - 1L], character(1))   # 2nd-to-last field
sex_ref  <- vapply(sp, function(z) z[1], character(1))
injury   <- vapply(sp, function(z) z[3], character(1))
time_ref <- vapply(sp, function(z) z[4], character(1))
rep_ref  <- vapply(sp, function(z) z[5], character(1))

keep <- injury == "Naive"
cat("   naive cells:", sum(keep), "\n")
ref  <- ref[, keep, drop = FALSE]
celltype <- celltype[keep]; sex_ref <- sex_ref[keep]
rep_ref <- rep_ref[keep];   time_ref <- time_ref[keep]

# drop cell types with too few naive cells (unstable signature)
ct_n <- table(celltype)
drop_ct <- names(ct_n)[ct_n < MIN_CELLS_PER_TYPE]
if (length(drop_ct)) {
  cat("   dropping rare cell types (<", MIN_CELLS_PER_TYPE, "cells):",
      paste(sprintf("%s(%d)", drop_ct, ct_n[drop_ct]), collapse = ", "), "\n")
  keep2 <- !celltype %in% drop_ct
  ref <- ref[, keep2, drop = FALSE]
  celltype <- celltype[keep2]; sex_ref <- sex_ref[keep2]; rep_ref <- rep_ref[keep2]
}

# dedupe reference gene symbols (keep highest-total row per symbol) -> sparse-safe
if (any(duplicated(rownames(ref)))) {
  cat("   deduping", sum(duplicated(rownames(ref))), "duplicate reference symbols\n")
  ord <- order(Matrix::rowSums(ref), decreasing = TRUE)
  ref <- ref[ord, , drop = FALSE]
  ref <- ref[!duplicated(rownames(ref)), , drop = FALSE]
}

# MuSiC "subject" grouping: biological replicate (sex x rep) for cross-subject variance
sampleID <- paste(sex_ref, rep_ref, sep = "_")
cat("   cell types kept:", length(unique(celltype)),
    "| subjects:", length(unique(sampleID)), "\n")
print(sort(table(celltype), decreasing = TRUE))

sce <- SingleCellExperiment(
  assays  = list(counts = ref),
  colData = DataFrame(cellType = celltype, sampleID = sampleID, sex = sex_ref)
)
celltypes <- sort(unique(celltype))

# =====================================================================
# 2. Bulk: raw gene counts, aggregated to unique SYMBOLS
# =====================================================================
cat(">> reading bulk counts:", BULK, "\n"); flush.console()
bt  <- read.delim(BULK, check.names = FALSE, stringsAsFactors = FALSE)
sym <- bt$gene_name
mat <- as.matrix(bt[, !(colnames(bt) %in% c("gene_id", "gene_name")), drop = FALSE])
storage.mode(mat) <- "double"
ok  <- !is.na(sym) & sym != ""
mat <- mat[ok, , drop = FALSE]; sym <- sym[ok]
bulk <- rowsum(mat, group = sym)                 # sum counts of duplicate symbols
bulk <- round(bulk)                              # MuSiC expects integer-ish counts
storage.mode(bulk) <- "integer"
cat("   bulk:", nrow(bulk), "symbols x", ncol(bulk), "samples\n")
cat("   samples:", paste(colnames(bulk), collapse = ", "), "\n")

common <- intersect(rownames(bulk), rownames(sce))
cat("   genes shared bulk<->reference:", length(common), "\n")
stopifnot(length(common) > 5000)

# =====================================================================
# 3. Run MuSiC
# =====================================================================
cat(">> running MuSiC::music_prop ...\n"); flush.console()
set_seed_ok <- TRUE
est <- music_prop(
  bulk.mtx  = bulk,
  sc.sce    = sce,
  clusters  = "cellType",
  samples   = "sampleID",
  select.ct = celltypes,
  verbose   = TRUE
)
prop <- est$Est.prop.weighted                    # samples x cell types
prop <- prop[colnames(bulk), , drop = FALSE]
cat("   proportions:", nrow(prop), "samples x", ncol(prop), "cell types\n")
cat("   row sums (should ~1):", paste(round(rowSums(prop), 3), collapse = " "), "\n")

# =====================================================================
# 4. Sample metadata + long proportions
# =====================================================================
samp <- rownames(prop)
cond <- sub("_.*$", "", samp)                    # NP / Sham / SNI
sx   <- sub("^.*_([FM]).*$", "\\1", samp)        # sex letter after the underscore
grp  <- paste0(cond, "_", sx)                    # NP_F / Sham_M / ... (matches GROUP_SEX)
meta <- data.frame(sample = samp, condition = cond, sex = sx, group = grp,
                   stringsAsFactors = FALSE)

long <- do.call(rbind, lapply(seq_along(samp), function(i) {
  data.frame(sample = samp[i], condition = cond[i], sex = sx[i], group = grp[i],
             cellType = colnames(prop), proportion = as.numeric(prop[i, ]),
             stringsAsFactors = FALSE)
}))

# =====================================================================
# 5. Per-cell-type stats for EVERY contrast
# =====================================================================
side_samples <- function(token) {
  # token is a group (NP_M) or a condition (NP/Sham); return matching bulk samples
  meta$sample[meta$group == token | meta$condition == token]
}
safe_wilcox <- function(a, b) tryCatch(
  suppressWarnings(wilcox.test(a, b, exact = FALSE)$p.value), error = function(e) NA_real_)
safe_t <- function(a, b) tryCatch(
  suppressWarnings(t.test(a, b)$p.value), error = function(e) NA_real_)

stats <- do.call(rbind, lapply(names(CONTRAST_GROUPS), function(cn) {
  sides <- strsplit(cn, "_vs_")[[1]]                      # e.g. NP_M , Sham_M
  A <- side_samples(sides[1]); B <- side_samples(sides[2])
  do.call(rbind, lapply(colnames(prop), function(ct) {
    a <- prop[A, ct]; b <- prop[B, ct]
    data.frame(
      contrast = cn, cellType = ct,
      mean_A = mean(a), mean_B = mean(b),
      diff   = mean(a) - mean(b),
      log2FC = log2((mean(a) + 1e-6) / (mean(b) + 1e-6)),
      p_wilcox = safe_wilcox(a, b), p_t = safe_t(a, b),
      n_A = length(A), n_B = length(B),
      stringsAsFactors = FALSE)
  }))
}))
# BH adjust across cell types within each contrast
stats <- do.call(rbind, lapply(split(stats, stats$contrast), function(d) {
  d$padj_wilcox <- p.adjust(d$p_wilcox, "BH")
  d$padj_t      <- p.adjust(d$p_t, "BH"); d
}))
rownames(stats) <- NULL

# =====================================================================
# 6. Save
# =====================================================================
ref_info <- list(
  source    = "Renthal et al. 2020, mouse DRG snRNA-seq",
  accession = "GEO GSE154659 (C57, Naive cells only)",
  n_cells   = ncol(sce), n_celltypes = length(celltypes),
  n_genes_shared = length(common), method = "MuSiC (weighted NNLS)",
  min_cells_per_type = MIN_CELLS_PER_TYPE)

saveRDS(list(prop_long = long, prop_wide = prop, meta = meta, stats = stats,
             celltypes = colnames(prop), contrast_groups = CONTRAST_GROUPS,
             ref_info = ref_info),
        OUT)
cat(">> wrote", OUT, "\n")
cat(">> DONE. cell types:", paste(colnames(prop), collapse = ", "), "\n")
