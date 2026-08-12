# =====================================================================
# 02_pseudobulk.R — aggregate raw single-cell UMI counts into pseudobulk.
#
# Two aggregations, both by SUMMING raw counts (never averaging normalised
# values — pseudobulk DE needs counts so edgeR/limma can model the mean-variance
# relationship properly):
#   data/pb_sample.rds           genes x sample                (whole DRG)
#   data/pb_sample_celltype.rds  genes x (sample x celltype)
#
# Aggregation uses one sparse matrix product (counts %*% indicator) rather than a
# per-group loop — seconds instead of minutes on 141k cells.
#
# Verification: total counts must be conserved exactly, and per-sample column
# sums must match colSums over the corresponding cells.
# =====================================================================
suppressMessages(library(Matrix))
source("/hpf/projects/msalter/sghazis/painseq/analysis/00_common.R")

# grouping vector -> sparse cells x groups 0/1 indicator
indicator <- function(g) {
  g <- factor(g)
  sparseMatrix(i = seq_along(g), j = as.integer(g), x = 1,
               dims = c(length(g), nlevels(g)), dimnames = list(NULL, levels(g)))
}

aggregate_counts <- function(counts, g) {
  M  <- indicator(g)
  pb <- counts %*% M                       # genes x groups, still sparse
  as.matrix(pb)                            # groups are few; dense is fine downstream
}

build <- function(counts_path, meta_path, prefix) {
  msg("=== ", prefix, " ===")
  x  <- read_counts(counts_path)
  md <- readRDS(meta_path)
  stopifnot(identical(colnames(x), md$cell))     # row order must match the matrix
  total_in <- sum(x)
  msg("  ", nrow(x), " genes x ", ncol(x), " cells | total UMI ", format(total_in, big.mark = ","))

  # ---- 1. whole-DRG pseudobulk, one column per sample ----------------------
  pb_s <- aggregate_counts(x, md$sample)
  n_s  <- as.integer(table(factor(md$sample))[colnames(pb_s)])
  meta_s <- unique(md[, c("sample", "sex", "geno", "inj", "time", "rep", "group")])
  meta_s <- meta_s[match(colnames(pb_s), meta_s$sample), ]
  meta_s$n_cells <- n_s
  rownames(meta_s) <- NULL
  stopifnot(sum(pb_s) == total_in)               # counts conserved
  msg("  pb_sample: ", nrow(pb_s), " genes x ", ncol(pb_s), " samples")

  # ---- 2. sample x celltype pseudobulk -------------------------------------
  key   <- paste(md$sample, md$celltype, sep = "|")
  pb_sc <- aggregate_counts(x, key)
  n_sc  <- as.integer(table(factor(key))[colnames(pb_sc)])
  sp    <- strsplit(colnames(pb_sc), "|", fixed = TRUE)
  meta_sc <- data.frame(unit     = colnames(pb_sc),
                        sample   = vapply(sp, `[`, "", 1L),
                        celltype = vapply(sp, `[`, "", 2L),
                        n_cells  = n_sc, stringsAsFactors = FALSE)
  j <- match(meta_sc$sample, meta_s$sample)
  meta_sc <- cbind(meta_sc, meta_s[j, c("sex", "geno", "inj", "time", "rep", "group")])
  rownames(meta_sc) <- NULL
  stopifnot(sum(pb_sc) == total_in)
  msg("  pb_sample_celltype: ", nrow(pb_sc), " genes x ", ncol(pb_sc), " units",
      " | units with >= ", MIN_CELLS, " cells: ", sum(meta_sc$n_cells >= MIN_CELLS))

  # ---- verification: spot-check three samples against a direct colSums ------
  set.seed(1)
  for (s in sample(meta_s$sample, min(3L, nrow(meta_s)))) {
    direct <- sum(x[, md$sample == s, drop = FALSE])
    stopifnot(identical(as.numeric(direct), as.numeric(sum(pb_s[, s]))))
  }
  msg("  verified: total counts conserved and 3 random samples match direct colSums")

  saveRDS(list(counts = pb_s,  meta = meta_s),  file.path(DATA, paste0(prefix, "_sample.rds")))
  saveRDS(list(counts = pb_sc, meta = meta_sc), file.path(DATA, paste0(prefix, "_sample_celltype.rds")))
  msg("  wrote data/", prefix, "_sample.rds and data/", prefix, "_sample_celltype.rds")
  rm(x); invisible(gc())
}

build(C57_RDS, file.path(DATA, "cell_meta.rds"),      "pb")
build(ATF3_GZ, file.path(DATA, "cell_meta_atf3.rds"), "pb_atf3")

msg("done")
