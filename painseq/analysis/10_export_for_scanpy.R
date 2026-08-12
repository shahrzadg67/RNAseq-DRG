# =====================================================================
# 10_export_for_scanpy.R — hand the count matrices to python.
#
# Writes the dgCMatrix CSC triplets as raw binary (indices/pointers/values) plus
# plain-text gene and cell metadata. Binary beats MatrixMarket here: the C57
# matrix has 181M non-zeros, which is ~3.5 GB as MM text and slow for scipy to
# parse, versus ~1.4 GB binary that np.fromfile reads in seconds.
#
# Counts are small integers (max 2282), so values go out as int32, not float64.
#
# Out (per dataset, under data/export_<name>/):
#   indices.i32  p.i32  values.i32   -> scipy.sparse.csc_matrix((v, i, p))
#   shape.txt  genes.txt  obs.csv
# =====================================================================
suppressMessages(library(Matrix))
source("/hpf/projects/msalter/sghazis/painseq/analysis/00_common.R")

export_one <- function(counts_path, meta_path, name) {
  msg("=== exporting ", name, " ===")
  x  <- read_counts(counts_path)
  md <- readRDS(meta_path)
  stopifnot(identical(colnames(x), md$cell))
  stopifnot(inherits(x, "dgCMatrix"))

  d <- file.path(DATA, paste0("export_", name)); dir.create(d, showWarnings = FALSE, recursive = TRUE)

  stopifnot(max(x@x) < .Machine$integer.max, all(x@x == floor(x@x)))
  writeBin(as.integer(x@i), file.path(d, "indices.i32"), size = 4L)
  writeBin(as.integer(x@p), file.path(d, "p.i32"),       size = 4L)
  writeBin(as.integer(x@x), file.path(d, "values.i32"),  size = 4L)
  writeLines(as.character(dim(x)), file.path(d, "shape.txt"))     # genes, cells
  writeLines(rownames(x), file.path(d, "genes.txt"))

  # GEO's cell names are NOT unique: 5,593 of the 141,093 C57 names collide,
  # because two nuclei from the same sample and cell type can draw the same
  # barcode suffix. All aggregation in this project is positional, so nothing
  # upstream is affected, but emit a guaranteed-unique id so anything downstream
  # can join safely instead of silently fanning out.
  md$cell_uid <- paste0(md$cell, "#", seq_len(nrow(md)))

  # QC metrics computed here (cheap in R, and keeps python from needing all genes)
  mt <- grep("^mt-", rownames(x), value = TRUE)
  md$total_umi <- Matrix::colSums(x)
  md$n_genes   <- Matrix::colSums(x > 0)
  md$pct_mito  <- if (length(mt)) 100 * Matrix::colSums(x[mt, , drop = FALSE]) / md$total_umi else 0
  write.csv(md, file.path(d, "obs.csv"), row.names = FALSE)

  msg("  ", nrow(x), " x ", ncol(x), " | nnz ", format(length(x@x), big.mark = ","),
      " | ", length(mt), " mito genes | wrote ", d)
  msg("  QC: median UMI ", median(md$total_umi), ", median genes ", median(md$n_genes),
      ", median %mito ", round(median(md$pct_mito), 3))
  rm(x); invisible(gc())
}

export_one(C57_RDS, file.path(DATA, "cell_meta.rds"),      "c57")
export_one(ATF3_GZ, file.path(DATA, "cell_meta_atf3.rds"), "atf3")

# the panel gene list, so python can score cells without re-deriving it
pg <- readRDS(file.path(DATA, "panel_genes.rds"))
writeLines(pg$matrix_sym[pg$present & pg$role == "panel"], file.path(DATA, "panel_symbols.txt"))
msg("wrote data/panel_symbols.txt (", sum(pg$present & pg$role == "panel"), " genes)")
msg("done")
