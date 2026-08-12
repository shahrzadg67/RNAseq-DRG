# =====================================================================
# 01_parse_meta.R — build the cell metadata table from the barcode strings.
#
# GSE154659 ships bare dgCMatrix count matrices with NO metadata slot; every
# annotation lives in the column name. This script decodes them, verifies the
# decode against the known design, and caches the result.
#
# Out: data/cell_meta.rds       (C57 atlas, 141,093 rows)
#      data/cell_meta_atf3.rds  (Atf3 WT/KO, 17,665 rows)
#      data/design_summary.txt  (human-readable design tables)
# =====================================================================
suppressMessages(library(Matrix))
source("/hpf/projects/msalter/sghazis/painseq/analysis/00_common.R")

sink_both <- file(file.path(DATA, "design_summary.txt"), open = "wt")

report <- function(md, label, expect_cells, expect_samples, expect_types) {
  md$sample <- paste(md$sex, md$geno, md$inj, md$time, md$rep, sep = "_")
  md$inj    <- factor(md$inj, levels = intersect(INJ_ORDER, unique(md$inj)))
  md$group  <- paste(md$inj, md$time, sep = "_")          # injury x time

  # render the whole block to a character vector once, then emit it twice
  as_text <- function(x) capture.output(print(x))
  # count distinct SAMPLES, not distinct rep labels — female and Mrgprd-cre
  # samples reuse rep1/rep2, so unique(inj,time,rep) would undercount them.
  r <- unique(md[, c("sample", "inj", "time")])
  blk <- c(
    sprintf("\n=========== %s ===========", label),
    sprintf("cells: %d | samples: %d | cell types: %d | injury x time groups: %d",
            nrow(md), length(unique(md$sample)), length(unique(md$celltype)),
            length(unique(md$group))),
    "\n-- injury x time (cells) --",                as_text(table(md$inj, md$time)),
    "\n-- replicate samples per injury x time --",  as_text(table(r$inj, r$time)),
    "\n-- cell types --",                           as_text(sort(table(md$celltype), decreasing = TRUE)),
    "\n-- sex x genotype --",                       as_text(table(md$sex, md$geno)))
  writeLines(blk)
  writeLines(blk, sink_both)

  # ---- verification: assert the decode matches the known design -------------
  stopifnot(nrow(md) == expect_cells)
  stopifnot(length(unique(md$sample)) == expect_samples)
  stopifnot(length(unique(md$celltype)) == expect_types)
  # the four underscore-containing cell types must survive intact
  underscore_types <- grep("_| ", unique(md$celltype), value = TRUE)
  msg("underscore/space cell types preserved: ", paste(sort(underscore_types), collapse = ", "))
  stopifnot(!any(md$celltype %in% c("M", "N", "cLTMR2", "schwann")))  # the 80_deconv.R mangling
  invisible(md)
}

# ---- C57 atlas -----------------------------------------------------------
msg("reading C57 counts: ", C57_RDS)
x <- read_counts(C57_RDS)
msg("  ", nrow(x), " genes x ", ncol(x), " cells")
md <- parse_barcodes(colnames(x))
md <- report(md, "C57 atlas", expect_cells = 141093L, expect_samples = 66L, expect_types = 20L)
stopifnot(setequal(unique(md$celltype), names(CELLTYPE_FULL)))
stopifnot(sum(md$celltype == "p_cLTMR2") == 3000L,
          sum(md$celltype == "Repair schwann_N") == 1824L)
saveRDS(md, file.path(DATA, "cell_meta.rds"))
msg("wrote data/cell_meta.rds")

# panel gene availability check, recorded once here so later scripts can assume it
pr <- panel_rows(rownames(x))
msg("panel genes present in matrix: ", sum(pr$present[pr$role == "panel"]), "/47",
    "  | missing: ", paste(pr$gene[!pr$present], collapse = ", "))
det <- Matrix::rowSums(x[pr$matrix_sym[pr$present], , drop = FALSE] > 0)
umi <- Matrix::rowSums(x[pr$matrix_sym[pr$present], , drop = FALSE])
pr$pct_cells <- NA_real_; pr$total_umi <- NA_integer_
pr$pct_cells[pr$present] <- round(100 * det / ncol(x), 3)
pr$total_umi[pr$present] <- as.integer(umi)
pr$low_detection <- !is.na(pr$pct_cells) & pr$pct_cells < 0.5
saveRDS(pr, file.path(DATA, "panel_genes.rds"))
msg("wrote data/panel_genes.rds — ", sum(pr$low_detection, na.rm = TRUE),
    " genes flagged low-detection (<0.5% of cells)")
rm(x); invisible(gc())

# ---- Atf3 WT/KO ----------------------------------------------------------
msg("reading Atf3 WT/KO counts: ", ATF3_GZ)
a <- read_counts(ATF3_GZ)
msg("  ", nrow(a), " genes x ", ncol(a), " cells")
mda <- parse_barcodes(colnames(a))
mda <- report(mda, "Atf3 WT/KO", expect_cells = 17665L, expect_samples = 12L, expect_types = 9L)
stopifnot(setequal(unique(mda$celltype), NEURONAL_TYPES))
saveRDS(mda, file.path(DATA, "cell_meta_atf3.rds"))
msg("wrote data/cell_meta_atf3.rds")

close(sink_both)
msg("done — design summary in data/design_summary.txt")
