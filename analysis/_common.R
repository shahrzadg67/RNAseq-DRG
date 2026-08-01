# ============================================================================
# _common.R  —  shared config, paths, contrasts, and helpers for TUY35595
# Sourced by every analysis script so gene- and transcript-level runs stay
# perfectly parallel and the 5 contrasts are defined in exactly one place.
# ============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(tibble)
})

## ---- paths -----------------------------------------------------------------
PROJ     <- "/hpf/projects/msalter/sghazis/rnaseq_TUY35595"
# Dev override: set TUY_DEV_DATA to a dir holding synthetic salmon.merged.*.tsv
# to develop/test before the real pipeline output lands.
RESULTS  <- if (nzchar(Sys.getenv("TUY_DEV_DATA"))) Sys.getenv("TUY_DEV_DATA") else file.path(PROJ, "results", "star_salmon")
META_CSV <- file.path(PROJ, "metadata.csv")
APP_DATA <- file.path(PROJ, "app", "data")     # ONLY app-bundled artifacts go here
CACHE    <- file.path(PROJ, "analysis", "cache") # heavy intermediates (NOT bundled)
dir.create(APP_DATA, showWarnings = FALSE, recursive = TRUE)
dir.create(CACHE, showWarnings = FALSE, recursive = TRUE)
save_cache <- function(obj, name) saveRDS(obj, file.path(CACHE, paste0(name, ".rds")))
read_cache <- function(name) readRDS(file.path(CACHE, paste0(name, ".rds")))

# Count matrices produced by nf-core star_salmon (both levels, in parallel).
COUNTS <- list(
  gene       = file.path(RESULTS, "salmon.merged.gene_counts.tsv"),
  transcript = file.path(RESULTS, "salmon.merged.transcript_counts.tsv")
)
TX2GENE <- file.path(RESULTS, "tx2gene.tsv")   # nf-core writes this in star_salmon/
LEVELS  <- c("gene", "transcript")

## ---- experimental design ---------------------------------------------------
# 6 groups exist (incl. SNI). SNI has n=1/sex (2 samples, replicate 1 only), so its
# contrasts are EXPLORATORY — fold-changes are estimable (dispersion is shared across
# all samples) but p-values/FDR are unreliable and confounded with sex. SNI is now
# INCLUDED in the DE fits so the six pooled-SNI contrasts below can be computed; the
# app flags every SNI view as exploratory.
DE_GROUPS <- c("NP_F", "NP_M", "Sham_F", "Sham_M", "SNI_F", "SNI_M")  # order defines contrast vectors

# Each contrast = a numeric vector over DE_GROUPS (a linear combination of group means).
# This same definition drives DESeq2 (contrast=) and edgeR (makeContrasts).
# NOTE: every vec MUST name all six groups (0 for uninvolved) or contrast_vec() aborts.
CONTRASTS <- list(
  NP_vs_Sham = list(
    label  = "NP (M+F) vs Sham (M+F)",
    vec    = c(NP_F = 0.5, NP_M = 0.5, Sham_F = -0.5, Sham_M = -0.5, SNI_F = 0, SNI_M = 0)),
  NP_M_vs_Sham_M = list(
    label  = "NP male vs Sham male",
    vec    = c(NP_F = 0, NP_M = 1, Sham_F = 0, Sham_M = -1, SNI_F = 0, SNI_M = 0)),
  NP_F_vs_Sham_F = list(
    label  = "NP female vs Sham female",
    vec    = c(NP_F = 1, NP_M = 0, Sham_F = -1, Sham_M = 0, SNI_F = 0, SNI_M = 0)),
  Sham_M_vs_Sham_F = list(
    label  = "Sham male vs Sham female",
    vec    = c(NP_F = 0, NP_M = 0, Sham_F = -1, Sham_M = 1, SNI_F = 0, SNI_M = 0)),
  NP_M_vs_NP_F = list(
    label  = "NP male vs NP female",
    vec    = c(NP_F = -1, NP_M = 1, Sham_F = 0, Sham_M = 0, SNI_F = 0, SNI_M = 0)),
  # --- SNI contrasts (pooled SNI_F + SNI_M on the left; exploratory, n=1/sex) ---
  SNI_vs_Sham = list(
    label  = "SNI (M+F) vs Sham (M+F)",
    vec    = c(NP_F = 0, NP_M = 0, Sham_F = -0.5, Sham_M = -0.5, SNI_F = 0.5, SNI_M = 0.5)),
  SNI_vs_NP = list(
    label  = "SNI (M+F) vs NP (M+F)",
    vec    = c(NP_F = -0.5, NP_M = -0.5, Sham_F = 0, Sham_M = 0, SNI_F = 0.5, SNI_M = 0.5)),
  SNI_vs_Sham_M = list(
    label  = "SNI (M+F) vs Sham male",
    vec    = c(NP_F = 0, NP_M = 0, Sham_F = 0, Sham_M = -1, SNI_F = 0.5, SNI_M = 0.5)),
  SNI_vs_Sham_F = list(
    label  = "SNI (M+F) vs Sham female",
    vec    = c(NP_F = 0, NP_M = 0, Sham_F = -1, Sham_M = 0, SNI_F = 0.5, SNI_M = 0.5)),
  SNI_vs_NP_M = list(
    label  = "SNI (M+F) vs NP male",
    vec    = c(NP_F = 0, NP_M = -1, Sham_F = 0, Sham_M = 0, SNI_F = 0.5, SNI_M = 0.5)),
  SNI_vs_NP_F = list(
    label  = "SNI (M+F) vs NP female",
    vec    = c(NP_F = -1, NP_M = 0, Sham_F = 0, Sham_M = 0, SNI_F = 0.5, SNI_M = 0.5))
)
# Ensure every contrast vector is ordered to match DE_GROUPS exactly.
contrast_vec <- function(name) {
  v <- CONTRASTS[[name]]$vec[DE_GROUPS]
  stopifnot(!any(is.na(v)))
  v
}

## ---- metadata --------------------------------------------------------------
load_metadata <- function() {
  m <- read.csv(META_CSV, stringsAsFactors = FALSE)
  rownames(m) <- m$sample_id
  m$condition <- factor(m$condition, levels = c("Sham", "NP", "SNI"))  # Sham = ref
  m$sex       <- factor(m$sex, levels = c("F", "M"))
  m$group     <- factor(m$group)
  m
}

## ---- count-matrix loader (robust to nf-core column layout) ------------------
# Returns list(mat = integer matrix [features x samples], annot = data.frame).
# gene file cols : gene_id, gene_name, <samples...>
# tx   file cols : tx (or tx_id), gene_id, <samples...>
load_counts <- function(level) {
  f <- COUNTS[[level]]
  if (!file.exists(f)) stop("Count file not found (pipeline may still be running): ", f)
  raw <- read.delim(f, header = TRUE, check.names = FALSE)

  id_cols   <- intersect(c("gene_id", "gene_name", "tx", "tx_id", "transcript_id"),
                         colnames(raw))
  # feature key: gene_id for gene level; the transcript-id column for tx level.
  tx_key    <- intersect(c("tx", "tx_id", "transcript_id"), colnames(raw))[1]
  feat_key  <- if (level == "gene") "gene_id" else tx_key
  if (is.na(feat_key)) stop("Could not find feature-id column in ", f)
  samp_cols <- setdiff(colnames(raw), id_cols)

  mat <- as.matrix(raw[, samp_cols, drop = FALSE])
  storage.mode(mat) <- "double"
  rownames(mat) <- raw[[feat_key]]
  mat <- round(mat); storage.mode(mat) <- "integer"

  annot <- raw[, id_cols, drop = FALSE]
  annot$feature_id <- raw[[feat_key]]
  list(mat = mat, annot = annot)
}

## ---- transcript -> gene map + symbol annotation ----------------------------
load_tx2gene <- function() {
  if (file.exists(TX2GENE)) {
    t2g <- read.delim(TX2GENE, header = FALSE, stringsAsFactors = FALSE)
    # nf-core tx2gene: V1=tx, V2=gene_id, V3=gene_name (cols may vary)
    nm <- c("tx_id", "gene_id", "gene_name")[seq_len(ncol(t2g))]
    colnames(t2g) <- nm
    return(t2g)
  }
  NULL
}

# Map any feature ids to gene SYMBOL + ENTREZID via org.Mm.eg.db (for GSEA/labels).
# Robust to NA / non-Ensembl / unknown ids: only valid keys are looked up, the
# rest come back as NA instead of erroring.
annotate_genes <- function(gene_ids) {
  suppressPackageStartupMessages(library(org.Mm.eg.db))
  ens <- sub("\\..*$", "", gene_ids)                 # ENSMUSG000...".7" -> base id
  valid <- keys(org.Mm.eg.db, keytype = "ENSEMBL")
  q <- unique(ens[!is.na(ens) & ens %in% valid])
  if (length(q)) {
    ann <- suppressWarnings(AnnotationDbi::select(
      org.Mm.eg.db, keys = q,
      columns = c("SYMBOL", "ENTREZID"), keytype = "ENSEMBL"))
    ann <- ann[!duplicated(ann$ENSEMBL), ]
  } else {
    ann <- data.frame(ENSEMBL = character(), SYMBOL = character(), ENTREZID = character())
  }
  data.frame(gene_id = gene_ids, ensembl = ens,
             symbol  = ann$SYMBOL[match(ens, ann$ENSEMBL)],
             entrez  = ann$ENTREZID[match(ens, ann$ENSEMBL)],
             stringsAsFactors = FALSE)
}

## ---- small IO helpers ------------------------------------------------------
save_artifact <- function(obj, name) {
  saveRDS(obj, file.path(APP_DATA, paste0(name, ".rds")))
  invisible(file.path(APP_DATA, paste0(name, ".rds")))
}
# Per-contrast CSVs are for download/inspection — keep them OUT of the app bundle.
write_table_artifact <- function(df, name) {
  d <- file.path(CACHE, "tables"); dir.create(d, showWarnings = FALSE, recursive = TRUE)
  utils::write.csv(df, file.path(d, paste0(name, ".csv")), row.names = FALSE)
}
