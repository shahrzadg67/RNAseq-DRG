# =====================================================================
# 40_extract_percell.R — per-cell log-normalised expression for a defined gene
# universe. This is the substrate for everything pseudobulk cannot show.
#
# Gene universe:
#   - the 49 panel rows (47 conserved injury genes + Sprr1a, Sox11)
#   - the 17 MINP-associated genes from the bulk study
#   - data-driven subtype markers: top N per cell type, derived from NAIVE cells
#     only, so the identity score in 43_identity.R is not built on hand-picked
#     markers chosen with the answer in mind
#
# Normalisation: CP10K + log1p (counts / total UMI * 1e4, then log1p) — the same
# transform scanpy applied for the panel score, so the two are comparable.
#
# Out: data/percell_expr.rds   full 141,093 cells (analysis)
#      data/percell_app.rds    20,000-cell stratified subsample (Shiny explorer)
#      data/gene_universe.rds  the gene table with provenance per gene
# =====================================================================
suppressMessages(library(Matrix))
source("/hpf/projects/msalter/sghazis/painseq/analysis/00_common.R")

N_MARKERS_PER_TYPE <- 10
N_APP_CELLS        <- 20000

msg("reading atlas")
x  <- read_counts(C57_RDS)
md <- readRDS(file.path(DATA, "cell_meta.rds"))
stopifnot(identical(colnames(x), md$cell))

# ---- normalisation factor (per cell) -------------------------------------
tot <- Matrix::colSums(x)
stopifnot(all(tot > 0))

# lognorm for a set of gene rows, returned as a dense matrix (genes x cells)
lognorm_rows <- function(genes) {
  g <- intersect(genes, rownames(x))
  sub <- x[g, , drop = FALSE]
  # scale each cell to 1e4 then log1p; sub is small (few hundred rows) so dense is fine
  m <- as.matrix(sub)
  m <- log1p(sweep(m, 2, tot / 1e4, "/"))
  m
}

# =====================================================================
# 1. Data-driven subtype markers from NAIVE cells only
# =====================================================================
msg("deriving subtype markers from naive cells (", sum(md$inj == "Naive"), " cells)")
nv <- md$inj == "Naive"
xn <- x[, nv, drop = FALSE]; tn <- tot[nv]; ctn <- md$celltype[nv]

# mean CP10K-log1p per cell type, computed on the sparse matrix without densifying:
# sum of log1p is not separable, so use mean normalised (not logged) expression for
# ranking — monotonic in the same direction and cheap.
norm_sum <- function(mat, cols) Matrix::rowSums(mat[, cols, drop = FALSE])
cts <- names(CELLTYPE_FULL)
cts <- cts[cts %in% unique(ctn)]
# fraction of cells in the type expressing each gene, and mean normalised level
mark <- list()
for (ct in cts) {
  k <- which(ctn == ct)
  if (length(k) < 20) { msg("  skip ", ct, " (", length(k), " naive cells)"); next }
  inn  <- Matrix::rowMeans(xn[, k, drop = FALSE] > 0)
  outn <- Matrix::rowMeans(xn[, -k, drop = FALSE] > 0)
  # rank by detection-rate difference: robust in shallow snRNA-seq, where mean
  # level is dominated by a handful of deep cells
  d <- inn - outn
  top <- names(sort(d, decreasing = TRUE))[seq_len(N_MARKERS_PER_TYPE)]
  mark[[ct]] <- top
  msg(sprintf("  %-18s %s", ct, paste(head(top, 5), collapse = ", ")))
}
saveRDS(mark, file.path(DATA, "subtype_markers.rds"))

# =====================================================================
# 2. Assemble the gene universe
# =====================================================================
pg <- readRDS(file.path(DATA, "panel_genes.rds"))
panel_syms <- pg$matrix_sym[pg$present]

minp <- readRDS(file.path(APP, "app/data/minp_specific.rds"))
minp_syms <- intersect(unique(as.character(minp$Gene)), rownames(x))

marker_syms <- unique(unlist(mark))

universe <- unique(c(panel_syms, minp_syms, marker_syms))
universe <- intersect(universe, rownames(x))

gu <- data.frame(
  gene    = universe,
  panel   = universe %in% panel_syms,
  minp    = universe %in% minp_syms,
  marker  = universe %in% marker_syms,
  stringsAsFactors = FALSE)
gu$marker_for <- vapply(gu$gene, function(g) {
  hits <- names(mark)[vapply(mark, function(v) g %in% v, logical(1))]
  if (length(hits)) paste(hits, collapse = ",") else NA_character_
}, character(1))
# Display label. Panel genes keep the * / ‡ annotation used on every other figure,
# and show their DISPLAY symbol rather than the 2020-era matrix symbol (so the row
# reads "Insyn2a", not "Fam196a"). panel_genes.rds has no `label` column — that is
# built in 20_build_app_artifact.R — so construct it the same way here.
pg_label <- ifelse(pg$role == "reference", paste0(pg$gene, " ‡"),
            ifelse(pg$low_detection %in% TRUE, paste0(pg$gene, " *"), pg$gene))
gu$label <- gu$gene
i <- match(gu$gene, pg$matrix_sym)
gu$label[!is.na(i)] <- pg_label[i[!is.na(i)]]

msg("gene universe: ", nrow(gu), " genes (",
    sum(gu$panel), " panel, ", sum(gu$minp), " MINP, ", sum(gu$marker), " subtype markers)")
saveRDS(gu, file.path(DATA, "gene_universe.rds"))

# =====================================================================
# 3. Extract per-cell expression
# =====================================================================
msg("extracting per-cell CP10K+log1p for ", nrow(gu), " genes x ", ncol(x), " cells")
E <- lognorm_rows(gu$gene)
E <- round(E, 2)                      # 2 dp is far finer than the noise floor here
E <- Matrix(E, sparse = TRUE)         # most entries are exactly 0
msg("  ", nrow(E), " x ", ncol(E), " | nnz ", format(length(E@x), big.mark = ","),
    " (", round(100 * length(E@x) / prod(dim(E)), 1), "% non-zero) | ",
    format(object.size(E), units = "MB"))

saveRDS(list(expr = E, genes = gu, meta = md[, c("cell", "sex", "geno", "inj", "time",
                                                 "rep", "celltype", "sample", "group")]),
        file.path(DATA, "percell_expr.rds"), compress = "xz")
msg("wrote data/percell_expr.rds (",
    round(file.size(file.path(DATA, "percell_expr.rds")) / 1e6, 2), " MB)")

# =====================================================================
# 4. App subsample — stratified by cell type, with a floor for rare types
# =====================================================================
set.seed(11)
frac <- N_APP_CELLS / ncol(E)
idx <- unlist(lapply(split(seq_len(nrow(md)), md$celltype), function(ii)
  if (length(ii) > 400) sample(ii, max(400, round(length(ii) * frac))) else ii))
idx <- sort(idx)
msg("app subsample: ", length(idx), " cells")

u <- read.csv(file.path(DATA, "umap_c57.csv"), stringsAsFactors = FALSE)
stopifnot(nrow(u) == nrow(md))        # scanpy removed no cells, so order is preserved
app <- list(
  expr = E[, idx, drop = FALSE],
  genes = gu,
  # keep the row indices into cell_meta.rds / umap_c57.csv so downstream can align
  # exactly by position. Cell NAMES are not unique (5,593 GEO collisions), so a
  # name join would silently fan out.
  idx = idx,
  meta = data.frame(
    x = round(u$UMAP1[idx], 2), y = round(u$UMAP2[idx], 2),
    celltype = factor(md$celltype[idx]),
    inj      = factor(as.character(md$inj)[idx], levels = INJ_ORDER),
    time     = as.integer(md$time[idx]),
    score    = round(u$panel_score[idx], 3),
    stringsAsFactors = FALSE))
saveRDS(app, file.path(DATA, "percell_app.rds"), compress = "xz")
msg("wrote data/percell_app.rds (",
    round(file.size(file.path(DATA, "percell_app.rds")) / 1e6, 2), " MB, ",
    format(object.size(app), units = "MB"), " resident)")
rm(x, E); invisible(gc())

# =====================================================================
# 5. The Atf3-WT / Atf3-KO arm, same gene universe
# ---------------------------------------------------------------------
# This object was previously only used at whole-neuron pseudobulk, where the
# genotype x injury interaction has n=3 samples per cell and almost nothing
# reaches significance. Per-cell is where this experiment actually has power:
# the FRACTION of neurons that switch a gene on is estimated from thousands of
# cells, not from three samples. Extract it so 45_atf3_percell.R can ask whether
# losing Atf3 blocks recruitment.
# =====================================================================
msg("extracting per-cell expression for the Atf3-WT/KO arm")
xa  <- read_counts(ATF3_GZ)
mda <- readRDS(file.path(DATA, "cell_meta_atf3.rds"))
stopifnot(identical(colnames(xa), mda$cell))
tota <- Matrix::colSums(xa)

ga <- intersect(gu$gene, rownames(xa))
msg("  ", length(ga), " of ", nrow(gu), " universe genes present in the Atf3 matrix")
Ea <- as.matrix(xa[ga, , drop = FALSE])
Ea <- round(log1p(sweep(Ea, 2, tota / 1e4, "/")), 2)
Ea <- Matrix(Ea, sparse = TRUE)
msg("  ", nrow(Ea), " x ", ncol(Ea), " | nnz ", format(length(Ea@x), big.mark = ","),
    " | ", format(object.size(Ea), units = "MB"))

saveRDS(list(expr = Ea, genes = gu[gu$gene %in% ga, ],
             meta = mda[, c("cell", "sex", "geno", "inj", "time", "rep", "celltype", "sample")]),
        file.path(DATA, "percell_expr_atf3.rds"), compress = "xz")
msg("wrote data/percell_expr_atf3.rds (",
    round(file.size(file.path(DATA, "percell_expr_atf3.rds")) / 1e6, 2), " MB)")
msg("done")
