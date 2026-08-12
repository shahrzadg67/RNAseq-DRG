# =====================================================================
# 20_build_app_artifact.R — pack the painseq single-cell results into ONE compact
# .rds for the TUY35595 Shiny app.
#
# Deliberately small: the app runs on shinyapps.io where a previous OOM was caused
# by large per-tab artifacts (see the memoised load_rds() in app/R/helpers.R). We
# ship only the summarised matrices the tab actually plots — never the 141k-cell
# count matrix or the per-cell metadata.
#
# Out: rnaseq_TUY35595/app/data/sc_injury.rds
# =====================================================================
source("/hpf/projects/msalter/sghazis/painseq/analysis/00_common.R")

tc    <- readRDS(file.path(DATA, "panel_timecourse.rds"))
fits  <- readRDS(file.path(DATA, "panel_celltype_fits.rds"))
pg    <- readRDS(file.path(DATA, "panel_genes.rds"))
peak  <- read.csv(file.path(DATA, "panel_celltype_peak.csv"), stringsAsFactors = FALSE)
cmeta <- tc$all$cmeta

rows <- pg[order(pg$role != "panel", -pg$mean_sni), ]
rows$label <- with(rows, ifelse(role == "reference", paste0(gene, " ‡"),
                        ifelse(low_detection %in% TRUE, paste0(gene, " *"), gene)))

# ---- time-course matrices (genes x injury x time) -------------------------
mat_of <- function(src, padj_src) {
  M <- matrix(NA_real_, nrow(rows), ncol(src), dimnames = list(rows$gene, colnames(src)))
  P <- M
  hit <- rows$matrix_sym %in% rownames(src)
  M[hit, ] <- src[rows$matrix_sym[hit], , drop = FALSE]
  if (!is.null(padj_src)) P[hit, ] <- padj_src[rows$matrix_sym[hit], , drop = FALSE]
  list(M = M, P = P)
}
A  <- mat_of(tc$all$l2fc,  tc$all$padj)
Bz <- mat_of(tc$all$gmean, NULL)
Bz$M <- t(scale(t(Bz$M))); Bz$M[!is.finite(Bz$M)] <- NA_real_
Bz$P <- A$P
Am <- mat_of(tc$male$l2fc, tc$male$padj)

cmeta$inj <- as.character(cmeta$inj)
cmeta$time_label <- ifelse(cmeta$time == 0, "0",
                    ifelse(cmeta$time < 168, paste0(cmeta$time, "h"), paste0(cmeta$time / 24, "d")))

# ---- cell-type matrices (genes x cell type, axotomy peak vs Naive) --------
PEAK_TIMES <- c(72L, 168L)
ct_names <- names(fits)
Mct <- matrix(NA_real_, nrow(rows), length(ct_names), dimnames = list(rows$gene, ct_names))
Pct <- Mct
for (ct in ct_names) {
  f  <- fits[[ct]]
  pc <- f$cmeta$col[f$cmeta$inj %in% AXOTOMY & f$cmeta$time %in% PEAK_TIMES]
  if (!length(pc)) next
  hit <- rows$matrix_sym %in% rownames(f$l2fc)
  Mct[hit, ct] <- rowMeans(f$l2fc[rows$matrix_sym[hit], pc, drop = FALSE], na.rm = TRUE)
  Pct[hit, ct] <- apply(f$padj[rows$matrix_sym[hit], pc, drop = FALSE], 1, min, na.rm = TRUE)
}
Pct[!is.finite(Pct)] <- NA_real_

# ---- per-cell-type time-course, for the faceted view ---------------------
ct_tc <- lapply(fits, function(f) {
  hit <- rows$matrix_sym %in% rownames(f$l2fc)
  M <- matrix(NA_real_, nrow(rows), nrow(f$cmeta), dimnames = list(rows$gene, f$cmeta$col))
  M[hit, ] <- f$l2fc[rows$matrix_sym[hit], f$cmeta$col, drop = FALSE]
  list(M = round(M, 3), cmeta = f$cmeta[, c("col", "inj", "time")], n_cells = f$n_cells)
})

# ---- headline numbers, recomputed here so the tab never hard-codes them ---
panel_syms <- rows$gene[rows$role == "panel"]
ax  <- cmeta$col[cmeta$inj %in% AXOTOMY & cmeta$time %in% c(24L, 72L, 168L)]
non <- cmeta$col[cmeta$inj %in% NON_AXOTOMY]
sumr <- function(cols) {
  V <- A$M[panel_syms, cols, drop = FALSE]; P <- A$P[panel_syms, cols, drop = FALSE]
  list(n_cols = length(cols), mean_l2fc = mean(V, na.rm = TRUE),
       pct_sig_up = 100 * mean(P < 0.05 & V > 1, na.rm = TRUE),
       n_sig_up = sum(P < 0.05 & V > 1, na.rm = TRUE))
}

# bulk cross-check, recomputed
de <- readRDS(file.path(APP, "app/data/de_long.rds"))
de <- de[de$engine == "DESeq2" & de$level == "gene" & de$contrast == "SNI_vs_Sham", ]
sp <- cmeta$col[cmeta$inj == "SpNT" & cmeta$time == 168L]
sc_v <- A$M[panel_syms, sp]; bulk <- de$log2FC[match(panel_syms, as.character(de$symbol))]
ok <- is.finite(sc_v) & is.finite(bulk)

stats <- list(
  n_cells = 141093L, n_samples = 66L, n_celltypes = 20L, n_groups = nrow(cmeta),
  panel_present = sum(rows$role == "panel" & rows$present),
  # "modelled" = survived filterByExpr, i.e. has a fitted value in a real contrast
  # (column 2; column 1 is Naive, which is 0 by construction for every gene)
  panel_modelled = sum(!is.na(A$M[panel_syms, 2])),
  axotomy = sumr(ax), non_axotomy = sumr(non),
  bulk_r = cor(sc_v[ok], bulk[ok]), bulk_rho = cor(sc_v[ok], bulk[ok], method = "spearman"),
  bulk_n = sum(ok), bulk_same_dir = mean(sign(sc_v[ok]) == sign(bulk[ok])),
  sensitivity_r = tc$sensitivity_r,
  no_baseline = setdiff(names(CELLTYPE_FULL), ct_names)
)

# =====================================================================
# Stage 2/3 artifacts — each optional, so the tab degrades gracefully if a
# stage has not been run yet.
# =====================================================================

# ---- UMAP: downsampled, stratified by cell type --------------------------
# 141,093 points would be far too heavy for plotly inside a shinyapps container
# (both memory and browser). Take ~20% stratified by cell type with a floor so the
# rare populations stay visible, and round the coordinates to 2 dp.
umap <- NULL
uf <- file.path(DATA, "umap_c57.csv")
if (file.exists(uf)) {
  u <- read.csv(uf, stringsAsFactors = FALSE)
  n_all <- nrow(u)
  set.seed(7)
  idx <- unlist(lapply(split(seq_len(nrow(u)), u$celltype), function(ii)
    if (length(ii) > 500) sample(ii, max(500, round(length(ii) * 0.2))) else ii))
  u <- u[sort(idx), ]
  umap <- data.frame(
    x = round(u$UMAP1, 2), y = round(u$UMAP2, 2),
    celltype = factor(u$celltype), inj = factor(u$inj, levels = INJ_ORDER),
    time = as.integer(u$time), score = round(u$panel_score, 3),
    stringsAsFactors = FALSE)
  msg("umap: ", nrow(umap), " of ", n_all, " cells kept (stratified ~20% with a floor for rare types)")
}

# ---- MINP cell-type concordance (Stage 3c) --------------------------------
minp_cc <- NULL
mf <- file.path(DATA, "minp_celltype_concordance.rds")
if (file.exists(mf)) {
  minp_cc <- readRDS(mf)
  minp_cc$r <- round(minp_cc$r, 4); minp_cc$rho <- round(minp_cc$rho, 4)
  msg("minp concordance: ", nrow(minp_cc), " rows")
}

# ---- Atf3 dependence (Stage 3d) -------------------------------------------
atf3 <- NULL
af <- file.path(DATA, "atf3_dependency.rds")
if (file.exists(af)) { atf3 <- readRDS(af)$table; msg("atf3 dependency: ", nrow(atf3), " genes") }

# ---- deconvolution comparison (Stage 3b) ----------------------------------
deconv <- NULL
df_ <- file.path(DATA, "deconv_compare.rds")
if (file.exists(df_)) {
  dc <- readRDS(df_)
  dc$long$prop <- round(dc$long$prop, 5)
  deconv <- list(long = dc$long, refs = dc$refs)
  msg("deconv comparison: ", nrow(dc$long), " rows")
}

# ---- Stage 4/5 summaries: recruitment, classifier, identity, Atf3, peak, MINP ----
opt <- function(f, fn = readRDS) if (file.exists(f)) fn(f) else NULL

recruit <- opt(file.path(DATA, "recruitment.rds"))
if (!is.null(recruit)) {
  # ship only the pooled (whole-DRG) view; the per-cell-type table is 24k rows
  recruit <- list(pooled = recruit$pooled[, c("gene", "group", "inj", "time", "n_cells",
                                              "f", "m", "recruitment", "level", "total")],
                  drivers = read.csv(file.path(DATA, "recruitment_drivers.csv"), stringsAsFactors = FALSE))
  recruit$pooled[, c("f", "m", "recruitment", "level", "total")] <-
    round(recruit$pooled[, c("f", "m", "recruitment", "level", "total")], 4)
  msg("recruitment: ", nrow(recruit$pooled), " rows")
}
classif <- opt(file.path(DATA, "classifier.rds"))
if (!is.null(classif)) { classif$roc$fpr <- round(classif$roc$fpr, 4); classif$roc$tpr <- round(classif$roc$tpr, 4) }
ident   <- opt(file.path(DATA, "identity.rds"));   if (!is.null(ident)) ident$cells <- NULL  # per-cell goes in sc_percell
atf3_pc <- opt(file.path(DATA, "atf3_percell.rds")); if (!is.null(atf3_pc)) atf3_pc$score <- NULL
peakviz <- opt(file.path(DATA, "peak_viz.rds"))
beyond  <- opt(file.path(DATA, "minp_beyond.rds"))
if (!is.null(beyond)) {
  beyond$sets <- NULL                       # the gene sets themselves are not plotted
  beyond$gsea$score <- round(beyond$gsea$score, 3)
  msg("minp beyond-injury: ", nrow(beyond$gsea), " rows, male-vs-female r = ",
      round(beyond$r_male_female, 3))
}
minp_sc <- opt(file.path(DATA, "minp_genes_sc.rds"))
if (!is.null(minp_sc)) minp_sc$minp <- minp_sc$minp[, c("Gene", "symbol", "named", "chr",
                                                        "MINP_l2fc", "MINP_padj", "SNI_l2fc", "SNI_padj",
                                                        "well_annotated", "MINP_up", "female_biased")]

obj <- list(
  rows = rows[, c("gene", "matrix_sym", "label", "role", "mean_sni", "n_sig", "minp_l2fc",
                  "present", "pct_cells", "total_umi", "low_detection")],
  umap = umap, minp_cc = minp_cc, atf3 = atf3, deconv = deconv,
  recruit = recruit, classif = classif, ident = ident, atf3_pc = atf3_pc,
  peakviz = peakviz, minp_sc = minp_sc, beyond = beyond,
  # the gene-explorer dropdown is populated from HERE, not from sc_percell.rds, so
  # that opening the app does not drag the per-cell matrix into memory
  gene_universe = opt(file.path(DATA, "gene_universe.rds")),
  l2fc = round(A$M, 3), padj = signif(A$P, 3),
  zscore = round(Bz$M, 3),
  l2fc_male = round(Am$M, 3), padj_male = signif(Am$P, 3),
  cmeta = cmeta[, c("col", "inj", "time", "time_label", "n_samples")],
  ct_l2fc = round(Mct, 3), ct_padj = signif(Pct, 3),
  ct_tc = ct_tc,
  celltype_full = CELLTYPE_FULL, neuronal = NEURONAL_TYPES,
  inj_full = INJ_FULL, inj_order = INJ_ORDER, axotomy = AXOTOMY, non_axotomy = NON_AXOTOMY,
  peak = peak, stats = stats,
  provenance = list(
    dataset = "GEO GSE154659 (Renthal et al. 2020, mouse DRG snRNA-seq)",
    built   = format(Sys.Date()),
    scripts = "painseq/analysis/01..05 + 20_build_app_artifact.R")
)

out <- file.path(APP, "app/data/sc_injury.rds")
saveRDS(obj, out, compress = "xz")
msg("wrote ", out, " — ", round(file.size(out) / 1e6, 2), " MB")
msg("  in-memory size: ", format(object.size(obj), units = "MB"))

# =====================================================================
# sc_percell.rds — the per-cell matrix for the interactive gene explorer.
#
# Kept in a SEPARATE file on purpose. load_rds() (app/R/helpers.R:16-29) is
# memoised per file, so a separate artifact is read only when the explorer panel
# is actually opened, instead of being pulled into memory at app start alongside
# sc_injury.rds. That matters: a previous shinyapps OOM was caused by large
# per-tab artifacts loading eagerly.
# =====================================================================
pcf <- file.path(DATA, "percell_app.rds")
if (file.exists(pcf)) {
  pcapp <- readRDS(pcf)
  # attach the per-cell derived scores computed on the FULL data
  idf <- file.path(DATA, "identity.rds")
  if (file.exists(idf) && !is.null(pcapp$idx)) {
    # 43_identity.R computes an identity score for every cell in cell_meta order;
    # pcapp$idx are the row positions of the subsample, so this aligns exactly.
    # (Joining on cell NAME would be wrong — 5,593 GEO names collide.)
    md_full  <- readRDS(file.path(DATA, "cell_meta.rds"))
    ident_all <- rep(NA_real_, nrow(md_full))
    ic <- readRDS(idf)$cells
    ident_all[ic$row] <- ic$identity
    pcapp$meta$identity <- round(ident_all[pcapp$idx], 3)
    msg("  identity attached for ", sum(!is.na(pcapp$meta$identity)), " of ",
        nrow(pcapp$meta), " subsampled cells")
  }
  cl <- opt(file.path(DATA, "classifier.rds"))
  if (!is.null(cl)) pcapp$meta$injured <- pcapp$meta$score > cl$threshold
  outp <- file.path(APP, "app/data/sc_percell.rds")
  saveRDS(pcapp, outp, compress = "xz")
  sz_disk <- file.size(outp) / 1e6; sz_mem <- as.numeric(object.size(pcapp)) / 1e6
  msg("wrote ", outp, " — ", round(sz_disk, 2), " MB on disk, ",
      round(sz_mem, 1), " MB resident (", nrow(pcapp$meta), " cells x ",
      nrow(pcapp$expr), " genes)")
  # budget assertions from the plan
  if (sz_disk > 5)  stop("sc_percell.rds exceeds the 5 MB disk budget: ", round(sz_disk, 2))
  if (sz_mem  > 15) stop("sc_percell.rds exceeds the 15 MB resident budget: ", round(sz_mem, 1))
}

tot <- sum(file.size(list.files(file.path(APP, "app/data"), pattern = "\\.rds$", full.names = TRUE)))
msg("app/data/*.rds total: ", round(tot / 1e6, 1), " MB")
msg("  axotomy mean log2FC ", round(stats$axotomy$mean_l2fc, 2),
    " (", round(stats$axotomy$pct_sig_up, 1), "% sig-up) vs non-axotomy ",
    round(stats$non_axotomy$mean_l2fc, 2), " (", stats$non_axotomy$n_sig_up, " genes sig-up)")
msg("  bulk cross-check r = ", round(stats$bulk_r, 3), " over n = ", stats$bulk_n)
