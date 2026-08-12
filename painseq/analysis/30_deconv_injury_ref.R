# =====================================================================
# 30_deconv_injury_ref.R — Stage 3b: does an INJURY-STATE deconvolution reference
# change the answer?
#
# The app's existing deconvolution (rnaseq_TUY35595/analysis/80_deconv.R) builds its
# MuSiC reference from NAIVE cells only. That is a defensible choice — naive
# signatures are not injury-reprogrammed — but it has a hard consequence: the two
# injury-induced populations, `Repair schwann_N` and `Repair fibroblast`, do not
# exist in naive tissue, so a naive-only reference is structurally incapable of
# reporting them, no matter how much of the bulk they explain.
#
# Here we build both references from the SAME atlas and deconvolve the same bulk:
#   naive_fixed : naive cells only, with the corrected cell-type labels
#   injury      : naive + injured cells, all 20 types incl. the Repair populations
#
# It also fixes the label bug inherited by the app: 80_deconv.R:60 takes the
# second-to-last underscore field as the cell type, turning p_cLTMR2 -> cLTMR2,
# Schwann_M -> M, Schwann_N -> N and `Repair schwann_N` -> schwann.
#
# Out: data/deconv_compare.rds, figs/deconv_naive_vs_injury.{png,pdf}
# =====================================================================
suppressMessages({ library(Matrix); library(SingleCellExperiment); library(MuSiC); library(ggplot2) })
source("/hpf/projects/msalter/sghazis/painseq/analysis/00_common.R")

MIN_CELLS_PER_TYPE <- 20
MAX_CELLS_PER_TYPE <- 1500   # cap per cell type: MuSiC cost scales with cells, and
                             # the injury reference would otherwise be 8x the naive one

BULKS <- list(
  mouse  = file.path(APP, "results/star_salmon/salmon.merged.gene_counts.tsv"),
  rat    = "/hpf/projects/msalter/sghazis/rnaseq_rat_DRG/results/star_salmon/salmon.merged.gene_counts.tsv",
  monkey = "/hpf/projects/msalter/sghazis/rnaseq_monkey_DRG/results/star_salmon/salmon.merged.gene_counts.tsv"
)
ORTHO <- list(
  rat    = list(f = "ortholog_rat_mouse.rds",    id = "rat_gene_id"),
  monkey = list(f = "ortholog_monkey_mouse.rds", id = "monkey_gene_id")
)

# ---------------------------------------------------------------------
# Reference construction
# ---------------------------------------------------------------------
build_ref <- function(counts, md, which = c("naive", "injury")) {
  which <- match.arg(which)
  keep <- if (which == "naive") md$inj == "Naive" else rep(TRUE, nrow(md))
  cnt <- counts[, keep, drop = FALSE]; m <- md[keep, ]

  ct_n <- table(m$celltype)
  drop <- names(ct_n)[ct_n < MIN_CELLS_PER_TYPE]
  if (length(drop)) {
    msg("  dropping rare types (<", MIN_CELLS_PER_TYPE, "): ",
        paste(sprintf("%s(%d)", drop, ct_n[drop]), collapse = ", "))
    k <- !m$celltype %in% drop; cnt <- cnt[, k, drop = FALSE]; m <- m[k, ]
  }
  # subsample within cell type for tractability, stratified so every subject stays represented
  set.seed(42)
  idx <- unlist(lapply(split(seq_len(nrow(m)), m$celltype), function(ii)
    if (length(ii) > MAX_CELLS_PER_TYPE) sample(ii, MAX_CELLS_PER_TYPE) else ii))
  idx <- sort(idx)
  cnt <- cnt[, idx, drop = FALSE]; m <- m[idx, ]

  if (any(duplicated(rownames(cnt)))) {
    ord <- order(Matrix::rowSums(cnt), decreasing = TRUE)
    cnt <- cnt[ord, , drop = FALSE]; cnt <- cnt[!duplicated(rownames(cnt)), , drop = FALSE]
  }
  # MuSiC "subject" = biological replicate. For the injury reference the subject
  # must include the injury x time group, else cells from different conditions are
  # pooled into one pseudo-subject and the cross-subject variance is meaningless.
  subj <- if (which == "naive") paste(m$sex, m$rep, sep = "_")
          else paste(m$sex, m$inj, m$time, m$rep, sep = "_")
  msg("  ", which, " reference: ", ncol(cnt), " cells, ",
      length(unique(m$celltype)), " types, ", length(unique(subj)), " subjects")
  list(sce = SingleCellExperiment(assays = list(counts = cnt),
         colData = DataFrame(cellType = m$celltype, sampleID = subj, sex = m$sex)),
       celltypes = sort(unique(m$celltype)))
}

# ---------------------------------------------------------------------
# Bulk loading, mapped into mouse SYMBOL space
# ---------------------------------------------------------------------
load_bulk <- function(path, species) {
  bt  <- read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
  idc <- intersect(c("gene_id", "gene_name"), colnames(bt))
  mat <- as.matrix(bt[, !(colnames(bt) %in% idc), drop = FALSE])
  storage.mode(mat) <- "double"
  if (species == "mouse") {
    sym <- bt$gene_name
  } else {
    o <- readRDS(file.path(APP, "app/data", ORTHO[[species]]$f))
    sym <- as.character(o$mouse_symbol[match(bt$gene_id, o[[ORTHO[[species]]$id]])])
  }
  ok <- !is.na(sym) & sym != ""
  b  <- rowsum(mat[ok, , drop = FALSE], group = sym[ok])
  storage.mode(b) <- "integer"
  msg("  bulk ", species, ": ", nrow(b), " mouse symbols x ", ncol(b), " samples")
  b
}

run_music <- function(bulk, ref) {
  common <- intersect(rownames(bulk), rownames(ref$sce))
  stopifnot(length(common) > 5000)
  est <- music_prop(bulk.mtx = bulk, sc.sce = ref$sce, clusters = "cellType",
                    samples = "sampleID", select.ct = ref$celltypes, verbose = FALSE)
  p <- est$Est.prop.weighted
  p[colnames(bulk), , drop = FALSE]
}

# ---------------------------------------------------------------------
msg("reading atlas")
counts <- read_counts(C57_RDS)
md     <- readRDS(file.path(DATA, "cell_meta.rds"))
stopifnot(identical(colnames(counts), md$cell))

refs <- list(naive_fixed = build_ref(counts, md, "naive"),
             injury      = build_ref(counts, md, "injury"))
rm(counts); invisible(gc())

res <- list()
for (sp in names(BULKS)) {
  if (!file.exists(BULKS[[sp]])) { msg("SKIP ", sp, " — no bulk matrix"); next }
  msg("=== ", sp, " ===")
  b <- load_bulk(BULKS[[sp]], sp)
  for (rn in names(refs)) {
    msg("  MuSiC: ", sp, " x ", rn, " reference ...")
    p <- tryCatch(run_music(b, refs[[rn]]),
                  error = function(e) { msg("    FAILED: ", conditionMessage(e)); NULL })
    if (!is.null(p)) res[[paste(sp, rn, sep = "|")]] <- p
  }
}

# ---------------------------------------------------------------------
# Compare: what does the injury reference reveal that naive-only cannot?
# ---------------------------------------------------------------------
long <- do.call(rbind, lapply(names(res), function(k) {
  sp <- sub("\\|.*", "", k); rf <- sub(".*\\|", "", k); p <- res[[k]]
  data.frame(species = sp, reference = rf, sample = rep(rownames(p), ncol(p)),
             celltype = rep(colnames(p), each = nrow(p)),
             prop = as.vector(p), stringsAsFactors = FALSE)
}))
saveRDS(list(prop = res, long = long,
             refs = lapply(refs, function(r) r$celltypes)),
        file.path(DATA, "deconv_compare.rds"))

msg("\n=== mean proportion by reference (mouse) ===")
mm <- long[long$species == "mouse", ]
if (nrow(mm)) {
  agg <- aggregate(prop ~ celltype + reference, data = mm, FUN = mean)
  w <- reshape(agg, idvar = "celltype", timevar = "reference", direction = "wide")
  names(w) <- sub("^prop\\.", "", names(w))
  w[is.na(w)] <- 0
  w <- w[order(-w$injury), ]
  print(data.frame(celltype = w$celltype,
                   naive_only = round(100 * w$naive_fixed, 2),
                   injury_ref = round(100 * w$injury, 2),
                   diff_pp = round(100 * (w$injury - w$naive_fixed), 2)), row.names = FALSE)
  inj_only <- setdiff(refs$injury$celltypes, refs$naive_fixed$celltypes)
  msg("\ncell types present ONLY in the injury reference: ",
      if (length(inj_only)) paste(inj_only, collapse = ", ") else "(none)")
  if (length(inj_only)) {
    sub <- mm[mm$reference == "injury" & mm$celltype %in% inj_only, ]
    a2 <- aggregate(prop ~ celltype, data = sub, FUN = function(v) c(mean = mean(v), max = max(v)))
    for (i in seq_len(nrow(a2)))
      msg(sprintf("  %-20s mean %.2f%% of bulk, max %.2f%% in any sample",
                  a2$celltype[i], 100 * a2$prop[i, "mean"], 100 * a2$prop[i, "max"]))
  }
}

# ---------------------------------------------------------------------
gg <- ggplot(long, aes(x = celltype, y = 100 * prop, fill = reference)) +
  geom_boxplot(outlier.size = .5, linewidth = .3, position = position_dodge(preserve = "single")) +
  facet_wrap(~ species, ncol = 1, scales = "free_y") +
  scale_fill_manual(values = c(naive_fixed = "#4E79A7", injury = "#E15759"), name = "Reference") +
  labs(title = "Deconvolution: naive-only vs injury-state reference",
       subtitle = paste(strwrap(paste(
         "Both references come from GSE154659. The naive-only reference cannot report",
         "`Repair schwann_N` or `Repair fibroblast` at all, because those populations do not",
         "exist in naive tissue. Cell-type labels are the corrected ones."), width = 120), collapse = "\n"),
       x = NULL, y = "Estimated proportion (%)") +
  theme_pub(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        plot.subtitle = element_text(colour = "grey30", size = 8.5),
        strip.background = element_rect(fill = "grey93", colour = "black"),
        strip.text = element_text(face = "bold"))
ggsave(file.path(FIGS, "deconv_naive_vs_injury.png"), gg, width = 12, height = 10, dpi = 300, bg = "white")
ggsave(file.path(FIGS, "deconv_naive_vs_injury.pdf"), gg, width = 12, height = 10, device = cairo_pdf)
msg("wrote figs/deconv_naive_vs_injury.{png,pdf} and data/deconv_compare.rds")
