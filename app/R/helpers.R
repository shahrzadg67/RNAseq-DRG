# Shared helpers for the app -------------------------------------------------
`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

# Read the build config written by publish.sh (target + per-stage status).
# Falls back to a staging build showing everything if absent (e.g. local dev).
read_build_config <- function() {
  f <- "data/build_config.json"
  if (file.exists(f)) {
    cfg <- jsonlite::fromJSON(f, simplifyVector = FALSE)
  } else {
    cfg <- list(target = "staging", project = list(), stages = list())
  }
  cfg
}

# Load an RDS artifact if present, else NULL (lets tabs show "coming soon").
# Memoised loader: each .rds is read from disk at most once per process and the
# object is shared (read-only, copy-on-write) across every module and session.
# This prevents heavy files (e.g. de_long.rds ~195 MB) from being duplicated in
# memory once per tab that loads them — the cause of shinyapps OOM crashes.
# Modules only ever subset/copy the returned object, never mutate it in place.
.rds_cache <- new.env(parent = emptyenv())
load_rds <- function(name) {
  if (!is.null(.rds_cache[[name]])) return(.rds_cache[[name]])
  f <- file.path("data", paste0(name, ".rds"))
  obj <- if (file.exists(f)) readRDS(f) else NULL
  if (!is.null(obj)) assign(name, obj, envir = .rds_cache)   # cache real objects only
  obj
}

# Standard "not released / not computed yet" placeholder panel.
coming_soon <- function(msg = "This section is being prepared and will appear here once released.") {
  div(class = "banner", icon("hourglass-half"), tags$b(" Coming soon. "), msg)
}

status_chip <- function(status) {
  span(class = paste0("status-chip s-", status), toupper(gsub("_", " ", status)))
}

# Gene/transcript annotation (symbol, full description, chromosome) built by
# analysis/60_annotate.R. Loaded once; used to enrich hover tooltips + tables.
ANNOT <- load_rds("gene_annot")

# Attach `description` and `chr` columns to a data.frame `d` that has feature_id,
# matched within the given level ("gene"/"transcript"). Safe if ANNOT is missing.
attach_annot <- function(d, level) {
  if (is.null(ANNOT) || !nrow(d)) { d$description <- NA_character_; d$chr <- NA_character_; return(d) }
  a <- ANNOT[ANNOT$level == level, ]
  i <- match(d$feature_id, a$feature_id)
  d$description <- a$description[i]
  d$chr <- a$chr[i]
  d
}

# Look up a single feature's (description, chr); returns a 1-row list.
annot_of <- function(feature_id, level) {
  if (is.null(ANNOT)) return(list(description = NA, chr = NA))
  a <- ANNOT[ANNOT$level == level & ANNOT$feature_id == feature_id, ]
  if (!nrow(a)) return(list(description = NA, chr = NA))
  list(description = a$description[1], chr = a$chr[1])
}

# ---- Lab figure conventions -------------------------------------------------
# NP is displayed everywhere as "MINP" (data keys stay NP internally).
# Sex is the primary split: Male = blue, Female = red.
# Condition: Sham = black, MINP = orange (+ SNI grey for the original with-SNI views).
# Shapes: Sham = circle, MINP = triangle, SNI = square.
SEX_COLORS  <- c(F = "#941200", M = "#021893")
COND_COLORS <- c(Sham = "#000000", NP = "#E67E22", SNI = "#7F8C8D")
COND_SHAPES <- c(Sham = "circle", NP = "triangle-up", SNI = "square")
GROUP_SEX   <- c(NP_F = "F", NP_M = "M", Sham_F = "F", Sham_M = "M", SNI_F = "F", SNI_M = "M")

# Display relabel: NP -> MINP (never double-applies; internal strings never contain MINP)
disp <- function(x) gsub("NP", "MINP", as.character(x))

# Ordered levels + matching colours for a colour-by variable ("sex"/"condition"/"group")
color_spec <- function(vals, by) {
  vals <- as.character(vals)
  if (by == "sex") { lv <- intersect(c("F","M"), vals); return(list(levels = lv, colors = unname(SEX_COLORS[lv]))) }
  if (by == "condition") { lv <- intersect(names(COND_COLORS), vals); return(list(levels = lv, colors = unname(COND_COLORS[lv]))) }
  lv <- intersect(names(GROUP_SEX), vals); list(levels = lv, colors = unname(SEX_COLORS[GROUP_SEX[lv]]))  # group -> sex colour
}
# Ordered levels + matching plotly symbols, keyed by condition
shape_spec <- function(vals) { lv <- intersect(names(COND_SHAPES), as.character(vals)); list(levels = lv, symbols = unname(COND_SHAPES[lv])) }

# ---- publication styling ----------------------------------------------------
# ggplot theme: rectangular border, no gridlines, larger fonts (used by the
# static volcano, Venn, heatmap, and analysis/93_pub_figs.R).
theme_pub <- function(base_size = 15) {
  ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      panel.border = ggplot2::element_rect(colour = "black", fill = NA, linewidth = 0.9),
      panel.grid   = ggplot2::element_blank(),
      axis.ticks   = ggplot2::element_line(colour = "black"),
      axis.text    = ggplot2::element_text(colour = "black"),
      plot.title   = ggplot2::element_text(face = "bold"),
      plot.subtitle= ggplot2::element_text(colour = "grey30"),
      legend.key   = ggplot2::element_blank())
}

# plotly axis styling to match: no grid, black mirrored border, bigger font.
# Uses a second layout() call — plotly merges these axis attrs with any titles
# already set by the plot (so titles are preserved).
ply_pub <- function(p, fontsize = 14) {
  ax <- list(showgrid = FALSE, zeroline = FALSE, showline = TRUE, mirror = TRUE,
             linecolor = "black", linewidth = 1.3, ticks = "outside", tickcolor = "black")
  plotly::layout(p, font = list(size = fontsize), xaxis = ax, yaxis = ax)
}

# Shared PCA scatter: one plotly trace per level of `grp`, each with a DISTINCT
# colour AND shape, so every group is separable on two channels and appears as a
# legend entry. df needs numeric columns xcol/ycol and a grouping column grpcol.
PCA_SHAPES <- c("circle","square","triangle-up","diamond","triangle-down","cross","x",
                "star","pentagon","hexagon","triangle-left","triangle-right","star-square","hourglass")
pca_scatter <- function(df, xcol, ycol, grpcol, title = "", xlab = "", ylab = "",
                        hovertext = NULL, palette = NULL, legend_title = "") {
  g <- factor(df[[grpcol]]); lv <- levels(g); n <- length(lv)
  cols <- if (!is.null(palette) && all(lv %in% names(palette))) palette
          else stats::setNames(scales::hue_pal()(n), lv)
  shp  <- stats::setNames(PCA_SHAPES[((seq_len(n) - 1) %% length(PCA_SHAPES)) + 1], lv)
  if (is.null(hovertext)) hovertext <- as.character(g)
  p <- plotly::plot_ly()
  for (l in lv) { idx <- which(g == l); if (!length(idx)) next
    p <- plotly::add_trace(p, x = df[[xcol]][idx], y = df[[ycol]][idx], type = "scatter", mode = "markers",
          name = l, marker = list(size = 11, color = unname(cols[l]), symbol = unname(shp[l]),
                                  line = list(width = .5, color = "#fff")),
          text = hovertext[idx], hoverinfo = "text") }
  plotly::layout(p, title = title, xaxis = list(title = xlab), yaxis = list(title = ylab),
                 legend = list(title = list(text = legend_title))) |> ply_pub()
}

# Hierarchical-clustering order for a heatmap matrix; returns row/col index orders.
# Falls back to identity order on any failure or when too few rows/cols.
cluster_order <- function(m, do_row = TRUE, do_col = FALSE) {
  ro <- seq_len(nrow(m)); co <- seq_len(ncol(m))
  safe <- function(x) { x[!is.finite(x)] <- 0; x }
  if (isTRUE(do_row) && nrow(m) > 2)
    ro <- tryCatch(stats::hclust(stats::dist(safe(m)))$order, error = function(e) ro)
  if (isTRUE(do_col) && ncol(m) > 2)
    co <- tryCatch(stats::hclust(stats::dist(t(safe(m))))$order, error = function(e) co)
  list(row = ro, col = co)
}

# Publication volcano (DESeq2): points + top-N gene-name labels (ggrepel),
# dotted cutoff lines, adjustable axes. `d` needs log2FC, padj, symbol, feature_id.
volcano_gg_app <- function(d, padj_cut, lfc_cut, nlab, xmax, ymax, title) {
  d <- d[!is.na(d$padj), ]
  d$sig <- d$padj < padj_cut & abs(d$log2FC) >= lfc_cut
  d$lab <- ifelse(is.na(d$symbol), as.character(d$feature_id), as.character(d$symbol))  # factor-safe
  ds <- d[d$sig, ]; ds <- ds[order(ds$padj), ]; lab_d <- utils::head(ds, nlab)
  ggplot2::ggplot(d, ggplot2::aes(log2FC, -log10(padj))) +
    ggplot2::geom_point(ggplot2::aes(colour = sig), size = 1, alpha = 0.55) +
    ggplot2::scale_colour_manual(values = c("FALSE" = "#B0B8BF", "TRUE" = "#E67E22"), guide = "none") +
    ggplot2::geom_vline(xintercept = c(-lfc_cut, lfc_cut), linetype = "dotted", colour = "grey40") +
    ggplot2::geom_hline(yintercept = -log10(padj_cut), linetype = "dotted", colour = "grey40") +
    ggrepel::geom_text_repel(data = lab_d, ggplot2::aes(label = lab), size = 4.2, colour = "black",
                             max.overlaps = Inf, min.segment.length = 0, box.padding = 0.45, seed = 1) +
    ggplot2::coord_cartesian(xlim = c(-xmax, xmax), ylim = c(0, ymax)) +
    ggplot2::labs(title = title,
                  subtitle = sprintf("%d significant (adj.p<%.3g, |log2FC|>=%g)", sum(d$sig), padj_cut, lfc_cut),
                  x = "log2 fold-change", y = "-log10 adj.p") +
    theme_pub()
}

# open a graphics device by format for server-side figure downloads
gg_device <- function(file, fmt, w = 8, h = 6) {
  if (fmt == "png")      ragg::agg_png(file, width = w, height = h, units = "in", res = 300)
  else if (fmt == "svg") grDevices::svg(file, width = w, height = h)
  else                   grDevices::cairo_pdf(file, width = w, height = h)
}

# Download links (PNG/SVG/PDF) for a pre-rendered figure stem under www/figs/.
# The files are produced by analysis/93_pub_figs.R (300 dpi PNG + vector SVG/PDF).
fig_dl <- function(stem, label = "Download figure") {
  base <- file.path("figs", stem)
  fmt <- function(ext) tags$a(href = paste0(base, ".", ext), download = NA,
                              class = "btn btn-outline-secondary btn-sm py-0", toupper(ext))
  div(class = "mt-2", tags$small(class = "text-muted", label, ": "),
      fmt("png"), " ", fmt("svg"), " ", fmt("pdf"))
}

# which groups each contrast actually compares (for the heatmap sample scope)
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

# pretty contrast labels (names = internal keys; display uses MINP)
CONTRAST_LABELS <- c(
  NP_vs_Sham        = "MINP (M+F) vs Sham (M+F)",
  NP_M_vs_Sham_M    = "MINP male vs Sham male",
  NP_F_vs_Sham_F    = "MINP female vs Sham female",
  Sham_M_vs_Sham_F  = "Sham male vs Sham female",
  NP_M_vs_NP_F      = "MINP male vs MINP female",
  SNI_vs_Sham       = "SNI (M+F) vs Sham (M+F)",
  SNI_vs_NP         = "SNI (M+F) vs MINP (M+F)",
  SNI_vs_Sham_M     = "SNI (M+F) vs Sham male",
  SNI_vs_Sham_F     = "SNI (M+F) vs Sham female",
  SNI_vs_NP_M       = "SNI (M+F) vs MINP male",
  SNI_vs_NP_F       = "SNI (M+F) vs MINP female"
)

# Rat cross-species contrasts (rat DRG, GRCr8; Naive = baseline). 12 types × DRG
# level {U,L} = 24 contrasts (upper vs lower DRG kept separate). Keys end in _U/_L.
.RAT_TYPES <- c(
  SNI_M_vs_SNI_F     = "SNI male vs SNI female",
  SL_M_vs_SL_F       = "Skin-lesion male vs Skin-lesion female",
  Naive_M_vs_Naive_F = "Naïve male vs Naïve female",
  SNI_M_vs_Naive_M   = "SNI male vs Naïve male",
  SNI_F_vs_Naive_F   = "SNI female vs Naïve female",
  SL_M_vs_Naive_M    = "Skin-lesion male vs Naïve male",
  SL_F_vs_Naive_F    = "Skin-lesion female vs Naïve female",
  SNI_M_vs_SL_M      = "SNI male vs Skin-lesion male",
  SNI_F_vs_SL_F      = "SNI female vs Skin-lesion female",
  SNI_vs_Naive       = "SNI (M+F) vs Naïve (M+F)",
  SNI_vs_SL          = "SNI (M+F) vs Skin-lesion (M+F)",
  SL_vs_Naive        = "Skin-lesion (M+F) vs Naïve (M+F)")
.RAT_TYPE_AB <- list(   # A/B condition_sex component(s) per type (level appended below)
  SNI_M_vs_SNI_F=list(A="SNI_M",B="SNI_F"), SL_M_vs_SL_F=list(A="SL_M",B="SL_F"),
  Naive_M_vs_Naive_F=list(A="Naive_M",B="Naive_F"), SNI_M_vs_Naive_M=list(A="SNI_M",B="Naive_M"),
  SNI_F_vs_Naive_F=list(A="SNI_F",B="Naive_F"), SL_M_vs_Naive_M=list(A="SL_M",B="Naive_M"),
  SL_F_vs_Naive_F=list(A="SL_F",B="Naive_F"), SNI_M_vs_SL_M=list(A="SNI_M",B="SL_M"),
  SNI_F_vs_SL_F=list(A="SNI_F",B="SL_F"),
  SNI_vs_Naive=list(A=c("SNI_M","SNI_F"),B=c("Naive_M","Naive_F")),
  SNI_vs_SL=list(A=c("SNI_M","SNI_F"),B=c("SL_M","SL_F")),
  SL_vs_Naive=list(A=c("SL_M","SL_F"),B=c("Naive_M","Naive_F")))
RAT_CONTRAST_LABELS <- unlist(lapply(c("U","L"), function(lv)
  setNames(paste0(.RAT_TYPES, " (", lv, ")"), paste0(names(.RAT_TYPES), "_", lv))))
RAT_CONTRAST_GROUPS <- do.call(c, lapply(c("U","L"), function(lv)
  setNames(lapply(.RAT_TYPE_AB, function(ab) list(A=paste0(ab$A,"_",lv), B=paste0(ab$B,"_",lv))),
           paste0(names(.RAT_TYPE_AB), "_", lv))))
# A/B sample-id sets for a rat contrast key (meta must have $group = condition_sex_level, $sample)
rat_sides <- function(cname, meta) {
  g <- RAT_CONTRAST_GROUPS[[cname]]
  list(A = meta$sample[meta$group %in% g$A], B = meta$sample[meta$group %in% g$B])
}

# Monkey cross-species contrasts (macaque DRG, Macaca_fascicularis_6.0; transection
# axotomy: ipsi = injured L4-6 side, contra = contralateral, control = naïve). 12
# contrasts, condition_sex groups (no DRG-level split). control = baseline.
MONKEY_CONTRAST_LABELS <- c(
  ipsi_M_vs_ipsi_F       = "Ipsi male vs Ipsi female",
  contra_M_vs_contra_F   = "Contra male vs Contra female",
  control_M_vs_control_F = "Control male vs Control female",
  ipsi_M_vs_control_M    = "Ipsi male vs Control male",
  ipsi_F_vs_control_F    = "Ipsi female vs Control female",
  contra_M_vs_control_M  = "Contra male vs Control male",
  contra_F_vs_control_F  = "Contra female vs Control female",
  ipsi_M_vs_contra_M     = "Ipsi male vs Contra male",
  ipsi_F_vs_contra_F     = "Ipsi female vs Contra female",
  ipsi_vs_control        = "Ipsi (M+F) vs Control (M+F)",
  ipsi_vs_contra         = "Ipsi (M+F) vs Contra (M+F)",
  contra_vs_control      = "Contra (M+F) vs Control (M+F)")
MONKEY_CONTRAST_GROUPS <- list(
  ipsi_M_vs_ipsi_F=list(A="ipsi_M",B="ipsi_F"), contra_M_vs_contra_F=list(A="contra_M",B="contra_F"),
  control_M_vs_control_F=list(A="control_M",B="control_F"), ipsi_M_vs_control_M=list(A="ipsi_M",B="control_M"),
  ipsi_F_vs_control_F=list(A="ipsi_F",B="control_F"), contra_M_vs_control_M=list(A="contra_M",B="control_M"),
  contra_F_vs_control_F=list(A="contra_F",B="control_F"), ipsi_M_vs_contra_M=list(A="ipsi_M",B="contra_M"),
  ipsi_F_vs_contra_F=list(A="ipsi_F",B="contra_F"),
  ipsi_vs_control=list(A=c("ipsi_M","ipsi_F"),B=c("control_M","control_F")),
  ipsi_vs_contra=list(A=c("ipsi_M","ipsi_F"),B=c("contra_M","contra_F")),
  contra_vs_control=list(A=c("contra_M","contra_F"),B=c("control_M","control_F")))
# A/B sample-id sets for a monkey contrast key (meta must have $group = condition_sex, $sample)
monkey_sides <- function(cname, meta) {
  g <- MONKEY_CONTRAST_GROUPS[[cname]]
  list(A = meta$sample[meta$group %in% g$A], B = meta$sample[meta$group %in% g$B])
}

# SNI contrasts are exploratory (SNI is n=1/sex, 2 samples) — flag them in the UI.
is_sni <- function(cname) grepl("^SNI", cname %||% "")
sni_note <- function(cname) if (isTRUE(is_sni(cname)))
  div(class = "banner", style = "border-left-color:#e74c3c; background:#fdecea;",
      tags$b("⚠ SNI is n = 1 per sex (2 samples). "),
      "These contrasts are ", tags$b("exploratory"), ": fold-changes are shown, but ",
      "p-values / FDR are unreliable and confounded with sex.") else NULL

# Full names for the deconvolution cell types (Renthal et al. 2020 DRG atlas
# nomenclature). M/N are the two Schwann-cell subtypes (confirmed by the lab).
CELLTYPE_FULL <- c(
  NP          = "Non-peptidergic nociceptor",
  PEP1        = "Peptidergic nociceptor 1",
  PEP2        = "Peptidergic nociceptor 2",
  NF1         = "Myelinated A-fibre LTMR / neurofilament 1",
  NF2         = "Myelinated A-fibre LTMR / neurofilament 2",
  NF3         = "Proprioceptor / neurofilament 3",
  cLTMR1      = "C-fibre low-threshold mechanoreceptor 1",
  cLTMR2      = "C-fibre low-threshold mechanoreceptor 2",
  SST         = "Somatostatin+ pruriceptor",
  Satglia     = "Satellite glia",
  M           = "Myelinating Schwann cell",
  N           = "Non-myelinating Schwann cell",
  Endothelial = "Endothelial cell",
  Pericyte    = "Pericyte",
  Fibroblast  = "Fibroblast",
  Macrophage  = "Macrophage",
  `B cell`    = "B cell",
  Neutrophil  = "Neutrophil"
)
# Sensory-neuron subtypes in the deconvolution reference (everything else is
# non-neuronal: satellite glia, Schwann, vascular, fibroblast, immune).
NEURONAL_TYPES <- c("NP","PEP1","PEP2","NF1","NF2","NF3","cLTMR1","cLTMR2","SST")

# "Full name (ABBR)" for a cell-type abbreviation; unknown or self-same -> the abbreviation
ct_full <- function(abbr) {
  a <- as.character(abbr); f <- unname(CELLTYPE_FULL[a])
  ifelse(is.na(f) | f == a, a, sprintf("%s (%s)", f, a))
}
