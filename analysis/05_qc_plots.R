# =====================================================================
# 05_qc_plots.R — count-level QC plots for the app's QC tab.
# Reads the pipeline gene counts, builds an edgeR DGEList, and renders
# three 300-dpi PNGs into app/www/qc/ that tab_qc.R embeds:
#   libsize.png       — per-sample library size
#   density.png       — log-CPM density, before vs after filterByExpr
#   tmm_boxplots.png  — TMM-normalised log-CPM boxplots
# Lab palette (sex blue/red, condition), NP shown as MINP, pub theme.
# =====================================================================
suppressMessages({ library(edgeR); library(ggplot2) })

PROJ <- "/hpf/projects/msalter/sghazis/rnaseq_TUY35595"
BULK <- file.path(PROJ, "results/star_salmon/salmon.merged.gene_counts.tsv")
OUT  <- file.path(PROJ, "app/www/qc"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

# ---- lab conventions (kept in sync with app/R/helpers.R) ----
SEX_COLORS <- c(F = "#941200", M = "#021893")
disp <- function(x) gsub("NP", "MINP", as.character(x))
theme_pub <- function(base_size = 14) {
  theme_bw(base_size = base_size) +
    theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.9),
          panel.grid = element_blank(), axis.ticks = element_line(colour = "black"),
          axis.text = element_text(colour = "black"), plot.title = element_text(face = "bold"),
          legend.key = element_blank())
}
GRP_ORDER <- c("NP_F","NP_M","Sham_F","Sham_M","SNI_F","SNI_M")

# ---- counts + sample metadata ----
bt  <- read.delim(BULK, check.names = FALSE, stringsAsFactors = FALSE)
mat <- round(as.matrix(bt[, !(colnames(bt) %in% c("gene_id","gene_name")), drop = FALSE]))
rownames(mat) <- bt$gene_id
samp <- colnames(mat)
cond <- sub("_.*$", "", samp)
sx   <- sub("^.*_([FM]).*$", "\\1", samp)
grp  <- paste0(cond, "_", sx)
ord  <- order(match(grp, GRP_ORDER), samp)
samp <- samp[ord]; cond <- cond[ord]; sx <- sx[ord]; grp <- grp[ord]; mat <- mat[, ord, drop = FALSE]
lab  <- factor(disp(samp), levels = disp(samp))                 # x-axis labels, MINP display, group-ordered

dge  <- DGEList(counts = mat)
keep <- filterByExpr(dge, group = grp)

save_png <- function(file, p, w = 10, h = 6) {
  ragg::agg_png(file, width = w, height = h, units = "in", res = 300); print(p); grDevices::dev.off()
  cat("wrote", file, "\n")
}

# ---- 1. library sizes ----
d1 <- data.frame(sample = lab, sex = sx, condition = disp(cond), libM = colSums(mat) / 1e6)
p1 <- ggplot(d1, aes(sample, libM, fill = sex)) +
  geom_col(width = 0.8, colour = "black", linewidth = 0.25) +
  scale_fill_manual(values = SEX_COLORS, name = "Sex") +
  labs(title = "Library size per sample", subtitle = "Total assigned gene counts",
       x = NULL, y = "Library size (millions of reads)") +
  theme_pub() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_png(file.path(OUT, "libsize.png"), p1)

# ---- 2. log-CPM density, before vs after filtering ----
lcpm_all  <- cpm(dge, log = TRUE)
lcpm_keep <- cpm(dge[keep, , keep.lib.sizes = FALSE], log = TRUE)
mk_long <- function(m, panel) {
  data.frame(logCPM = as.vector(m),
             sample = rep(disp(samp), each = nrow(m)),
             sex    = rep(sx, each = nrow(m)),
             panel  = panel, stringsAsFactors = FALSE)
}
d2 <- rbind(mk_long(lcpm_all,  sprintf("Before filtering (%d genes)", nrow(lcpm_all))),
            mk_long(lcpm_keep, sprintf("After filterByExpr (%d genes)", nrow(lcpm_keep))))
d2$panel <- factor(d2$panel, levels = unique(d2$panel))
p2 <- ggplot(d2, aes(logCPM, group = sample, colour = sex)) +
  geom_density(linewidth = 0.5, alpha = 0.8) +
  facet_wrap(~ panel, nrow = 1) +
  scale_colour_manual(values = SEX_COLORS, name = "Sex") +
  labs(title = "Expression distribution (log-CPM)",
       subtitle = "One curve per sample; low-count genes removed by edgeR::filterByExpr",
       x = "log2 CPM", y = "Density") +
  theme_pub() + theme(strip.background = element_rect(fill = "grey92", colour = "black"))
save_png(file.path(OUT, "density.png"), p2, w = 11, h = 5)

# ---- 3. TMM-normalised log-CPM boxplots ----
dge2 <- calcNormFactors(dge[keep, , keep.lib.sizes = FALSE], method = "TMM")
lcpm <- cpm(dge2, log = TRUE)
d3 <- data.frame(logCPM = as.vector(lcpm),
                 sample = factor(rep(disp(samp), each = nrow(lcpm)), levels = levels(lab)),
                 sex    = rep(sx, each = nrow(lcpm)))
p3 <- ggplot(d3, aes(sample, logCPM, fill = sex)) +
  geom_boxplot(outlier.size = 0.3, linewidth = 0.3) +
  scale_fill_manual(values = SEX_COLORS, name = "Sex") +
  labs(title = "TMM-normalised expression per sample",
       subtitle = "log-CPM after TMM normalisation (filtered genes)",
       x = NULL, y = "log2 CPM (TMM)") +
  theme_pub() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_png(file.path(OUT, "tmm_boxplots.png"), p3)

cat(">> QC plots done:", paste(list.files(OUT), collapse = ", "), "\n")
