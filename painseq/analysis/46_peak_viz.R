# =====================================================================
# 46_peak_viz.R — visualise the peak table.
#
# data/panel_celltype_peak.csv records, for every panel gene, WHERE (cell type),
# WHEN (injury model + timepoint) and HOW STRONGLY it reaches its maximum
# induction. As a table that is 49 rows of numbers; as a picture it is a
# temporal wave — immediate-early genes peaking within a day, effector genes
# peaking at a week — which is the part worth seeing.
#
# Out: figs/peak_landscape.{png,pdf}, data/peak_viz.rds (for the Shiny tab)
# =====================================================================
suppressMessages({ library(ggplot2); library(patchwork) })
source("/hpf/projects/msalter/sghazis/painseq/analysis/00_common.R")

pk <- read.csv(file.path(DATA, "panel_celltype_peak.csv"), stringsAsFactors = FALSE)
pk <- pk[!is.na(pk$log2FC), ]
msg("peak table: ", nrow(pk), " genes with a fitted peak")

TIME_LV <- c(6, 12, 24, 36, 48, 72, 168, 336, 672, 1440, 2160)
fmt_time <- function(h) ifelse(h == 0, "0", ifelse(h < 168, paste0(h, "h"), paste0(h / 24, "d")))

pk$tlab <- factor(fmt_time(pk$peak_time_h), levels = fmt_time(TIME_LV))
pk$peak_celltype <- factor(pk$peak_celltype,
                           levels = intersect(names(CELLTYPE_FULL), unique(pk$peak_celltype)))
pk$peak_injury <- factor(pk$peak_injury, levels = intersect(INJ_ORDER, unique(pk$peak_injury)))
pk$neuronal <- pk$peak_celltype %in% NEURONAL_TYPES
pk$lab <- ifelse(pk$role == "reference", paste0(pk$Gene, " ‡"), pk$Gene)

# time-class: how fast does each gene reach its maximum
pk$wave <- cut(pk$peak_time_h, breaks = c(-1, 24, 72, 168, Inf),
               labels = c("Early (≤24h)", "Mid (36-72h)", "1 week", "Late (>1 week)"))
msg("\n--- when do panel genes peak? ---")
print(table(pk$wave))
msg("\n--- where do they peak? ---")
print(sort(table(as.character(pk$peak_celltype)), decreasing = TRUE))
msg("\n--- in which injury model? ---")
print(sort(table(as.character(pk$peak_injury)), decreasing = TRUE))

PAL <- setNames(c("#4E79A7","#F28E2B","#E15759","#76B7B2","#59A14F","#EDC948","#B07AA1",
                  "#FF9DA7","#9C755F","#BAB0AC","#1F77B4","#D62728","#2CA02C","#9467BD",
                  "#8C564B","#E377C2","#7F7F7F","#BCBD22")[seq_len(nlevels(pk$peak_celltype))],
                levels(pk$peak_celltype))

# ---- main panel: when vs how strongly, coloured by where --------------------
main <- ggplot(pk, aes(tlab, log2FC)) +
  geom_hline(yintercept = 0, linewidth = .25, colour = "grey75") +
  geom_point(aes(colour = peak_celltype, shape = peak_injury), size = 3.1, alpha = .92) +
  ggrepel::geom_text_repel(aes(label = lab), size = 2.8, max.overlaps = 30, seed = 4,
                           colour = "grey15", segment.colour = "grey70", segment.size = .25) +
  scale_colour_manual(values = PAL, name = "Peak cell type") +
  scale_shape_manual(values = c(Crush = 16, ScNT = 17, SpNT = 15, CFA = 3, Paclitaxel = 4),
                     name = "Peak injury model") +
  labs(title = "Peak landscape of the conserved injury panel",
       subtitle = paste(strwrap(paste(
         "Each gene at the cell type, injury model and timepoint where it reaches maximum induction.",
         "The x axis is a temporal wave: transcription factors and immediate-early genes crest within",
         "a day, secreted effectors and structural genes a week later. ‡ = RAG reference gene,",
         "not one of the 47."), width = 118), collapse = "\n"),
       x = "Time at which the gene peaks", y = "Peak log2FC vs Naive") +
  theme_pub(base_size = 12) +
  theme(plot.subtitle = element_text(colour = "grey30", size = 9),
        legend.box = "vertical", legend.key.height = unit(12, "pt"))

# ---- marginals -------------------------------------------------------------
ct_n <- as.data.frame(table(peak_celltype = pk$peak_celltype))
ct_n <- ct_n[ct_n$Freq > 0, ]
p_ct <- ggplot(ct_n, aes(reorder(peak_celltype, Freq), Freq, fill = peak_celltype)) +
  geom_col(width = .72, show.legend = FALSE) +
  geom_text(aes(label = Freq), hjust = -0.3, size = 3, colour = "grey25") +
  scale_fill_manual(values = PAL) + coord_flip() +
  scale_y_continuous(expand = expansion(mult = c(0, .18))) +
  labs(title = "Genes by peak cell type", x = NULL, y = "genes") +
  theme_pub(base_size = 11) + theme(plot.title = element_text(size = 11))

t_n <- as.data.frame(table(tlab = pk$tlab)); t_n <- t_n[t_n$Freq > 0, ]
p_t <- ggplot(t_n, aes(tlab, Freq)) +
  geom_col(width = .72, fill = "#B2182B", alpha = .85) +
  geom_text(aes(label = Freq), vjust = -0.35, size = 3, colour = "grey25") +
  scale_y_continuous(expand = expansion(mult = c(0, .18))) +
  labs(title = "Genes by peak timepoint", x = NULL, y = "genes") +
  theme_pub(base_size = 11) + theme(plot.title = element_text(size = 11))

out <- main / (p_ct | p_t) + patchwork::plot_layout(heights = c(2.5, 1))
ggsave(file.path(FIGS, "peak_landscape.png"), out, width = 13, height = 13, dpi = 300, bg = "white")
ggsave(file.path(FIGS, "peak_landscape.pdf"), out, width = 13, height = 13, device = cairo_pdf)
msg("wrote figs/peak_landscape.{png,pdf}")

saveRDS(pk[, c("Gene", "lab", "role", "peak_celltype", "peak_celltype_full", "peak_injury",
               "peak_time_h", "tlab", "wave", "log2FC", "adj_p", "pct_cells_detected")],
        file.path(DATA, "peak_viz.rds"))
msg("wrote data/peak_viz.rds")
