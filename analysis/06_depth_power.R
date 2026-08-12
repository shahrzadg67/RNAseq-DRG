# =====================================================================
# 06_depth_power.R — depth-normalisation + sensitivity/specificity check
# for the Sham F vs M contrast (the sex comparison with unequal depth).
# Produces app/www/qc/depth_concordance.png for the QC "Depth & power" panel,
# and prints the headline numbers used in the QC / Methods / Summary text.
# =====================================================================
suppressMessages({ library(DESeq2); library(edgeR); library(ggplot2) })
set.seed(1)
PROJ <- "/hpf/projects/msalter/sghazis/rnaseq_TUY35595"
OUT  <- file.path(PROJ, "app/www/qc"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
PADJ <- 0.05; LFC <- 1

SEX_COLORS <- c(F = "#941200", M = "#021893")
theme_pub <- function(bs = 14) theme_bw(base_size = bs) +
  theme(panel.border = element_rect(colour="black", fill=NA, linewidth=0.9),
        panel.grid = element_blank(), axis.text = element_text(colour="black"),
        plot.title = element_text(face="bold"), legend.key = element_blank())

bt  <- read.delim(file.path(PROJ,"results/star_salmon/salmon.merged.gene_counts.tsv"), check.names=FALSE)
ann <- readRDS(file.path(PROJ,"app/data/gene_annot.rds")); ann <- ann[ann$level=="gene",]
mat <- round(as.matrix(bt[,-(1:2)])); rownames(mat) <- bt[,1]; symv <- bt[,2]
sham <- grep("^Sham_", colnames(mat), value=TRUE)
m <- mat[, sham]; sex <- factor(sub("^Sham_([FM]).*","\\1", sham), levels=c("M","F"))
keep <- rowSums(m >= 10) >= 3; m <- m[keep,]; sym <- symv[keep]
chr <- ann$chr[match(rownames(m), ann$feature_id)]

run_de <- function(counts, grp) {
  dds <- DESeqDataSetFromMatrix(counts, DataFrame(g=grp), ~ g); dds <- DESeq(dds, quiet=TRUE)
  lv <- levels(grp); as.data.frame(results(dds, contrast=c("g", lv[2], lv[1])))
}
is_deg <- function(r) !is.na(r$padj) & r$padj < PADJ & abs(r$log2FoldChange) > LFC

obs <- run_de(m, sex); deg <- is_deg(obs); obs_ids <- rownames(m)[deg]
pos <- unique(c(rownames(m)[!is.na(chr) & chr=="Y"], rownames(m)[sym=="Xist"]))
sens <- mean(pos %in% obs_ids)

# permutation specificity
combs <- combn(6,3); keepc <- which(combs[1,]==1); true_set <- sort(which(sex=="F"))
nullc <- c()
for (j in keepc) { gA <- combs[,j]
  if (identical(sort(gA),true_set) || identical(sort(setdiff(1:6,gA)),true_set)) next
  g <- factor(ifelse(1:6 %in% gA,"A","B"), levels=c("B","A"))
  r <- tryCatch(run_de(m,g), error=function(e) NULL); if(!is.null(r)) nullc <- c(nullc, sum(is_deg(r)))
}

# equal-depth downsampling
target <- min(colSums(m)); thin <- edgeR::thinCounts(m, target.size=target)
ds <- run_de(thin, sex); ds_ids <- rownames(thin)[is_deg(ds)]
recov <- mean(obs_ids %in% ds_ids)
common <- intersect(rownames(obs), rownames(ds))
xf <- obs$log2FoldChange[match(common,rownames(obs))]; xe <- ds$log2FoldChange[match(common,rownames(ds))]
rr <- cor(xf, xe, use="complete.obs")

# ---- concordance figure: log2FC full vs equal-depth ----
df <- data.frame(full=xf, equal=xe, id=common, sym=sym[match(common,rownames(m))],
                 gold=common %in% pos)
lim <- range(c(df$full, df$equal), na.rm=TRUE)
p <- ggplot(df, aes(full, equal)) +
  geom_abline(slope=1, intercept=0, linetype="dotted", colour="grey40") +
  geom_point(data=subset(df,!gold), size=0.5, alpha=0.25, colour="grey55") +
  geom_point(data=subset(df, gold), colour="#E67E22", size=2.4) +
  ggrepel::geom_text_repel(data=subset(df, gold), aes(label=sym), size=3.6, colour="black",
                           max.overlaps=Inf, min.segment.length=0, seed=1) +
  coord_equal(xlim=lim, ylim=lim) +
  labs(title=sprintf("Sham F vs M effect sizes are unchanged by equal depth (r = %.3f)", rr),
       subtitle=sprintf("Each point a gene · all 6 libraries down-sampled to %.1fM reads · orange = gold-standard sex genes", target/1e6),
       x="log2 fold-change (full depth)", y="log2 fold-change (equal depth)") +
  theme_pub()
ragg::agg_png(file.path(OUT,"depth_concordance.png"), width=8, height=7, units="in", res=300)
print(p); grDevices::dev.off()

cat("=== HEADLINE NUMBERS ===\n")
cat(sprintf("observed DEGs (F vs M): %d\n", length(obs_ids)))
cat(sprintf("sensitivity (gold-standard %d genes): %.0f%%\n", length(pos), 100*sens))
cat(sprintf("null DEG counts (9 random splits): %s | median %g | mean %.1f\n",
    paste(sort(nullc),collapse=","), median(nullc), mean(nullc)))
cat(sprintf("specificity: median FPR %.3f%% (spec %.2f%%) | mean FPR %.3f%% (spec %.2f%%)\n",
    100*median(nullc)/nrow(m), 100*(1-median(nullc)/nrow(m)),
    100*mean(nullc)/nrow(m),   100*(1-mean(nullc)/nrow(m))))
cat(sprintf("equal-depth: recovered %d/%d = %.0f%% | log2FC r = %.3f | pos-control sens %.0f%%\n",
    sum(obs_ids %in% ds_ids), length(obs_ids), 100*recov, rr, 100*mean(pos %in% ds_ids)))
cat(">> wrote", file.path(OUT,"depth_concordance.png"), "\n")
