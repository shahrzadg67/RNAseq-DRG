#!/usr/bin/env Rscript
# ============================================================================
# 51_gsea_bc.R  —  BATCH-CORRECTED GSEA. Identical to 50_gsea.R but reads the
# batch-corrected DESeq2 cache (deseq2_bc) and writes to SEPARATE output paths
# so nothing existing is overwritten:
#   result tables -> app/data/gsea_bc.rds
#   running-score PNGs -> app/www/gsea_bc/<level>/<contrast>/<category>/<id>.png
#   pathview KEGG maps -> app/www/pathview_bc/<contrast>/
# Ranking metric = DESeq2 Wald stat from the batch-adjusted model.
# ============================================================================
source("/hpf/projects/msalter/sghazis/rnaseq_TUY35595/analysis/_common.R")
suppressPackageStartupMessages({
  library(clusterProfiler); library(org.Mm.eg.db); library(enrichplot)
  library(ggplot2)
})
options(bitmapType = "cairo")          # headless node: PNG via cairo, not X11
PATHVIEW_OK <- requireNamespace("pathview", quietly = TRUE)

# robust PNG save that never halts the run on a graphics-device hiccup
save_png <- function(file, plot, w = 7, h = 4, dpi = 85) {
  ok <- tryCatch({ ggplot2::ggsave(file, plot, width = w, height = h, dpi = dpi); TRUE },
                 error = function(e) tryCatch({
                   grDevices::png(file, width = w*dpi, height = h*dpi, res = dpi, type = "cairo")
                   print(plot); grDevices::dev.off(); TRUE }, error = function(e2) FALSE))
  ok
}

de <- read_cache("deseq2_bc")
WWW      <- file.path(PROJ, "app", "www")
GSEA_DIR <- file.path(WWW, "gsea_bc"); dir.create(GSEA_DIR, recursive = TRUE, showWarnings = FALSE)
PV_DIR   <- file.path(WWW, "pathview_bc"); dir.create(PV_DIR, recursive = TRUE, showWarnings = FALSE)
TOP_PLOTS <- 12L     # running-score PNGs per (level,contrast,category) — capped to keep the app bundle small

# Build an ENTREZ-keyed ranked vector from a DE table (collapse tx->gene).
ranked_list <- function(df) {
  d <- df[!is.na(df$entrez) & !is.na(df$stat), c("entrez","stat")]
  d <- d[order(-abs(d$stat)), ]
  d <- d[!duplicated(d$entrez), ]                 # one stat per gene (best)
  v <- sort(setNames(d$stat, d$entrez), decreasing = TRUE)
  v
}

run_one <- function(geneList, category) {
  set.seed(1)
  if (category == "KEGG") {
    gseKEGG(geneList, organism = "mmu", minGSSize = 15, maxGSSize = 500,
            pvalueCutoff = 0.25, verbose = FALSE)
  } else {
    gseGO(geneList, OrgDb = org.Mm.eg.db, ont = category, keyType = "ENTREZID",
          minGSSize = 15, maxGSSize = 500, pvalueCutoff = 0.25, verbose = FALSE)
  }
}

CATS <- c("BP","CC","MF","KEGG")
gsea <- list(); plot_counts <- list()

for (lv in LEVELS) for (cn in names(CONTRASTS)) {
  gl <- ranked_list(de[[lv]]$results[[cn]])
  if (length(gl) < 25) { message("skip ", lv,"/",cn," (too few ranked genes)"); next }
  for (cat in CATS) {
    key <- paste(lv, cn, cat, sep = "::")
    message("GSEA(bc) :: ", key, "  (", length(gl), " genes)")
    res <- tryCatch(run_one(gl, cat), error = function(e) { message("  ", e$message); NULL })
    if (is.null(res) || nrow(as.data.frame(res)) == 0) next
    rdf <- as.data.frame(res)
    rdf$level <- lv; rdf$contrast <- cn; rdf$category <- cat
    gsea[[key]] <- rdf

    # precompute running-score PNGs for the top sets (by p.adjust)
    ord <- order(rdf$p.adjust)[seq_len(min(TOP_PLOTS, nrow(rdf)))]
    odir <- file.path(GSEA_DIR, lv, cn, cat); dir.create(odir, recursive = TRUE, showWarnings = FALSE)
    for (i in ord) {
      sid <- rdf$ID[i]
      p <- tryCatch(gseaplot2(res, geneSetID = sid,
                              title = paste0(sid, " — ", rdf$Description[i])),
                    error = function(e) NULL)
      if (!is.null(p)) save_png(file.path(odir, paste0(gsub("[:/]", "_", sid), ".png")), p)
    }
    plot_counts[[key]] <- length(ord)

    # pathview KEGG maps for top enriched pathways (gene-level only to avoid dup work)
    if (cat == "KEGG" && PATHVIEW_OK && lv == "gene") {
      fc <- setNames(de[[lv]]$results[[cn]]$log2FC, de[[lv]]$results[[cn]]$entrez)
      fc <- fc[!is.na(names(fc))]
      pvd <- file.path(PV_DIR, cn); dir.create(pvd, recursive = TRUE, showWarnings = FALSE)
      old <- setwd(pvd); on.exit(setwd(old), add = TRUE)
      for (pid in rdf$ID[order(rdf$p.adjust)][seq_len(min(5, nrow(rdf)))]) {
        try(pathview::pathview(gene.data = fc, pathway.id = sub("^mmu","",pid),
                               species = "mmu", limit = list(gene = 3)), silent = TRUE)
      }
      setwd(old)
    }
  }
}

saveRDS(gsea, file.path(APP_DATA, "gsea_bc.rds"))
message("\nGSEA (batch-corrected) done. result sets: ", length(gsea),
        " | running-score PNGs: ", sum(unlist(plot_counts)),
        " (capped at ", TOP_PLOTS, "/category)")
