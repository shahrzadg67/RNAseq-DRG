#!/usr/bin/env Rscript
# gsea_core_noxy.R  —  GSEA for the sex-chromosome-removed pipeline. Sourced by
# run_all_noxy.R (uses the in-memory `de` = deseq2_noxy results). Writes to
# app/data/gsea_noxy.rds and app/www/{gsea_noxy,pathview_noxy}/ so nothing else
# is overwritten. Mirrors 51_gsea_bc.R.
suppressPackageStartupMessages({
  library(clusterProfiler); library(org.Mm.eg.db); library(enrichplot); library(ggplot2)
})
options(bitmapType = "cairo")
PATHVIEW_OK <- requireNamespace("pathview", quietly = TRUE)
save_png <- function(file, plot, w=7, h=4, dpi=85) tryCatch(
  { ggplot2::ggsave(file, plot, width=w, height=h, dpi=dpi); TRUE },
  error=function(e) tryCatch({ grDevices::png(file, width=w*dpi, height=h*dpi, res=dpi, type="cairo")
    print(plot); grDevices::dev.off(); TRUE }, error=function(e2) FALSE))

WWW <- file.path(PROJ,"app","www")
GDIR <- file.path(WWW,"gsea_noxy"); dir.create(GDIR, recursive=TRUE, showWarnings=FALSE)
PDIR <- file.path(WWW,"pathview_noxy"); dir.create(PDIR, recursive=TRUE, showWarnings=FALSE)
TOP_PLOTS <- 12L

ranked_list <- function(df){ d <- df[!is.na(df$entrez) & !is.na(df$stat), c("entrez","stat")]
  d <- d[order(-abs(d$stat)), ]; d <- d[!duplicated(d$entrez), ]
  sort(setNames(d$stat, d$entrez), decreasing=TRUE) }
run_one <- function(gl, cat){ set.seed(1)
  if (cat=="KEGG") gseKEGG(gl, organism="mmu", minGSSize=15, maxGSSize=500, pvalueCutoff=0.25, verbose=FALSE)
  else gseGO(gl, OrgDb=org.Mm.eg.db, ont=cat, keyType="ENTREZID", minGSSize=15, maxGSSize=500, pvalueCutoff=0.25, verbose=FALSE) }

CATS <- c("BP","CC","MF","KEGG"); gsea <- list()
for (lv in LEVELS) for (cn in names(CONTRASTS)) {
  gl <- ranked_list(de[[lv]]$results[[cn]])
  if (length(gl) < 25) { message("skip ",lv,"/",cn); next }
  for (cat in CATS) {
    key <- paste(lv, cn, cat, sep="::"); message("GSEA(noXY) :: ", key, " (", length(gl), ")")
    res <- tryCatch(run_one(gl, cat), error=function(e){ message("  ", e$message); NULL })
    if (is.null(res) || nrow(as.data.frame(res))==0) next
    rdf <- as.data.frame(res); rdf$level <- lv; rdf$contrast <- cn; rdf$category <- cat
    gsea[[key]] <- rdf
    ord <- order(rdf$p.adjust)[seq_len(min(TOP_PLOTS, nrow(rdf)))]
    odir <- file.path(GDIR, lv, cn, cat); dir.create(odir, recursive=TRUE, showWarnings=FALSE)
    for (i in ord) { sid <- rdf$ID[i]
      p <- tryCatch(gseaplot2(res, geneSetID=sid, title=paste0(sid," — ",rdf$Description[i])), error=function(e) NULL)
      if (!is.null(p)) save_png(file.path(odir, paste0(gsub("[:/]","_",sid),".png")), p) }
    if (cat=="KEGG" && PATHVIEW_OK && lv=="gene") {
      fc <- setNames(de[[lv]]$results[[cn]]$log2FC, de[[lv]]$results[[cn]]$entrez); fc <- fc[!is.na(names(fc))]
      pvd <- file.path(PDIR, cn); dir.create(pvd, recursive=TRUE, showWarnings=FALSE)
      old <- setwd(pvd); on.exit(setwd(old), add=TRUE)
      for (pid in rdf$ID[order(rdf$p.adjust)][seq_len(min(5, nrow(rdf)))])
        try(pathview::pathview(gene.data=fc, pathway.id=sub("^mmu","",pid), species="mmu", limit=list(gene=3)), silent=TRUE)
      setwd(old)
    }
  }
}
saveRDS(gsea, file.path(APP_DATA,"gsea_noxy.rds"))
message("GSEA (noXY) done. sets: ", length(gsea))
