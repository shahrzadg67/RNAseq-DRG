#!/usr/bin/env Rscript
# ============================================================================
# 91_interactive_html.R  —  self-contained INTERACTIVE plotly HTML for BOTH the
# original (uncorrected, with SNI) and the batch-corrected analyses. No Shiny
# server, no websocket, no WebAssembly: plain static HTML that works over any
# HTTP proxy and can be opened directly in VS Code. Widgets share one lib/ dir.
# Output: batch_corrected_results/interactive/ with a two-section index.html.
# ============================================================================
source("/hpf/projects/msalter/sghazis/rnaseq_TUY35595/analysis/_common.R")
suppressPackageStartupMessages({ library(plotly); library(htmlwidgets) })

OUT <- file.path(PROJ, "batch_corrected_results", "interactive")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
# clean out any prior run so the index reflects exactly what we build now
old <- list.files(OUT, pattern = "\\.html$", full.names = TRUE); file.remove(old)

save_widget <- function(w, file) {
  htmlwidgets::saveWidget(w, file.path(OUT, file), selfcontained = FALSE,
                          libdir = "lib", title = sub("\\.html$", "", file))
}

## ---- reusable builders -----------------------------------------------------
pca_widget <- function(co, ve, cb, title) {
  co$.col <- co[[cb]]
  plot_ly(co, x = ~PC1, y = ~PC2, color = ~.col, type = "scatter", mode = "markers",
          marker = list(size = 13, line = list(width = 1, color = "#fff")),
          text = ~paste0("<b>", sample_id, "</b><br>group: ", group,
                         "<br>condition: ", condition, "<br>sex: ", sex,
                         "<br>replicate: ", replicate),
          hoverinfo = "text") |>
    layout(title = title,
           xaxis = list(title = sprintf("PC1 (%.1f%%)", ve$pct[1])),
           yaxis = list(title = sprintf("PC2 (%.1f%%)", ve$pct[2])),
           legend = list(title = list(text = cb)))
}

volcano_widget <- function(d, title) {
  d <- d[!is.na(d$padj), ]
  d$sig <- d$padj < 0.05 & abs(d$log2FC) >= 1
  d$lab <- ifelse(is.na(d$symbol), d$feature_id, d$symbol)
  d$sigcat <- ifelse(d$sig, "significant", "NS")
  ns <- sum(d$sig)
  if (sum(!d$sig) > 15000) {                 # keep all sig; thin the NS cloud (memory/browser)
    keep_ns <- sample(which(!d$sig), 15000)
    d <- d[sort(c(which(d$sig), keep_ns)), ]
  }
  plot_ly(d, x = ~log2FC, y = ~ -log10(padj), color = ~sigcat,
          colors = c(NS = "#b0b8bf", significant = "#e74c3c"),
          type = "scattergl", mode = "markers",
          marker = list(size = 5, opacity = 0.55),
          text = ~paste0("<b>", lab, "</b><br>", feature_id,
                         "<br>log2FC: ", round(log2FC, 2), "<br>adj.p: ", signif(padj, 3)),
          hoverinfo = "text") |>
    layout(title = sprintf("%s · %d sig (padj<0.05, |log2FC|>=1)", title, ns),
           xaxis = list(title = "log2 fold-change"), yaxis = list(title = "-log10 adj.p"),
           shapes = list(
             list(type="line", x0=1,  x1=1,  y0=0, y1=1, yref="paper", line=list(dash="dot", color="#888")),
             list(type="line", x0=-1, x1=-1, y0=0, y1=1, yref="paper", line=list(dash="dot", color="#888"))))
}

## ---- ORIGINAL (uncorrected, with SNI) --------------------------------------
pca0 <- readRDS(file.path(APP_DATA, "pca.rds"))
de0  <- readRDS(file.path(APP_DATA, "de_long.rds"))
for (key in names(pca0)) {                    # gene_withSNI, gene_noSNI, transcript_*
  d <- pca0[[key]]
  for (cb in c("group","condition","sex")) {
    w <- pca_widget(d$coords, d$var_explained, cb,
                    sprintf("PCA (original) — %s, coloured by %s", key, cb))
    tryCatch(save_widget(w, sprintf("orig_PCA_%s_by_%s.html", key, cb)),
             error = function(e) message("save fail: ", conditionMessage(e)))
    rm(w); gc(FALSE)
  }
}
for (en in c("DESeq2","edgeR")) for (lv in c("gene","transcript")) for (cn in names(CONTRASTS)) {
  d <- de0[de0$engine==en & de0$level==lv & de0$contrast==cn, ]
  if (!nrow(d)) next
  w <- volcano_widget(d, sprintf("%s — %s — %s (original)", en, lv, cn))
  tryCatch(save_widget(w, sprintf("orig_volcano_%s_%s_%s.html", en, lv, cn)),
           error = function(e) message("save fail: ", conditionMessage(e)))
  rm(w, d); gc(FALSE)
}
message("original HTML written")

## ---- BATCH-CORRECTED (SNI excluded, replicate removed) ---------------------
pcaB <- readRDS(file.path(APP_DATA, "pca_bc.rds"))
deB  <- readRDS(file.path(APP_DATA, "de_long_bc.rds"))
for (lv in names(pcaB)) {                     # gene, transcript
  d <- pcaB[[lv]]
  for (cb in c("group","condition","sex")) {
    w <- pca_widget(d$coords, d$var_explained, cb,
                    sprintf("PCA (batch-corrected) — %s level, coloured by %s", lv, cb))
    tryCatch(save_widget(w, sprintf("bc_PCA_%s_by_%s.html", lv, cb)),
             error = function(e) message("save fail: ", conditionMessage(e)))
    rm(w); gc(FALSE)
  }
}
for (en in c("DESeq2","edgeR")) for (lv in c("gene","transcript")) for (cn in names(CONTRASTS)) {
  d <- deB[deB$engine==en & deB$level==lv & deB$contrast==cn, ]
  if (!nrow(d)) next
  w <- volcano_widget(d, sprintf("%s — %s — %s (batch-corrected)", en, lv, cn))
  tryCatch(save_widget(w, sprintf("bc_volcano_%s_%s_%s.html", en, lv, cn)),
           error = function(e) message("save fail: ", conditionMessage(e)))
  rm(w, d); gc(FALSE)
}
message("batch-corrected HTML written")

## ---- index.html (two sections, built from files on disk) -------------------
li <- function(pat) {
  fs <- sort(list.files(OUT, pattern = pat))
  if (!length(fs)) return("<li><i>none</i></li>")
  paste0("<li><a href='", fs, "'>", sub("\\.html$", "", fs), "</a></li>", collapse = "\n")
}
html <- paste0(
"<!doctype html><html><head><meta charset='utf-8'><title>TUY35595 — interactive plots</title>",
"<style>body{font-family:system-ui,Arial,sans-serif;max-width:960px;margin:2rem auto;padding:0 1rem;color:#222}",
"h1{color:#18bc9c}h2{margin-top:1.8rem;background:#f4f6f8;border-left:4px solid #18bc9c;padding:.4rem .8rem}",
"h3{margin:1rem 0 .3rem;color:#555}ul{line-height:1.85;columns:2}a{color:#2c7fb8;text-decoration:none}",
"a:hover{text-decoration:underline}.note{background:#fff8e1;border-left:4px solid #f39c12;padding:.6rem 1rem;border-radius:4px}</style></head><body>",
"<h1>TUY35595 — interactive plots</h1>",
"<p class='note'>Interactive plotly figures (hover for gene names, zoom, pan). Two versions are provided.</p>",

"<h2>Original (uncorrected, includes SNI)</h2>",
"<h3>PCA</h3><ul>", li("^orig_PCA_.*\\.html$"), "</ul>",
"<h3>Volcano</h3><ul>", li("^orig_volcano_.*\\.html$"), "</ul>",

"<h2>Batch-corrected (SNI excluded, replicate effect removed)</h2>",
"<h3>PCA</h3><ul>", li("^bc_PCA_.*\\.html$"), "</ul>",
"<h3>Volcano</h3><ul>", li("^bc_volcano_.*\\.html$"), "</ul>",

"<p style='margin-top:1.5rem'><b>Static PNGs & DE tables:</b> <code>batch_corrected_results/{PCA_plots,volcano_plots,GSEA_plots,DE_tables}/</code></p>",
"</body></html>")
writeLines(html, file.path(OUT, "index.html"))
n <- length(list.files(OUT, pattern = "\\.html$"))
message("\nInteractive report -> ", file.path(OUT, "index.html"), "  (", n, " html total)")
