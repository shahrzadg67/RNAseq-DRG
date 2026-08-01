# =====================================================================
# 07_gsea_leading_edge.R — add human-readable leading-edge (core-enrichment)
# gene SYMBOLS to the precomputed GSEA results, so the app can show "which
# genes drive this term" WITHOUT needing org.Mm.eg.db at runtime.
# Adds columns: core_symbols (comma-separated) and core_n (count).
# =====================================================================
suppressMessages({ library(org.Mm.eg.db); library(AnnotationDbi) })
PROJ <- "/hpf/projects/msalter/sghazis/rnaseq_TUY35595"

add_symbols <- function(f) {
  path <- file.path(PROJ, "app/data", f)
  if (!file.exists(path)) { cat("skip (absent):", f, "\n"); return(invisible()) }
  g <- readRDS(path)
  # collect every Entrez id used across all sets, map once
  all_ids <- unique(unlist(lapply(g, function(d)
    if (is.data.frame(d) && "core_enrichment" %in% names(d))
      unlist(strsplit(d$core_enrichment, "/")) else character(0))))
  all_ids <- all_ids[!is.na(all_ids) & all_ids != ""]
  map <- AnnotationDbi::mapIds(org.Mm.eg.db, keys = all_ids, column = "SYMBOL",
                               keytype = "ENTREZID", multiVals = "first")
  g <- lapply(g, function(d) {
    if (!is.data.frame(d) || !"core_enrichment" %in% names(d)) return(d)
    parts <- strsplit(d$core_enrichment, "/")
    d$core_symbols <- vapply(parts, function(ids) {
      s <- map[ids]; s[is.na(s)] <- ids[is.na(s)]        # fall back to Entrez if unmapped
      paste(s, collapse = ", ")
    }, character(1))
    d$core_n <- vapply(parts, function(ids) sum(ids != "" & !is.na(ids)), integer(1))
    d
  })
  saveRDS(g, path)
  cat("updated", f, "-", length(g), "sets; example genes:\n  ",
      substr(g[[1]]$core_symbols[1], 1, 120), "...\n")
}

for (f in c("gsea.rds", "gsea_noxy.rds")) add_symbols(f)
cat(">> done.\n")
