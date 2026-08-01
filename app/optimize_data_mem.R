# Slim the deployed .rds data to cut shinyapps memory (OOM fix). SAFE transforms only:
#  - factor-encode repeated character columns EXCEPT `feature_id` (used as a match()/
#    merge() key in attach_annot/venn — kept character to avoid factor-match pitfalls).
#  - drop constant single-value columns that no visible tab references (rat/monkey level,engine).
# Numeric log2FC/padj and feature_id are left as-is. No rows dropped (all features kept).
setwd("/hpf/projects/msalter/sghazis/rnaseq_TUY35595/app")
fact <- function(df, cols) { for (c in intersect(cols, names(df))) if (!is.factor(df[[c]])) df[[c]] <- factor(df[[c]]); df }
sz   <- function(x) as.numeric(object.size(x))/1e6
report <- function(nm, before, after) cat(sprintf("  %-24s %6.1f -> %6.1f MB\n", nm, before, after))

opt <- function(name, factor_cols, drop_cols = character(0)) {
  f <- file.path("data", paste0(name, ".rds"))
  if (!file.exists(f)) { cat("  (skip, missing)", name, "\n"); return(invisible()) }
  d <- readRDS(f); b <- sz(d)
  d <- d[, setdiff(names(d), drop_cols), drop = FALSE]
  d <- fact(d, factor_cols)
  saveRDS(d, f, compress = "xz")
  report(name, b, sz(d))
}

cat("Optimising deployed data files:\n")
# NOTE: de_long / de_long_noxy are intentionally left untouched — their `symbol`
# column feeds Venn/cross-species SET operations (unique/union -> ggVennDiagram),
# where a factor could misbehave. Their contrast/level/engine are already factors.
#
# gene annotation lookup: safe to factor ALL char cols incl feature_id, because it
# is only ever a match()/== target (match() coerces factors to character).
opt("gene_annot", c("level","feature_id","gene_id","symbol","description","chr"))
# rat/monkey DE: `symbol` here is only shown in tables/hovers (never a set op),
# so factoring is safe. Drop constant unused level+engine columns.
opt("rat_de_long",      c("symbol","contrast","rat_gene_id"), drop_cols = c("level","engine"))
opt("rat_de_long_noxy", c("symbol","contrast","rat_gene_id"), drop_cols = c("level","engine"))
opt("monkey_de_long",      c("symbol","contrast","monkey_gene_id"), drop_cols = c("level","engine"))
opt("monkey_de_long_noxy", c("symbol","contrast","monkey_gene_id"), drop_cols = c("level","engine"))

cat("\nNew total in-memory (deployed non-bc files):\n")
fs <- list.files("data", pattern="\\.rds$", full.names=TRUE); fs <- fs[!grepl("_bc\\.rds$", fs)]
tot <- sum(vapply(fs, function(f) sz(readRDS(f)), numeric(1)))
cat(sprintf("  TOTAL: %.0f MB\n", tot))
