#!/usr/bin/env Rscript
# Install all R packages needed for the TUY35595 downstream analysis + app build.
# Installs into a PROJECT-LOCAL library (home is only 10 GB on this cluster).
#
#   module load R/4.5.1
#   export R_LIBS_USER=/hpf/projects/msalter/sghazis/rnaseq_TUY35595/.Rlib
#   Rscript analysis/install_packages.R
#
# Heavy build (DESeq2 / clusterProfiler etc.). Run in the background and tail the log.

lib <- Sys.getenv("R_LIBS_USER")
if (nzchar(lib)) dir.create(lib, showWarnings = FALSE, recursive = TRUE)
.libPaths(c(lib, .libPaths()))
options(Ncpus = max(1, parallel::detectCores() - 2),
        repos = c(CRAN = "https://cloud.r-project.org"))

message("Library target: ", .libPaths()[1])

have <- function(p) requireNamespace(p, quietly = TRUE)
ensure_cran <- function(pkgs) {
  todo <- pkgs[!vapply(pkgs, have, logical(1))]
  if (length(todo)) install.packages(todo)
}

ensure_cran(c("BiocManager", "remotes",
              "ggplot2", "dplyr", "tidyr", "tibble", "readr", "stringr", "purrr",
              "plotly", "DT", "shiny", "bslib", "cowplot", "scales",
              "yaml", "jsonlite", "matrixStats", "R.utils", "htmlwidgets"))

library(BiocManager)
# Pin to the Bioconductor release matching R 4.5.x
bioc_pkgs <- c("DESeq2", "edgeR", "limma", "tximport", "tximeta",
               "clusterProfiler", "org.Mm.eg.db", "enrichplot", "pathview",
               "EnhancedVolcano", "apeglm", "AnnotationDbi", "fgsea")
todo <- bioc_pkgs[!vapply(bioc_pkgs, have, logical(1))]
if (length(todo)) BiocManager::install(todo, update = FALSE, ask = FALSE)

# shinylive (CRAN) — compiles the app to a static WebAssembly site
ensure_cran("shinylive")

# Final report
all_pkgs <- c("DESeq2","edgeR","limma","tximport","tximeta","clusterProfiler",
              "org.Mm.eg.db","enrichplot","pathview","EnhancedVolcano","apeglm",
              "fgsea","plotly","DT","shiny","bslib","shinylive","yaml","jsonlite",
              "ggplot2","dplyr","tidyr","cowplot","matrixStats")
status <- data.frame(pkg = all_pkgs,
                     installed = vapply(all_pkgs, have, logical(1)))
print(status)
missing <- status$pkg[!status$installed]
if (length(missing)) {
  message("STILL MISSING: ", paste(missing, collapse = ", "))
  quit(status = 1)
}
message("All packages installed OK.")
