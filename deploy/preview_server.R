#!/usr/bin/env Rscript
# Robust local preview for the shinylive site: serves deploy/site_staging with
# the Cross-Origin-Isolation headers shinylive/WebR need (so it works behind the
# VS Code port-forward proxy without relying on the service worker).
#   module load R/4.5.1 libarchive/3.8.1
#   Rscript deploy/preview_server.R [port]
library(httpuv)
args <- commandArgs(trailingOnly = TRUE)
port <- if (length(args)) as.integer(args[1]) else 8008L
site <- "/hpf/projects/msalter/sghazis/rnaseq_TUY35595/deploy/site_staging"

sp <- staticPath(
  site,
  indexhtml = TRUE,
  fallthrough = FALSE,
  headers = list(
    "Cross-Origin-Opener-Policy"   = "same-origin",
    "Cross-Origin-Embedder-Policy" = "require-corp",
    "Cross-Origin-Resource-Policy" = "cross-origin",
    "Cache-Control"                = "no-cache"
  )
)
cat(sprintf("Serving %s on 0.0.0.0:%d  (COOP/COEP enabled)\n", site, port))
runServer("0.0.0.0", port, list(staticPaths = list("/" = sp)))
