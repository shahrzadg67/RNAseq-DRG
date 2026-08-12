#!/usr/bin/env Rscript
# Deploy the app to shinyapps.io. Run via deploy/deploy_shinyapps.sh (which sets
# the module env + OpenSSL shim so uploads work from this HPC node).
#
# ONE-TIME account setup first (copy the 3-arg line from your shinyapps.io
# dashboard → Account → Tokens → Show):
#   Rscript -e 'rsconnect::setAccountInfo(name="<acct>", token="<token>", secret="<secret>")'
#
# Then: bash deploy/deploy_shinyapps.sh
suppressPackageStartupMessages(library(rsconnect))
APP <- "/hpf/projects/msalter/sghazis/rnaseq_TUY35595/deploy/app_deploy"
rsconnect::deployApp(
  appDir       = APP,
  appName      = "tuy35595-rnaseq",
  appTitle     = "RNA-seq Differential Expression — TUY35595",
  account      = if (length(rsconnect::accounts()$name)) rsconnect::accounts()$name[1] else NULL,
  forceUpdate  = TRUE,
  launch.browser = FALSE,
  logLevel     = "normal")
