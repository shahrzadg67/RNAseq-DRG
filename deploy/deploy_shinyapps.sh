#!/usr/bin/env bash
# Deploy the app to shinyapps.io from this HPC node (sets the module env +
# OpenSSL 1.0.2 shim so R's https uploads work here).
#   1) create a free account at https://www.shinyapps.io
#   2) Account -> Tokens -> Show -> copy the setAccountInfo(...) line, run once:
#        module load R/4.5.2 ; export LD_LIBRARY_PATH=<proj>/.extralib:$LD_LIBRARY_PATH
#        Rscript -e 'rsconnect::setAccountInfo(name="..", token="..", secret="..")'
#   3) bash deploy/deploy_shinyapps.sh
set -euo pipefail
PROJ="/hpf/projects/msalter/sghazis/rnaseq_TUY35595"
module load R/4.5.2 libarchive/3.8.1 2>/dev/null || true
export R_LIBS_USER="$PROJ/.Rlib"
export LD_LIBRARY_PATH="$PROJ/.extralib:${LD_LIBRARY_PATH:-}"   # libssl.so.10 for https
cd "$PROJ"
Rscript deploy/shinyapps.R
