#!/usr/bin/env bash
# Launch the live Shiny review app on port 8009, fully detached so it survives
# the launching shell. Writes its PID to app/shiny.pid. To stop: kill $(cat app/shiny.pid)
set -euo pipefail
ROOT=/hpf/projects/msalter/sghazis/rnaseq_TUY35595
cd "$ROOT/app"
module load R/4.5.2 libarchive/3.8.1 2>/dev/null || true
export R_LIBS_USER="$ROOT/.Rlib"
PORT="${1:-8009}"
setsid Rscript -e "shiny::runApp(getwd(), host='0.0.0.0', port=${PORT}L, launch.browser=FALSE)" \
  > "$ROOT/logs/live_shiny.log" 2>&1 < /dev/null &
echo $! > "$ROOT/app/shiny.pid"
echo "launched shiny pid $(cat "$ROOT/app/shiny.pid") on port ${PORT}"
