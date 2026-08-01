#!/usr/bin/env bash
# Start BOTH review servers fully detached (setsid) so they survive the launching
# shell AND the launching session ending. They only die if the whole
# OnDemand SLURM session ends (node change). Re-run this any time they're down.
#   bash start_servers.sh
ROOT=/hpf/projects/msalter/sghazis/rnaseq_TUY35595
LOGS="$ROOT/logs"
export R_LIBS_USER="$ROOT/.Rlib"
module load R/4.5.2 libarchive/3.8.1 2>/dev/null || true

# stop any existing instances by PID (never pkill -f, which self-matches)
for p in 8009 8010; do
  pid=$(ss -ltnp 2>/dev/null | grep ":$p" | grep -oE 'pid=[0-9]+' | head -1 | cut -d= -f2)
  [ -n "$pid" ] && kill "$pid" 2>/dev/null || true
done
sleep 2

# 8009 — live Shiny app
setsid bash -c "cd '$ROOT/app'; export R_LIBS_USER='$ROOT/.Rlib'; \
  exec Rscript -e 'shiny::runApp(getwd(), host=\"0.0.0.0\", port=8009L, launch.browser=FALSE)'" \
  > "$LOGS/live_shiny.log" 2>&1 < /dev/null &
disown 2>/dev/null || true

# 8010 — static interactive plots (both original + batch-corrected)
setsid bash -c "cd '$ROOT/batch_corrected_results/interactive'; \
  exec python3 -m http.server 8010 --bind 0.0.0.0" \
  > "$LOGS/static_8010.log" 2>&1 < /dev/null &
disown 2>/dev/null || true

echo "launched (detached). give the Shiny app ~15s to load."
