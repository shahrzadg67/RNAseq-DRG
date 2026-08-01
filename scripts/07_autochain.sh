#!/bin/bash
#
# Self-driving launcher: wait for container pre-staging to finish, then submit
# the smoke test and chain the real run to it (real runs only if smoke succeeds).
# Run fully detached so it survives any disconnect:
#   cd /hpf/projects/msalter/sghazis/rnaseq_TUY35595
#   setsid bash scripts/07_autochain.sh >logs/autochain.log 2>&1 </dev/null &
#
# Everything it submits is a Slurm job, so the pipeline runs to completion with
# nobody connected. Job IDs are written to logs/autochain_jobids.txt.

set -uo pipefail
BASE=/hpf/projects/msalter/sghazis/rnaseq_TUY35595
cd "$BASE"

echo "[$(date)] autochain start"

# 1. Wait for the login-node container pre-stage to finish (cap at ~90 min;
#    if the process already died, pgrep returns nothing and we proceed -- the
#    compute nodes have internet so the run can pull anything still missing).
waited=0
while pgrep -f 'nf-core pipelines download' >/dev/null 2>&1; do
    sleep 30
    waited=$((waited+30))
    [[ $waited -ge 5400 ]] && { echo "[$(date)] prestage wait cap reached, proceeding"; break; }
done
echo "[$(date)] prestage finished/absent; singularity_cache = $(du -sh singularity_cache 2>/dev/null | cut -f1)"

# 2. Sanity: all 28 fastq.gz present and non-empty.
n=$(find fastq -name '*.fastq.gz' -size +0c 2>/dev/null | wc -l)
echo "[$(date)] fastq.gz ready: $n/28"
if [[ "$n" -ne 28 ]]; then
    echo "[$(date)] ERROR: expected 28 fastq.gz, found $n -- aborting chain." >&2
    exit 1
fi

# 3. Submit smoke test, then real run gated on smoke success (afterok).
SMOKE=$(sbatch --parsable scripts/05_smoke_test.sh)
echo "[$(date)] submitted smoke test: job $SMOKE"
REAL=$(sbatch --parsable --dependency=afterok:"$SMOKE" run_rnaseq.sh)
echo "[$(date)] submitted real run: job $REAL (runs after $SMOKE succeeds)"

printf 'smoke=%s\nreal=%s\n' "$SMOKE" "$REAL" > logs/autochain_jobids.txt
echo "[$(date)] autochain done. See logs/autochain_jobids.txt"
