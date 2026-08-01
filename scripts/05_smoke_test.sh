#!/bin/bash
#SBATCH --job-name=rnaseq_test
#SBATCH --account=msalter
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=06:00:00
#SBATCH --output=/hpf/projects/msalter/sghazis/rnaseq_TUY35595/logs/rnaseq_test_%j.out
#SBATCH --error=/hpf/projects/msalter/sghazis/rnaseq_TUY35595/logs/rnaseq_test_%j.err
#
# Smoke test: runs nf-core/rnaseq on its tiny built-in dataset to confirm Slurm
# submission, Singularity, and cache redirection all work BEFORE the real run.
#   sbatch scripts/05_smoke_test.sh
# Watch logs/rnaseq_test_<jobid>.out. If it completes, the real run is good to go.
# If it fails on a QOS/account error, fix clusterOptions in sickkids_slurm.config.

set -euo pipefail
BASE=/hpf/projects/msalter/sghazis/rnaseq_TUY35595
cd "$BASE"

export NXF_HOME="$BASE/nf_cache"
export NXF_WORK="$BASE/work_test"
export NXF_SINGULARITY_CACHEDIR="$BASE/singularity_cache"
export SINGULARITY_CACHEDIR="$BASE/singularity_cache/.singularity"
export APPTAINER_CACHEDIR="$SINGULARITY_CACHEDIR"
export NXF_OPTS='-Xms1g -Xmx4g'
mkdir -p "$NXF_HOME" "$NXF_WORK" "$NXF_SINGULARITY_CACHEDIR" "$SINGULARITY_CACHEDIR"

module load nf-core/4.2.0
RNASEQ_VER=$(cat "$BASE/.rnaseq_ver" 2>/dev/null || echo 3.19.0)

nextflow run nf-core/rnaseq \
    -r "$RNASEQ_VER" \
    -profile test,singularity \
    -c "$BASE/sickkids_slurm.config" \
    --outdir "$BASE/results_test" \
    -work-dir "$NXF_WORK"

echo "Smoke test OK. Clean up with: rm -rf $BASE/work_test $BASE/results_test"
