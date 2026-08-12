#!/usr/bin/env bash
#SBATCH --job-name=render_png
#SBATCH --account=msalter
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=00:40:00
#SBATCH --output=/hpf/projects/msalter/sghazis/rnaseq_TUY35595/logs/render_png_%j.out
#SBATCH --error=/hpf/projects/msalter/sghazis/rnaseq_TUY35595/logs/render_png_%j.err
# Render ALL PNGs (original + batch-corrected) on a SEPARATE node with its own
# 16G memory, so it can never OOM-kill the 2G OnDemand code-server session.
set -euo pipefail
ROOT=/hpf/projects/msalter/sghazis/rnaseq_TUY35595
module load R/4.5.2 libarchive/3.8.1 2>/dev/null || true
export R_LIBS_USER="$ROOT/.Rlib"
cd "$ROOT"
Rscript analysis/92_all_png.R
echo "RENDER_PNG_DONE"
