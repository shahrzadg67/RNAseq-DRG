#!/usr/bin/env bash
#SBATCH --job-name=noxy_pipe
#SBATCH --account=msalter
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=03:00:00
#SBATCH --output=/hpf/projects/msalter/sghazis/rnaseq_TUY35595/logs/noxy_%j.out
#SBATCH --error=/hpf/projects/msalter/sghazis/rnaseq_TUY35595/logs/noxy_%j.err
# Batch-corrected + sex-chromosome-removed pipeline (DE + PCA + GSEA), on a
# separate node with 32G so it can never OOM-kill the OnDemand code-server session.
set -euo pipefail
ROOT=/hpf/projects/msalter/sghazis/rnaseq_TUY35595
module load R/4.5.2 libarchive/3.8.1 2>/dev/null || true
export R_LIBS_USER="$ROOT/.Rlib"
# OpenSSL 1.0.2 shim so org.Mm.eg.db / RSQLite load on any node (libssl.so.10)
export LD_LIBRARY_PATH="$ROOT/.extralib:${LD_LIBRARY_PATH:-}"
cd "$ROOT"
echo "START $(date) on $(hostname)"
Rscript analysis/run_all_noxy.R
echo "NOXY_PIPELINE_DONE $(date)"
