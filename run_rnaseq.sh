#!/bin/bash
#SBATCH --job-name=rnaseq_ctl
#SBATCH --account=msalter
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=72:00:00
#SBATCH --output=/hpf/projects/msalter/sghazis/rnaseq_TUY35595/logs/rnaseq_ctl_%j.out
#SBATCH --error=/hpf/projects/msalter/sghazis/rnaseq_TUY35595/logs/rnaseq_ctl_%j.err
#
# nf-core/rnaseq controller. This job is just the Nextflow orchestrator (light:
# 2 CPU / 8 GB) -- the heavy STAR/Salmon work is submitted by Nextflow as its own
# Slurm jobs. Survives logout because it IS a Slurm job (Option A).
#
#   OPTION A (recommended): sbatch run_rnaseq.sh
#   OPTION B (tmux):        tmux new -s rnaseq ; then run the body below by hand,
#                           detach with Ctrl-b d, reattach with tmux attach -t rnaseq
#
# First do a smoke test (tiny built-in dataset) -- see scripts/05_smoke_test.sh.

set -euo pipefail

BASE=/hpf/projects/msalter/sghazis/rnaseq_TUY35595
cd "$BASE"

# --- keep ALL caches off the 10 GB home ------------------------------------
export NXF_HOME="$BASE/nf_cache"
export NXF_WORK="$BASE/work"
export NXF_SINGULARITY_CACHEDIR="$BASE/singularity_cache"
export SINGULARITY_CACHEDIR="$BASE/singularity_cache/.singularity"
export APPTAINER_CACHEDIR="$SINGULARITY_CACHEDIR"
export NXF_OPTS='-Xms1g -Xmx4g'
mkdir -p "$NXF_HOME" "$NXF_WORK" "$NXF_SINGULARITY_CACHEDIR" "$SINGULARITY_CACHEDIR"

module load nf-core/4.2.0    # nextflow + Singularity + python

# --- pinned versions (written by the pre-stage / reference scripts) --------
RNASEQ_VER=$(cat "$BASE/.rnaseq_ver" 2>/dev/null || echo 3.19.0)
REL=$(cat "$BASE/refs/.ensembl_release" 2>/dev/null || echo "")
FASTA="$BASE/refs/Mus_musculus.GRCm39.dna_sm.primary_assembly.fa.gz"
GTF="$BASE/refs/Mus_musculus.GRCm39.${REL}.gtf.gz"

echo "rnaseq=$RNASEQ_VER  ensembl=$REL"
echo "FASTA=$FASTA"
echo "GTF=$GTF"
[[ -s "$FASTA" ]] || { echo "Missing FASTA -- run 02_get_reference.sh first"; exit 1; }
[[ -s "$GTF"   ]] || { echo "Missing GTF -- run 02_get_reference.sh first";   exit 1; }

nextflow run nf-core/rnaseq \
    -r "$RNASEQ_VER" \
    -profile singularity \
    -c "$BASE/sickkids_slurm.config" \
    --input    "$BASE/samplesheet.csv" \
    --outdir   "$BASE/results" \
    --fasta    "$FASTA" \
    --gtf      "$GTF" \
    --aligner  star_salmon \
    --save_reference \
    -work-dir  "$NXF_WORK" \
    -resume

echo
echo "=== Run finished. Two count matrices: ==="
echo "  gene:       $BASE/results/star_salmon/salmon.merged.gene_counts.tsv"
echo "  transcript: $BASE/results/star_salmon/salmon.merged.transcript_counts.tsv"
echo
echo "If the run succeeded, reclaim space (minimize footprint) with:"
echo "  cd $BASE && nextflow clean -f -before $(date +%s) ; rm -rf $NXF_WORK"
