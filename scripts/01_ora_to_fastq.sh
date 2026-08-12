#!/bin/bash
#SBATCH --job-name=ora2fq
#SBATCH --account=msalter
#SBATCH --array=0-27%8
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=/hpf/projects/msalter/sghazis/rnaseq_TUY35595/logs/ora2fq_%A_%a.out
#SBATCH --error=/hpf/projects/msalter/sghazis/rnaseq_TUY35595/logs/ora2fq_%A_%a.err
#
# Decompress Illumina ORA (.fastq.ora) -> .fastq.gz for nf-core/rnaseq.
# One array task per .ora file (28 files). %8 = at most 8 running at once.
# orad finds its reference via ORA_REF_PATH (set by the orad module).

set -euo pipefail

BASE=/hpf/projects/msalter/sghazis/rnaseq_TUY35595
SRC=/hpf/projects/msalter/TUY35595.20260429/20260424_LH00403_0176_B23TKL3LT4
OUT="$BASE/fastq"

module load orad/2.7.0

# Stable, sorted list of the .ora files; pick this task's file by array index.
mapfile -t FILES < <(ls -1 "$SRC"/*.fastq.ora | sort)
FILE="${FILES[$SLURM_ARRAY_TASK_ID]}"
NAME=$(basename "$FILE" .fastq.ora)

echo "[$(date)] task $SLURM_ARRAY_TASK_ID -> $FILE"

# Skip if already done (lets the array be re-run safely).
if [[ -s "$OUT/${NAME}.fastq.gz" ]]; then
    echo "Already present: $OUT/${NAME}.fastq.gz -- skipping."
    exit 0
fi

# -P writes the .gz into $OUT; -t threads; -f overwrite partials; --gz default lvl 5.
orad --gz -t "${SLURM_CPUS_PER_TASK}" -P "$OUT" -f "$FILE"

echo "[$(date)] done -> $OUT/${NAME}.fastq.gz"
ls -lh "$OUT/${NAME}.fastq.gz"
