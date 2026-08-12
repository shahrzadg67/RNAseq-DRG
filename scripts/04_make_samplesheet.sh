#!/bin/bash
#
# Build nf-core/rnaseq samplesheet.csv from the (to-be) decompressed fastq.gz files.
# Sample name = filename up to the _S### token (e.g. NP_F1, Sham_M2, SNI_M).
# strandedness=auto lets Salmon infer it. One row per sample (single lane L001).
# Safe to run before decompression: paths are deterministic from the .ora names.

set -euo pipefail

BASE=/hpf/projects/msalter/sghazis/rnaseq_TUY35595
SRC=/hpf/projects/msalter/TUY35595.20260429/20260424_LH00403_0176_B23TKL3LT4
FQ="$BASE/fastq"
OUT="$BASE/samplesheet.csv"

echo "sample,fastq_1,fastq_2,strandedness" > "$OUT"

# Iterate R1 ora files; derive sample name and matching R2.
for r1 in $(ls -1 "$SRC"/*_R1_001.fastq.ora | sort); do
    bn=$(basename "$r1" .fastq.ora)              # e.g. NP_F1_S130_L001_R1_001
    sample=$(echo "$bn" | sed -E 's/_S[0-9]+_L[0-9]+_R1_001$//')   # -> NP_F1
    r1_gz="$FQ/${bn}.fastq.gz"
    r2_gz="$FQ/${bn/_R1_001/_R2_001}.fastq.gz"
    echo "${sample},${r1_gz},${r2_gz},auto" >> "$OUT"
done

echo "Wrote $OUT :"
cat "$OUT"
echo
echo "Sample count: $(($(wc -l < "$OUT") - 1))"
