#!/bin/bash
#
# Reclaim space AFTER a successful run (minimize-footprint choice).
# Keeps: results/, fastq/, refs/ (incl. saved indices). Removes Nextflow work dirs.
# Run only once you've confirmed results/star_salmon/*.merged.*_counts.tsv exist.

set -euo pipefail
BASE=/hpf/projects/msalter/sghazis/rnaseq_TUY35595
cd "$BASE"

GENE="$BASE/results/star_salmon/salmon.merged.gene_counts.tsv"
TX="$BASE/results/star_salmon/salmon.merged.transcript_counts.tsv"
if [[ ! -s "$GENE" || ! -s "$TX" ]]; then
    echo "Refusing to clean: expected count matrices not found."
    echo "  $GENE"
    echo "  $TX"
    exit 1
fi

echo "Both count matrices present. Removing Nextflow work dirs ..."
module load nf-core/4.2.0 2>/dev/null || true
nextflow clean -f 2>/dev/null || true
rm -rf "$BASE/work" "$BASE/work_test" "$BASE/results_test"

echo "Done. Remaining footprint:"
du -sh "$BASE"/{fastq,refs,results} 2>/dev/null
df -h /hpf/projects/msalter
