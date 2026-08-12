#!/bin/bash
#
# Download mouse GRCm39 genome + annotation from Ensembl.
# This SickKids login node has outbound internet to Ensembl, so it can run here
# directly (no data-mover SSH hop needed). It also works on data.ccm.sickkids.ca.
#     bash /hpf/projects/msalter/sghazis/rnaseq_TUY35595/scripts/02_get_reference.sh
#
# Release is read from Ensembl's VERSION file (the current_* dir isn't browsable
# over HTTPS) so the GTF name, which embeds the release number, is always correct.

set -euo pipefail

REFS=/hpf/projects/msalter/sghazis/rnaseq_TUY35595/refs
cd "$REFS"

echo "Detecting current Ensembl release ..."
REL="${ENSEMBL_REL:-$(curl -fsSL https://ftp.ensembl.org/pub/VERSION)}"
if ! [[ "$REL" =~ ^[0-9]+$ ]]; then
    echo "ERROR: could not detect Ensembl release (got '$REL')." >&2
    exit 1
fi
echo "Ensembl release: $REL"

BASEURL="https://ftp.ensembl.org/pub/release-${REL}"
FASTA_URL="$BASEURL/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna_sm.primary_assembly.fa.gz"
GTF_URL="$BASEURL/gtf/mus_musculus/Mus_musculus.GRCm39.${REL}.gtf.gz"

echo "Downloading genome FASTA ..."
wget -c -q --show-progress "$FASTA_URL"
echo "Downloading GTF ..."
wget -c -q --show-progress "$GTF_URL"

# Also grab Ensembl's md5 checksums for the record.
wget -qO CHECKSUMS.fasta "$BASEURL/fasta/mus_musculus/dna/CHECKSUMS" || true
wget -qO CHECKSUMS.gtf    "$BASEURL/gtf/mus_musculus/CHECKSUMS"        || true

echo
echo "Downloaded into $REFS :"
ls -lh "$REFS"/*.gz
echo
echo "Record this for the run command:"
echo "  FASTA = $REFS/Mus_musculus.GRCm39.dna_sm.primary_assembly.fa.gz"
echo "  GTF   = $REFS/Mus_musculus.GRCm39.${REL}.gtf.gz"
echo "$REL" > "$REFS/.ensembl_release"
