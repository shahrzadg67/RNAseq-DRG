#!/bin/bash
#
# Pre-pull the nf-core/rnaseq pipeline code + all Singularity containers into the
# project-space cache, so the real run never needs the container registry from a
# compute node. Run where BOTH the nf-core module AND internet are available:
# try the data-mover node first; if `module` isn't there, run on a login node.
#     bash /hpf/projects/msalter/sghazis/rnaseq_TUY35595/scripts/03_prestage_pipeline.sh
#
# Pin the release with: RNASEQ_VER=3.19.0 bash 03_prestage_pipeline.sh
# (Confirm the latest stable tag at https://github.com/nf-core/rnaseq/releases)

set -euo pipefail

BASE=/hpf/projects/msalter/sghazis/rnaseq_TUY35595
RNASEQ_VER="${RNASEQ_VER:-3.19.0}"

# Keep every cache OUT of the 10 GB home.
export NXF_HOME="$BASE/nf_cache"
export NXF_SINGULARITY_CACHEDIR="$BASE/singularity_cache"
export SINGULARITY_CACHEDIR="$BASE/singularity_cache/.singularity"
export APPTAINER_CACHEDIR="$SINGULARITY_CACHEDIR"
export NXF_OPTS='-Xms1g -Xmx4g'
mkdir -p "$NXF_HOME" "$NXF_SINGULARITY_CACHEDIR" "$SINGULARITY_CACHEDIR"

module load nf-core/4.2.0

echo "Pinned nf-core/rnaseq version: $RNASEQ_VER"
echo "$RNASEQ_VER" > "$BASE/.rnaseq_ver"

# nf-core download fetches pipeline code + all referenced Singularity images into
# the cache, ready for offline use on compute nodes.
nf-core pipelines download rnaseq \
    --revision "$RNASEQ_VER" \
    --compress none \
    --container-system singularity \
    --container-cache-utilisation amend \
    --outdir "$BASE/nfcore_rnaseq_${RNASEQ_VER}" \
    --download-configuration yes

echo
echo "Pipeline + containers staged."
echo "  NXF_SINGULARITY_CACHEDIR = $NXF_SINGULARITY_CACHEDIR"
du -sh "$NXF_SINGULARITY_CACHEDIR" 2>/dev/null || true
