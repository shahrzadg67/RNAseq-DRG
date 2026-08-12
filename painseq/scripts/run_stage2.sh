#!/usr/bin/env bash
# Stage 2 — export the count matrices from R, then QC + UMAP in scanpy.
#   bash scripts/run_stage2.sh
#
# scanpy 1.11.1 exists ONLY in the alma8 python tree; `module load python/3.13.1`
# resolves to alma9, which has no scanpy. The alma8 binary also needs libffi 3.2.1
# on LD_LIBRARY_PATH when run from an alma9 node, or it dies with
#   ImportError: libffi.so.6: cannot open shared object file
set -euo pipefail
PROJ=/hpf/projects/msalter/sghazis/painseq
APP=/hpf/projects/msalter/sghazis/rnaseq_TUY35595
# Venv built with --system-site-packages off the alma8 python, so scanpy/anndata
# are inherited and only the two missing pieces were installed into it:
#   scikit-misc  -> required by highly_variable_genes(flavor="seurat_v3")
#   harmonypy    -> batch integration across samples
# Rebuild with:
#   /hpf/tools/alma8/python/3.13.1/bin/python3 -m venv --system-site-packages .venv-sc
#   ./.venv-sc/bin/python -m pip install scikit-misc harmonypy
PY="$PROJ/.venv-sc/bin/python"
cd "$PROJ"

# skip the export if it is already on disk (it takes ~1 min and never changes)
if [ -f data/export_c57/values.i32 ] && [ -f data/export_atf3/values.i32 ] && [ "${1:-}" != "force" ]; then
  echo "=== export already present, skipping 10_export_for_scanpy.R (pass 'force' to redo) ==="
else
echo "=== 10_export_for_scanpy.R ==="
module load R/4.5.2 libarchive/3.8.1 2>/dev/null || true
export R_LIBS_USER="$APP/.Rlib"
export LD_LIBRARY_PATH="$APP/.extralib:${LD_LIBRARY_PATH:-}"
Rscript analysis/10_export_for_scanpy.R 2>&1 | tee logs/10_export_for_scanpy.log
fi

echo "=== 11_qc_umap.py ==="
module load libffi/3.2.1 2>/dev/null || true
export LD_LIBRARY_PATH=/hpf/tools/alma8/python/3.13.1/lib:/hpf/tools/alma9/libffi/3.2.1/lib64:$LD_LIBRARY_PATH
# keep caches off /home — it is at ~8.6G against a ~10G cap
export XDG_CACHE_HOME="$PROJ/.cache"
export MPLCONFIGDIR="$PROJ/.cache/mpl"
export NUMBA_CACHE_DIR="$PROJ/.cache/numba"
mkdir -p "$XDG_CACHE_HOME" "$MPLCONFIGDIR" "$NUMBA_CACHE_DIR"
"$PY" analysis/11_qc_umap.py 2>&1 | tee logs/11_qc_umap.log

echo "=== stage 2 complete ==="
