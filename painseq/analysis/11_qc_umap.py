#!/usr/bin/env python3
"""
11_qc_umap.py — QC and UMAP embedding for the GSE154659 objects.

Reads the binary CSC export written by 10_export_for_scanpy.R, builds an AnnData,
runs standard QC + normalisation + HVG + PCA + neighbours + UMAP, scores each cell
on the 47-gene panel, and writes the per-cell results back out as CSV for R to plot
in the app's house style.

NOTE on QC: this is single-NUCLEUS data. Median %mito is ~0.4 and the max is <10,
so a mito filter is not meaningful here. The authors also pre-filtered — the observed
floors are 566 UMI / 501 genes with UMI capped at 14,999 — so our QC is confirmatory
and is expected to remove very few cells.

Clustering is NOT run: the atlas ships author-assigned cell-type labels, so we embed
and colour by those rather than inventing new ones.

Out: data/umap_<name>.csv, data/qc_summary_<name>.csv, data/<name>.h5ad
"""
import os, sys, time
import numpy as np
import pandas as pd
import scipy.sparse as sp
import scanpy as sc
import anndata as ad

PROJ = "/hpf/projects/msalter/sghazis/painseq"
DATA = os.path.join(PROJ, "data")

# QC thresholds. Deliberately permissive — the authors already filtered, and we do
# not want to silently discard the injury-induced populations, which are rarer and
# shallower than the neurons.
MIN_GENES, MIN_CELLS_PER_GENE, MAX_PCT_MITO = 200, 3, 20.0
N_HVG, N_PCS = 2000, 50


def log(*a):
    print(f"[{time.strftime('%H:%M:%S')}]", *a, flush=True)


def load(name):
    d = os.path.join(DATA, f"export_{name}")
    n_genes, n_cells = [int(v) for v in open(os.path.join(d, "shape.txt")).read().split()]
    i = np.fromfile(os.path.join(d, "indices.i32"), dtype=np.int32)
    p = np.fromfile(os.path.join(d, "p.i32"), dtype=np.int32)
    v = np.fromfile(os.path.join(d, "values.i32"), dtype=np.int32).astype(np.float32)
    # R gives genes x cells in CSC; AnnData wants cells x genes, so transpose to CSR
    X = sp.csc_matrix((v, i, p), shape=(n_genes, n_cells)).T.tocsr()
    genes = [l.rstrip("\n") for l in open(os.path.join(d, "genes.txt"))]
    obs = pd.read_csv(os.path.join(d, "obs.csv"))
    obs.index = obs["cell"].astype(str)
    a = ad.AnnData(X=X, obs=obs, var=pd.DataFrame(index=pd.Index(genes, name="gene")))
    log(f"{name}: {a.n_obs:,} cells x {a.n_vars:,} genes, nnz {X.nnz:,}")
    return a


def run(name):
    log(f"===== {name} =====")
    a = load(name)
    a.layers["counts"] = a.X.copy()

    # ---- QC ---------------------------------------------------------------
    a.var["mt"] = a.var_names.str.startswith("mt-")
    sc.pp.calculate_qc_metrics(a, qc_vars=["mt"], inplace=True, log1p=False, percent_top=None)
    before = a.n_obs
    qc_rows = [{
        "stage": "raw", "n_cells": a.n_obs, "n_genes": a.n_vars,
        "median_umi": float(np.median(a.obs["total_counts"])),
        "median_genes": float(np.median(a.obs["n_genes_by_counts"])),
        "median_pct_mt": float(np.median(a.obs["pct_counts_mt"])),
        "max_pct_mt": float(np.max(a.obs["pct_counts_mt"])),
    }]
    sc.pp.filter_cells(a, min_genes=MIN_GENES)
    a = a[a.obs["pct_counts_mt"] < MAX_PCT_MITO].copy()
    sc.pp.filter_genes(a, min_cells=MIN_CELLS_PER_GENE)
    log(f"  QC removed {before - a.n_obs} of {before} cells ({100*(before-a.n_obs)/before:.2f}%)"
        f" -> {a.n_obs:,} cells x {a.n_vars:,} genes")
    qc_rows.append({
        "stage": "filtered", "n_cells": a.n_obs, "n_genes": a.n_vars,
        "median_umi": float(np.median(a.obs["total_counts"])),
        "median_genes": float(np.median(a.obs["n_genes_by_counts"])),
        "median_pct_mt": float(np.median(a.obs["pct_counts_mt"])),
        "max_pct_mt": float(np.max(a.obs["pct_counts_mt"])),
    })
    pd.DataFrame(qc_rows).to_csv(os.path.join(DATA, f"qc_summary_{name}.csv"), index=False)

    # ---- normalise / HVG / PCA -------------------------------------------
    # seurat_v3 runs on RAW counts and must precede log1p. It fits a loess per
    # gene, which fails ("reciprocal condition number") when a batch is tiny —
    # several samples here have only ~280 cells — so it is attempted without a
    # batch_key and falls back to the log-space "seurat" flavor if it still dies.
    hvg_done = False
    try:
        sc.pp.highly_variable_genes(a, n_top_genes=N_HVG, flavor="seurat_v3")
        hvg_done = True
        log("  HVG: seurat_v3 on raw counts")
    except Exception as e:
        log(f"  HVG: seurat_v3 failed ({type(e).__name__}: {e}) -> falling back to flavor='seurat'")

    sc.pp.normalize_total(a, target_sum=1e4)
    sc.pp.log1p(a)
    a.raw = a
    if not hvg_done:
        sc.pp.highly_variable_genes(a, n_top_genes=N_HVG, flavor="seurat")
        log("  HVG: seurat flavor on log-normalised data")
    log(f"  {int(a.var['highly_variable'].sum())} highly variable genes")

    # ---- panel score, on log-normalised data, before scaling --------------
    panel = [g.strip() for g in open(os.path.join(DATA, "panel_symbols.txt")) if g.strip()]
    panel = [g for g in panel if g in a.var_names]
    sc.tl.score_genes(a, panel, score_name="panel_score")
    log(f"  scored {len(panel)} panel genes -> obs['panel_score']")

    a = a[:, a.var["highly_variable"]].copy()
    sc.pp.scale(a, max_value=10)
    sc.tl.pca(a, n_comps=N_PCS, svd_solver="arpack")

    # ---- batch integration: deliberately NOT applied -----------------------
    # `sample` is perfectly nested within injury x time — every sample belongs to
    # exactly one condition. Integrating on it would regress out the condition
    # differences along with any technical batch effect, i.e. delete the biology
    # this analysis exists to show. GEO supplies no independent technical batch
    # variable to correct on instead, so the embedding is left un-integrated and
    # the figures say so. (Harmony is installed and does converge on this data;
    # it is omitted on statistical grounds, not because it is unavailable.)
    rep = "X_pca"
    a.uns["integration"] = "none"
    log("  no batch integration: 'sample' is nested within injury x time, so correcting")
    log("    on it would remove the condition effect. Embedding is un-integrated by design.")

    # ---- neighbours + UMAP ------------------------------------------------
    log("  neighbors ...")
    sc.pp.neighbors(a, n_neighbors=15, n_pcs=N_PCS, use_rep=rep)
    log("  umap ...")
    sc.tl.umap(a, min_dist=0.3)

    out = a.obs.copy()
    out["UMAP1"] = a.obsm["X_umap"][:, 0]
    out["UMAP2"] = a.obsm["X_umap"][:, 1]
    out["integration"] = a.uns["integration"]
    keep = [c for c in ["cell", "sex", "geno", "inj", "time", "rep", "celltype", "sample",
                        "total_counts", "n_genes_by_counts", "pct_counts_mt",
                        "panel_score", "UMAP1", "UMAP2", "integration"] if c in out.columns]
    out[keep].to_csv(os.path.join(DATA, f"umap_{name}.csv"), index=False)
    log(f"  wrote data/umap_{name}.csv ({a.n_obs:,} cells)")

    # anndata refuses to write when the obs index NAME matches a column whose
    # values differ. Our index is the (non-unique) cell name and there is also a
    # "cell" column, so clear the index name before writing.
    a.obs.index.name = None
    a.write_h5ad(os.path.join(DATA, f"{name}.h5ad"), compression="gzip")
    log(f"  wrote data/{name}.h5ad")


if __name__ == "__main__":
    for nm in (sys.argv[1:] or ["c57", "atf3"]):
        run(nm)
    log("done")
