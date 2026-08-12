# painseq — the 47-gene conserved injury panel in the Renthal DRG injury atlas

Testing the TUY35595 **47-gene conserved injury panel** (SNI-positive / MINP-negative) against
single-cell data: which cells express it, when it peaks, and whether it is specific to axotomy.

**Dataset:** Renthal et al. 2020, mouse DRG snRNA-seq (GEO **GSE154659**).

| Object | Dim | Design |
|---|---|---|
| C57 atlas | 25,105 genes × 141,093 cells, 66 samples | 6 injury models × 12 timepoints, 20 cell types |
| Atf3 WT/KO | 24,597 × 17,665, 12 samples | Atf3-WT vs Atf3-KO, Naive / Crush 36h,168h, neurons only |

## Headline result

The panel is **axotomy-specific**, exactly as the bulk cross-species work predicted.

| | axotomy (Crush/ScNT/SpNT, 24h–7d) | non-axotomy (CFA 48h, Paclitaxel 7d) |
|---|---|---|
| mean log2FC | **+1.82** | +0.18 |
| panel genes significantly up | **59.5%** | **0 of 45** |

Single-cell SpNT-7d log2FC correlates with the app's bulk mouse `SNI_vs_Sham` at **r = 0.815**
(Spearman 0.733), 95% same-direction — independent confirmation of the bulk signature.

Two further results close out the story:

- **MINP is not hiding in a cell type.** Correlating the bulk log2FC vectors against every
  single-cell injury signature (384 cell type × injury × timepoint combinations, all shared genes):
  bulk SNI reaches **r = 0.499** with 144/384 combinations above r = 0.2; bulk MINP never exceeds
  **r = 0.152** and clears 0.2 **nowhere**. The bulk null result is not a dilution artefact.
- **The injury-state deconvolution reference is worse, not better** — a negative result.
  A reference built from injured cells reports **0.000%** repair cells in every mouse sample,
  including SNI. The naive-only reference puts `Repair schwann_N` at **1.56% in SNI** against
  0.31% in Sham and 0.40% in MINP — a ~5× elevation in the one condition that is a true axotomy.
  MuSiC's weighted NNLS drives the rare, extreme injury signatures to exactly zero while other
  types absorb the signal (NF3 2.09→6.68%, PEP1 2.88→0.03%). **Keep the naive reference.** The
  caveat to carry is that absolute proportions are sensitive to reference construction and should
  be read as relative comparisons across samples, not as absolute tissue composition.

## Running

```bash
bash scripts/run_analysis.sh all     # Stage 1: 01 | 02 | 03 | 04 | 05
bash scripts/run_stage2.sh           # Stage 2: export + scanpy QC/UMAP (add 'force' to re-export)
```

`run_analysis.sh` loads `R/4.5.2 libarchive/3.8.1`, points `R_LIBS_USER` at the TUY35595 `.Rlib`, and
adds the `.extralib` OpenSSL shim (R's `curl` needs `libssl.so.10` here). Everything Stage 1 needs is
already installed.

Stage 2 uses `.venv-sc`, built with `--system-site-packages` off the **alma8** python 3.13.1 so
scanpy is inherited; only `scikit-misc` (for `seurat_v3` HVG) and `harmonypy` were installed into it.
Rebuild with:

```bash
/hpf/tools/alma8/python/3.13.1/bin/python3 -m venv --system-site-packages .venv-sc
./.venv-sc/bin/python -m pip install scikit-misc harmonypy
```

## Scripts

| Script | Output |
|---|---|
| `00_common.R` | paths, barcode parser, gene panel, house styling — sourced by the rest |
| `01_parse_meta.R` | `data/cell_meta{,_atf3}.rds`, `data/panel_genes.rds`, `data/design_summary.txt` |
| `02_pseudobulk.R` | `data/pb{,_atf3}_sample.rds`, `data/pb{,_atf3}_sample_celltype.rds` |
| `03_panel_timecourse.R` | Panels A/B + male-C57 sensitivity, `data/panel_timecourse.rds` |
| `04_panel_celltype.R` | Figs C/D/E, `data/panel_celltype_{fits.rds,peak.csv}` |
| `05_verify.R` | `data/verification_report.txt` |
| `10_export_for_scanpy.R` | `data/export_{c57,atf3}/` binary CSC triplets + `obs.csv` |
| `11_qc_umap.py` | `data/umap_{c57,atf3}.csv`, `data/qc_summary_*.csv`, `data/{c57,atf3}.h5ad` |
| `13_umap_figs.R` | `figs/umap_*`, `figs/qc_*` |
| `20_build_app_artifact.R` | `rnaseq_TUY35595/app/data/sc_injury.rds` (the Shiny tab's data) |
| `30_deconv_injury_ref.R` | `data/deconv_compare.rds`, `figs/deconv_naive_vs_injury.*` |
| `31_minp_celltype.R` | `data/minp_celltype_concordance.*`, `figs/minp_celltype_concordance.*` |
| `32_atf3_dependency.R` | `data/atf3_dependency.*`, `figs/atf3_dependency.*` |
| `33_deconv_report.R` | `data/deconv_report.txt`, `figs/deconv_repair_by_group.*` (re-runnable without redoing MuSiC) |

## In the Shiny app

Stage 16 `sc_injury` — **"Single-cell injury panel"** tab in `rnaseq_TUY35595`
(`app/R/tab_sc_injury.R`; data `app/data/sc_injury.rds` 0.32 MB on disk / 4.4 MB resident, plus the
lazily-loaded `app/data/sc_percell.rds` 1.01 MB / 12.4 MB — see Stage 4 below). Sub-panels:
Injury × time, By cell type, Cell type × time, Gene explorer, Distributions, Recruitment,
Injured state, UMAP, MINP specificity, Atf3 dependence, Deconvolution reference, Peak landscape,
Methods & caveats.
Every heatmap has independent **row and column clustering** toggles; clustering the columns
suppresses the injury-model block separators, since the columns are then no longer in
biological order.

Rebuild the tab's data after re-running any analysis:

```bash
Rscript analysis/20_build_app_artifact.R
# then, before deploying, sync the STAGED copy — deploy/deploy_shinyapps.sh ships
# deploy/app_deploy/, not app/, and nothing syncs them automatically:
cd /hpf/projects/msalter/sghazis/rnaseq_TUY35595
Rscript deploy/make_build_config.R prod
rsync -a --delete --exclude 'rsconnect/' --exclude '*_bc.rds' \
      --exclude 'www/gsea_bc/' --exclude 'www/pathview_bc/' app/ deploy/app_deploy/
```

## Figures

| File | What |
|---|---|
| `figs/panelA_log2fc_timecourse` | **primary** — 49 genes × 25 injury×time, log2FC vs Naive |
| `figs/panelB_zscore_timecourse` | same layout, row z-scored logCPM |
| `figs/panelA_log2fc_timecourse_maleC57` | sensitivity: male C57 only (r = 0.996 vs all cells) |
| `figs/figC_panel_by_celltype` | 49 genes × 18 cell types, axotomy peak vs Naive |
| `figs/figD_panel_by_neuron_subtype` | the 9 sensory-neuron subtypes |
| `figs/figE_panel_celltype_timecourse` | small multiples, genes × time per cell type (SpNT) |
| `figs/umap_c57_{celltype,injury,time,genotype,panelscore}` | UMAP of all 141,093 cells |
| `figs/umap_atf3_*` | UMAP of the Atf3-WT/KO arm |
| `figs/qc_{c57,atf3}` | per-sample QC violins |
| `figs/minp_celltype_concordance` | MINP vs SNI concordance across cell type × injury × time |
| `figs/atf3_dependency` | crush induction in Atf3-WT vs Atf3-KO |
| `figs/deconv_naive_vs_injury` | naive-only vs injury-state deconvolution reference |

Styling mirrors the TUY35595 app (`app/R/helpers.R`): `#2166AC`→`#F7F7F7`→`#B2182B` diverging scale
centred at 0, `theme_pub()` / `ply_pub()`.

## Things that bite

- **The `.RDS.gz` files in this directory are double-gzipped** (gzip around an already-compressed
  RDS). Plain `readRDS` fails with `unknown input format`; use `readRDS(gzcon(gzfile(f, "rb")))`, or
  the `read_counts()` helper in `00_common.R`. The C57 file is also already **decompressed** at
  `rnaseq_TUY35595/refs/scref/GSE154659_C57_Raw_counts.RDS` — read that one. Only the Atf3 WT/KO
  file is genuinely new; the C57 `.gz` here duplicates 396 MB that already existed.
- **Cell-type parsing.** The type is fields 6..n-1 joined by `_`, **not** the second-to-last field.
  Four types contain a separator: `p_cLTMR2`, `Schwann_M`, `Schwann_N`, `Repair schwann_N`. The
  second-to-last-field shortcut in `rnaseq_TUY35595/analysis/80_deconv.R:60` truncates them to
  `cLTMR2` / `M` / `N` / `schwann`, and that mangled nomenclature is what ended up in the app's
  `CELLTYPE_FULL` and `NEURONAL_TYPES` (`app/R/helpers.R:288-310`). `00_common.R` carries the
  corrected 20-type table.
- **`Repair schwann_N` and `Repair fibroblast` have no Naive baseline** — they are injury-induced and
  do not exist in naive tissue, so no "vs Naive" contrast is definable and they are absent from
  Figs C/D/E. This is precisely why a naive-only deconvolution reference cannot see them.
- **Gene symbols are 2020-era.** `Insyn2a` is present only as `Fam196a`; `Crybg1` (old name `Aim1`)
  is absent entirely → 46/47 resolvable, 45 surviving the expression filter.
- **Ambient RNA.** Neuron-restricted transcripts (Npy, Gal) show apparent induction in every cell
  type including B cells and neutrophils. Real neuronal enrichment exists (mean panel log2FC 2.14
  neuronal vs 1.41 non-neuronal; 30/33 genes higher in neurons), but non-neuronal values are
  inflated by ambient leakage. **Read Fig C for relative ranking across cell types, not as absolute
  per-type induction.**
- **Shallow depth.** Median 846 genes/cell. Ten panel genes are detected in <0.5% of cells and are
  marked `*` on every figure; `Pappa2` (29 UMI total), `Ano7` (69) and `Crybg1` are not usable.
- **`Sprr1a` and `Sox11` are not in the 47** — Sprr1a only because it lacks a macaque ortholog
  (`app/R/tab_cross_pca.R:156-157`). Both are carried as `‡`-marked reference rows.
- **Cell names are not unique.** 5,593 of the 141,093 C57 names collide, because two nuclei from the
  same sample and cell type can draw the same barcode suffix. All aggregation here is positional, so
  nothing is affected, but never join on the cell name — `10_export_for_scanpy.R` emits `cell_uid`
  for that.
- **The UMAP is not batch-integrated, deliberately.** `sample` is perfectly nested within
  injury × time, so integrating on it would regress out the condition effect along with any technical
  batch effect. Harmony is installed and converges on this data; it is omitted on statistical grounds,
  not availability. GEO provides no independent technical batch variable to correct on instead.
- **The Atf3 knockout cannot be verified from these data.** Crush induces `Atf3` mRNA *more* in the
  KO (+5.99) than in the WT (+4.73) — the allele still produces a transcript that 3'-end
  single-nucleus sequencing counts. Any Atf3-dependence result is conditional on the published
  genotype labels. That arm also has only 3 pseudobulk samples per genotype per state, so the
  interaction test is badly underpowered: read the retention percentages, not the p-values.
- **scanpy lives only in the alma8 tree.** `module load python/3.13.1` resolves to alma9, which has
  no scanpy, and the alma8 binary needs `libffi/3.2.1` on `LD_LIBRARY_PATH` when run from an alma9
  node. `scripts/run_stage2.sh` handles both.
- **`seurat_v3` HVG fails with a `batch_key`** here (`reciprocal condition number` from the per-batch
  loess — several samples have only ~280 cells). It is run without one, with a fallback to the
  log-space `seurat` flavor.

---

## Stage 4 — per-cell analyses (beyond pseudobulk)

Pseudobulk collapses each sample to a mean and discards the distribution across cells. That
distribution carries the mechanism.

### Injury recruits cells, it does not upregulate them

With `f` = fraction of cells expressing a gene and `m` = mean among those cells, the mean over all
cells is `f × m`, so the log2 fold-change splits **exactly** into a recruitment term and a level term.

```
Atf3 in NF1 neurons, SpNT
  time    % cells expressing    mean among expressing
  naive          0.2%                   2.50
  24h           53.0%                   2.29
  48h           63.1%                   2.28
  7d            42.9%                   1.82
```

**32 of 47 panel genes are recruitment-driven; of the 32 induced more than two-fold, _none_ is
level-driven.** Median 81% of the fold-change is recruitment. Pseudobulk reports this as
"log2FC ≈ +5" and cannot distinguish it from every cell rising uniformly.

### The panel is a single-cell classifier of the injured state

Thresholded at the 99th percentile of naive, the per-cell panel score separates naive from injured
sensory neurons with **AUC 0.978** (SpNT), 0.752 (ScNT), 0.729 (Crush) — and **near chance for the
non-axotomy insults**, 0.568 (CFA) and 0.542 (paclitaxel). Peak fraction of neurons in the injured
state: SpNT 98.3%, ScNT 70.5%, Crush 63.9% versus CFA 1.9% and paclitaxel 12.5%, on a naive baseline
of 0.8%. The Stage 1 specificity result holds cell by cell, not just in group means.

Recruitment order is consistent: **PEP1 activates by 6 h, NF1/NF2 lag to 12–24 h, satellite glia
never exceed 39% and Schwann cells 17%** — independent confirmation that the non-neuronal signal in
Fig C is ambient RNA.

### The panel marks loss of subtype identity

Scoring each neuron on **its own** subtype's naive-derived markers: **r = −0.512 after axotomy**
versus −0.125 in naive, with mean identity dropping **0.60 SD**. Panel-high cells are exactly the
cells that have drifted furthest from what they were — the convergence Renthal et al. describe.
Crush *recovers* identity by 60–90 days; the transection models do not.

### Atf3 loss blocks entry into the injured state

`32_atf3_dependency.R` analysed this arm as whole-neuron pseudobulk and found it underpowered —
3 samples per genotype per state, almost nothing significant. That was the **wrong test**. Injury
induction is recruitment, which is a per-cell quantity estimated from 17,665 neurons:

| | Atf3-WT | Atf3-KO | difference |
|---|---|---|---|
| naive | 1.0% | 0.3% | — |
| crush 36 h | **49.8%** | **23.1%** | −26.7 pp |
| crush 7 d | **48.7%** | **11.9%** | −36.9 pp |

Knockout neurons partially initiate the response but cannot sustain it (83% of WT recruitment
retained at 36 h, 59% at 7 d). Most Atf3-blocked: Sprr1a, Npy, Cryba2, Sez6l, Gal, Mmp16, Csf1, Fosl1.

**Two caveats stand.** Atf3 mRNA is *not* reduced in the KO (crush induces it in 34% of KO cells vs
38% of WT) — the allele still yields a transcript that 3'-end snRNA-seq counts, so the genotype
cannot be verified from these data. And a proportion test over thousands of cells from 1–3 mice
overstates significance; read the effect sizes.

### The MINP-specific gene set is fragile

Only **5 of the 17** MINP-associated rows carry a gene symbol; the rest are unannotated Ensembl loci
carrying the most extreme fold-changes with non-significant pooled p-values. Of the two the app calls
"core": **Myh7 is detected in 0.185% of DRG nuclei (270 UMI across all 141,093 cells)** and Slc15a2
in 1.5%. Neither is injury-responsive. Treat MINP-specific calls as provisional pending a targeted assay.

### Peak landscape

The peak table plotted against time is a temporal wave: 7 genes peak within 24 h (Socs3, Acp7, Jun,
Gadd45a, Itga7), 23 between 36 and 72 h (Atf3, Sema6a, Flrt3, Adcyap1, Gal, Ecel1), 13 at one week
(Npy, Sprr1a, Stmn4, Cdkn1a, Gap43). 33 of 44 peak in **SpNT** — the model with no regenerative outlet.

### Scripts and figures

| Script | Output |
|---|---|
| `40_extract_percell.R` | `data/percell_expr{,_atf3}.rds`, `percell_app.rds`, `gene_universe.rds`, `subtype_markers.rds` |
| `41_recruitment.R` | `data/recruitment*.{rds,csv}`, `figs/recruitment_{heatmap,scatter}.*` |
| `42_classifier.R` | `data/classifier*.{rds,csv}`, `figs/classifier_{roc,curves}.*` |
| `43_identity.R` | `data/identity*.{rds,csv}`, `figs/identity_{scatter,trajectory}.*` |
| `44_minp_genes.R` | `data/minp_genes_sc.{rds,csv}`, `figs/minp_genes_sc.*` |
| `45_atf3_percell.R` | `data/atf3_percell*.{rds,csv}`, `figs/atf3_{entry,recruitment}.*` |
| `46_peak_viz.R` | `data/peak_viz.rds`, `figs/peak_landscape.*` |

### In the app

New sub-panels: **Gene explorer** (UMAP feature plot + ridge/violin/box distribution, grouped by cell
type / injury / injury×time), **Distributions** (panel-score ridges over time), **Recruitment**,
**Injured state** (classifier + subtype identity), **Peak landscape** (plot / counts / table), plus
MINP-genes and Atf3-per-cell sub-panels.

Per-cell data lives in a **separate** `app/data/sc_percell.rds` (1.01 MB on disk, 12.4 MB resident,
21,209 cells × 201 genes). `load_rds()` is memoised per file, so it is read only when a per-cell plot
first renders — verified: with no output rendered, only `sc_injury.rds` is in the cache. Dropdowns are
populated from `sc_injury.rds` precisely so that opening the app does not drag the matrix into memory.
`app/data/*.rds` totals 90.5 MB (was 89.4).

**Ridges are drawn over expressing cells only**, with % expressing printed beside each group label —
panel genes are so zero-inflated that a density over all cells is one spike at zero. Violin and box
views show the zeros directly. The in-app ridges are hand-built from `toself`-filled plotly polygons
rather than `ggridges`, because the app also publishes as a shinylive/WebR static site where an extra
package would need a wasm build; `ggridges` is available on the cluster for static figures.

---

## Stage 5 — is MINP doing something the injury axis cannot see?

Every earlier test asked "does MINP look like nerve injury?" and answered no. That leaves a fair
objection: MINP might be doing something real that simply isn't an injury programme — a glial
reaction, an immune infiltrate, or a change confined to a cell type too rare to move whole-DRG bulk.

`47_minp_beyond_injury.R` tests exactly that, asking by competitive gene-set test whether the bulk
MINP contrast is enriched for any **cell type's marker programme** (50 markers per type, derived from
naive atlas cells).

**Taken at face value the answer looks like yes** — 43 cell-type sets reach adj.p < 0.05 across the
MINP contrasts, with all non-neuronal markers moving one way and all neuronal markers the other.

**It isn't real.** The male and female MINP contrasts are almost perfectly **anti-correlated**
(**r = −0.940**): in males the non-neuronal markers fall and neuronal rise; in females the exact
reverse. A genuine MINP effect cannot be anti-correlated with itself. This is compositional wobble —
how much nerve and connective tissue each dissection happened to capture alongside the ganglion,
differing by chance between groups of three animals. A competitive test aggregates ~50 markers at
once, so it picks that up easily even when no individual gene is significant. It is also consistent
with the app's Deconvolution tab, which finds no significant MINP-vs-Sham shift in any cell type.

**The method works — the positive control proves it.** In SNI both injury reference programmes come
up strongly: the atlas-derived injury programme at p ≈ 2×10⁻³¹ and the conserved 47-gene panel at
p ≈ 2×10⁻¹². So the absence of a coherent MINP signal is a real absence, not a failure to look.

**Method note.** Ranking uses the DESeq2 **Wald statistic** from the full-precision `DE_results/*.csv`
tables, not `app/data/de_long.rds` — that table retains only ~1,684 distinct fold-change values across
20,830 genes, leaving ~92% of the ranking tied, which fgsea correctly refused to trust. The test is
`limma::cameraPR` (competitive, tie-tolerant, corrects for inter-gene correlation, no permutation —
fgsea's BiocParallel backend also cannot open a port on this cluster). SNI has no full-precision
table, so its ranking stays coarse; that costs power only, so its positive result stands.

| Script | Output |
|---|---|
| `47_minp_beyond_injury.R` | `data/minp_beyond.{rds,csv}`, `data/celltype_gene_sets.rds`, `figs/minp_beyond_{injury,consistency}.*` |

In the app: **MINP specificity → Beyond the injury axis**.

## Recruitment tab — reading REC and LVL

The tab now explains the two quantities inline rather than assuming them:

- **REC (recruitment)** — log2 change in the **fraction of cells** with any detectable transcript
- **LVL (level)** — log2 change in the **mean amount per expressing cell** (zeros excluded)

Since the mean over all cells is `f × m`, **REC + LVL = the pseudobulk log2 fold-change, exactly**.
The tab carries the worked Atf3/NF1 example, a "one gene in detail" view plotting the two components
as separate time-courses, and the dropout rebuttal: Npy in NF1 at 7 days averages 10.3 transcripts
per cell — which predicts ~100% detection — yet 41% of those neurons have *exactly zero* while the
rest carry ~18 copies each. That is genuine on/off bimodality, not a detection artefact.
