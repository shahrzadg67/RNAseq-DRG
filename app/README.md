# TUY35595 RNA-seq — Client App (shinylive)

Interactive, fully static report compiled with **shinylive** (runs in the browser via
WebAssembly; no R server). All heavy analysis is precomputed on the HPC; the app only
visualises the artifacts in `app/data/` and images in `app/www/`.

## Tabs
1. **Methods & Pipeline** — what was done + live per-stage release status (always visible).
2. **Quality Control** — MultiQC report + library size / density / TMM plots.
3. **PCA** — gene & transcript level, with/without SNI, elbow (scree), free PC-pair selector.
4. **Differential Expression** — DESeq2 vs edgeR side-by-side volcano (hover = gene/transcript
   name, log2FC, adj-p), adjustable cutoffs, concordance summary, results table, and a
   gene/transcript expression boxplot lookup with stats alongside.
5. **GSEA** — GO BP/CC/MF + KEGG; interactive dotplot + clickable table → running-score
   plot (and KEGG pathway map).

## Staged, gated release (staging vs prod)
Visibility is controlled by **`config/stage_status.yaml`**. Each stage has:
- `status`: `planned` | `in_progress` | `complete` (shown transparently in the Methods tab)
- `released`: `true|false` (whether the tab is visible on the **public/prod** build)

The **staging** build shows every tab (for your review) and is **local-only**. The **prod**
build hides any tab whose stage is not `released`, and is the only one published to GitHub.

### To release a stage to the client
1. Edit `config/stage_status.yaml`: set that stage `status: in_progress` and `released: true`.
2. `bash deploy/publish.sh prod --push`  → rebuilds prod and pushes to the public repo.

Nothing is shown to the client until you do this — data is released **on your command**, not
when the analysis happens to finish.

### Build / preview staging locally (your private review copy)
```bash
bash deploy/publish.sh staging      # full build -> deploy/site_staging
Rscript -e "httpuv::runStaticServer('deploy/site_staging', port=8008)"
# then forward port 8008 in VS Code to preview in your browser
```

### Publish prod to GitHub (PAT auth, authored as you)
```bash
export PROD_REPO=https://github.com/<USER>/tuy35595-rnaseq.git
export GITHUB_USER=<USER>
export GITHUB_TOKEN=<your fine-scoped PAT>      # never printed/committed
# (optional, defaults to your name/email)
export GIT_AUTHOR_NAME="Shahrzad Ghazisaeidi"
export GIT_AUTHOR_EMAIL="shahrzad67@gmail.com"
bash deploy/publish.sh prod --push
```
Commits are authored with `GIT_AUTHOR_NAME/EMAIL`; the PAT is used only for the push URL and
is masked in logs. Enable GitHub Pages (gh-pages branch) once in the repo's Settings.

## Regenerating artifacts (after pipeline finishes)
```bash
module load R/4.5.1
export R_LIBS_USER=/hpf/projects/msalter/sghazis/rnaseq_TUY35595/.Rlib
Rscript analysis/run_all.R          # DESeq2, edgeR, PCA, assemble, GSEA -> app/data, app/www
```
