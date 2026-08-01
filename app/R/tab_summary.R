# Summary of Results tab — methodology, packages, the batch story, and the
# main findings with interpretation. Static content (no reactives).
tab_summary_ui <- function(id) {
  ns <- NS(id)
  nav_panel(
    title = "Summary of Results", icon = icon("clipboard-list"),
    div(class = "container py-3", style = "max-width: 950px;",

      h3("Summary of results"),
      div(class = "banner",
          "Bulk RNA-seq of mouse DRG/nerve comparing ", tags$b("MINP"), " (nucleus pulposus) and ",
          tags$b("Sham"), " (plus SNI for QC), in males and females. This page summarises the ",
          "methodology, why a batch effect was corrected, and what the results mean."),

      # ---- Methodology ----
      h4("Methodology & software"),
      tags$ul(
        tags$li(tags$b("Primary pipeline: "), "nf-core/rnaseq v3.19.0 — Trim Galore (QC/trimming) → ",
                "STAR (alignment) → Salmon (quantification); mouse ", tags$b("GRCm39, Ensembl 116"),
                ". Gene- and transcript-level count matrices analysed in parallel."),
        tags$li(tags$b("Differential expression: "), tags$b("DESeq2"), " and ",
                tags$b("edgeR / limma-voom"), " run side-by-side (TMM + voom for edgeR); five contrasts ",
                "across the MINP/Sham × sex design; Sham as reference."),
        tags$li(tags$b("Batch correction: "), "replicate (processing cohort) modelled as an additive ",
                "covariate in DE (", tags$code("~ 0 + group + replicate"), "); for PCA and expression ",
                "boxplots the replicate effect is removed with ", tags$b("limma::removeBatchEffect"),
                " (condition and sex preserved)."),
        tags$li(tags$b("Sex chromosomes: "), "chrX/chrY genes removed in the corrected analysis so sex ",
                "does not dominate; genes annotated (symbol, description, chromosome) via ",
                tags$b("org.Mm.eg.db"), "."),
        tags$li(tags$b("Enrichment: "), tags$b("clusterProfiler"), " GSEA over GO (BP/CC/MF) and KEGG ",
                "(", tags$code("gseGO"), " / ", tags$code("gseKEGG"), "), ranked by the DESeq2 statistic."),
        tags$li(tags$b("Visualisation: "), tags$b("pheatmap"), " (clustered heatmaps), ",
                tags$b("ggVennDiagram / UpSetR"), " (DEG overlaps), ", tags$b("ggplot2 + ggrepel"),
                " (labelled volcanoes); interactive app in ",
                tags$b("shiny / bslib / plotly / DT"), "."),
        tags$li(tags$b("Cell-type deconvolution: "), tags$b("MuSiC"), " estimates each bulk sample's ",
                "cell-type composition from a mouse DRG single-cell reference (Renthal et al. 2020, ",
                "GSE154659; naive cells, 18 cell types); proportions and per-contrast comparisons ",
                "are precomputed for the app.")
      ),

      # ---- Batch effect ----
      h4("The batch effect — what it was and why we removed it"),
      div(class = "banner",
          "The three biological replicates are three ", tags$b("processing cohorts"),
          " operated and extracted on different dates (cohort 3 ~3 weeks after cohorts 1–2), ",
          "documented in the study's sequencing-planning sheet (see the Methods tab)."),
      tags$ul(
        tags$li("In the uncorrected PCA, ", tags$b("cohort/replicate explained ~60% of PC1"),
                " — the single largest axis of variance — while the MINP-vs-Sham condition explained ~0%. ",
                "That is a batch effect, not biology."),
        tags$li("Removing chrX/Y confirmed it: the ", tags$b("cohort split on PC1 stayed strong (R² 0.60 → 0.68)"),
                " while the sex axis collapsed (PC2 sex R² 0.75 → 0.06). So the male cohort-2/3 clustering ",
                "was the processing batch (plus sex), not a new effect."),
        tags$li(tags$b("Correction: "), "replicate is adjusted for in the DE models and removed from the ",
                "PCA/boxplot visualisations; ", tags$b("SNI is excluded"), " from DE (a separate operator/",
                "injury model with n=1 per sex).")
      ),

      # ---- Main findings ----
      h4("Main findings & what they mean"),
      tags$ol(
        tags$li(tags$b("Sex is the strongest biological signal. "),
                "Sham-male vs Sham-female is a clean, reproducible difference (all three replicates per ",
                "sex separate with no overlap) — ~256 gene-level DEGs of genuine ", tags$b("autosomal sexual dimorphism"),
                " (e.g. Sst ↑ in females; Reln, Sptbn5, Sema5a ↑ in males). It is ", tags$b("not driven by an outlier"),
                " and serves as a positive control that the pipeline detects real differences."),
        tags$li(tags$b("The sex effect is validated — and is not a sequencing-depth artefact. "),
                "Sham females were sequenced ~1.67× deeper than males, but depth is removed by size-factor / ",
                "TMM normalisation (an offset, not a covariate). On the Sham F-vs-M comparison we measured ",
                tags$b("100% sensitivity"), " on gold-standard sex genes (chrY genes + ", tags$i("Xist"), "), ",
                tags$b("~99.98% specificity"), " by label permutation, and ", tags$b("94% DEG recovery"),
                " with log2FC concordance ", tags$b("r = 0.97"), " after down-sampling all libraries to equal ",
                "depth — so the differences are real biology. Details in the QC ", tags$b("“Depth & power”"), " panel."),
        tags$li(tags$b("MINP vs Sham is a modest, gene-specific effect. "),
                "The two conditions do not separate strongly in global PCA, but real DEGs emerge once the ",
                "cohort batch and sex are accounted for. ", tags$b("Batch correction increased power"),
                " — DESeq2 recovered roughly ", tags$b("2–3× more significant transcripts"),
                " after modelling replicate; edgeR is more conservative (it spends residual degrees of ",
                "freedom on the batch term)."),
        tags$li(tags$b("Interpretation. "),
                "The MINP-vs-Sham biology of interest is ", tags$b("subtle relative to sex and cohort batch"),
                ", so handling those two nuisances is essential to see it. Where DESeq2 and edgeR disagree ",
                "(transcript level), trust the genes both engines call. All thresholds and gene-set ",
                "selections here are ", tags$b("exploratory and pending experimental validation"), "."),
        tags$li(tags$b("Cell-type composition — what the tissue is made of, and why it doesn't explain the MINP effect. "),
                "Averaged across samples, MuSiC estimates ", tags$b("~54% neurons"), " — dominated by ",
                "neurofilament / myelinated A-fibre subtypes (", tags$b("NF2 ~19%, NF1 ~17%"), "), with ",
                "C-fibre low-threshold mechanoreceptors (~7%), non-peptidergic (~4%) and peptidergic (~3%) ",
                "nociceptors and somatostatin⁺ neurons (~1%) — alongside ", tags$b("satellite glia ~12%"), ", ",
                tags$b("Schwann cells ~13%"), " (myelinating + non-myelinating), ", tags$b("vascular ~11%"),
                " (endothelial + pericyte), ", tags$b("fibroblasts ~7%"), " and ", tags$b("immune ~4%"),
                " (macrophage, B cell, neutrophil) — a plausible DRG/nerve make-up. Aside from these proportions, ",
                tags$b("MINP vs Sham shows no significant shift in any cell type"), " (Wilcoxon adj.p ≈ 1; the ",
                "largest, non-significant differences are only a few percent), so the MINP-vs-Sham DEGs reflect ",
                tags$b("within-cell-type expression changes"), " rather than a change in cell proportions. With ",
                "n=3 per group the pooled contrast is best-powered; read differences alongside the effect sizes. ",
                "The Deconvolution tab shows the full per-sample and per-group composition.")
      ),

      h4(class = "mt-4", "Cross-species validation (rat & macaque DRG)"),
      div(class = "banner",
          "The mouse study was benchmarked against well-powered rat DRG (SNI / skin-lesion / naïve × sex × DRG level; ",
          "Ensembl GRCr8) and macaque DRG (ipsi / contra / control nerve transection × sex; Macaca_fascicularis_6.0), ",
          "each mapped to mouse via one-to-one orthologs. See the ",
          tags$b("Cross-species (rat / macaque)"), " tabs and the ",
          tags$b("Cross-species meta-analysis"), " tab, which holds the joint PCA, the cell-type ",
          "composition shifts, the conserved injury-marker panel and the ", tags$b("Gene lookup"),
          " sub-tab."),
      tags$ol(
        tags$li(tags$b("Nerve injury has a strongly conserved transcriptional program. "),
                "The classic regeneration-associated gene (RAG) response — ", tags$b("Atf3, Sox11, Sprr1a, Gal, Npy, Vgf, Gap43, Ecel1, Fosl1, Jun, Socs3"),
                " — is induced in SNI across mouse, rat and macaque (", tags$b("~279 genes concordant in all three species"),
                "; rat SNI 1,813 DEGs, macaque ipsi 3,892). Because the mouse SNI (n=1/sex) fold-changes track the ",
                "well-powered rat and macaque in direction, the underpowered mouse nerve-injury signal is ",
                tags$b("cross-species validated"), "."),
        tags$li(tags$b("MINP is distinct from nerve injury — “Sprr1a-negative”. "),
                "MINP-vs-Sham shows only ~9 DEGs (≈200× fewer than nerve injury, on the scale of a skin lesion), does ",
                tags$b("not"), " engage the RAG program (", tags$i("Sprr1a"), " −1.1, ", tags$i("Atf3"), " ≈ 0, both n.s.), and its ",
                "fold-changes are uncorrelated with the conserved injury signature (r ≈ 0). A ", tags$b("47-gene panel"),
                " (SNI-positive / MINP-negative; e.g. ", tags$i("Sprr1a, Atf3, Npy, Gal, Ecel1"),
                ") cleanly separates true nerve injury from MINP. MINP is also female-biased (16 DEGs in females, 0 in males)."),
        tags$li(tags$b("The sex difference is conserved but sex-chromosome-driven. "),
                "Baseline male-vs-female DRG differs in every species, but the ", tags$b("conserved"), " component is the ",
                "sex-chromosome set (", tags$i("Xist"), " in females; ", tags$i("Ddx3y, Uty, Kdm5d, Eif2s3y"),
                " in males). There is ", tags$b("no conserved autosomal sex signature and no conserved sex-specific injury response"),
                " — after SNI the male/female difference is essentially the same sex-chromosome genes as baseline."),
        tags$li(tags$b("A conserved cellular echo of injury. "),
                "Deconvolution (same 18 Renthal cell types in all species) shows a conserved compositional shift: ",
                tags$b("non-myelinating (repair) Schwann cells rise"), " — the classic Schwann repair/proliferation ",
                "response to axotomy — while ", tags$b("myelinated (NF2) neurons and satellite glia fall"),
                ", in the same direction across all three species (significant in the well-powered rat). ",
                "Overall the response is mainly a ", tags$b("within-cell-type program"),
                " (the RAG genes) rather than a large shift in cell proportions; in rat, injury is concentrated in the ",
                tags$b("lower (L4/5) DRG"), ".")
      ),

      div(class = "text-muted small mt-3",
          "Explore the tabs above for the uncorrected vs corrected PCA, differential expression, GSEA, ",
          "DEG overlaps (Venn) and expression heatmaps; every figure and table is downloadable.")
    )
  )
}

tab_summary_server <- function(id) {
  moduleServer(id, function(input, output, session) { })   # static content only
}
