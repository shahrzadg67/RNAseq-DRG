# QC tab — MultiQC report + custom QC plots ----------------------------------

# "Depth & power" panel — documents the sex/depth imbalance, the normalisation
# that removes it, and the sensitivity/specificity validation (Sham F vs M).
tab_qc_depth_ui <- function(ns) {
  sfrow <- function(s, sx, raw, sf, tmm)
    tags$tr(tags$td(s), tags$td(sx), tags$td(raw), tags$td(sf), tags$td(tmm))
  div(class = "py-2", style = "max-width: 900px;",
    div(class = "banner",
        "In the Sham group, females were sequenced ", tags$b("~1.67× deeper"), " than males ",
        "(F ≈ 53M, M ≈ 32M reads). Sequencing depth is removed by ", tags$b("size-factor / TMM"),
        " normalisation before any test — as a per-sample offset, not a design covariate — and we ",
        "verified the sex differences are real biology, not a depth artefact."),
    h5("Depth and normalisation (Sham samples)"),
    tags$table(class = "table table-sm table-bordered",
      tags$thead(tags$tr(lapply(c("Sample","Sex","Raw depth (M reads)","DESeq2 size factor","edgeR TMM factor"), tags$th))),
      tags$tbody(
        sfrow("Sham_F1","F","54.6","1.40","1.02"), sfrow("Sham_F2","F","51.8","1.32","1.02"),
        sfrow("Sham_F3","F","53.9","1.33","0.98"), sfrow("Sham_M1","M","30.8","0.79","1.05"),
        sfrow("Sham_M2","M","34.2","0.81","0.95"), sfrow("Sham_M3","M","30.9","0.76","0.99"))),
    tags$p(class = "text-muted small",
      "Size factors absorb the depth (F ≈ 1.35 vs M ≈ 0.79 ≈ the 1.67× ratio): after normalisation the ",
      "total counts F/M ≈ 0.98 (depth removed). TMM factors ≈ 1.0 mean there is ",
      tags$b("no composition bias"), " — the condition under which this normalisation is valid."),
    h5("Validation of the Sham F-vs-M sex effect"),
    tags$ul(
      tags$li(tags$b("Sensitivity = 100%. "), "All 6 gold-standard sex genes are detected as DEGs — ",
              "chrY ", tags$i("Ddx3y, Uty, Eif2s3y, Kdm5d, Uba1y"), " and ", tags$i("Xist"),
              " (huge effect sizes; ", tags$i("Xist"), " adj.p ≈ 1e-152)."),
      tags$li(tags$b("Specificity ≈ 99.98%. "), "Random 3-vs-3 relabellings of the six samples yield a ",
              tags$b("median of 4"), " spurious “DEGs” (~0.02% of genes) versus ", tags$b("31"),
              " for the true sex split. The one relabelling that produced many (303) is the split that ",
              "aligns with the ", tags$b("processing cohort"), " — real batch structure, not noise ",
              "(exactly what the batch correction handles)."),
      tags$li(tags$b("Not a depth artefact. "), "Down-sampling every library to equal depth (30.7M reads) ",
              "recovers ", tags$b("94% of DEGs"), " (29/31), with log2FC concordance ", tags$b("r = 0.97"),
              " and 100% gold-standard sensitivity retained — the sex signal is essentially unchanged by depth.")
    ),
    uiOutput(ns("depthfig")),
    tags$p(class = "text-muted small mt-2",
      "Note: unequal depth affects statistical ", tags$i("power"), " (deeper samples detect low-count genes ",
      "slightly better), which is modelled by the negative-binomial dispersion and independent filtering; ",
      "it does not create false positives.")
  )
}

tab_qc_ui <- function(id) {
  ns <- NS(id)
  nav_panel(
    title = "Quality Control", icon = icon("clipboard-check"),
    div(class = "container-fluid py-3",
      h3("Quality Control"),
      div(class = "banner", "Per-sample sequencing and alignment QC. The full ",
          tags$b("MultiQC"), " report aggregates FastQC, trimming, STAR alignment, ",
          "Salmon quantification, duplication and strand metrics for all samples."),
      navset_tab(
        nav_panel("MultiQC report", uiOutput(ns("multiqc"))),
        nav_panel("Library sizes",  uiOutput(ns("libsize"))),
        nav_panel("Expression density", uiOutput(ns("density"))),
        nav_panel("Normalisation (TMM)", uiOutput(ns("tmm"))),
        nav_panel("Depth & power", tab_qc_depth_ui(ns))
      )
    )
  )
}

tab_qc_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    img_or_soon <- function(rel, alt) {
      if (file.exists(file.path("www", rel)))
        tags$img(src = rel, style = "max-width:100%;border:1px solid #eee;border-radius:6px")
      else coming_soon(alt)
    }
    output$multiqc <- renderUI({
      if (file.exists("www/multiqc/multiqc_report.html"))
        tags$iframe(src = "multiqc/multiqc_report.html",
                    style = "width:100%;height:80vh;border:1px solid #ddd;border-radius:6px")
      else coming_soon("The MultiQC report will be embedded here once the pipeline finishes.")
    })
    output$libsize <- renderUI(img_or_soon("qc/libsize.png", "Library-size bar chart pending."))
    output$density <- renderUI(img_or_soon("qc/density.png", "log-CPM density (pre/post filter) pending."))
    output$tmm     <- renderUI(img_or_soon("qc/tmm_boxplots.png", "TMM normalisation boxplots pending."))
    output$depthfig<- renderUI(img_or_soon("qc/depth_concordance.png", "Depth-concordance figure pending."))
  })
}
