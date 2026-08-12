# Single-cell validation of the 47-gene conserved injury panel.
#
# Data: app/data/sc_injury.rds, built by painseq/analysis/20_build_app_artifact.R
# from GEO GSE154659 (Renthal et al. 2020, mouse DRG snRNA-seq — the same atlas
# this app already uses as the naive deconvolution reference, but here the full
# injury time-course). Everything shown is pseudobulk: raw UMI summed within
# sample (or sample x cell type), then edgeR TMM + limma.
#
# Every heatmap here offers independent row and column clustering. When columns
# are clustered the injury-model block separators are suppressed, because the
# columns are no longer in biological order.
SC_DIVERGE <- list(list(0, "#2166AC"), list(0.5, "#F7F7F7"), list(1, "#B2182B"))
# sequential scale for the feature plot: pale grey (not expressed) -> deep red
SC_SEQ <- list(list(0, "#E9ECEF"), list(0.25, "#FADFD6"), list(0.6, "#DE7B65"), list(1, "#8B1A16"))
SC_INJ_COL <- c(Naive = "#7F8C8D", Crush = "#4E79A7", ScNT = "#59A14F",
                SpNT = "#B2182B", CFA = "#F28E2B", Paclitaxel = "#9C755F")

# ---------------------------------------------------------------------
# Ridgeline built by hand from filled plotly polygons.
#
# plotly has no ridgeline trace, and the app also publishes as a shinylive/WebR
# static site where an extra R package would need a wasm build — so ggridges is
# used only for the static figures and the interactive version is assembled here
# from `toself`-filled scatter traces offset on the y axis.
#
# `vals` is a list of numeric vectors, one per ridge, named by group.
# `pct` (optional) is a per-ridge annotation, used to print "% expressing" when the
# distribution is computed over expressing cells only.
# ---------------------------------------------------------------------
#
# `levels` fixes the y ordering across several ridge plots so they can be combined
# into facets with a shared axis; groups absent from one facet leave a gap rather
# than shifting every other row.
sc_ridges <- function(vals, title, xlab, pct = NULL, cols = NULL, overlap = 1.6,
                      levels = NULL, show_ticks = TRUE) {
  keep <- vapply(vals, function(v) length(v) >= 5 && stats::sd(v) > 0, logical(1))
  vals <- vals[keep]
  lv <- levels %||% names(vals)
  if (!length(vals)) return(NULL)
  n <- length(lv)
  dens <- lapply(vals, function(v) stats::density(v, adjust = 1.1))
  hmax <- max(vapply(dens, function(d) max(d$y), 0))
  p <- plotly::plot_ly()
  for (g in names(vals)) {
    i <- match(g, lv); if (is.na(i)) next
    d <- dens[[g]]; off <- (n - i)                      # first level at the top
    yy <- off + overlap * d$y / hmax
    col <- if (!is.null(cols) && g %in% names(cols)) unname(cols[g]) else "#B2182B"
    p <- plotly::add_trace(p, x = c(d$x, rev(d$x)), y = c(yy, rep(off, length(yy))),
          type = "scatter", mode = "lines", fill = "toself",
          fillcolor = plotly::toRGB(col, 0.62),
          line = list(color = col, width = 1.1), name = g, showlegend = FALSE,
          hoverinfo = "text",
          text = sprintf("<b>%s</b><br>n = %d<br>median %.2f%s", g, length(vals[[g]]),
                         stats::median(vals[[g]]),
                         if (!is.null(pct) && g %in% names(pct))
                           sprintf("<br>%.1f%% expressing", pct[[g]]) else ""))
  }
  ticktext <- if (!is.null(pct)) vapply(lv, function(k)
      if (k %in% names(pct)) sprintf("%s  (%.0f%%)", k, pct[[k]]) else k, character(1)) else lv
  plotly::layout(p, title = title,
    xaxis = list(title = xlab),
    yaxis = list(tickmode = "array", tickvals = rev(seq_len(n)) - 1,
                 ticktext = if (show_ticks) ticktext else rep("", n),
                 title = "", tickfont = list(size = 9.5))) |> ply_pub()
}

# violin alternative over the same grouped values
sc_violin <- function(vals, title, ylab, kind = c("violin", "box"), cols = NULL,
                      levels = NULL, show_ticks = TRUE) {
  kind <- match.arg(kind)
  vals <- vals[vapply(vals, length, 0L) > 0]
  if (!length(vals)) return(NULL)
  lv <- levels %||% names(vals)
  p <- plotly::plot_ly()
  for (k in intersect(lv, names(vals))) {
    col <- if (!is.null(cols) && k %in% names(cols)) unname(cols[k]) else "#B2182B"
    p <- plotly::add_trace(p, x = rep(k, length(vals[[k]])), y = vals[[k]], type = kind, name = k,
          fillcolor = plotly::toRGB(col, 0.55), line = list(color = col, width = 1.1),
          marker = list(color = col, size = 3, opacity = .4),
          box = if (kind == "violin") list(visible = TRUE) else NULL,
          meanline = if (kind == "violin") list(visible = TRUE) else NULL,
          points = FALSE, showlegend = FALSE)
  }
  plotly::layout(p, title = title,
    yaxis = list(title = ylab),
    xaxis = list(title = "", tickangle = -40, categoryorder = "array", categoryarray = lv,
                 tickfont = list(size = 9), showticklabels = show_ticks)) |> ply_pub()
}

tab_sc_injury_ui <- function(id) {
  ns <- NS(id)
  nav_panel(
    title = "Single-cell injury panel", icon = icon("dna"),
    div(class = "container-fluid py-3",
      h3("The 47-gene conserved injury panel in single-cell DRG"),
      div(class = "banner",
          "The ", tags$b("SNI-positive / MINP-negative"), " panel was derived from bulk RNA-seq across ",
          "mouse, rat and macaque (see the ", tags$b("Cross-species meta-analysis"), " tab). Here it is ",
          "tested against ", tags$b("GSE154659"), " — a mouse DRG single-nucleus atlas of 141,093 cells ",
          "covering six injury models across twelve timepoints, with author-assigned cell-type labels. ",
          "This answers what the bulk data could not: ", tags$b("which cells"), " express the panel, ",
          tags$b("when"), " it peaks, and whether it is specific to ", tags$b("axotomy"),
          " rather than to injury in general."),
      uiOutput(ns("headline")),
      navset_tab(
        # ---------------- injury x time ----------------
        nav_panel("Injury × time",
          layout_sidebar(
            sidebar = sidebar(width = 330,
              radioButtons(ns("value"), "Colour by",
                c("log2FC vs Naive" = "l2fc", "z-score of mean logCPM" = "z")),
              radioButtons(ns("cells"), "Cells included",
                c("All 141,093 cells (sex + genotype as covariates)" = "all",
                  "Male C57 only (sensitivity check)" = "male")),
              tags$hr(),
              checkboxInput(ns("crow1"), "Cluster genes (rows)", FALSE),
              checkboxInput(ns("ccol1"), "Cluster columns", FALSE),
              tags$small(class = "text-muted",
                "Unclustered, genes are ordered by mean SNI log2FC from the bulk panel and columns ",
                "run Naive → Crush → ScNT → SpNT → CFA → Paclitaxel with time ascending. ",
                "Clustering columns removes the injury-model separators."),
              tags$hr(),
              sliderInput(ns("clip1"), "Colour clip (±)", min = 1, max = 6, value = 3, step = 0.5)
            ),
            div(plotly::plotlyOutput(ns("heat_time"), height = "760px"),
                uiOutput(ns("legend_note")))
          )
        ),
        # ---------------- by cell type ----------------
        nav_panel("By cell type",
          layout_sidebar(
            sidebar = sidebar(width = 330,
              radioButtons(ns("ctset"), "Cell types",
                c("All types" = "all", "Sensory neurons only" = "neuronal")),
              tags$hr(),
              checkboxInput(ns("crow2"), "Cluster genes (rows)", FALSE),
              checkboxInput(ns("ccol2"), "Cluster cell types (columns)", FALSE),
              tags$hr(),
              sliderInput(ns("clip2"), "Colour clip (±)", min = 1, max = 6, value = 3, step = 0.5),
              tags$small(class = "text-muted",
                "Each cell type is fitted independently against its own Naive baseline, so library ",
                "sizes, expression filters and dispersions are not shared across cell types.")
            ),
            div(plotly::plotlyOutput(ns("heat_ct"), height = "760px"),
                uiOutput(ns("ambient_note")))
          )
        ),
        # ---------------- one cell type over time ----------------
        nav_panel("Cell type × time",
          layout_sidebar(
            sidebar = sidebar(width = 330,
              selectInput(ns("ct_one"), "Cell type", choices = NULL),
              tags$hr(),
              checkboxInput(ns("crow3"), "Cluster genes (rows)", FALSE),
              checkboxInput(ns("ccol3"), "Cluster columns", FALSE),
              sliderInput(ns("clip3"), "Colour clip (±)", min = 1, max = 6, value = 3, step = 0.5),
              uiOutput(ns("ct_info"))
            ),
            plotly::plotlyOutput(ns("heat_ct_time"), height = "760px")
          )
        ),
        # ---------------- gene explorer ----------------
        nav_panel("Gene explorer",
          layout_sidebar(
            sidebar = sidebar(width = 340,
              selectizeInput(ns("ge_gene"), "Gene", choices = NULL,
                             options = list(maxOptions = 60, placeholder = "type a gene…")),
              radioButtons(ns("ge_kind"), "Distribution as",
                c("Ridges" = "ridge", "Violin" = "violin"), inline = TRUE),
              radioButtons(ns("ge_group"), "Rows within each panel",
                c("Cell type" = "celltype", "Time after injury" = "time")),
              checkboxGroupInput(ns("ge_injs"), "Injury models to show (panels)", choices = NULL),
              checkboxInput(ns("ge_nonzero"), "Ridges over expressing cells only", TRUE),
              tags$hr(),
              uiOutput(ns("ge_info"))
            ),
            div(
              plotly::plotlyOutput(ns("ge_feature"), height = "440px"),
              tags$hr(),
              plotly::plotlyOutput(ns("ge_dist"), height = "430px"),
              uiOutput(ns("ge_note"))
            )
          )
        ),
        # ---------------- panel score distributions ----------------
        nav_panel("Distributions",
          layout_sidebar(
            sidebar = sidebar(width = 330,
              selectInput(ns("ds_inj"), "Injury model", choices = NULL),
              radioButtons(ns("ds_kind"), "Show as",
                c("Ridges" = "ridge", "Violin" = "violin", "Box" = "box"), inline = TRUE),
              selectInput(ns("ds_ct"), "Restrict to cell type",
                          choices = c("All sensory neurons" = "__neurons__"), selected = "__neurons__"),
              tags$small(class = "text-muted",
                "The per-cell panel score is continuous, so unlike individual genes it has no ",
                "zero-inflation problem — this is the cleanest view of the injured state emerging ",
                "over time.")
            ),
            div(plotly::plotlyOutput(ns("ds_plot"), height = "700px"),
                uiOutput(ns("ds_note")))
          )
        ),
        # ---------------- recruitment ----------------
        nav_panel("Recruitment",
          div(class = "container-fluid",
            h4(class = "mt-2", "What actually changes when a gene is \"induced\"?"),
            div(class = "banner",
              "A bulk fold-change of +5 has two completely different cellular explanations, and ",
              "pseudobulk cannot tell them apart because both raise the group mean identically:",
              tags$ul(class = "mb-1 mt-2",
                tags$li(tags$b("RECRUITMENT"), " — almost no cell expressed the gene before, and now ",
                        "many do, each at an ordinary level. ", tags$i("A switch being thrown.")),
                tags$li(tags$b("LEVEL"), " — the same cells that already expressed it now make much ",
                        "more of it. ", tags$i("A dial being turned up."))),
              "These are different biology: a switch is a per-cell state decision, a dial is a graded ",
              "response. Single-cell data separates them exactly."),
            div(class = "row g-3 my-1",
              div(class = "col-md-6",
                div(class = "card h-100", div(class = "card-body py-2",
                  tags$h6(class = "card-title mb-1", "The two quantities plotted here"),
                  tags$table(class = "table table-sm mb-1",
                    tags$tbody(
                      tags$tr(tags$td(tags$b("REC"), " (recruitment)"),
                              tags$td("log2 change in the ", tags$b("fraction of cells"),
                                      " with any detectable transcript")),
                      tags$tr(tags$td(tags$b("LVL"), " (level)"),
                              tags$td("log2 change in the ", tags$b("mean amount per expressing cell"),
                                      " — cells with zero are excluded")))),
                  tags$small(class = "text-muted",
                    "With f = fraction expressing and m = mean among those cells, the mean over all ",
                    "cells is f × m. So ", tags$b("REC + LVL = the pseudobulk log2 fold-change"),
                    ", exactly. A gene sitting far right and near zero vertically is pure recruitment.")))),
              div(class = "col-md-6",
                div(class = "card h-100", div(class = "card-body py-2",
                  tags$h6(class = "card-title mb-1", "Worked example — Atf3 in NF1 neurons, SpNT"),
                  tags$table(class = "table table-sm table-borderless mb-1",
                    tags$thead(tags$tr(tags$th(""), tags$th("% cells expressing"),
                                       tags$th("mean among expressing"))),
                    tags$tbody(
                      tags$tr(tags$td("naive"), tags$td("0.2%"), tags$td("2.50")),
                      tags$tr(tags$td("24h"),   tags$td("53.0%"), tags$td("2.29")),
                      tags$tr(tags$td("48h"),   tags$td(tags$b("63.1%")), tags$td("2.28")),
                      tags$tr(tags$td("7d"),    tags$td("42.9%"), tags$td("1.82")))),
                  tags$small(class = "text-muted",
                    "The fraction moves 300-fold; the amount per expressing cell does not move at all. ",
                    "Bulk calls this \"log2FC ≈ +5\", but no single neuron is in a 5-fold state — cells ",
                    "are simply on or off."))))
            )
          ),
          uiOutput(ns("rc_note")),
          layout_sidebar(
            sidebar = sidebar(width = 330,
              radioButtons(ns("rc_view"), "View",
                c("One gene in detail" = "gene",
                  "Heatmap — all genes" = "heat",
                  "Scatter at axotomy peak" = "scatter")),
              conditionalPanel(sprintf("input['%s'] == 'gene'", ns("rc_view")),
                selectInput(ns("rc_gene"), "Gene", choices = NULL)),
              conditionalPanel(sprintf("input['%s'] == 'heat'", ns("rc_view")),
                checkboxInput(ns("crow5"), "Cluster genes (rows)", FALSE),
                checkboxInput(ns("ccol5"), "Cluster columns", FALSE),
                sliderInput(ns("clip5"), "Colour clip (±)", min = 1, max = 6, value = 4, step = .5)),
              tags$hr(),
              tags$small(class = "text-muted",
                tags$b("Is this just dropout? "), "No. If a gene were expressed uniformly at a low ",
                "level, the fraction of cells with at least one transcript would match the Poisson ",
                "expectation. It does not: Npy in NF1 neurons at 7 days averages 10.3 transcripts per ",
                "cell, which predicts ~100% detection, yet 41% of those neurons have ",
                tags$b("exactly zero"), " while the rest carry ~18 copies each. That is genuine ",
                "on/off bimodality, not a detection artefact.")
            ),
            plotly::plotlyOutput(ns("rc_plot"), height = "700px")
          )
        ),
        # ---------------- injured-state classifier + identity ----------------
        nav_panel("Injured state",
          navset_pill(
            nav_panel("Classifier",
              uiOutput(ns("cl_note")),
              plotly::plotlyOutput(ns("cl_roc"), height = "440px"),
              tags$hr(),
              layout_sidebar(
                sidebar = sidebar(width = 300,
                  selectInput(ns("cl_ct"), "Cell types", choices = NULL, multiple = TRUE)),
                plotly::plotlyOutput(ns("cl_curves"), height = "460px"))
            ),
            nav_panel("Subtype identity",
              uiOutput(ns("id_note")),
              plotly::plotlyOutput(ns("id_traj"), height = "430px"),
              tags$hr(),
              plotly::plotlyOutput(ns("id_scatter"), height = "470px"))
          )
        ),
        # ---------------- UMAP ----------------
        nav_panel("UMAP",
          layout_sidebar(
            sidebar = sidebar(width = 330,
              radioButtons(ns("umap_by"), "Colour by",
                c("Cell type" = "celltype", "Injury model" = "inj",
                  "Time after injury" = "time", "Panel score" = "score")),
              tags$small(class = "text-muted",
                "Embedding computed in scanpy: 2,000 highly variable genes, 50 PCs, ",
                "UMAP on a 15-nearest-neighbour graph. Cell-type labels are the authors', not ",
                "re-clustered here. Points are a 20% stratified subsample so the plot stays ",
                "responsive; rare populations are kept in full."),
              tags$hr(),
              tags$small(class = "text-muted",
                tags$b("Not batch-integrated, deliberately. "),
                "Every sample belongs to exactly one injury × time group, so integrating on ",
                "sample would regress out the condition effect along with any technical batch ",
                "effect — it would delete the biology this analysis exists to show.")
            ),
            plotly::plotlyOutput(ns("umap"), height = "760px")
          )
        ),
        # ---------------- MINP specificity ----------------
        nav_panel("MINP specificity",
          navset_pill(
            nav_panel("Concordance",
              layout_sidebar(
                sidebar = sidebar(width = 330,
                  radioButtons(ns("cc_contrast"), "Bulk contrast",
                    c("MINP vs Sham" = "MINP", "SNI vs Sham (positive control)" = "SNI")),
                  checkboxInput(ns("crow4"), "Cluster cell types (rows)", FALSE),
                  checkboxInput(ns("ccol4"), "Cluster columns", FALSE),
                  tags$small(class = "text-muted",
                    "Correlation of the mouse bulk log2FC vector against each single-cell injury ",
                    "signature, over all shared genes — not just the 47, which would beg the question.")
                ),
                div(plotly::plotlyOutput(ns("minp_cc"), height = "700px"),
                    uiOutput(ns("minp_note")))
              )
            ),
            nav_panel("MINP genes in single cell",
              uiOutput(ns("mg_note")),
              plotly::plotlyOutput(ns("mg_plot"), height = "560px"),
              tags$hr(), DT::DTOutput(ns("mg_tab"))),
            nav_panel("Beyond the injury axis",
              div(class = "banner",
                  "Every test so far asked ", tags$b("\"does MINP look like nerve injury?\""),
                  " and answered no. That leaves a fair objection open: MINP might be doing ",
                  "something real that simply is not an injury programme — a glial reaction, an ",
                  "immune infiltrate, or a change confined to a cell type too rare to move ",
                  "whole-DRG bulk. This tests for exactly that, by asking whether the bulk MINP ",
                  "contrast is enriched for any ", tags$b("cell type's marker programme"),
                  " (50 markers per type, derived from naive atlas cells)."),
              uiOutput(ns("by_note")),
              plotly::plotlyOutput(ns("by_plot"), height = "620px"),
              tags$hr(),
              plotly::plotlyOutput(ns("by_consistency"), height = "540px"),
              uiOutput(ns("by_foot")))
          )
        ),
        # ---------------- Atf3 dependence ----------------
        nav_panel("Atf3 dependence",
          div(class = "banner",
              "Is the panel downstream of ", tags$b("Atf3"), ", the master regeneration ",
              "transcription factor? The Atf3-WT / Atf3-KO arm of GSE154659 (crush, neurons only, ",
              "17,665 cells) is a direct loss-of-function test."),
          navset_pill(
            nav_panel("Per-cell (recommended)",
              uiOutput(ns("ap_note")),
              plotly::plotlyOutput(ns("ap_entry"), height = "500px"),
              tags$hr(),
              plotly::plotlyOutput(ns("ap_rec"), height = "560px")),
            nav_panel("Pseudobulk",
              uiOutput(ns("atf3_note")),
              plotly::plotlyOutput(ns("atf3_plot"), height = "620px"),
              tags$hr(), DT::DTOutput(ns("atf3_tab")))
          )
        ),
        # ---------------- deconvolution reference ----------------
        nav_panel("Deconvolution reference",
          div(class = "banner",
              "The ", tags$b("Deconvolution"), " tab builds its MuSiC reference from ",
              tags$b("naive"), " cells only. That reference cannot represent ",
              tags$i("Repair fibroblast"), " at all — the whole atlas contains just six naive ones. ",
              "Would a reference built from injured cells do better?"),
          uiOutput(ns("deconv_note")),
          plotly::plotlyOutput(ns("deconv_plot"), height = "520px")
        ),
        # ---------------- peak landscape ----------------
        nav_panel("Peak landscape",
          div(class = "banner",
              "For every panel gene: the cell type, injury model and timepoint at which it reaches ",
              "its maximum induction, searched across all fitted cell types and contrasts. Plotted ",
              "against time this becomes a ", tags$b("temporal wave"),
              " — immediate-early genes cresting within a day, secreted effectors a week later."),
          navset_pill(
            nav_panel("Landscape",
              layout_sidebar(
                sidebar = sidebar(width = 300,
                  radioButtons(ns("pk_colour"), "Colour by",
                    c("Peak cell type" = "peak_celltype", "Peak injury model" = "peak_injury",
                      "Timing wave" = "wave")),
                  checkboxInput(ns("pk_labels"), "Show gene labels", TRUE),
                  tags$small(class = "text-muted",
                    "x = when the gene peaks, y = how strongly. Hover any point for the full record.")
                ),
                div(plotly::plotlyOutput(ns("pk_plot"), height = "620px"),
                    uiOutput(ns("pk_note")))
              )
            ),
            nav_panel("Counts",
              plotly::plotlyOutput(ns("pk_bars"), height = "560px")),
            nav_panel("Table",
              DT::DTOutput(ns("peak")),
              downloadButton(ns("dl_peak"), "Download CSV", class = "btn-sm btn-outline-secondary mt-2"))
          )
        ),
        # ---------------- documentation ----------------
        nav_panel("Methods & caveats", uiOutput(ns("methods")))
      )
    )
  )
}

tab_sc_injury_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    SC <- load_rds("sc_injury")

    # ---- shared plotting helper -------------------------------------------
    # M: numeric matrix (genes x columns). Row/col clustering optional; block
    # separators only drawn when columns are in their biological order.
    draw_heat <- function(M, P, collab, blocks, title, zlab, clip,
                          do_row, do_col, rowlab, hover_extra = NULL) {
      keep <- rowSums(is.finite(M)) > 0
      if (!any(keep)) return(NULL)
      ord <- cluster_order(M, do_row, do_col)
      M <- M[ord$row, ord$col, drop = FALSE]; P <- P[ord$row, ord$col, drop = FALSE]
      rl <- rowlab[ord$row]; cl <- collab[ord$col]
      bl <- if (do_col) NULL else blocks[ord$col]

      hov <- matrix("", nrow(M), ncol(M))
      for (i in seq_len(nrow(M))) for (j in seq_len(ncol(M)))
        hov[i, j] <- sprintf("<b>%s</b><br>%s<br>%s: %s%s",
          rownames(M)[i], cl[j], zlab,
          if (is.na(M[i, j])) "n/a" else sprintf("%+.2f", M[i, j]),
          if (!is.null(P) && !is.na(P[i, j])) sprintf("<br>adj.p: %s", format.pval(P[i, j], digits = 2)) else "")

      p <- plotly::plot_ly(x = seq_len(ncol(M)), y = rl,
            z = pmax(pmin(M, clip), -clip), type = "heatmap",
            colorscale = SC_DIVERGE, zmid = 0, zmin = -clip, zmax = clip,
            text = hov, hoverinfo = "text",
            colorbar = list(title = list(text = gsub(" ", "<br>", zlab))))
      shp <- list()
      if (!is.null(bl)) {
        b <- cumsum(rle(as.character(bl))$lengths); b <- b[-length(b)]
        shp <- lapply(b, function(k) list(type = "line", x0 = k + .5, x1 = k + .5,
                 y0 = -.5, y1 = nrow(M) - .5, line = list(color = "black", width = 2.5)))
      }
      plotly::layout(p, title = title, shapes = shp,
        margin = list(b = 130),
        xaxis = list(tickmode = "array", tickvals = seq_len(ncol(M)), ticktext = cl,
                     tickangle = -45, tickfont = list(size = 9)),
        yaxis = list(autorange = "reversed", categoryorder = "array", categoryarray = rev(rl),
                     tickfont = list(size = 8.5))) |> ply_pub()
    }

    rowlab_of <- function(SC) setNames(SC$rows$label, SC$rows$gene)

    # ---- headline banner ---------------------------------------------------
    output$headline <- renderUI({
      if (is.null(SC)) return(NULL)
      s <- SC$stats
      div(class = "banner", style = "border-left-color:#B2182B",
        tags$b("Result: the panel is axotomy-specific. "),
        sprintf("Across nerve injury (Crush / ScNT / SpNT at 24h–7d) the panel shows a mean log2FC of %+.2f, with %.0f%% of gene × column combinations significantly up. ",
                s$axotomy$mean_l2fc, s$axotomy$pct_sig_up),
        sprintf("Across the two non-axotomy insults (CFA inflammation, paclitaxel neuropathy) the mean log2FC is %+.2f and %s panel genes are significantly induced. ",
                s$non_axotomy$mean_l2fc,
                if (s$non_axotomy$n_sig_up == 0) tags$b("zero") else s$non_axotomy$n_sig_up),
        sprintf("Single-cell SpNT-7d fold-changes correlate with this study's bulk mouse SNI-vs-Sham at r = %.2f (%.0f%% same-direction, n = %d).",
                s$bulk_r, 100 * s$bulk_same_dir, s$bulk_n))
    })

    # ---- 1. injury x time --------------------------------------------------
    output$heat_time <- plotly::renderPlotly({
      req(SC)
      cm <- SC$cmeta
      M <- if (input$value == "z") SC$zscore else if (input$cells == "male") SC$l2fc_male else SC$l2fc
      P <- if (input$cells == "male") SC$padj_male else SC$padj
      M <- M[, cm$col, drop = FALSE]; P <- P[, cm$col, drop = FALSE]
      zlab <- if (input$value == "z") "z-score" else "log2FC vs Naive"
      ttl <- sprintf("Conserved injury panel — %s%s", zlab,
                     if (input$cells == "male") " · male C57 only" else "")
      draw_heat(M, P, collab = paste(cm$inj, cm$time_label), blocks = cm$inj,
                title = ttl, zlab = zlab, clip = input$clip1,
                do_row = isTRUE(input$crow1), do_col = isTRUE(input$ccol1),
                rowlab = rowlab_of(SC))
    })

    output$legend_note <- renderUI({
      req(SC); s <- SC$stats
      div(class = "small text-muted mt-2",
        tags$b("Reading this: "), "each column is one injury model at one timepoint. ",
        "Naive is the reference, so its column is 0 by construction. ",
        tags$b("*"), " marks genes detected in fewer than 0.5% of cells — shallow single-nucleus ",
        "coverage makes those estimates unreliable. ", tags$b("‡"), " marks ", tags$i("Sprr1a"),
        " and ", tags$i("Sox11"), ", the classic regeneration-associated genes, shown for reference: ",
        "they are not part of the 47 (", tags$i("Sprr1a"), " has no macaque ortholog, so it could not ",
        "survive the three-species join). Grey cells are genes absent from the atlas or below its ",
        "expression filter. ",
        sprintf("The male-C57-only sensitivity fit agrees with the all-cells fit at r = %.3f.", s$sensitivity_r))
    })

    # ---- 2. by cell type ---------------------------------------------------
    output$heat_ct <- plotly::renderPlotly({
      req(SC)
      cts <- colnames(SC$ct_l2fc)
      if (input$ctset == "neuronal") cts <- intersect(SC$neuronal, cts)
      M <- SC$ct_l2fc[, cts, drop = FALSE]; P <- SC$ct_padj[, cts, drop = FALSE]
      draw_heat(M, P, collab = cts, blocks = NULL,
                title = "Panel by cell type — axotomy peak (72h + 7d) vs Naive",
                zlab = "log2FC vs Naive", clip = input$clip2,
                do_row = isTRUE(input$crow2), do_col = isTRUE(input$ccol2),
                rowlab = rowlab_of(SC))
    })

    output$ambient_note <- renderUI({
      req(SC)
      div(class = "small text-muted mt-2",
        tags$b("Important caveat — ambient RNA. "),
        "Neuron-restricted transcripts such as ", tags$i("Npy"), " and ", tags$i("Gal"),
        " appear induced in every cell type here, including B cells and neutrophils. Real neuronal ",
        "enrichment does exist (mean panel log2FC 2.14 in neurons vs 1.41 in non-neuronal types, with ",
        "30 of 33 genes higher in neurons), but in single-nucleus data a transcript induced to very ",
        "high levels in one population leaks into the ambient pool and inflates apparent induction ",
        "everywhere. ", tags$b("Read this heatmap as a relative ranking across cell types, not as ",
        "absolute per-type induction."),
        tags$br(), tags$br(),
        tags$b("Two cell types are absent: "), paste(SC$stats$no_baseline, collapse = " and "),
        ". They are injury-induced populations that do not exist in naive tissue, so no ",
        tags$i("versus Naive"), " contrast is definable for them. That is precisely why the ",
        "naive-only reference used on the ", tags$b("Deconvolution"), " tab cannot detect them.")
    })

    # ---- 3. one cell type over time ---------------------------------------
    observe({
      req(SC)
      cts <- names(SC$ct_tc)
      lab <- setNames(cts, sprintf("%s — %s", cts, SC$celltype_full[cts]))
      updateSelectInput(session, "ct_one", choices = lab,
                        selected = if ("NF1" %in% cts) "NF1" else cts[1])
    })

    output$ct_info <- renderUI({
      req(SC, input$ct_one)
      f <- SC$ct_tc[[input$ct_one]]; req(f)
      div(class = "small text-muted mt-2",
          tags$b(SC$celltype_full[[input$ct_one]]), tags$br(),
          sprintf("%s cells across %d contrasts.", format(f$n_cells, big.mark = ","), nrow(f$cmeta)))
    })

    output$heat_ct_time <- plotly::renderPlotly({
      req(SC, input$ct_one)
      f <- SC$ct_tc[[input$ct_one]]; req(f)
      cm <- f$cmeta
      cm <- cm[order(match(cm$inj, SC$inj_order), cm$time), ]
      M <- f$M[, cm$col, drop = FALSE]
      tl <- ifelse(cm$time == 0, "0", ifelse(cm$time < 168, paste0(cm$time, "h"), paste0(cm$time / 24, "d")))
      draw_heat(M, NULL, collab = paste(cm$inj, tl), blocks = cm$inj,
                title = sprintf("%s — panel across injury and time", input$ct_one),
                zlab = "log2FC vs Naive", clip = input$clip3,
                do_row = isTRUE(input$crow3), do_col = isTRUE(input$ccol3),
                rowlab = rowlab_of(SC))
    })

    # =====================================================================
    # Gene explorer — feature plot + distribution
    # ---------------------------------------------------------------------
    # sc_percell.rds is a SEPARATE artifact loaded lazily: load_rds() is memoised
    # per file, so opening this panel is what pulls the per-cell matrix into
    # memory, not app start.
    # =====================================================================
    PC <- reactive(load_rds("sc_percell"))

    # Populate from sc_injury's copy of the gene table, NOT from sc_percell — an
    # observer runs at server init, so touching PC() here would load the per-cell
    # matrix on app start and defeat the lazy split.
    observe({
      g <- SC$gene_universe; req(g)
      lab <- setNames(g$gene, ifelse(g$panel, paste0(g$label, "  — panel"),
                             ifelse(g$minp, paste0(g$gene, "  — MINP-associated"),
                             ifelse(!is.na(g$marker_for), paste0(g$gene, "  — ", g$marker_for, " marker"),
                                    g$gene))))
      updateSelectizeInput(session, "ge_gene", choices = lab,
                           selected = if ("Atf3" %in% g$gene) "Atf3" else g$gene[1], server = TRUE)
    })

    observe({
      req(SC)
      updateCheckboxGroupInput(session, "ge_injs", choices = SC$inj_order,
        selected = intersect(c("Naive", "Crush", "ScNT", "SpNT"), SC$inj_order))
    })

    # values split by injury (the facet) and, within each, by the chosen row variable
    ge_vals <- reactive({
      p <- PC(); req(p, input$ge_gene, input$ge_gene %in% rownames(p$expr))
      req(length(input$ge_injs) > 0)
      v <- as.numeric(p$expr[input$ge_gene, ])
      m <- p$meta
      tlab <- ifelse(m$time == 0, "naive",
               ifelse(m$time < 168, paste0(m$time, "h"), paste0(m$time / 24, "d")))
      key <- if (input$ge_group == "celltype") as.character(m$celltype) else tlab
      # a single shared row ordering across every panel, so the facets line up
      lv <- if (input$ge_group == "celltype") intersect(names(SC$celltype_full), unique(key))
            else c("naive", unique(tlab[m$time > 0][order(m$time[m$time > 0])]))
      lv <- intersect(lv, unique(key))
      injs <- intersect(SC$inj_order, input$ge_injs)
      list(v = v, key = key, levels = lv, inj = as.character(m$inj), injs = injs, meta = m)
    })

    output$ge_feature <- plotly::renderPlotly({
      d <- ge_vals(); req(d)
      m <- d$meta; v <- d$v
      o <- order(v)                                     # expressing cells drawn last
      p <- plotly::plot_ly(x = m$x[o], y = m$y[o], type = "scattergl", mode = "markers",
            marker = list(size = 3.2, opacity = .75, color = v[o], colorscale = SC_SEQ,
                          cmin = 0, cmax = max(v, 0.1),
                          colorbar = list(title = list(text = "log1p<br>CP10K"))),
            text = sprintf("%s<br>%s %sh<br>%s = %.2f", m$celltype[o], m$inj[o], m$time[o],
                           input$ge_gene, v[o]),
            hoverinfo = "text")
      plotly::layout(p, title = sprintf("%s — expression across the DRG atlas", input$ge_gene),
        xaxis = list(title = "UMAP 1", showticklabels = FALSE),
        yaxis = list(title = "UMAP 2", showticklabels = FALSE)) |> ply_pub()
    })

    # one sub-plot per injury model, sharing the row ordering and the value axis
    output$ge_dist <- plotly::renderPlotly({
      d <- ge_vals(); req(d)
      rng <- range(d$v[d$v > 0], 0, na.rm = TRUE)
      panels <- lapply(seq_along(d$injs), function(i) {
        m <- d$injs[i]
        k <- d$inj == m
        sp <- split(d$v[k], factor(d$key[k], levels = d$levels))
        sp <- sp[vapply(sp, length, 0L) > 0]
        if (!length(sp)) return(NULL)
        pct <- lapply(sp, function(z) 100 * mean(z > 0))
        first <- i == 1L
        if (input$ge_kind == "ridge") {
          vv <- if (isTRUE(input$ge_nonzero)) lapply(sp, function(z) z[z > 0]) else sp
          vv <- vv[vapply(vv, length, 0L) >= 5]
          if (!length(vv)) return(NULL)
          p <- sc_ridges(vv, title = "", xlab = "", pct = pct, levels = d$levels,
                         show_ticks = first, cols = setNames(rep(unname(SC_INJ_COL[m]),
                                                                 length(d$levels)), d$levels))
          plotly::layout(p, xaxis = list(range = rng))
        } else {
          p <- sc_violin(sp, title = "", ylab = "", kind = "violin", levels = d$levels,
                         show_ticks = TRUE,
                         cols = setNames(rep(unname(SC_INJ_COL[m]), length(d$levels)), d$levels))
          plotly::layout(p, yaxis = list(range = rng, showticklabels = first))
        }
      })
      names(panels) <- d$injs
      keep <- !vapply(panels, is.null, logical(1))
      req(any(keep))
      panels <- panels[keep]; labs <- d$injs[keep]

      pp <- plotly::subplot(panels, nrows = 1,
                            shareY = input$ge_kind == "ridge",
                            shareX = input$ge_kind == "violin",
                            margin = 0.006, titleX = FALSE, titleY = FALSE)
      # column headers, positioned over each sub-plot
      xs <- (seq_along(labs) - 0.5) / length(labs)
      ann <- lapply(seq_along(labs), function(i) list(
        x = xs[i], y = 1.045, xref = "paper", yref = "paper", text = paste0("<b>", labs[i], "</b>"),
        showarrow = FALSE, font = list(size = 12, color = unname(SC_INJ_COL[labs[i]]))))
      plotly::layout(pp, annotations = ann, margin = list(t = 72),
        title = sprintf("%s — %s by %s, split by injury model", input$ge_gene,
                        if (input$ge_kind == "ridge") "distribution" else "violin",
                        if (input$ge_group == "celltype") "cell type" else "time"))
    })

    output$ge_info <- renderUI({
      p <- PC(); req(p, input$ge_gene)
      g <- p$genes[p$genes$gene == input$ge_gene, ]
      if (!nrow(g)) return(NULL)
      v <- as.numeric(p$expr[input$ge_gene, ])
      div(class = "small",
        tags$b(input$ge_gene), tags$br(),
        if (isTRUE(g$panel)) tags$span(class = "badge bg-danger", "conserved injury panel"),
        if (isTRUE(g$minp)) tags$span(class = "badge bg-warning text-dark", "MINP-associated"),
        if (!is.na(g$marker_for)) tags$span(class = "badge bg-secondary",
                                            paste("marker:", g$marker_for)),
        tags$br(), tags$br(),
        sprintf("Detected in %.1f%% of the %s displayed cells; mean %.2f among expressing cells.",
                100 * mean(v > 0), format(length(v), big.mark = ","),
                if (any(v > 0)) mean(v[v > 0]) else 0))
    })

    output$ge_note <- renderUI({
      div(class = "small text-muted mt-2",
        tags$b("Why ridges default to expressing cells only. "),
        "Panel gene expression is extremely zero-inflated — Atf3 is detected in 0.2% of naive NF1 ",
        "neurons — so a density over all cells is one spike at zero and shows nothing. Ridges are ",
        "therefore computed over expressing cells, with the ", tags$b("% expressing"),
        " printed beside each row label in the leftmost panel and in every hover tooltip, so the ",
        "zeros are never hidden. Switch to violin to see the full distribution including zeros. ",
        tags$br(), tags$br(),
        "Each panel is one injury model and all panels share the same rows and the same value axis, ",
        "so they can be read across. Cells are a stratified subsample of the atlas.")
    })

    # =====================================================================
    # Panel score distributions
    # =====================================================================
    # likewise sourced from sc_injury, so the picker exists before the per-cell
    # matrix is ever read
    observe({
      req(SC)
      updateSelectInput(session, "ds_inj",
        choices = setdiff(SC$inj_order, "Naive"), selected = "SpNT")
      updateSelectInput(session, "ds_ct",
        choices = c("All sensory neurons" = "__neurons__",
                    setNames(names(SC$celltype_full), names(SC$celltype_full))),
        selected = "__neurons__")
    })

    output$ds_plot <- plotly::renderPlotly({
      p <- PC(); req(p, input$ds_inj)
      m <- p$meta
      keep <- if (identical(input$ds_ct, "__neurons__")) m$celltype %in% SC$neuronal
              else m$celltype == input$ds_ct
      s <- m[keep & m$inj %in% c("Naive", input$ds_inj), ]
      req(nrow(s) > 50)
      s$tl <- ifelse(s$time == 0, "naive",
               ifelse(s$time < 168, paste0(s$time, "h"), paste0(s$time / 24, "d")))
      lv <- c("naive", unique(s$tl[s$time > 0][order(s$time[s$time > 0])]))
      sp <- split(s$score, factor(s$tl, levels = lv))
      sp <- sp[vapply(sp, length, 0L) >= 20]
      ttl <- sprintf("Panel score over time — %s, %s", input$ds_inj,
                     if (identical(input$ds_ct, "__neurons__")) "sensory neurons" else input$ds_ct)
      if (input$ds_kind == "ridge") sc_ridges(sp, ttl, "per-cell panel score")
      else sc_violin(sp, ttl, "per-cell panel score", kind = input$ds_kind)
    })

    output$ds_note <- renderUI({
      req(SC, SC$classif)
      div(class = "small text-muted mt-2",
        sprintf("The injured-state threshold is the 99th percentile of naive cells (%.3f). ",
                SC$classif$threshold),
        "Watch the whole distribution shift right and then, in Crush, drift back — the ",
        "regeneration-permissive model is the only one where neurons return toward baseline.")
    })

    # =====================================================================
    # Recruitment
    # =====================================================================
    output$rc_note <- renderUI({
      req(SC, SC$recruit)
      d <- SC$recruit$drivers
      up <- d[d$total > 1, ]
      div(class = "banner", style = "border-left-color:#B2182B",
        tags$b("Injury switches genes on, it does not turn them up. "),
        sprintf("Across the %d panel genes, the median share of the fold-change attributable to recruitment is %.0f%%. ",
                nrow(d), median(d$pct_recruitment, na.rm = TRUE)),
        sprintf("Of the %d genes induced more than two-fold at the axotomy peak, %d are recruitment-driven, %d mixed and %s level-driven. ",
                nrow(up), sum(up$driver == "recruitment"), sum(up$driver == "mixed"),
                if (sum(up$driver == "level") == 0) tags$b("none") else sum(up$driver == "level")),
        "In NF1 neurons after spinal nerve transection, Atf3 goes from 0.2% to 63% of cells ",
        "expressing it while the amount per expressing cell does not move.")
    })

    observe({
      req(SC, SC$recruit)
      g <- SC$recruit$drivers
      g <- g[order(-g$total), ]
      updateSelectInput(session, "rc_gene",
        choices = setNames(g$gene, sprintf("%s  (%.0f%% recruitment)", g$gene, g$pct_recruitment)),
        selected = if ("Atf3" %in% g$gene) "Atf3" else g$gene[1])
    })

    output$rc_plot <- plotly::renderPlotly({
      req(SC, SC$recruit)
      d <- SC$recruit$pooled

      # ---- one gene in detail: the two components as separate curves ----------
      if (identical(input$rc_view, "gene")) {
        req(input$rc_gene)
        s <- d[d$gene == input$rc_gene, ]
        req(nrow(s) > 0)
        s$inj <- factor(s$inj, levels = SC$inj_order)
        s <- s[order(match(s$inj, SC$inj_order), s$time), ]
        base <- s[s$inj == "Naive", ]
        p1 <- plotly::plot_ly(); p2 <- plotly::plot_ly()
        for (m in setdiff(levels(droplevels(s$inj)), "Naive")) {
          z <- rbind(base, s[s$inj == m, ]); z <- z[order(z$time), ]
          col <- unname(SC_INJ_COL[m])
          p1 <- plotly::add_trace(p1, x = pmax(z$time, 3), y = 100 * z$f, type = "scatter",
                 mode = "lines+markers", name = m, legendgroup = m,
                 line = list(color = col, width = 2), marker = list(color = col, size = 6),
                 hovertext = sprintf("%s · %sh<br>%.1f%% of cells expressing", m, z$time, 100 * z$f),
                 hoverinfo = "text")
          p2 <- plotly::add_trace(p2, x = pmax(z$time, 3), y = z$m, type = "scatter",
                 mode = "lines+markers", name = m, legendgroup = m, showlegend = FALSE,
                 line = list(color = col, width = 2), marker = list(color = col, size = 6),
                 hovertext = sprintf("%s · %sh<br>mean %.2f among expressing cells", m, z$time, z$m),
                 hoverinfo = "text")
        }
        ax <- list(type = "log", tickvals = c(6, 24, 72, 168, 672, 2160),
                   ticktext = c("6h", "24h", "72h", "7d", "28d", "90d"),
                   title = "Time after injury")
        drv <- SC$recruit$drivers[SC$recruit$drivers$gene == input$rc_gene, ]
        return(plotly::subplot(
            plotly::layout(p1, yaxis = list(title = "RECRUITMENT<br>% of cells expressing"), xaxis = ax),
            plotly::layout(p2, yaxis = list(title = "LEVEL<br>mean among expressing"), xaxis = ax),
            nrows = 2, shareX = TRUE, titleY = TRUE, margin = .07) |>
          plotly::layout(title = sprintf(
            "%s — recruitment (top) versus level (bottom)%s", input$rc_gene,
            if (nrow(drv)) sprintf("  ·  %.0f%% of its fold-change is recruitment · %s-driven",
                                   drv$pct_recruitment[1], drv$driver[1]) else "")) |> ply_pub())
      }

      if (identical(input$rc_view, "scatter")) {
        a <- SC$recruit$drivers
        cols <- c(recruitment = "#B2182B", mixed = "#E8A33D", level = "#2166AC")
        p <- plotly::plot_ly()
        for (k in unique(a$driver)) {
          s <- a[a$driver == k, ]
          p <- plotly::add_trace(p, x = s$recruitment, y = s$level, type = "scatter",
                mode = "markers+text", name = k,
                marker = list(size = 9, color = unname(cols[k]), line = list(width = .5, color = "#fff")),
                text = s$gene, textposition = "top center", textfont = list(size = 8),
                hovertext = sprintf("<b>%s</b><br>recruitment %+.2f<br>level %+.2f<br>total %+.2f<br>%.0f%% recruitment",
                                    s$gene, s$recruitment, s$level, s$total, s$pct_recruitment),
                hoverinfo = "text")
        }
        lim <- max(abs(c(a$recruitment, a$level)), na.rm = TRUE)
        return(plotly::layout(p, title = "Recruitment versus level at the axotomy peak",
          shapes = list(list(type = "line", x0 = -lim, x1 = lim, y0 = -lim, y1 = lim,
                             line = list(color = "grey60", dash = "dot", width = 1))),
          xaxis = list(title = "Recruitment — log2 change in fraction expressing", range = c(-lim, lim)),
          yaxis = list(title = "Level — log2 change among expressing cells", range = c(-lim, lim)),
          legend = list(orientation = "h", y = -0.16)) |> ply_pub())
      }
      cm <- SC$cmeta
      mk <- function(val) {
        M <- matrix(NA_real_, length(unique(d$gene)), nrow(cm),
                    dimnames = list(unique(d$gene), cm$col))
        gi <- match(d$gene, rownames(M)); ci <- match(d$group, cm$col)
        ok <- !is.na(gi) & !is.na(ci)
        M[cbind(gi[ok], ci[ok])] <- d[[val]][ok]
        M
      }
      Mr <- mk("recruitment"); Ml <- mk("level")
      ord <- SC$recruit$drivers$gene[order(-SC$recruit$drivers$total)]
      ord <- intersect(ord, rownames(Mr))
      M <- cbind(Mr[ord, , drop = FALSE], Ml[ord, , drop = FALSE])
      # cmeta carries `time_label`, not `label` — using the wrong name silently
      # collapses colnames to length 2 and the matrix assignment then fails.
      # Label columns with injury + time and carry the REC/LVL split in the block
      # separator instead, so the axis reads as biology rather than as a code.
      colnames(M) <- c(paste0(cm$inj, " ", cm$time_label, " · REC"),
                       paste0(cm$inj, " ", cm$time_label, " · LVL"))
      draw_heat(M, NULL, collab = colnames(M),
        blocks = rep(c("← RECRUITMENT: log2 change in % of cells expressing",
                       "LEVEL: log2 change in amount per expressing cell →"), each = nrow(cm)),
        title = paste("Left block = RECRUITMENT (how many cells express it) ·",
                      "Right block = LEVEL (how much each expressing cell makes)"),
        zlab = "log2 vs Naive", clip = input$clip5,
        do_row = isTRUE(input$crow5), do_col = isTRUE(input$ccol5),
        rowlab = setNames(rownames(M), rownames(M)))
    })

    # =====================================================================
    # Injured-state classifier
    # =====================================================================
    output$cl_note <- renderUI({
      req(SC, SC$classif)
      a <- SC$classif$auc
      n <- a[a$compartment == "Neurons", ]
      div(class = "banner", style = "border-left-color:#B2182B",
        tags$b("The panel classifies single cells, not just samples. "),
        sprintf("In sensory neurons the per-cell panel score separates naive from injured with AUC %.3f for spinal nerve transection, %.3f for sciatic transection and %.3f for crush. ",
                n$auc[n$model == "SpNT"], n$auc[n$model == "ScNT"], n$auc[n$model == "Crush"]),
        sprintf("For the two non-axotomy insults it is near chance — %.3f for CFA and %.3f for paclitaxel. ",
                n$auc[n$model == "CFA"], n$auc[n$model == "Paclitaxel"]),
        "The specificity established in pseudobulk holds at single-cell resolution.")
    })

    output$cl_roc <- plotly::renderPlotly({
      req(SC, SC$classif)
      r <- SC$classif$roc; r <- r[r$compartment == "Neurons", ]
      p <- plotly::plot_ly()
      for (m in unique(r$model)) {
        s <- r[r$model == m, ]
        au <- SC$classif$auc$auc[SC$classif$auc$compartment == "Neurons" & SC$classif$auc$model == m]
        p <- plotly::add_trace(p, x = s$fpr, y = s$tpr, type = "scatter", mode = "lines",
              name = sprintf("%s (AUC %.3f)", m, au),
              line = list(color = unname(SC_INJ_COL[m]), width = 2.4))
      }
      plotly::layout(p, title = "ROC — naive versus injured sensory neurons (≥24h)",
        shapes = list(list(type = "line", x0 = 0, x1 = 1, y0 = 0, y1 = 1,
                           line = list(color = "grey60", dash = "dash", width = 1))),
        xaxis = list(title = "False positive rate", range = c(0, 1)),
        yaxis = list(title = "True positive rate", range = c(0, 1)),
        legend = list(x = .55, y = .12)) |> ply_pub()
    })

    observe({
      req(SC, SC$classif)
      ct <- intersect(names(SC$celltype_full), unique(SC$classif$curves$celltype))
      updateSelectInput(session, "cl_ct", choices = ct,
                        selected = intersect(c("NF1", "NF2", "NP", "PEP1", "Satglia", "Schwann_M"), ct))
    })

    output$cl_curves <- plotly::renderPlotly({
      req(SC, SC$classif, input$cl_ct)
      d <- SC$classif$curves
      d <- d[d$celltype %in% input$cl_ct, ]
      req(nrow(d) > 0)
      p <- plotly::plot_ly()
      for (ct in unique(d$celltype)) for (m in unique(d$inj[d$celltype == ct])) {
        s <- d[d$celltype == ct & d$inj == m, ]; s <- s[order(s$time), ]
        p <- plotly::add_trace(p, x = pmax(s$time, 3), y = s$pct, type = "scatter",
              mode = "lines+markers", name = paste(ct, m),
              line = list(color = unname(SC_INJ_COL[m]), width = 1.8),
              marker = list(size = 5, color = unname(SC_INJ_COL[m])),
              legendgroup = m,
              hovertext = sprintf("<b>%s</b><br>%s<br>%s<br>%.1f%% injured (n=%d)", ct, m,
                                  ifelse(s$time < 168, paste0(s$time, "h"), paste0(s$time / 24, "d")),
                                  s$pct, s$n),
              hoverinfo = "text")
      }
      plotly::layout(p, title = "Fraction of cells entering the injured state",
        xaxis = list(title = "Time after injury (log scale)", type = "log",
                     tickvals = c(6, 24, 72, 168, 672, 2160),
                     ticktext = c("6h", "24h", "72h", "7d", "28d", "90d")),
        yaxis = list(title = "% of cells in the injured state")) |> ply_pub()
    })

    # =====================================================================
    # Subtype identity
    # =====================================================================
    output$id_note <- renderUI({
      req(SC, SC$ident)
      div(class = "banner", style = "border-left-color:#B2182B",
        tags$b("The panel marks loss of subtype identity. "),
        "Renthal et al.'s central claim is that injured neurons lose their subtype identity and ",
        "converge on a common injured state. Scoring every neuron on ", tags$b("its own"),
        " subtype's naive markers and correlating against the panel score gives ",
        sprintf("r = %+.3f in naive cells but r = %+.3f after axotomy (≥24h). ",
                SC$ident$r_naive, SC$ident$r_axotomy),
        "Panel-high cells are exactly the cells that have drifted furthest from what they were. ",
        tags$br(), tags$br(),
        tags$b("Markers were derived from naive cells only"), " — picking them using injured cells ",
        "would have built the answer into the question.")
    })

    output$id_traj <- plotly::renderPlotly({
      req(SC, SC$ident)
      d <- SC$ident$by_group; d <- d[d$inj != "Naive", ]
      p <- plotly::plot_ly()
      for (m in unique(d$inj)) {
        s <- d[d$inj == m, ]; s <- s[order(s$time), ]
        p <- plotly::add_trace(p, x = pmax(s$time, 3), y = s$mean_identity, type = "scatter",
              mode = "lines+markers", name = m,
              line = list(color = unname(SC_INJ_COL[m]), width = 2.2),
              marker = list(size = 6, color = unname(SC_INJ_COL[m])),
              hovertext = sprintf("<b>%s</b><br>%s<br>identity %+.2f SD<br>r = %+.2f (n=%d)", m,
                                  ifelse(s$time < 168, paste0(s$time, "h"), paste0(s$time / 24, "d")),
                                  s$mean_identity, s$r, s$n),
              hoverinfo = "text")
      }
      plotly::layout(p, title = "Subtype identity over time",
        shapes = list(list(type = "line", x0 = 3, x1 = 2200, y0 = 0, y1 = 0,
                           line = list(color = "grey60", dash = "dot", width = 1))),
        xaxis = list(title = "Time after injury (log scale)", type = "log",
                     tickvals = c(6, 24, 72, 168, 672, 2160),
                     ticktext = c("6h", "24h", "72h", "7d", "28d", "90d")),
        yaxis = list(title = "Mean identity (SD from that subtype's naive mean)")) |> ply_pub()
    })

    output$id_scatter <- plotly::renderPlotly({
      p <- PC(); req(p)
      m <- p$meta
      s <- m[!is.na(m$identity) & m$celltype %in% SC$neuronal, ]
      req(nrow(s) > 100)
      s$grp <- ifelse(s$inj == "Naive", "Naive",
               ifelse(s$inj %in% SC$axotomy & s$time >= 24, "Axotomy ≥24h", "other"))
      s <- s[s$grp != "other", ]
      pp <- plotly::plot_ly()
      for (g in c("Naive", "Axotomy ≥24h")) {
        z <- s[s$grp == g, ]; if (!nrow(z)) next
        pp <- plotly::add_trace(pp, x = z$identity, y = z$score, type = "scattergl", mode = "markers",
              name = g, marker = list(size = 3, opacity = .38,
                                      color = if (g == "Naive") "#7F8C8D" else "#B2182B"),
              hoverinfo = "skip")
      }
      plotly::layout(pp, title = "Subtype identity versus panel score, per neuron",
        xaxis = list(title = "Subtype identity (SD from naive)"),
        yaxis = list(title = "Panel score"),
        legend = list(orientation = "h", y = -0.16)) |> ply_pub()
    })

    # =====================================================================
    # Peak landscape
    # =====================================================================
    output$pk_plot <- plotly::renderPlotly({
      req(SC, SC$peakviz)
      d <- SC$peakviz
      d$tnum <- match(as.character(d$tlab), levels(d$tlab))
      key <- factor(as.character(d[[input$pk_colour]]))
      pal <- if (input$pk_colour == "peak_injury") SC_INJ_COL else
             setNames(rep(c("#4E79A7","#F28E2B","#E15759","#76B7B2","#59A14F","#EDC948","#B07AA1",
                            "#FF9DA7","#9C755F","#BAB0AC","#1F77B4","#D62728"),
                          length.out = nlevels(key)), levels(key))
      p <- plotly::plot_ly()
      for (k in levels(key)) {
        s <- d[key == k, ]; if (!nrow(s)) next
        p <- plotly::add_trace(p, x = s$tnum, y = s$log2FC, type = "scatter",
              mode = if (isTRUE(input$pk_labels)) "markers+text" else "markers",
              name = k, marker = list(size = 10, color = unname(pal[k]),
                                      line = list(width = .6, color = "#fff")),
              text = s$lab, textposition = "top center", textfont = list(size = 8.5),
              hovertext = sprintf("<b>%s</b><br>peaks in %s<br>%s at %s<br>log2FC %+.2f (adj.p %.2g)<br>detected in %.2f%% of cells",
                                  s$Gene, s$peak_celltype_full, s$peak_injury, s$tlab,
                                  s$log2FC, s$adj_p, s$pct_cells_detected),
              hoverinfo = "text")
      }
      plotly::layout(p, title = "Peak landscape of the conserved injury panel",
        xaxis = list(title = "Time at which the gene peaks", tickmode = "array",
                     tickvals = seq_along(levels(d$tlab)), ticktext = levels(d$tlab)),
        yaxis = list(title = "Peak log2FC vs Naive"),
        legend = list(orientation = "h", y = -0.18)) |> ply_pub()
    })

    output$pk_note <- renderUI({
      req(SC, SC$peakviz)
      d <- SC$peakviz
      w <- table(d$wave); ct <- sort(table(as.character(d$peak_celltype)), decreasing = TRUE)
      inj <- sort(table(as.character(d$peak_injury)), decreasing = TRUE)
      div(class = "small text-muted mt-2",
        tags$b("Reading the wave. "),
        sprintf("%d genes peak within 24 hours, %d between 36 and 72 hours, %d at one week. ",
                w[["Early (≤24h)"]], w[["Mid (36-72h)"]], w[["1 week"]]),
        sprintf("The earliest are transcription factors and immediate-early genes (Socs3, Jun, Gadd45a, Atf3); the latest are secreted and structural effectors (Npy, Sprr1a, Gap43, Stmn4). "),
        sprintf("Most genes peak in %s (%d genes) and in %s (%d of %d) — the model with no regenerative outlet.",
                names(ct)[1], ct[[1]], names(inj)[1], inj[[1]], nrow(d)))
    })

    output$pk_bars <- plotly::renderPlotly({
      req(SC, SC$peakviz)
      d <- SC$peakviz
      ct <- sort(table(as.character(d$peak_celltype)))
      tt <- table(d$tlab); tt <- tt[tt > 0]
      p1 <- plotly::plot_ly(x = as.numeric(ct), y = names(ct), type = "bar", orientation = "h",
              marker = list(color = "#4E79A7"), name = "cell type", showlegend = FALSE)
      p2 <- plotly::plot_ly(x = names(tt), y = as.numeric(tt), type = "bar",
              marker = list(color = "#B2182B"), name = "timepoint", showlegend = FALSE)
      plotly::subplot(plotly::layout(p1, xaxis = list(title = "genes"), yaxis = list(title = "")),
                      plotly::layout(p2, xaxis = list(title = ""), yaxis = list(title = "genes")),
                      nrows = 2, heights = c(.55, .45), margin = .1, titleX = TRUE, titleY = TRUE) |>
        plotly::layout(title = "Genes by peak cell type (top) and peak timepoint (bottom)") |> ply_pub()
    })

    # =====================================================================
    # MINP genes in single cell
    # =====================================================================
    output$mg_note <- renderUI({
      req(SC, SC$minp_sc)
      m <- SC$minp_sc; det <- m$detection
      core <- det[det$gene %in% c("Myh7", "Slc15a2"), ]
      div(class = "banner", style = "border-left-color:#f39c12",
        tags$b("Treat the MINP-specific gene set as provisional. "),
        sprintf("Only %d of the %d MINP-associated rows carry a gene symbol at all — the rest are unannotated Ensembl loci, and they are also the rows with the most extreme fold-changes and non-significant pooled p-values, which is the signature of low-count noise. ",
                sum(m$minp$named), nrow(m$minp)),
        if (nrow(core)) sprintf("Of the two genes the app calls 'core', %s is detected in %.2f%% of DRG nuclei (%s UMI across all 141,093 cells) and %s in %.2f%%. ",
                core$gene[core$gene == "Myh7"], core$pct_cells[core$gene == "Myh7"],
                format(core$total_umi[core$gene == "Myh7"], big.mark = ","),
                core$gene[core$gene == "Slc15a2"], core$pct_cells[core$gene == "Slc15a2"]) else "",
        "None of them is injury-responsive. A bulk call of MINP-specific on a gene that is near-undetectable ",
        "across the whole DRG should be confirmed by a targeted assay before anything is built on it.")
    })

    output$mg_plot <- plotly::renderPlotly({
      req(SC, SC$minp_sc)
      b <- SC$minp_sc$by_celltype
      cts <- intersect(names(SC$celltype_full), unique(b$celltype))
      gs <- unique(b$gene)
      M <- matrix(NA_real_, length(gs), length(cts), dimnames = list(gs, cts))
      M[cbind(match(b$gene, gs), match(b$celltype, cts))] <- b$pct
      p <- plotly::plot_ly(x = cts, y = gs, z = M, type = "heatmap",
            colorscale = list(list(0, "#F7F7F7"), list(1, "#B2182B")),
            text = matrix(sprintf("%.2f%% of cells", M), nrow(M)), hoverinfo = "x+y+text",
            colorbar = list(title = list(text = "% cells<br>expressing")))
      plotly::layout(p, title = "MINP-associated genes across DRG cell types",
        xaxis = list(tickangle = -40), yaxis = list(title = "")) |> ply_pub()
    })

    output$mg_tab <- DT::renderDT({
      req(SC, SC$minp_sc)
      DT::datatable(SC$minp_sc$detection, rownames = FALSE, extensions = "Buttons",
        options = list(pageLength = 10, dom = "Bfrtip", scrollX = TRUE,
          buttons = list(list(extend = "excel", text = "⤓ Excel", filename = "minp_genes_singlecell"))))
    })

    # =====================================================================
    # MINP beyond the injury axis
    # =====================================================================
    output$by_note <- renderUI({
      req(SC, SC$beyond)
      b <- SC$beyond; g <- b$gsea
      nsig <- sum(g$contrast != "SNI vs Sham (control)" & !is.na(g$padj) & g$padj < 0.05)
      div(class = "banner", style = "border-left-color:#f39c12",
        tags$b("A cautionary result — read the second plot before the first. "),
        sprintf("Taken at face value, %d cell-type marker sets are significantly shifted in the MINP contrasts, which looks like a real compositional or glial signal. ", nsig),
        "It is not. ",
        sprintf("The male and female MINP contrasts are almost perfectly ANTI-correlated across the cell-type sets (r = %+.3f): ", b$r_male_female),
        "in males the non-neuronal markers fall and neuronal markers rise, in females the exact ",
        "reverse. A real MINP effect cannot be anti-correlated with itself. ",
        tags$br(), tags$br(),
        tags$b("What this most likely is: "), "how much nerve and connective tissue each dissection ",
        "happened to capture alongside the ganglion, differing by chance between groups of three ",
        "animals. A competitive gene-set test aggregates ~50 markers at once, so it detects that ",
        "wobble easily even when no individual gene is significant. ",
        tags$br(), tags$br(),
        tags$b("The method works — the positive control confirms it. "),
        "In SNI, both injury reference programmes come up strongly (the atlas injury programme at ",
        "p ≈ 2×10⁻³¹, the conserved 47-gene panel at p ≈ 2×10⁻¹²). So the absence of a coherent ",
        "MINP signal is a real absence, not a failure to look.")
    })

    output$by_plot <- plotly::renderPlotly({
      req(SC, SC$beyond)
      g <- SC$beyond$gsea
      g$lab <- sub("^CT: ", "", sub("^REF: ", "", g$pathway))
      ref <- g$pathway[grepl("^REF:", g$pathway)]
      ord <- g$lab[g$contrast == "SNI vs Sham (control)"][order(g$score[g$contrast == "SNI vs Sham (control)"])]
      lv <- unique(ord)
      p <- plotly::plot_ly()
      for (cn in unique(g$contrast)) {
        s <- g[g$contrast == cn, ]; s <- s[match(lv, s$lab), ]
        p <- plotly::add_trace(p, x = s$score, y = factor(lv, levels = lv), type = "bar",
              orientation = "h", name = cn,
              marker = list(color = if (grepl("SNI", cn)) "#2166AC" else
                                    if (grepl("female", cn)) "#B2182B" else
                                    if (grepl("male", cn)) "#E8A33D" else "#7F8C8D"),
              hovertext = sprintf("<b>%s</b><br>%s<br>%s, p = %.2g (adj %.2g)", s$lab, cn,
                                  s$Direction, s$pval, s$padj),
              hoverinfo = "text")
      }
      plotly::layout(p, barmode = "group",
        title = "Cell-type marker enrichment by contrast (signed -log10 p)",
        xaxis = list(title = "← set shifted DOWN      |      set shifted UP →"),
        yaxis = list(title = ""),
        legend = list(orientation = "h", y = -0.12)) |> ply_pub()
    })

    output$by_consistency <- plotly::renderPlotly({
      req(SC, SC$beyond, SC$beyond$wide)
      w <- SC$beyond$wide
      w$lab <- sub("^CT: ", "", w$pathway)
      lim <- max(abs(c(w$`MINP male`, w$`MINP female`)), na.rm = TRUE)
      p <- plotly::plot_ly()
      for (tp in unique(w$type)) {
        s <- w[w$type == tp, ]
        p <- plotly::add_trace(p, x = s$`MINP male`, y = s$`MINP female`,
              type = "scatter", mode = "markers+text", name = tp,
              marker = list(size = 10, color = if (tp == "neuronal") "#B2182B" else "#2166AC",
                            line = list(width = .5, color = "#fff")),
              text = s$lab, textposition = "top center", textfont = list(size = 8.5),
              hovertext = sprintf("<b>%s</b><br>male %+.2f<br>female %+.2f", s$lab,
                                  s$`MINP male`, s$`MINP female`),
              hoverinfo = "text")
      }
      plotly::layout(p,
        title = sprintf("Do the sexes agree about MINP?  No — r = %+.3f", SC$beyond$r_male_female),
        shapes = list(
          list(type = "line", x0 = -lim, x1 = lim, y0 = -lim, y1 = lim,
               line = list(color = "#59A14F", dash = "dash", width = 1.4)),
          list(type = "line", x0 = -lim, x1 = lim, y0 = lim, y1 = -lim,
               line = list(color = "grey65", dash = "dot", width = 1))),
        xaxis = list(title = "MINP male — signed -log10 p", range = c(-lim, lim)),
        yaxis = list(title = "MINP female — signed -log10 p", range = c(-lim, lim)),
        legend = list(orientation = "h", y = -0.14)) |> ply_pub()
    })

    output$by_foot <- renderUI({
      div(class = "small text-muted mt-2",
        "The green dashed line is where the points would lie if the two sexes agreed; the grey ",
        "dotted line is perfect disagreement. The points follow the grey line. ",
        tags$br(), tags$br(),
        tags$b("Method note. "), "Ranking uses the DESeq2 Wald statistic from the full-precision ",
        "DE_results tables, not the app's de_long.rds — that table retains only about 1,684 ",
        "distinct fold-change values across 20,830 genes, which leaves ~92% of the ranking tied and ",
        "is too coarse for a rank-based test. The test is limma's cameraPR, which is competitive, ",
        "tolerant of ties and corrects for inter-gene correlation. SNI has no full-precision table, ",
        "so its ranking is the coarse one — that costs power only, so its positive result stands.")
    })

    # =====================================================================
    # Atf3 — per-cell
    # =====================================================================
    output$ap_note <- renderUI({
      req(SC, SC$atf3_pc)
      p <- SC$atf3_pc$pooled
      w168 <- p[p$state == "Crush 168h", ]; w36 <- p[p$state == "Crush 36h", ]
      div(class = "banner", style = "border-left-color:#B2182B",
        tags$b("Losing Atf3 largely blocks entry into the injured state. "),
        sprintf("At 7 days after crush, %.1f%% of wild-type neurons have crossed the injured-state threshold against %.1f%% of Atf3-KO neurons — a %.1f point difference. At 36 hours it is %.1f%% versus %.1f%%. ",
                w168$pct[w168$geno == "Atf3-WT"], w168$pct[w168$geno == "Atf3-KO"],
                w168$pct[w168$geno == "Atf3-WT"] - w168$pct[w168$geno == "Atf3-KO"],
                w36$pct[w36$geno == "Atf3-WT"], w36$pct[w36$geno == "Atf3-KO"]),
        "Knockout neurons partially initiate the response but cannot sustain it.",
        tags$br(), tags$br(),
        tags$b("Why this differs from the pseudobulk tab. "),
        "The sample-level interaction test has three samples per genotype per state and finds almost ",
        "nothing. But injury induction is overwhelmingly ", tags$b("recruitment"),
        " — a per-cell quantity estimated here from 17,665 neurons. The pseudobulk analysis was not ",
        "wrong, it was simply the wrong test for this dataset.",
        tags$br(), tags$br(),
        tags$b("Two caveats stand. "), "Atf3 mRNA is not reduced in the knockout (the allele still ",
        "yields a transcript that 3'-end sequencing counts), so the genotype cannot be verified from ",
        "these data; and a proportion test over thousands of cells from 1–3 mice overstates ",
        "significance, so read the effect sizes rather than the p-values.")
    })

    output$ap_entry <- plotly::renderPlotly({
      req(SC, SC$atf3_pc)
      e <- SC$atf3_pc$entry
      cts <- intersect(names(SC$celltype_full), unique(e$celltype))
      e$celltype <- factor(e$celltype, levels = cts)
      p <- plotly::plot_ly()
      for (g in c("Atf3-WT", "Atf3-KO")) {
        s <- e[e$geno == g, ]
        p <- plotly::add_trace(p, x = list(as.character(s$celltype), as.character(s$state)),
              y = s$pct, type = "bar", name = g,
              marker = list(color = if (g == "Atf3-WT") "#4E79A7" else "#E15759"),
              hovertext = sprintf("<b>%s</b><br>%s · %s<br>%.1f%% injured (%d of %d cells)",
                                  g, s$celltype, s$state, s$pct, s$n_injured, s$n),
              hoverinfo = "text")
      }
      plotly::layout(p, barmode = "group",
        title = "Neurons entering the injured state — Atf3-WT vs Atf3-KO",
        xaxis = list(title = "", tickangle = -35),
        yaxis = list(title = "% of neurons in the injured state"),
        legend = list(orientation = "h", y = -0.25)) |> ply_pub()
    })

    output$ap_rec <- plotly::renderPlotly({
      req(SC, SC$atf3_pc)
      r <- SC$atf3_pc$recruitment; r <- r[r$state == "Crush 168h", ]
      req(nrow(r) > 0)
      r$blocked_f <- r$blocked > 1
      lim <- range(c(r$wt_rec, r$ko_rec), na.rm = TRUE)
      p <- plotly::plot_ly()
      for (b in c(FALSE, TRUE)) {
        s <- r[r$blocked_f == b, ]; if (!nrow(s)) next
        p <- plotly::add_trace(p, x = s$wt_rec, y = s$ko_rec, type = "scatter", mode = "markers+text",
              name = if (b) "recruitment blocked without Atf3" else "recruited normally",
              marker = list(size = 9, color = if (b) "#B2182B" else "#4E79A7",
                            line = list(width = .5, color = "#fff")),
              text = ifelse(s$wt_rec > 1 | s$blocked > 1, s$gene, ""),
              textposition = "top center", textfont = list(size = 8),
              hovertext = sprintf("<b>%s</b><br>WT %+.2f (%.1f%% of cells)<br>KO %+.2f (%.1f%% of cells)<br>blocked %+.2f",
                                  s$gene, s$wt_rec, s$wt_pct, s$ko_rec, s$ko_pct, s$blocked),
              hoverinfo = "text")
      }
      plotly::layout(p, title = "Recruitment of panel genes at 7 days: WT versus KO",
        shapes = list(list(type = "line", x0 = lim[1], x1 = lim[2], y0 = lim[1], y1 = lim[2],
                           line = list(color = "grey55", dash = "dash", width = 1.2))),
        xaxis = list(title = "Recruitment in Atf3-WT (log2 fraction expressing)"),
        yaxis = list(title = "Recruitment in Atf3-KO"),
        legend = list(orientation = "h", y = -0.16)) |> ply_pub()
    })

    # ---- 3b. UMAP ----------------------------------------------------------
    output$umap <- plotly::renderPlotly({
      req(SC, SC$umap)
      d <- SC$umap
      # scattergl (WebGL), not scatter — an SVG trace with ~28k points would crawl
      if (input$umap_by == "score") {
        p <- plotly::plot_ly(d, x = ~x, y = ~y, type = "scattergl", mode = "markers",
              marker = list(size = 3, opacity = .55, color = ~score, colorscale = SC_DIVERGE,
                            cmid = 0, colorbar = list(title = list(text = "panel<br>score"))),
              text = ~sprintf("%s<br>%s %sh<br>score %+.2f", celltype, inj, time, score),
              hoverinfo = "text")
      } else {
        d$.g <- if (input$umap_by == "time") factor(d$time, levels = sort(unique(d$time)))
                else factor(d[[input$umap_by]])
        p <- plotly::plot_ly()
        for (l in levels(d$.g)) {
          s <- d[d$.g == l, ]
          p <- plotly::add_trace(p, x = s$x, y = s$y, type = "scattergl", mode = "markers",
                name = as.character(l), marker = list(size = 3, opacity = .6),
                text = sprintf("%s<br>%s %sh", s$celltype, s$inj, s$time), hoverinfo = "text")
        }
      }
      plotly::layout(p, title = "GSE154659 DRG atlas — UMAP",
        xaxis = list(title = "UMAP 1", showticklabels = FALSE),
        yaxis = list(title = "UMAP 2", showticklabels = FALSE),
        legend = list(itemsizing = "constant")) |> ply_pub()
    })

    # ---- 3c. MINP specificity ---------------------------------------------
    cc_matrix <- reactive({
      req(SC, SC$minp_cc)
      s <- SC$minp_cc[SC$minp_cc$contrast == input$cc_contrast, ]
      s$col <- paste(s$injury, ifelse(s$time < 168, paste0(s$time, "h"), paste0(s$time / 24, "d")))
      cols <- unique(s[order(match(s$injury, SC$inj_order), s$time), c("col", "injury")])
      cts  <- intersect(names(SC$celltype_full), unique(s$celltype))
      M <- matrix(NA_real_, length(cts), nrow(cols), dimnames = list(cts, cols$col))
      M[cbind(match(s$celltype, cts), match(s$col, cols$col))] <- s$r
      list(M = M, blocks = cols$injury)
    })

    output$minp_cc <- plotly::renderPlotly({
      m <- cc_matrix(); req(m)
      lim <- max(abs(SC$minp_cc$r), na.rm = TRUE)
      draw_heat(m$M, NULL, collab = colnames(m$M), blocks = m$blocks,
        title = sprintf("Bulk %s vs single-cell injury signatures", input$cc_contrast),
        zlab = "Pearson r", clip = round(lim, 2),
        do_row = isTRUE(input$crow4), do_col = isTRUE(input$ccol4),
        rowlab = setNames(rownames(m$M), rownames(m$M)))
    })

    output$minp_note <- renderUI({
      req(SC, SC$minp_cc)
      s <- SC$minp_cc
      mx <- tapply(s$r, s$contrast, max); md <- tapply(s$r, s$contrast, median)
      n_over <- tapply(s$r, s$contrast, function(v) sum(v > 0.2))
      n_tot  <- tapply(s$r, s$contrast, length)
      div(class = "small text-muted mt-2",
        tags$b("Why this matters. "),
        "The bulk result — that MINP is uncorrelated with the conserved injury signature — has an ",
        "obvious escape hatch: a real MINP response confined to one cell type could be diluted below ",
        "detection in whole-DRG bulk. This closes it. ",
        sprintf("Across %d cell type × injury × timepoint combinations, bulk SNI reaches a maximum r of %+.3f (median %+.3f, %d combinations above r = 0.2), while MINP never exceeds %+.3f (median %+.3f, %d above 0.2). ",
                n_tot[["SNI"]], mx[["SNI"]], md[["SNI"]], n_over[["SNI"]],
                mx[["MINP"]], md[["MINP"]], n_over[["MINP"]]),
        tags$b("MINP does not resemble the injury program in any cell type, at any timepoint, in any of the six injury models."))
    })

    # ---- 3d. Atf3 dependence ----------------------------------------------
    output$atf3_note <- renderUI({
      req(SC, SC$atf3)
      a <- SC$atf3
      ind <- !is.na(a$WT_crush_l2fc) & a$WT_crush_l2fc > 1 & a$WT_crush_padj < 0.05
      wa <- a$WT_crush_l2fc[a$Gene == "Atf3"]; ka <- a$KO_crush_l2fc[a$Gene == "Atf3"]
      div(class = "banner", style = "border-left-color:#f39c12",
        tags$b("Read this with two caveats. "),
        sprintf("First, Atf3 mRNA is not reduced in the knockout — crush induces it %+.2f in WT and %+.2f in the KO. ", wa, ka),
        "The allele still produces a transcript that 3'-end single-nucleus sequencing counts, so the ",
        "knockout cannot be verified from these data; everything here is conditional on the published ",
        "genotype labels being correct. ",
        tags$b("Second, this arm has only three pseudobulk samples per genotype per state, "),
        "so the interaction test is badly underpowered. The retention percentages are the readable ",
        "signal; the near-absence of significant interaction p-values reflects sample size, not ",
        "evidence that the panel is Atf3-independent. ",
        sprintf("Of %d panel genes induced by crush in WT, the median retention of that induction in the KO is %.0f%%.",
                sum(ind), median(a$pct_retained_in_KO[ind], na.rm = TRUE)))
    })

    output$atf3_plot <- plotly::renderPlotly({
      req(SC, SC$atf3)
      a <- SC$atf3[!is.na(SC$atf3$WT_crush_l2fc), ]
      a$dep <- !is.na(a$dependence_padj) & a$dependence_padj < 0.05 & a$Atf3_dependence > 0
      lim <- range(c(a$WT_crush_l2fc, a$KO_crush_l2fc), na.rm = TRUE)
      p <- plotly::plot_ly()
      for (lv in c(FALSE, TRUE)) {
        s <- a[a$dep == lv, ]; if (!nrow(s)) next
        p <- plotly::add_trace(p, x = s$WT_crush_l2fc, y = s$KO_crush_l2fc,
              type = "scatter", mode = "markers+text",
              name = if (lv) "Atf3-dependent (adj.p<0.05)" else "not significant",
              marker = list(size = 9, color = if (lv) "#B2182B" else "#4E79A7",
                            line = list(width = .5, color = "#fff")),
              text = s$Gene, textposition = "top center", textfont = list(size = 8),
              hovertext = sprintf("<b>%s</b><br>WT %+.2f<br>KO %+.2f<br>dependence %+.2f",
                                  s$Gene, s$WT_crush_l2fc, s$KO_crush_l2fc, s$Atf3_dependence),
              hoverinfo = "text")
      }
      plotly::layout(p,
        title = "Crush-induced fold-change: Atf3-WT vs Atf3-KO",
        shapes = list(list(type = "line", x0 = lim[1], x1 = lim[2], y0 = lim[1], y1 = lim[2],
                           line = list(color = "grey45", dash = "dash", width = 1.2))),
        xaxis = list(title = "log2FC crush vs naive — Atf3-WT"),
        yaxis = list(title = "log2FC crush vs naive — Atf3-KO"),
        legend = list(orientation = "h", y = -0.15)) |> ply_pub()
    })

    output$atf3_tab <- DT::renderDT({
      req(SC, SC$atf3)
      DT::datatable(SC$atf3, rownames = FALSE, filter = "top", extensions = "Buttons",
        options = list(pageLength = 12, scrollX = TRUE, dom = "Bfrtip",
          order = list(list(6, "desc")),
          buttons = list(list(extend = "excel", text = "⤓ Excel", filename = "atf3_dependency"))))
    })

    # ---- 3e. deconvolution reference --------------------------------------
    deconv_mouse <- reactive({
      req(SC, SC$deconv)
      l <- SC$deconv$long
      l <- l[l$species == "mouse" & l$celltype %in% SC$stats$no_baseline, ]
      req(nrow(l) > 0)
      l$group <- sub("_.*", "", l$sample)
      aggregate(prop ~ reference + celltype + group, data = l, FUN = mean)
    })

    output$deconv_note <- renderUI({
      a <- deconv_mouse(); req(a)
      n <- a[a$reference == "naive_fixed" & a$celltype == "Repair schwann_N", ]
      v <- setNames(100 * n$prop, n$group)
      div(class = "banner", style = "border-left-color:#B2182B",
        tags$b("No — and the answer is the reverse of what you would expect. "),
        "The injury-state reference reports ", tags$b("0.000%"), " repair cells in every mouse ",
        "sample, including the two SNI samples. The naive-only reference does better: it puts ",
        tags$i("Repair schwann_N"),
        sprintf(" at %.2f%% in SNI against %.2f%% in Sham and %.2f%% in MINP — a %.0f-fold elevation, in the one condition where axotomy actually happened. ",
                v[["SNI"]], v[["Sham"]], v[["NP"]], v[["SNI"]] / max(v[["Sham"]], 1e-6)),
        tags$br(), tags$br(),
        tags$b("Why. "), "MuSiC solves a weighted non-negative least squares. Signatures built from ",
        "injured cells are far more extreme, and with 20 competing cell types the solver drives the ",
        "rare repair components to exactly zero while other types absorb the signal — NF3 rises from ",
        "2.09% to 6.68% and PEP1 collapses from 2.88% to 0.03%. ",
        tags$br(), tags$br(),
        tags$b("What this means for the Deconvolution tab. "), "Keep the naive-only reference; it is ",
        "the better choice. The caveat worth carrying is that absolute proportions are sensitive to ",
        "how the reference is built, so they should be read as relative comparisons across samples, ",
        "never as absolute tissue composition.")
    })

    output$deconv_plot <- plotly::renderPlotly({
      a <- deconv_mouse(); req(a)
      a$lab <- ifelse(a$group == "NP", "MINP", a$group)
      p <- plotly::plot_ly()
      for (rf in unique(a$reference)) {
        s <- a[a$reference == rf, ]
        p <- plotly::add_trace(p, x = paste(s$celltype, s$lab), y = 100 * s$prop,
              type = "bar", name = if (rf == "injury") "injury-state reference" else "naive-only reference",
              marker = list(color = if (rf == "injury") "#E15759" else "#4E79A7"),
              hovertext = sprintf("%s<br>%s<br>%.3f%%", s$celltype, s$lab, 100 * s$prop),
              hoverinfo = "text")
      }
      plotly::layout(p, barmode = "group",
        title = "Injury-induced cell types estimated in the mouse bulk",
        xaxis = list(title = "", tickangle = -30),
        yaxis = list(title = "Estimated proportion of bulk (%)"),
        legend = list(orientation = "h", y = -0.28)) |> ply_pub()
    })

    # ---- 4. peak table -----------------------------------------------------
    output$peak <- DT::renderDT({
      req(SC)
      d <- SC$peak
      d <- d[, intersect(c("Gene", "role", "peak_celltype", "peak_celltype_full", "peak_injury",
                           "peak_time", "log2FC", "adj_p", "pct_cells_detected"), names(d))]
      names(d) <- c("Gene", "Role", "Peak cell type", "Cell type (full)", "Injury",
                    "Time", "log2FC", "adj.p", "% cells detected")[seq_len(ncol(d))]
      DT::datatable(d, rownames = FALSE, filter = "top", extensions = "Buttons",
        options = list(pageLength = 15, scrollX = TRUE, dom = "Bfrtip",
          order = list(list(6, "desc")),
          buttons = list(list(extend = "excel", text = "⤓ Excel",
                              filename = "sc_injury_panel_peaks"))))
    })
    output$dl_peak <- downloadHandler(
      filename = function() "sc_injury_panel_peaks.csv",
      content  = function(f) utils::write.csv(SC$peak, f, row.names = FALSE))

    # ---- 5. documentation --------------------------------------------------
    output$methods <- renderUI({
      req(SC); s <- SC$stats
      div(class = "py-2",
        h4("What this analysis is"),
        tags$p("The bulk cross-species work identified 47 genes induced by nerve injury in at least ",
               "two of three species (mouse SNI, rat SNI, macaque ipsilateral) while remaining flat ",
               "in MINP. That result says nothing about which cells respond, or how the response ",
               "evolves. This tab tests the panel against a single-cell atlas that has both."),
        h4("Data"),
        tags$ul(
          tags$li(tags$b("GEO GSE154659"), " — Renthal et al. 2020, mouse DRG single-nucleus RNA-seq. ",
                  sprintf("%s cells, %d samples, %d author-assigned cell types.",
                          format(s$n_cells, big.mark = ","), s$n_samples, s$n_celltypes)),
          tags$li(tags$b("Six injury models"), ": sciatic nerve crush, sciatic nerve transection (ScNT), ",
                  "spinal nerve transection (SpNT), CFA inflammation, paclitaxel neuropathy, and naive. ",
                  sprintf("%d injury × timepoint combinations, spanning 6 hours to 90 days.", s$n_groups)),
          tags$li("This is the same atlas the ", tags$b("Deconvolution"), " tab uses as its reference, ",
                  "but that tab uses only the naive cells; here the full injury time-course is used.")
        ),
        h4("Method"),
        tags$ul(
          tags$li(tags$b("Pseudobulk, not per-cell tests."), " Raw UMI counts are summed within each ",
                  "sample (and within each sample × cell type), then normalised with edgeR TMM and ",
                  "modelled with limma. Summing counts before testing treats the biological replicate ",
                  "as the unit of replication, which per-cell tests do not."),
          tags$li(tags$b("Design."), " ~0 + injury×time group + sex + genotype, with every group ",
                  "contrasted against Naive. Moderated variance means the three single-replicate ",
                  "columns (CFA 48h, ScNT 6h, ScNT 24h) still yield estimates, but they carry less weight."),
          tags$li(tags$b("Cell-type fits are independent."), " Each cell type gets its own library sizes, ",
                  "expression filter and dispersion rather than sharing one model."),
          tags$li(tags$b("Minimum 10 cells"), " per sample × cell type unit; smaller units are dropped ",
                  "rather than contributing an unstable estimate.")
        ),
        h4("Verification"),
        tags$ul(
          tags$li(tags$b("Positive control. "), "Every canonical regeneration-associated gene is induced ",
                  "by axotomy: Sprr1a +8.07, Npy +5.98, Atf3 +5.05, Sox11 +3.51, Gal +3.33, Gap43 +1.84."),
          tags$li(tags$b("Negative control. "), sprintf(
                  "Zero of %d modelled panel genes are significantly induced in CFA or paclitaxel, ",
                  s$panel_modelled),
                  sprintf("against %.0f%% in axotomy — a %.0f-fold difference in mean fold-change.",
                          s$axotomy$pct_sig_up, s$axotomy$mean_l2fc / max(s$non_axotomy$mean_l2fc, 1e-6))),
          tags$li(tags$b("Independent bulk agreement. "), sprintf(
                  "Single-cell SpNT-7d versus this study's bulk mouse SNI-vs-Sham: r = %.3f (Spearman %.3f), %.0f%% same-direction across %d genes.",
                  s$bulk_r, s$bulk_rho, 100 * s$bulk_same_dir, s$bulk_n)),
          tags$li(tags$b("Sensitivity. "), sprintf(
                  "The full model uses all cells with sex and genotype as covariates. Refitting on male C57 cells alone gives r = %.3f, so those covariates are not absorbing real signal.",
                  s$sensitivity_r))
        ),
        h4("Embedding (UMAP tab)"),
        tags$ul(
          tags$li("scanpy: counts normalised to 10,000 per nucleus, log1p, 2,000 highly variable ",
                  "genes (seurat_v3 on raw counts), 50 principal components, UMAP on a ",
                  "15-nearest-neighbour graph."),
          tags$li(tags$b("Cell types are the authors' labels"), ", not re-clustered here. The atlas ",
                  "ships them, and inventing new ones would only make this analysis harder to ",
                  "compare against the published work."),
          tags$li(tags$b("No batch integration, by choice. "), "Each sample belongs to exactly one ",
                  "injury × time group, so `sample` is perfectly nested within condition. ",
                  "Integrating on it would regress out the condition effect together with any ",
                  "technical batch effect. GEO supplies no independent technical batch variable to ",
                  "correct on instead, so the embedding is left un-integrated and labelled as such."),
          tags$li(tags$b("QC removed no cells. "), "The authors pre-filtered — observed floors are ",
                  "566 UMI and 501 genes, with UMI capped at 14,999. Median mitochondrial content is ",
                  "0.4% because these are nuclei, so a mitochondrial filter is not informative. ",
                  "Our QC is confirmatory.")
        ),
        h4("Follow-on analyses"),
        tags$ul(
          tags$li(tags$b("MINP specificity. "), "The bulk finding that MINP is uncorrelated with the ",
                  "injury signature could in principle be a dilution artefact — a response confined to ",
                  "one cell type would be diluted in whole-DRG bulk. Correlating the bulk log2FC ",
                  "vectors against every cell type × injury × timepoint signature, over all shared ",
                  "genes rather than just the 47, rules that out."),
          tags$li(tags$b("Atf3 dependence. "), "A loss-of-function test using the Atf3-WT / Atf3-KO ",
                  "arm. Two hard limits apply: the KO transcript is still produced and counted, so the ",
                  "knockout is unverifiable from these data, and with three samples per genotype per ",
                  "state the interaction test is badly underpowered."),
          tags$li(tags$b("Injury-state deconvolution reference. "), "The ", tags$b("Deconvolution"),
                  " tab builds its MuSiC reference from naive cells only. That reference is ",
                  "structurally unable to report ", tags$i("Repair schwann_N"), " or ",
                  tags$i("Repair fibroblast"), ", which do not exist in naive tissue — the naive atlas ",
                  "contains only 6 repair fibroblasts, below any usable threshold. A reference built ",
                  "from injured cells can report them.")
        ),
        h4("Per-cell analyses"),
        tags$p("Everything above this point is pseudobulk: each sample collapsed to a mean. That is ",
               "the right unit for differential expression, but it discards the distribution across ",
               "cells, and the distribution turns out to carry the mechanism."),
        tags$ul(
          tags$li(tags$b("Recruitment versus level. "), "With f = the fraction of cells expressing a ",
                  "gene and m = its mean among those cells, the mean over all cells is f × m, so the ",
                  "log2 fold-change splits exactly into a recruitment term and a level term. ",
                  "Almost all of the panel's induction is recruitment: genes switch on in cells that ",
                  "were silent. In NF1 neurons after spinal nerve transection, Atf3 goes from 0.2% to ",
                  "63% of cells expressing it while the amount per expressing cell does not move. ",
                  "Pseudobulk reports that as log2FC ≈ +5 and cannot distinguish it from every cell ",
                  "rising uniformly — a materially different biological claim."),
          tags$li(tags$b("Injured-state classifier. "), "The per-cell panel score, thresholded at the ",
                  "99th percentile of naive, separates naive from injured sensory neurons with ",
                  "AUC 0.978 for spinal nerve transection. CFA and paclitaxel sit near chance, so the ",
                  "specificity result holds at single-cell resolution and not only in group means. ",
                  "The recruitment curves also show a clear order: peptidergic nociceptors activate ",
                  "within 6 hours, myelinated neurons lag to 12–24 hours, and glial and immune ",
                  "populations barely move — independent confirmation that their apparent signal in ",
                  "the cell-type heatmap is ambient RNA."),
          tags$li(tags$b("Subtype identity. "), "Each neuron is scored on its own subtype's markers, ",
                  "derived from naive cells only by detection-rate difference so that the markers are ",
                  "not chosen with the answer in view. Identity is standardised within subtype against ",
                  "that subtype's naive distribution."),
          tags$li(tags$b("Gene explorer. "), "A 201-gene universe — the 49 panel rows, the ",
                  "MINP-associated genes, and the top 10 naive markers per cell type — with per-cell ",
                  "CP10K + log1p values, matching the transform used for the panel score. Points are a ",
                  "stratified subsample of about 21,000 cells so the plots stay responsive; all ",
                  "statistics elsewhere are computed on the full 141,093."),
          tags$li(tags$b("Zero inflation. "), "Panel genes are detected in a small minority of cells, ",
                  "so a density over all cells would be a single spike at zero. Ridges are therefore ",
                  "drawn over expressing cells with the percentage expressing printed beside each ",
                  "group; violin and box views show the zeros directly.")
        ),
        h4("A label bug this work corrected"),
        tags$p("The deconvolution script derives each cell type from the second-to-last ",
               "underscore-delimited field of the barcode. Four cell types contain a separator, so ",
               "that shortcut truncates ", tags$i("p_cLTMR2"), " to ", tags$i("cLTMR2"), ", ",
               tags$i("Schwann_M"), " and ", tags$i("Schwann_N"), " to ", tags$i("M"), " and ",
               tags$i("N"), ", and ", tags$i("Repair schwann_N"), " to ", tags$i("schwann"),
               ". Those truncated names propagated into the app's cell-type table. The effect on the ",
               "published deconvolution numbers is small — it used naive cells only, with a 20-cell ",
               "minimum — but the labels were wrong. This tab uses the corrected 20-type nomenclature."),
        h4("Caveats"),
        tags$ul(
          tags$li(tags$b("Ambient RNA. "), "Apparent induction of neuron-restricted genes in non-neuronal ",
                  "cell types is ambient leakage, not expression. Cell-type comparisons are relative."),
          tags$li(tags$b("Shallow coverage. "), "Median 846 genes per nucleus. Ten panel genes are ",
                  "detected in under 0.5% of cells and are marked ", tags$b("*"), "; ",
                  tags$i("Pappa2"), " has 29 total UMI across all 141,093 cells and should be read as absent."),
          tags$li(tags$b("Symbol drift. "), "The atlas uses 2020-era symbols. ", tags$i("Insyn2a"),
                  " appears as ", tags$i("Fam196a"), "; ", tags$i("Crybg1"), " (formerly ", tags$i("Aim1"),
                  ") is absent entirely, so 46 of the 47 are resolvable and ",
                  sprintf("%d survive the expression filter.", s$panel_modelled)),
          tags$li(tags$b("Sex imbalance. "), "The atlas is 95% male; female cells appear in only two of ",
                  "the six injury models. Sex effects cannot be assessed here."),
          tags$li(tags$b("Different study. "), "GSE154659 is an independent experiment on different ",
                  "animals. Agreement with the bulk data is corroboration, not a technical replicate.")
        ),
        h4("Provenance"),
        tags$p(class = "small text-muted",
               sprintf("%s · built %s · %s", SC$provenance$dataset, SC$provenance$built,
                       SC$provenance$scripts))
      )
    })
  })
}
