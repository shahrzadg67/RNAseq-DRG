# GSEA tab (BATCH-CORRECTED) — parallel to tab_gsea.R, reads gsea_bc.rds and
# the gsea_bc/ + pathview_bc/ image trees written by 51_gsea_bc.R.
tab_gsea_bc_ui <- function(id) {
  ns <- NS(id)
  nav_panel(
    title = "GSEA (batch-corrected)", icon = icon("diagram-project"),
    div(class = "container-fluid py-3",
      h3("Gene Set Enrichment Analysis — batch-corrected"),
      div(class = "banner",
          "GSEA over GO (BP/CC/MF) and KEGG, ranked by the batch-adjusted DESeq2 statistic ",
          "(replicate modelled as a covariate; SNI excluded). Positive NES = enriched among ",
          "up-regulated genes; negative = down-regulated. Click a row in the table to open its ",
          "running-score plot. Gene-set selections are exploratory."),
      layout_sidebar(
        sidebar = sidebar(width = 320,
          radioButtons(ns("level"), "Level", c("Gene"="gene","Transcript"="transcript"), inline = TRUE),
          selectInput(ns("contrast"), "Contrast", choices = CONTRAST_LABELS),
          radioButtons(ns("cat"), "Gene-set collection",
                       c("GO: Biological Process"="BP","GO: Cellular Component"="CC",
                         "GO: Molecular Function"="MF","KEGG pathways"="KEGG"))
        ),
        navset_tab(
          nav_panel("Enrichment dotplot", uiOutput(ns("dot_wrap"))),
          nav_panel("Table + running score",
            div(class="row",
              div(class="col-md-7", uiOutput(ns("tbl_wrap"))),
              div(class="col-md-5", uiOutput(ns("plot_panel")))))
        )
      )
    )
  )
}

tab_gsea_bc_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    gsea <- load_rds("gsea_bc")
    cname <- reactive(names(CONTRAST_LABELS)[match(input$contrast, CONTRAST_LABELS)])
    key   <- reactive(paste(input$level, cname(), input$cat, sep = "::"))
    res   <- reactive({ req(gsea); gsea[[key()]] })

    output$dot_wrap <- renderUI({
      if (is.null(gsea)) return(coming_soon("Batch-corrected GSEA results pending."))
      if (is.null(res())) return(div(class="banner","No enriched sets for this selection."))
      plotly::plotlyOutput(ns("dot"), height = "560px")
    })
    output$dot <- plotly::renderPlotly({
      d <- res(); req(!is.null(d))
      d <- d[order(d$p.adjust), ]; d <- head(d, 25)
      d$Description <- factor(d$Description, levels = rev(d$Description))
      plotly::plot_ly(d, x = ~NES, y = ~Description,
        size = ~setSize, color = ~ -log10(p.adjust),
        type = "scatter", mode = "markers",
        marker = list(sizemode = "diameter", sizeref = max(d$setSize)/30),
        text = ~paste0("<b>", Description, "</b><br>", ID,
                       "<br>NES: ", round(NES,2), "<br>adj.p: ", signif(p.adjust,3),
                       "<br>size: ", setSize),
        hoverinfo = "text") |>
        plotly::layout(title = paste0(input$cat, " — ", input$contrast, " (batch-corrected)"),
                       xaxis = list(title = "NES"),
                       yaxis = list(title = ""),
                       margin = list(l = 260))
    })

    output$tbl_wrap <- renderUI({
      if (is.null(gsea) || is.null(res())) return(coming_soon("GSEA table pending."))
      DT::DTOutput(ns("tbl"))
    })
    output$tbl <- DT::renderDT({
      d <- res(); req(!is.null(d))
      d <- d[order(d$p.adjust), c("ID","Description","setSize","NES","pvalue","p.adjust")]
      DT::datatable(d, rownames = FALSE, selection = "single", filter = "top",
                    options = list(pageLength = 12, scrollX = TRUE)) |>
        DT::formatSignif(c("NES","pvalue","p.adjust"), 3)
    })

    output$plot_panel <- renderUI({
      d <- res(); sel <- input$tbl_rows_selected
      if (is.null(d)) return(NULL)
      if (is.null(sel)) return(div(class="banner","Select a row to view its running-score plot."))
      d <- d[order(d$p.adjust), ]
      sid <- d$ID[sel]
      png <- file.path("gsea_bc", input$level, cname(), input$cat,
                       paste0(gsub("[:/]","_", sid), ".png"))
      ui <- tagList(h5(d$Description[sel]),
                    if (file.exists(file.path("www", png)))
                      tags$img(src = png, style="max-width:100%;border:1px solid #eee;border-radius:6px")
                    else div(class="text-muted","Running-score plot not precomputed for this set."))
      if (input$cat == "KEGG") {
        kegg_png <- file.path("pathview_bc", cname(), paste0(sid, ".pathview.png"))
        if (file.exists(file.path("www", kegg_png)))
          ui <- tagList(ui, h6("KEGG pathway map"),
                        tags$img(src = kegg_png, style="max-width:100%;border:1px solid #eee;border-radius:6px"))
      }
      ui
    })
  })
}
