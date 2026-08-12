# GSEA tab (BATCH-CORRECTED, SEX CHROMOSOMES REMOVED) ------------------------
# Same as tab_gsea_bc but reads gsea_noxy.rds and the gsea_noxy/ + pathview_noxy/
# image trees written by the noxy pipeline.
tab_gsea_noxy_ui <- function(id) {
  ns <- NS(id)
  nav_panel(
    title = "GSEA (corrected)", icon = icon("diagram-project"),
    div(class = "container-fluid py-3",
      h3("Gene Set Enrichment — batch-corrected, sex chromosomes removed"),
      div(class = "banner",
          "GSEA over GO (BP/CC/MF) and KEGG, ranked by the batch-adjusted DESeq2 statistic with ",
          "chrX/chrY genes removed. Click a table row to open its running-score plot."),
      uiOutput(ns("sni_note")),
      layout_sidebar(
        sidebar = sidebar(width = 320,
          radioButtons(ns("level"), "Level", c("Gene"="gene","Transcript"="transcript"), inline = TRUE),
          selectInput(ns("contrast"), "Contrast", choices = unname(CONTRAST_LABELS)),
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

tab_gsea_noxy_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    gsea <- load_rds("gsea_noxy")
    cname <- reactive(names(CONTRAST_LABELS)[match(input$contrast, CONTRAST_LABELS)])
    output$sni_note <- renderUI(sni_note(cname()))
    key   <- reactive(paste(input$level, cname(), input$cat, sep = "::"))
    res   <- reactive({ req(gsea); gsea[[key()]] })

    output$dot_wrap <- renderUI({
      if (is.null(gsea)) return(coming_soon("Appears once the sex-chr-removed pipeline finishes."))
      if (is.null(res())) return(div(class="banner","No enriched sets for this selection."))
      plotly::plotlyOutput(ns("dot"), height = "560px")
    })
    output$dot <- plotly::renderPlotly({
      d <- res(); req(!is.null(d)); d <- d[order(d$p.adjust), ]; d <- head(d, 25)
      d$Description <- factor(d$Description, levels = rev(d$Description))
      d$msize <- scales::rescale(d$setSize, to = c(7, 22))
      plotly::plot_ly(d, x = ~NES, y = ~Description, type = "scatter", mode = "markers",
        marker = list(size = ~msize, color = ~ -log10(p.adjust), colorscale = "Viridis",
                      showscale = TRUE, colorbar = list(title = "-log10 adj.p"),
                      line = list(width = 0.5, color = "#333")),
        text = ~paste0("<b>", Description, "</b><br>", ID, "<br>NES: ", round(NES,2),
                       "<br>adj.p: ", signif(p.adjust,3), "<br>size: ", setSize), hoverinfo = "text") |>
        plotly::layout(title = paste0(input$cat, " — ", input$contrast, " (sex-chr removed)"),
                       xaxis = list(title = "NES", showgrid = FALSE, zeroline = FALSE),
                       yaxis = list(title = ""), font = list(size = 14), margin = list(l = 260))
    })

    output$tbl_wrap <- renderUI({
      if (is.null(gsea) || is.null(res())) return(coming_soon("GSEA table pending."))
      DT::DTOutput(ns("tbl"))
    })
    output$tbl <- DT::renderDT({
      d <- res(); req(!is.null(d))
      d <- d[order(d$p.adjust), c("ID","Description","setSize","NES","pvalue","p.adjust")]
      DT::datatable(d, rownames = FALSE, selection = "single", filter = "top", extensions = "Buttons",
        colnames = c("ID"="ID","Description"="Description","Set size"="setSize","NES"="NES",
                     "p-value"="pvalue","adj.p"="p.adjust"),
        options = list(pageLength = 12, scrollX = TRUE, dom = "Bfrtip",
          buttons = list(list(extend = "excel", text = "⤓ Excel",
                              filename = paste0("GSEA_corrected_", input$level, "_", disp(cname()), "_", input$cat))))) |>
        DT::formatSignif(4:6, 3)   # NES, p-value, adj.p (by position — headers are renamed)
    })

    output$plot_panel <- renderUI({
      d <- res(); sel <- input$tbl_rows_selected
      if (is.null(d)) return(NULL)
      if (is.null(sel)) return(div(class="banner","Select a row to view its running-score plot."))
      d <- d[order(d$p.adjust), ]; sid <- d$ID[sel]
      png <- file.path("gsea_noxy", input$level, cname(), input$cat, paste0(gsub("[:/]","_", sid), ".png"))
      ui <- tagList(h5(d$Description[sel]),
                    if (file.exists(file.path("www", png)))
                      tags$img(src = png, style="max-width:100%;border:1px solid #eee;border-radius:6px")
                    else div(class="text-muted","Running-score plot not precomputed for this set."),
                    h6(sprintf("Leading-edge genes (%d) — the genes driving this term",
                               d$core_n[sel] %||% 0L)),
                    div(style="max-height:220px;overflow:auto;font-size:.85rem;line-height:1.5;border:1px solid #eee;border-radius:6px;padding:.5rem;background:#fafafa;",
                        d$core_symbols[sel] %||% "not available"))
      if (input$cat == "KEGG") {
        kegg_png <- file.path("pathview_noxy", cname(), paste0(sid, ".pathview.png"))
        if (file.exists(file.path("www", kegg_png)))
          ui <- tagList(ui, h6("KEGG pathway map"),
                        tags$img(src = kegg_png, style="max-width:100%;border:1px solid #eee;border-radius:6px"))
      }
      ui
    })
  })
}
