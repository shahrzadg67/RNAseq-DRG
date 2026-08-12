# PCA tab (BATCH-CORRECTED) — gene/transcript, NP/Sham only, replicate removed
# Parallel to tab_pca.R but reads pca_bc.rds. SNI is excluded entirely, so the
# SNI toggle is dropped and the variant key is just the level ("gene"/"transcript").
tab_pca_bc_ui <- function(id) {
  ns <- NS(id)
  nav_panel(
    title = "PCA (batch-corrected)", icon = icon("braille"),
    div(class = "container-fluid py-3",
      h3("Principal Component Analysis — batch-corrected"),
      div(class = "banner",
          "SNI samples excluded. The replicate batch effect has been removed with ",
          tags$code("limma::removeBatchEffect"), " (preserving condition × sex) before PCA, ",
          "so the components now reflect biology rather than processing round. Compare with the ",
          "uncorrected ", tags$b("PCA"), " tab."),
      layout_sidebar(
        sidebar = sidebar(width = 300,
          radioButtons(ns("level"), "Quantification level",
                       c("Gene" = "gene", "Transcript" = "transcript"), inline = TRUE),
          selectInput(ns("pcx"), "X axis", paste0("PC", 1:10), "PC1"),
          selectInput(ns("pcy"), "Y axis", paste0("PC", 1:10), "PC2"),
          selectInput(ns("color"), "Colour by",
                      c("Group" = "group", "Condition" = "condition", "Sex" = "sex"))
        ),
        navset_tab(
          nav_panel("PCA scatter", uiOutput(ns("scatter_wrap"))),
          nav_panel("Elbow / scree", uiOutput(ns("scree_wrap")))
        )
      )
    )
  )
}

tab_pca_bc_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    pca <- load_rds("pca_bc")
    dat <- reactive({ req(pca); pca[[input$level]] })

    output$scatter_wrap <- renderUI({
      if (is.null(pca)) return(coming_soon("Batch-corrected PCA will appear once the analysis has run."))
      plotly::plotlyOutput(session$ns("scatter"), height = "560px")
    })
    output$scree_wrap <- renderUI({
      if (is.null(pca)) return(coming_soon("Scree plot will appear once the analysis has run."))
      plotly::plotlyOutput(session$ns("scree"), height = "460px")
    })

    output$scatter <- plotly::renderPlotly({
      d <- dat(); req(d)
      co <- d$coords; ve <- d$var_explained
      vx <- ve$pct[match(input$pcx, ve$PC)]; vy <- ve$pct[match(input$pcy, ve$PC)]
      co$.col <- co[[input$color]]
      plotly::plot_ly(
        co, x = ~get(input$pcx), y = ~get(input$pcy), color = ~.col,
        type = "scatter", mode = "markers",
        text = ~paste0("<b>", sample_id, "</b><br>group: ", group,
                       "<br>condition: ", condition, "<br>sex: ", sex,
                       "<br>replicate: ", replicate),
        hoverinfo = "text", marker = list(size = 12, line = list(width = 1, color = "#fff"))
      ) |>
        plotly::layout(
          title = sprintf("%s level — batch-corrected (replicate removed)",
                          tools::toTitleCase(input$level)),
          xaxis = list(title = sprintf("%s (%.1f%%)", input$pcx, vx %||% NA)),
          yaxis = list(title = sprintf("%s (%.1f%%)", input$pcy, vy %||% NA)),
          legend = list(title = list(text = input$color)))
    })

    output$scree <- plotly::renderPlotly({
      d <- dat(); req(d); ve <- d$var_explained
      plotly::plot_ly(ve, x = ~factor(PC, levels = PC)) |>
        plotly::add_bars(y = ~pct, name = "% variance",
                         text = ~paste0(pct, "%"), hoverinfo = "text") |>
        plotly::add_lines(y = ~cumpct, name = "cumulative %", yaxis = "y2",
                          line = list(color = "#e74c3c")) |>
        plotly::layout(
          title = "Variance explained per PC (elbow)",
          xaxis = list(title = "Principal component"),
          yaxis = list(title = "% variance"),
          yaxis2 = list(title = "cumulative %", overlaying = "y", side = "right",
                        range = c(0, 100)),
          legend = list(orientation = "h"))
    })
  })
}
