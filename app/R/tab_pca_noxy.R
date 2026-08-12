# PCA tab (BATCH-CORRECTED, SEX CHROMOSOMES REMOVED) -------------------------
# Same as tab_pca_bc but reads pca_noxy.rds (chrX/chrY features dropped from the
# count matrix, replicate batch removed). Diagnostic: does NP-vs-Sham structure
# emerge once both sex and batch are out of the way?
tab_pca_noxy_ui <- function(id) {
  ns <- NS(id)
  nav_panel(
    title = "PCA (corrected)", icon = icon("braille"),
    div(class = "container-fluid py-3",
      h3("PCA — batch-corrected, sex chromosomes removed"),
      div(class = "banner",
          "SNI excluded, replicate batch removed (", tags$code("limma::removeBatchEffect"),
          "), and all ", tags$b("chrX / chrY genes dropped from the count matrix"),
          " before PCA. This strips the two dominant nuisances (batch + sex) so any ",
          "NP-vs-Sham structure can show. (Global separation is not required for valid DE.)"),
      layout_sidebar(
        sidebar = sidebar(width = 300,
          radioButtons(ns("level"), "Quantification level",
                       c("Gene" = "gene"), inline = TRUE),  # gene-only (transcript skipped)
          selectInput(ns("pcx"), "X axis", paste0("PC", 1:10), "PC1"),
          selectInput(ns("pcy"), "Y axis", paste0("PC", 1:10), "PC2"),
          selectInput(ns("color"), "Colour by",
                      c("Sex" = "sex", "Condition" = "condition", "Group" = "group"))
        ),
        navset_tab(
          nav_panel("PCA scatter", uiOutput(ns("scatter_wrap")), uiOutput(ns("pca_dl"))),
          nav_panel("Elbow / scree", uiOutput(ns("scree_wrap")))
        )
      )
    )
  )
}

tab_pca_noxy_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    pca <- load_rds("pca_noxy")
    dat <- reactive({ req(pca); pca[[input$level]] })

    output$scatter_wrap <- renderUI({
      if (is.null(pca)) return(coming_soon("This PCA appears once the sex-chr-removed pipeline finishes."))
      plotly::plotlyOutput(session$ns("scatter"), height = "560px")
    })
    output$scree_wrap <- renderUI({
      if (is.null(pca)) return(coming_soon("Scree plot pending."))
      plotly::plotlyOutput(session$ns("scree"), height = "460px")
    })

    output$pca_dl <- renderUI({
      req(pca)
      fig_dl(sprintf("PCA_corrected_%s_by_%s", input$level, input$color), "Download this PCA")
    })

    output$scatter <- plotly::renderPlotly({
      d <- dat(); req(d); co <- d$coords; ve <- d$var_explained
      vx <- ve$pct[match(input$pcx, ve$PC)]; vy <- ve$pct[match(input$pcy, ve$PC)]
      cs  <- color_spec(co[[input$color]], input$color)
      co$.colf <- factor(disp(co[[input$color]]), levels = disp(cs$levels))
      shp <- shape_spec(co$condition)
      co$.shpf <- factor(disp(co$condition), levels = disp(shp$levels))
      plotly::plot_ly(co, x = ~get(input$pcx), y = ~get(input$pcy),
        color = ~.colf, colors = cs$colors, symbol = ~.shpf, symbols = shp$symbols,
        type = "scatter", mode = "markers",
        text = ~paste0("<b>", disp(sample_id), "</b><br>group: ", disp(group),
                       "<br>condition: ", disp(condition), "<br>sex: ", sex,
                       "<br>replicate: ", replicate),
        hoverinfo = "text", marker = list(size = 13, line = list(width = 1, color = "#fff"))) |>
        plotly::layout(
          title = sprintf("%s level — batch-corrected, sex chromosomes removed",
                          tools::toTitleCase(input$level)),
          xaxis = list(title = sprintf("%s (%.1f%%)", input$pcx, vx %||% NA)),
          yaxis = list(title = sprintf("%s (%.1f%%)", input$pcy, vy %||% NA)),
          legend = list(title = list(text = tools::toTitleCase(input$color)))) |>
        ply_pub()
    })

    output$scree <- plotly::renderPlotly({
      d <- dat(); req(d); ve <- d$var_explained
      plotly::plot_ly(ve, x = ~factor(PC, levels = PC)) |>
        plotly::add_bars(y = ~pct, name = "% variance", text = ~paste0(pct, "%"), hoverinfo = "text") |>
        plotly::add_lines(y = ~cumpct, name = "cumulative %", yaxis = "y2", line = list(color = "#e74c3c")) |>
        plotly::layout(title = "Variance explained per PC (elbow)",
          xaxis = list(title = "Principal component"), yaxis = list(title = "% variance"),
          yaxis2 = list(title = "cumulative %", overlaying = "y", side = "right", range = c(0, 100)),
          legend = list(orientation = "h"))
    })
  })
}
