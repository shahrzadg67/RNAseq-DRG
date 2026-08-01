# PCA tab — gene/transcript, with/without SNI, elbow, free PC selector --------
tab_pca_ui <- function(id) {
  ns <- NS(id)
  nav_panel(
    title = "PCA (uncorrected)", icon = icon("braille"),
    div(class = "container-fluid py-3",
      h3("Principal Component Analysis"),
      div(class = "banner",
          "Sample structure on the top variable features. Use the elbow (scree) plot to ",
          "judge how many components are informative, then pick any two PCs to view. ",
          "Toggle SNI in/out, and switch between gene and transcript level."),
      layout_sidebar(
        sidebar = sidebar(width = 300,
          radioButtons(ns("level"), "Quantification level",
                       c("Gene" = "gene", "Transcript" = "transcript"), inline = TRUE),
          radioButtons(ns("sni"), "SNI samples",
                       c("Include SNI" = "withSNI", "Exclude SNI" = "noSNI")),
          radioButtons(ns("xy"), "Sex chromosomes",
                       c("Include" = "", "Remove chrX/Y" = "_noXY")),
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

tab_pca_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    pca <- load_rds("pca")
    key <- reactive(sprintf("%s_%s%s", input$level, input$sni, input$xy))
    dat <- reactive({ req(pca); pca[[key()]] })

    output$scatter_wrap <- renderUI({
      if (is.null(pca)) return(coming_soon("PCA will appear once the analysis has run."))
      plotly::plotlyOutput(session$ns("scatter"), height = "560px")
    })
    output$scree_wrap <- renderUI({
      if (is.null(pca)) return(coming_soon("Scree plot will appear once the analysis has run."))
      plotly::plotlyOutput(session$ns("scree"), height = "460px")
    })

    output$pca_dl <- renderUI({
      req(pca)
      fig_dl(sprintf("PCA_uncorrected_%s_%s%s_by_%s", input$level, input$sni, input$xy, input$color),
             "Download this PCA")
    })

    output$scatter <- plotly::renderPlotly({
      d <- dat(); req(d)
      co <- d$coords; ve <- d$var_explained
      vx <- ve$pct[match(input$pcx, ve$PC)]; vy <- ve$pct[match(input$pcy, ve$PC)]
      cs  <- color_spec(co[[input$color]], input$color)
      co$.colf <- factor(disp(co[[input$color]]), levels = disp(cs$levels))
      shp <- shape_spec(co$condition)
      co$.shpf <- factor(disp(co$condition), levels = disp(shp$levels))
      plotly::plot_ly(
        co, x = ~get(input$pcx), y = ~get(input$pcy),
        color = ~.colf, colors = cs$colors, symbol = ~.shpf, symbols = shp$symbols,
        type = "scatter", mode = "markers",
        text = ~paste0("<b>", disp(sample_id), "</b><br>group: ", disp(group),
                       "<br>condition: ", disp(condition), "<br>sex: ", sex),
        hoverinfo = "text", marker = list(size = 13, line = list(width = 1, color = "#fff"))
      ) |>
        plotly::layout(
          title = sprintf("%s level — %s%s", tools::toTitleCase(input$level),
                          if (input$sni == "withSNI") "with SNI" else "without SNI",
                          if (identical(input$xy, "_noXY")) ", chrX/Y removed" else ""),
          xaxis = list(title = sprintf("%s (%.1f%%)", input$pcx, vx %||% NA)),
          yaxis = list(title = sprintf("%s (%.1f%%)", input$pcy, vy %||% NA)),
          legend = list(title = list(text = tools::toTitleCase(input$color)))) |>
        ply_pub()
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
