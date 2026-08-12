# Differential Expression tab (BATCH-CORRECTED, SEX CHROMOSOMES REMOVED) ------
# Same as tab_de_bc but reads the *_noxy artifacts (chrX/chrY features removed
# from the count matrix; replicate modelled as a batch covariate). Hover/tables
# include gene description + chromosome.
tab_de_noxy_ui <- function(id) {
  ns <- NS(id)
  nav_panel(
    title = "Differential Expression (corrected)", icon = icon("chart-column"),
    div(class = "container-fluid py-3",
      h3("Differential Expression — batch-corrected, sex chromosomes removed"),
      div(class = "banner",
          "chrX / chrY genes removed from the count matrix; SNI excluded; replicate modelled as an ",
          "additive batch covariate (", tags$code("~ 0 + group + replicate"), "). Hover any point for ",
          "gene name, description, chromosome, fold change and adjusted p-value."),
      uiOutput(ns("sni_note")),
      layout_sidebar(
        sidebar = sidebar(width = 320,
          radioButtons(ns("level"), "Level", c("Gene"="gene"), inline = TRUE),  # gene-only (transcript skipped)
          selectInput(ns("contrast"), "Contrast", choices = unname(CONTRAST_LABELS)),
          sliderInput(ns("padj"), "Adj. p-value cutoff", min = 0.01, max = 0.25, value = 0.05, step = 0.01),
          sliderInput(ns("lfc"), "|log2 fold-change| cutoff", min = 0, max = 3, value = 1, step = 0.1),
          uiOutput(ns("concordance"))
        ),
        navset_tab(
          nav_panel("Volcano",
            div(class = "row",
              div(class = "col-lg-3",
                sliderInput(ns("nlab"), "Label top-N genes", min = 0, max = 60, value = 20, step = 5),
                numericInput(ns("xmax"), "x-axis |log2FC| max", value = 6, min = 0.5, step = 0.5),
                numericInput(ns("ymax"), "y-axis −log10(adj.p) max", value = 30, min = 1, step = 1),
                radioButtons(ns("vfmt"), "Save as", c("PNG"="png","SVG"="svg","PDF"="pdf"), inline = TRUE),
                downloadButton(ns("dl_volc"), "Download volcano", class = "btn-sm btn-outline-secondary")),
              div(class = "col-lg-9", uiOutput(ns("volcano_wrap"))))),
          nav_panel("Results table", uiOutput(ns("table_wrap"))),
          nav_panel("Gene / transcript lookup",
            div(class = "row",
              div(class = "col-md-4",
                selectizeInput(ns("feat"), "Search gene symbol or feature ID",
                               choices = NULL, options = list(maxOptions = 50,
                                 placeholder = "type a symbol or ID…")),
                uiOutput(ns("feat_stats"))),
              div(class = "col-md-8", uiOutput(ns("box_wrap"))))
          )
        )
      )
    )
  )
}

tab_de_noxy_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    de_long <- load_rds("de_long_noxy")
    expr    <- load_rds("expr_noxy")
    idx     <- load_rds("search_index_noxy")

    cname <- reactive(names(CONTRAST_LABELS)[match(input$contrast, CONTRAST_LABELS)])
    output$sni_note <- renderUI(sni_note(cname()))
    sub <- function(engine) {
      d <- de_long
      d[d$level == input$level & d$contrast == cname() & d$engine == engine, ]
    }
    sig <- function(d) which(!is.na(d$padj) & d$padj < input$padj & abs(d$log2FC) >= input$lfc)

    vdat <- reactive({ req(de_long); attach_annot(sub("DESeq2"), input$level) })
    observeEvent(list(input$level, input$contrast), {
      d <- vdat(); d <- d[!is.na(d$padj), ]; req(nrow(d))
      updateNumericInput(session, "xmax", value = ceiling(max(abs(d$log2FC), na.rm = TRUE)))
      updateNumericInput(session, "ymax", value = max(5, ceiling(max(-log10(d$padj), na.rm = TRUE))))
    })
    output$volcano_wrap <- renderUI({
      if (is.null(de_long)) return(coming_soon("Appears once the sex-chr-removed pipeline finishes."))
      plotOutput(ns("volcano"), height = "600px")
    })
    volc_plot <- reactive({
      volcano_gg_app(vdat(), input$padj, input$lfc, input$nlab, input$xmax, input$ymax, disp(input$contrast))
    })
    output$volcano <- renderPlot({ req(de_long); print(volc_plot()) })
    output$dl_volc <- downloadHandler(
      filename = function() sprintf("volcano_corrected_%s_%s.%s", disp(cname()), input$level, input$vfmt),
      content  = function(file) { gg_device(file, input$vfmt, 8, 6.5); print(volc_plot()); grDevices::dev.off() })

    output$concordance <- renderUI({
      if (is.null(de_long)) return(NULL)
      dd <- sub("DESeq2"); de <- sub("edgeR")
      if (!nrow(dd) || !nrow(de)) return(NULL)
      s1 <- dd$feature_id[sig(dd)]; s2 <- de$feature_id[sig(de)]
      div(class = "mt-3", tags$b("Significant features (current cutoff)"),
        tags$ul(tags$li(sprintf("DESeq2: %d", length(s1))), tags$li(sprintf("edgeR: %d", length(s2))),
          tags$li(sprintf("Shared: %d", length(intersect(s1, s2)))),
          tags$li(sprintf("DESeq2-only: %d · edgeR-only: %d", length(setdiff(s1,s2)), length(setdiff(s2,s1))))))
    })

    output$table_wrap <- renderUI({
      if (is.null(de_long)) return(coming_soon("Results table pending."))
      DT::DTOutput(ns("table"))
    })
    output$table <- DT::renderDT({
      dd <- sub("DESeq2"); de <- sub("edgeR")
      m <- merge(dd[, c("feature_id","symbol","log2FC","padj")], de[, c("feature_id","log2FC","padj")],
                 by = "feature_id", suffixes = c("_DESeq2","_edgeR"))
      m <- attach_annot(m, input$level); m <- m[order(m$padj_DESeq2), ]
      m[, c("symbol","chr","description","feature_id","log2FC_DESeq2","padj_DESeq2","log2FC_edgeR","padj_edgeR")] |>
        DT::datatable(rownames = FALSE, filter = "top", extensions = "Buttons",
          colnames = c("Symbol"="symbol","Chromosome"="chr","Description"="description",
                       "Feature ID"="feature_id","log2FC (DESeq2)"="log2FC_DESeq2","adj.p (DESeq2)"="padj_DESeq2",
                       "log2FC (edgeR)"="log2FC_edgeR","adj.p (edgeR)"="padj_edgeR"),
          options = list(pageLength = 15, scrollX = TRUE, dom = "Bfrtip",
            buttons = list(list(extend = "excel", text = "⤓ Excel",
                                filename = paste0("DE_corrected_", disp(cname()), "_", input$level))))) |>
        DT::formatSignif(5:8, 3)   # the 4 numeric columns (by position — headers are renamed)
    })

    observe({ req(idx); ch <- idx[idx$level == input$level, ]
      updateSelectizeInput(session, "feat", choices = setNames(ch$feature_id, ch$label), server = TRUE) })
    output$box_wrap <- renderUI({
      if (is.null(expr)) return(coming_soon("Expression boxplots pending."))
      plotly::plotlyOutput(ns("box"), height = "460px")
    })
    output$box <- plotly::renderPlotly({
      req(input$feat, expr); E <- expr[[input$level]]; req(input$feat %in% rownames(E$vst))
      df <- data.frame(expr = E$vst[input$feat, ], E$meta)
      df$.x <- factor(disp(df$group), levels = disp(sort(unique(as.character(df$group)))))
      sc <- color_spec(df$sex, "sex"); df$.sexf <- factor(df$sex, levels = sc$levels)
      plotly::plot_ly(df, x = ~.x, y = ~expr, color = ~.sexf, colors = sc$colors, type = "box",
                      boxpoints = "all", jitter = 0.5, pointpos = 0,
                      text = ~paste0(disp(sample_id), " (", disp(condition), ", ", sex, ")"), hoverinfo = "text+y") |>
        plotly::layout(title = paste0("Expression (batch-corrected) — ", input$feat),
                       xaxis = list(title = ""), yaxis = list(title = "VST expression (replicate removed)"),
                       legend = list(title = list(text = "Sex"))) |>
        ply_pub()
    })
    output$feat_stats <- renderUI({
      req(input$feat, de_long)
      rows <- de_long[de_long$level==input$level & de_long$contrast==cname() & de_long$feature_id==input$feat, ]
      if (!nrow(rows)) return(div(class="text-muted","No DE stats for this feature/contrast."))
      sym <- rows$symbol[1]; an <- annot_of(input$feat, input$level)
      div(class = "banner mt-2",
        tags$b(if (is.na(sym)) input$feat else paste0(sym, " (", input$feat, ")")), tags$br(),
        tags$small(if (!is.na(an$chr)) paste0("chromosome ", an$chr) else ""),
        if (!is.na(an$description)) div(tags$small(em(an$description))),
        tags$small(em(input$contrast)), tags$hr(),
        lapply(split(rows, rows$engine), function(r)
          div(tags$b(r$engine[1]), ": log2FC ", round(r$log2FC[1], 2), " · adj.p ", signif(r$padj[1], 3))))
    })
  })
}
