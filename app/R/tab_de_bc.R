# Differential Expression tab (BATCH-CORRECTED) -----------------------------
# Parallel to tab_de.R but reads the *_bc artifacts. DESeq2 and edgeR were fit
# with replicate as an additive batch covariate (~ 0 + group + replicate);
# boxplots show batch-corrected (replicate-removed) VST expression.
tab_de_bc_ui <- function(id) {
  ns <- NS(id)
  nav_panel(
    title = "Differential Expression (batch-corrected)", icon = icon("chart-column"),
    div(class = "container-fluid py-3",
      h3("Differential Expression — batch-corrected"),
      div(class = "banner",
          "SNI excluded. Both engines modelled ", tags$b("replicate"), " as an additive batch ",
          "covariate (", tags$code("~ 0 + group + replicate"), "), so each contrast is adjusted ",
          "for processing round. Boxplots show replicate-removed expression. Hover any point for ",
          "the gene/transcript name, fold change and adjusted p-value."),
      layout_sidebar(
        sidebar = sidebar(width = 320,
          radioButtons(ns("level"), "Level", c("Gene"="gene","Transcript"="transcript"), inline = TRUE),
          selectInput(ns("contrast"), "Contrast", choices = CONTRAST_LABELS),
          sliderInput(ns("padj"), "Adj. p-value cutoff", min = 0.01, max = 0.25, value = 0.05, step = 0.01),
          sliderInput(ns("lfc"), "|log2 fold-change| cutoff", min = 0, max = 3, value = 1, step = 0.1),
          uiOutput(ns("concordance"))
        ),
        navset_tab(
          nav_panel("Volcano (DESeq2 vs edgeR)", uiOutput(ns("volcano_wrap"))),
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

tab_de_bc_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    de_long <- load_rds("de_long_bc")
    expr    <- load_rds("expr_bc")
    idx     <- load_rds("search_index_bc")

    cname <- reactive(names(CONTRAST_LABELS)[match(input$contrast, CONTRAST_LABELS)])

    sub <- function(engine) {
      d <- de_long
      d[d$level == input$level & d$contrast == cname() & d$engine == engine, ]
    }
    sig <- function(d) which(!is.na(d$padj) & d$padj < input$padj & abs(d$log2FC) >= input$lfc)

    # ---- volcano (side-by-side) ----
    output$volcano_wrap <- renderUI({
      if (is.null(de_long)) return(coming_soon("Batch-corrected differential expression is being computed."))
      plotly::plotlyOutput(ns("volcano"), height = "600px")
    })
    one_volcano <- function(engine) {
      d <- sub(engine); req(nrow(d) > 0)
      d <- attach_annot(d, input$level)
      d$sigcat <- "NS"
      d$sigcat[sig(d)] <- "sig"
      d$lab <- ifelse(is.na(d$symbol), d$feature_id, d$symbol)
      plotly::plot_ly(d, x = ~log2FC, y = ~ -log10(padj),
        color = ~sigcat, colors = c(NS = "#b0b8bf", sig = "#e74c3c"),
        type = "scatter", mode = "markers",
        marker = list(size = 5, opacity = 0.6),
        text = ~paste0("<b>", lab, "</b>  (chr ", ifelse(is.na(chr), "?", chr), ")<br>",
                       feature_id,
                       ifelse(is.na(description), "", paste0("<br>", description)),
                       "<br>log2FC: ", round(log2FC, 2),
                       "<br>adj.p: ", signif(padj, 3)),
        hoverinfo = "text") |>
        plotly::layout(title = engine,
          xaxis = list(title = "log2 fold-change"),
          yaxis = list(title = "-log10 adj.p"),
          shapes = list(
            list(type="line", x0=input$lfc, x1=input$lfc, y0=0, y1=1, yref="paper",
                 line=list(dash="dot", color="#888")),
            list(type="line", x0=-input$lfc, x1=-input$lfc, y0=0, y1=1, yref="paper",
                 line=list(dash="dot", color="#888"))))
    }
    output$volcano <- plotly::renderPlotly({
      plotly::subplot(one_volcano("DESeq2"), one_volcano("edgeR"),
                      nrows = 1, shareY = TRUE, titleX = TRUE) |>
        plotly::layout(showlegend = FALSE,
          annotations = list(
            list(x=0.2, y=1.04, text="<b>DESeq2</b>", showarrow=FALSE, xref="paper", yref="paper"),
            list(x=0.8, y=1.04, text="<b>edgeR</b>",  showarrow=FALSE, xref="paper", yref="paper")))
    })

    # ---- concordance summary ----
    output$concordance <- renderUI({
      if (is.null(de_long)) return(NULL)
      dd <- sub("DESeq2"); de <- sub("edgeR")
      if (!nrow(dd) || !nrow(de)) return(NULL)
      s1 <- dd$feature_id[sig(dd)]; s2 <- de$feature_id[sig(de)]
      shared <- length(intersect(s1, s2))
      div(class = "mt-3",
        tags$b("Significant features (current cutoff)"),
        tags$ul(
          tags$li(sprintf("DESeq2: %d", length(s1))),
          tags$li(sprintf("edgeR: %d", length(s2))),
          tags$li(sprintf("Shared: %d", shared)),
          tags$li(sprintf("DESeq2-only: %d · edgeR-only: %d",
                          length(setdiff(s1, s2)), length(setdiff(s2, s1))))))
    })

    # ---- results table (intersection-friendly: DESeq2 stats + edgeR padj) ----
    output$table_wrap <- renderUI({
      if (is.null(de_long)) return(coming_soon("Results table pending."))
      DT::DTOutput(ns("table"))
    })
    output$table <- DT::renderDT({
      dd <- sub("DESeq2"); de <- sub("edgeR")
      m <- merge(dd[, c("feature_id","symbol","log2FC","padj")],
                 de[, c("feature_id","log2FC","padj")],
                 by = "feature_id", suffixes = c("_DESeq2","_edgeR"))
      m <- attach_annot(m, input$level)
      m <- m[order(m$padj_DESeq2), ]
      m[, c("symbol","chr","description","feature_id",
            "log2FC_DESeq2","padj_DESeq2","log2FC_edgeR","padj_edgeR")] |>
        DT::datatable(rownames = FALSE, filter = "top",
                      colnames = c(Symbol="symbol", Chr="chr", Description="description"),
                      options = list(pageLength = 15, scrollX = TRUE)) |>
        DT::formatSignif(c("log2FC_DESeq2","padj_DESeq2","log2FC_edgeR","padj_edgeR"), 3)
    })

    # ---- gene/transcript expression lookup ----
    observe({
      req(idx)
      ch <- idx[idx$level == input$level, ]
      updateSelectizeInput(session, "feat",
        choices = setNames(ch$feature_id, ch$label), server = TRUE)
    })
    output$box_wrap <- renderUI({
      if (is.null(expr)) return(coming_soon("Expression boxplots pending."))
      plotly::plotlyOutput(ns("box"), height = "460px")
    })
    output$box <- plotly::renderPlotly({
      req(input$feat, expr)
      E <- expr[[input$level]]; req(input$feat %in% rownames(E$vst))
      df <- data.frame(expr = E$vst[input$feat, ], E$meta)
      plotly::plot_ly(df, x = ~group, y = ~expr, color = ~group, type = "box",
                      boxpoints = "all", jitter = 0.5, pointpos = 0,
                      text = ~sample_id, hoverinfo = "text+y") |>
        plotly::layout(title = paste0("Expression (batch-corrected) — ", input$feat),
                       yaxis = list(title = "VST expression (replicate removed)"), showlegend = FALSE)
    })
    output$feat_stats <- renderUI({
      req(input$feat, de_long)
      rows <- de_long[de_long$level==input$level & de_long$contrast==cname() &
                      de_long$feature_id==input$feat, ]
      if (!nrow(rows)) return(div(class="text-muted","No DE stats for this feature/contrast."))
      sym <- rows$symbol[1]; an <- annot_of(input$feat, input$level)
      div(class = "banner mt-2",
        tags$b(if (is.na(sym)) input$feat else paste0(sym, " (", input$feat, ")")), tags$br(),
        tags$small(if (!is.na(an$chr)) paste0("chromosome ", an$chr) else ""),
        if (!is.na(an$description)) div(tags$small(em(an$description))),
        tags$small(em(input$contrast)), tags$hr(),
        lapply(split(rows, rows$engine), function(r)
          div(tags$b(r$engine[1]), ": log2FC ", round(r$log2FC[1], 2),
              " · adj.p ", signif(r$padj[1], 3))))
    })
  })
}
