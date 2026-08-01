# Venn / overlap tab — pick 2–4 comparisons at a DEG cutoff, see the overlap,
# the per-region gene tables, and save. Corrected analysis, DESeq2. (>4 -> UpSet.)
tab_venn_ui <- function(id) {
  ns <- NS(id)
  nav_panel(
    title = "Venn / overlap", icon = icon("circle-nodes"),
    div(class = "container-fluid py-3",
      h3("DEG overlap across comparisons"),
      div(class = "banner",
          "Corrected analysis (batch + sex chromosomes), DESeq2. Pick 2–4 comparisons and a DEG ",
          "cutoff; the Venn shows how their gene sets overlap. The table lists the genes in each ",
          "region (downloadable). Selecting 5 shows an UpSet plot instead."),
      uiOutput(ns("sni_note")),
      layout_sidebar(
        sidebar = sidebar(width = 340,
          radioButtons(ns("level"), "Level", c("Gene"="gene"), inline = TRUE),  # gene-only (transcript skipped)
          checkboxGroupInput(ns("contrasts"), "Comparisons (pick 2–4)",
                             choices = unname(CONTRAST_LABELS),
                             selected = unname(CONTRAST_LABELS[c("NP_vs_Sham","NP_M_vs_Sham_M","NP_F_vs_Sham_F")])),
          sliderInput(ns("padj"), "Adj. p cutoff", min = 0.01, max = 0.25, value = 0.05, step = 0.01),
          sliderInput(ns("lfc"), "|log2FC| cutoff", min = 0, max = 3, value = 1, step = 0.1),
          radioButtons(ns("dir"), "Direction", c("Both"="both","Up only"="up","Down only"="down"), inline = TRUE),
          radioButtons(ns("vfmt"), "Save as", c("PNG"="png","SVG"="svg","PDF"="pdf"), inline = TRUE),
          downloadButton(ns("dl_venn"), "Download Venn", class = "btn-sm btn-outline-secondary")
        ),
        navset_tab(
          nav_panel("Diagram", uiOutput(ns("venn_wrap"))),
          nav_panel("Genes per region", uiOutput(ns("tbl_wrap"))),
          nav_panel("log2FC correlation",
            div(class = "row g-2 mt-1 mb-1",
              div(class = "col-md-5",
                  selectInput(ns("cx"), "Contrast X (horizontal)", choices = unname(CONTRAST_LABELS),
                              selected = unname(CONTRAST_LABELS["NP_M_vs_Sham_M"]))),
              div(class = "col-md-5",
                  selectInput(ns("cy"), "Contrast Y (vertical)", choices = unname(CONTRAST_LABELS),
                              selected = unname(CONTRAST_LABELS["NP_F_vs_Sham_F"]))),
              div(class = "col-md-2",
                  radioButtons(ns("cset"), "Genes",
                               c("DEG in either" = "deg", "All tested" = "all")))),
            tags$small(class = "text-muted",
                       "\"DEG in either\" uses the adj.p / |log2FC| cutoffs in the sidebar. Hover a point to see the gene."),
            uiOutput(ns("corr_stats")),
            plotly::plotlyOutput(ns("corr"), height = "520px"),
            tags$hr(), DT::DTOutput(ns("corr_tbl")))
        )
      )
    )
  )
}

tab_venn_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    de_long <- load_rds("de_long_noxy")

    # warn whenever any selected comparison (overlap checkboxes or correlation X/Y) is SNI
    output$sni_note <- renderUI({
      keys <- names(CONTRAST_LABELS)[match(c(input$contrasts, input$cx, input$cy), CONTRAST_LABELS)]
      if (any(is_sni(keys), na.rm = TRUE)) sni_note("SNI") else NULL
    })

    # named list of DEG feature-id sets, one per selected comparison
    sets <- reactive({
      req(de_long, length(input$contrasts) >= 2)
      cn <- names(CONTRAST_LABELS)[match(input$contrasts, CONTRAST_LABELS)]
      s <- lapply(cn, function(c1) {
        d <- de_long[de_long$engine=="DESeq2" & de_long$level==input$level & de_long$contrast==c1, ]
        keep <- !is.na(d$padj) & d$padj < input$padj & abs(d$log2FC) >= input$lfc
        if (input$dir == "up")   keep <- keep & d$log2FC > 0
        if (input$dir == "down") keep <- keep & d$log2FC < 0
        d$feature_id[keep]
      })
      setNames(s, disp(cn))
    })

    venn_plot <- reactive({
      s <- sets(); req(length(s) >= 2)
      if (length(s) <= 4) {
        ggVennDiagram::ggVennDiagram(s, label = "count", label_alpha = 0, set_size = 4.5) +
          ggplot2::scale_fill_gradient(low = "#f7fbff", high = "#E67E22", name = "genes") +
          ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = 0.25)) +
          ggplot2::labs(title = sprintf("DEG overlap (adj.p<%.3g, |log2FC|>=%g)", input$padj, input$lfc))
      } else {
        function() print(UpSetR::upset(UpSetR::fromList(s), nsets = length(s), order.by = "freq",
                                       main.bar.color = "#E67E22", sets.bar.color = "#021893"))
      }
    })

    output$venn_wrap <- renderUI({
      if (is.null(de_long)) return(coming_soon("Corrected DE results pending."))
      if (length(input$contrasts) < 2) return(div(class="banner","Pick at least 2 comparisons."))
      plotOutput(ns("venn"), height = "560px")
    })
    output$venn <- renderPlot({
      p <- venn_plot(); if (is.function(p)) p() else print(p)
    })
    output$dl_venn <- downloadHandler(
      filename = function() sprintf("venn_%s_%d-way.%s", input$level, length(sets()), input$vfmt),
      content  = function(file) {
        gg_device(file, input$vfmt, 8, 7); p <- venn_plot(); if (is.function(p)) p() else print(p); grDevices::dev.off()
      })

    # ---- genes per region (membership matrix) ----
    membership <- reactive({
      s <- sets(); req(length(s) >= 2)
      u <- sort(unique(unlist(s))); m <- data.frame(feature_id = u, stringsAsFactors = FALSE)
      for (nm in names(s)) m[[nm]] <- m$feature_id %in% s[[nm]]
      m$region <- apply(m[names(s)], 1, function(r) paste(names(s)[r], collapse = " & "))
      m <- attach_annot(m, input$level)
      a <- ANNOT[ANNOT$level == input$level, ]
      m$symbol <- a$symbol[match(m$feature_id, a$feature_id)]
      m[, c("region","symbol","chr","description","feature_id", names(s))]
    })
    output$tbl_wrap <- renderUI({
      if (is.null(de_long) || length(input$contrasts) < 2) return(coming_soon("Pick at least 2 comparisons."))
      DT::DTOutput(ns("tbl"))
    })
    output$tbl <- DT::renderDT({
      m <- membership(); req(nrow(m))
      DT::datatable(m, rownames = FALSE, filter = "top", extensions = "Buttons",
        colnames = c("Region"="region","Symbol"="symbol","Chromosome"="chr","Description"="description","Feature ID"="feature_id"),
        options = list(pageLength = 15, scrollX = TRUE, dom = "Bfrtip",
          buttons = list(list(extend = "excel", text = "⤓ Excel", filename = paste0("venn_genes_", input$level)))))
    })

    # ---- log2FC correlation between two chosen contrasts ----
    corr_data <- reactive({
      req(de_long, input$level, input$cx, input$cy, input$cset)
      cx <- names(CONTRAST_LABELS)[match(input$cx, CONTRAST_LABELS)]
      cy <- names(CONTRAST_LABELS)[match(input$cy, CONTRAST_LABELS)]
      req(length(cx) == 1L, length(cy) == 1L, !is.na(cx), !is.na(cy))
      base <- de_long[de_long$engine == "DESeq2" & de_long$level == input$level, ]
      dx <- base[base$contrast == cx, c("feature_id","symbol","log2FC","padj")]
      dy <- base[base$contrast == cy, c("feature_id","log2FC","padj")]
      m <- merge(dx, dy, by = "feature_id", suffixes = c("_x","_y"))
      m <- m[!is.na(m$log2FC_x) & !is.na(m$log2FC_y), ]
      m$symbol <- ifelse(is.na(m$symbol) | m$symbol == "", as.character(m$feature_id), as.character(m$symbol))  # factor-safe
      sigx <- !is.na(m$padj_x) & m$padj_x < input$padj & abs(m$log2FC_x) >= input$lfc
      sigy <- !is.na(m$padj_y) & m$padj_y < input$padj & abs(m$log2FC_y) >= input$lfc
      m$class <- ifelse(sigx & sigy, "DEG in both", ifelse(sigx | sigy, "DEG in one", "n.s."))
      if (input$cset == "deg") m <- m[sigx | sigy, ]
      m
    })

    output$corr_stats <- renderUI({
      if (is.null(de_long)) return(coming_soon("Corrected DE results pending."))
      m <- corr_data(); req(nrow(m) >= 3)
      pe <- suppressWarnings(cor.test(m$log2FC_x, m$log2FC_y, method = "pearson"))
      sp <- suppressWarnings(cor.test(m$log2FC_x, m$log2FC_y, method = "spearman"))
      div(class = "banner",
          sprintf("n = %d genes · Pearson r = %.3f (p = %.2g) · Spearman ρ = %.3f (p = %.2g)",
                  nrow(m), unname(pe$estimate), pe$p.value, unname(sp$estimate), sp$p.value))
    })

    output$corr <- plotly::renderPlotly({
      m <- corr_data(); req(nrow(m) >= 3)
      cols <- c("DEG in both" = "#E67E22", "DEG in one" = "#3477eb", "n.s." = "#B0B8BF")
      m$class <- factor(m$class, levels = names(cols))
      m$hover <- paste0("<b>", m$symbol, "</b><br>", input$cx, ": ", round(m$log2FC_x, 2),
                        "<br>", input$cy, ": ", round(m$log2FC_y, 2))
      rng <- range(m$log2FC_x, na.rm = TRUE)
      fit <- stats::lm(log2FC_y ~ log2FC_x, m)
      fy  <- as.numeric(predict(fit, data.frame(log2FC_x = rng)))
      plotly::plot_ly(m, x = ~log2FC_x, y = ~log2FC_y, type = "scatter", mode = "markers",
        color = ~class, colors = cols, text = ~hover, hoverinfo = "text",
        marker = list(size = 6, line = list(width = 0.3, color = "#fff"))) |>
        plotly::layout(title = sprintf("log2FC correlation — %s vs %s", input$cx, input$cy),
          xaxis = list(title = paste0("log2FC · ", input$cx)),
          yaxis = list(title = paste0("log2FC · ", input$cy)),
          shapes = list(
            list(type = "line", x0 = rng[1], x1 = rng[2], y0 = rng[1], y1 = rng[2],
                 line = list(color = "grey70", dash = "dot")),
            list(type = "line", x0 = rng[1], x1 = rng[2], y0 = fy[1], y1 = fy[2],
                 line = list(color = "black", dash = "dash")))) |>
        ply_pub()
    })

    output$corr_tbl <- DT::renderDT({
      m <- corr_data(); req(nrow(m))
      out <- data.frame(Symbol = m$symbol,
                        `log2FC (X)` = round(m$log2FC_x, 3), `log2FC (Y)` = round(m$log2FC_y, 3),
                        `adj.p (X)` = signif(m$padj_x, 3), `adj.p (Y)` = signif(m$padj_y, 3),
                        Class = m$class, check.names = FALSE, stringsAsFactors = FALSE)
      out <- out[order(-(abs(out[["log2FC (X)"]]) + abs(out[["log2FC (Y)"]]))), ]
      DT::datatable(out, rownames = FALSE, filter = "top", extensions = "Buttons",
        options = list(pageLength = 12, scrollX = TRUE, dom = "Bfrtip",
          buttons = list(list(extend = "excel", text = "⤓ Excel", filename = "log2FC_correlation"))))
    })
  })
}
