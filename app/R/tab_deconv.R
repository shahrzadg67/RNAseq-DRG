# Deconvolution tab — cell-type proportions of every bulk sample estimated from a
# mouse DRG single-cell reference (MuSiC). Stacked composition bars + per-contrast
# comparison boxplots (all 5 study contrasts) + downloadable table/figures.
# Shows a "processing" placeholder until app/data/deconv.rds exists.
tab_deconv_ui <- function(id) {
  ns <- NS(id)
  nav_panel(
    title = "Deconvolution", icon = icon("chart-pie"),
    div(class = "container-fluid py-3",
      h3("Cell-type deconvolution"),
      div(class = "banner",
          "Estimated proportion of each DRG cell type in every bulk sample, inferred from a mouse ",
          "DRG single-cell reference with ", tags$b("MuSiC"), ". Use it to see whether MINP-vs-Sham ",
          "(or the sex) differences are partly driven by shifts in cell-type composition."),
      uiOutput(ns("prov")),
      uiOutput(ns("sni_note")),
      layout_sidebar(
        sidebar = sidebar(width = 340,
          radioButtons(ns("compmode"), "Composition — group by",
                       c("Sample" = "sample", "Group (sex×condition)" = "group", "Condition" = "condition"),
                       inline = TRUE),
          radioButtons(ns("compres"), "Layout",
                       c("All 18 types (one plot)" = "all",
                         "Split: neuronal + non-neuronal (two plots)" = "split")),
          selectInput(ns("contrast"), "Comparison (Compare tab)", choices = unname(CONTRAST_LABELS)),
          tags$small(class = "text-muted",
                     "Composition can be shown per sample or averaged per group; the comparison box-plots and stats use the selected contrast."),
          tags$hr(),
          radioButtons(ns("fmt"), "Save figures as", c("PNG"="png","SVG"="svg","PDF"="pdf"), inline = TRUE),
          downloadButton(ns("dl_comp"), "Download composition", class = "btn-sm btn-outline-secondary"),
          div(class = "mt-1"),
          downloadButton(ns("dl_cmp"), "Download comparison", class = "btn-sm btn-outline-secondary")
        ),
        navset_tab(
          nav_panel("Composition",
            uiOutput(ns("comp_wrap")),
            tags$hr(), h5("Proportions table (per sample)"), DT::DTOutput(ns("comp_tbl"))),
          nav_panel("Compare by contrast",
            uiOutput(ns("cmp_wrap")),
            tags$hr(), h5("Per-cell-type statistics"), DT::DTOutput(ns("stat_tbl")))
        )
      )
    )
  )
}

tab_deconv_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    ns  <- session$ns
    dec <- load_rds("deconv")

    processing <- function()
      div(class = "banner", style = "border-left-color:#f39c12;",
          icon("gear"), tags$b(" Processing. "),
          "The single-cell reference is being downloaded and the deconvolution is running. ",
          "Cell-type proportions and the contrast comparisons will appear here once complete.")

    # ---- reference provenance caption ----
    output$prov <- renderUI({
      if (is.null(dec)) return(NULL)
      r <- dec$ref_info
      div(class = "text-muted small mb-2",
          sprintf("Reference: %s — %s. %s cells, %d cell types, %s genes shared with the bulk. Method: %s. Values are the fraction of each cell type per sample (columns sum to ~1).",
                  r$source, r$accession, format(r$n_cells, big.mark=","),
                  r$n_celltypes, format(r$n_genes_shared, big.mark=","), r$method))
    })

    # ---- ordering helpers ----
    grp_order <- c("NP_F","NP_M","Sham_F","Sham_M","SNI_F","SNI_M")
    samp_order <- reactive({
      m <- dec$meta; m[order(match(m$group, grp_order), m$sample), "sample"]
    })
    ct_order <- reactive({                      # cell types by descending mean proportion
      p <- dec$prop_wide; names(sort(colMeans(p), decreasing = TRUE))
    })

    # full 18-type palette (consistent colours whether combined or split)
    ct_pal <- reactive({ cts <- ct_order(); stats::setNames(scales::hue_pal()(length(cts)), cts) })
    comp_cols <- function(which) {
      cts <- ct_order()
      if (which == "neuro") intersect(cts, NEURONAL_TYPES)
      else if (which == "nonneuro") setdiff(cts, NEURONAL_TYPES)
      else cts
    }
    mode_lab <- reactive(switch(input$compmode %||% "sample",
                                group = "group", condition = "condition", "sample"))

    # aggregated matrix (rows per sample/group/condition) for a set of columns
    comp_matrix <- function(cols) {
      p <- dec$prop_wide[samp_order(), cols, drop = FALSE]
      meta <- dec$meta[match(rownames(p), dec$meta$sample), ]
      agg <- switch(input$compmode %||% "sample",
        group     = list(key = meta$group,     ord = intersect(grp_order, unique(meta$group))),
        condition = list(key = meta$condition, ord = intersect(c("NP","Sham","SNI"), unique(meta$condition))),
        NULL)
      if (!is.null(agg)) {
        mat <- t(vapply(agg$ord, function(k) colMeans(p[agg$key == k, , drop = FALSE]),
                        numeric(length(cols))))
        rownames(mat) <- agg$ord; colnames(mat) <- cols
        list(mat = mat, x = disp(agg$ord))
      } else list(mat = p, x = disp(rownames(p)))
    }

    # one interactive stacked-bar for a compartment ("all"/"neuro"/"nonneuro")
    comp_plotly <- function(which, title) {
      cols <- comp_cols(which); cm <- comp_matrix(cols); pal <- ct_pal()
      pw <- cm$mat[, cols, drop = FALSE]; xcat <- cm$x
      p <- plotly::plot_ly()
      for (ct in cols) {
        y <- as.numeric(pw[, ct])
        p <- plotly::add_bars(p, x = xcat, y = y, name = ct,
          marker = list(color = pal[[ct]], line = list(color = "white", width = 0.4)),
          text = sprintf("<b>%s</b><br>%s<br>proportion: %.3f", ct_full(ct), xcat, y),
          hoverinfo = "text")
      }
      plotly::layout(p, barmode = "stack", title = title,
        xaxis = list(title = "", categoryorder = "array", categoryarray = xcat, tickangle = -45),
        yaxis = list(title = "Estimated proportion"),
        legend = list(title = list(text = "Cell type"))) |> ply_pub()
    }

    # download figure: all types, or the two compartments as facets
    comp_plot <- function() {
      cts <- ct_order(); pal <- ct_pal(); xlv <- comp_matrix(cts)$x
      build_long <- function(which, comp) {
        cols <- comp_cols(which); cm <- comp_matrix(cols)
        do.call(rbind, lapply(seq_len(nrow(cm$mat)), function(i)
          data.frame(x = cm$x[i], cellType = cols, proportion = as.numeric(cm$mat[i, cols]),
                     compartment = comp, stringsAsFactors = FALSE)))
      }
      split <- identical(input$compres, "split")
      d <- if (split) rbind(build_long("neuro", "Neuronal"), build_long("nonneuro", "Non-neuronal"))
           else build_long("all", "All")
      d$x <- factor(d$x, levels = xlv)
      d$cellType <- factor(d$cellType, levels = rev(cts))
      d$compartment <- factor(d$compartment, levels = c("Neuronal", "Non-neuronal"))
      g <- ggplot2::ggplot(d, ggplot2::aes(x, proportion, fill = cellType)) +
        ggplot2::geom_col(width = 0.82, colour = "white", linewidth = 0.15) +
        ggplot2::scale_fill_manual(values = pal, breaks = cts, name = "Cell type") +
        ggplot2::labs(title = sprintf("Estimated cell-type composition per %s", mode_lab()),
                      x = NULL, y = "Estimated proportion") +
        theme_pub(14) +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                       legend.key.size = ggplot2::unit(0.8, "lines"))
      if (split) g <- g + ggplot2::facet_wrap(~ compartment, ncol = 1, scales = "free_y")
      g
    }

    output$comp_wrap <- renderUI({
      if (is.null(dec)) return(processing())
      if (identical(input$compres, "split"))
        tagList(
          h6(sprintf("Neuronal subtypes — per %s", mode_lab())),
          plotly::plotlyOutput(ns("comp_neuro"), height = "340px"),
          tags$hr(),
          h6(sprintf("Non-neuronal cell types — per %s", mode_lab())),
          plotly::plotlyOutput(ns("comp_nonneuro"), height = "340px"))
      else plotly::plotlyOutput(ns("comp"), height = "520px")
    })
    output$comp <- plotly::renderPlotly({ req(!is.null(dec))
      comp_plotly("all", sprintf("Estimated cell-type composition per %s", mode_lab())) })
    output$comp_neuro <- plotly::renderPlotly({ req(!is.null(dec))
      comp_plotly("neuro", sprintf("Neuronal subtypes — per %s", mode_lab())) })
    output$comp_nonneuro <- plotly::renderPlotly({ req(!is.null(dec))
      comp_plotly("nonneuro", sprintf("Non-neuronal cell types — per %s", mode_lab())) })
    output$dl_comp <- downloadHandler(
      filename = function() paste0("deconv_composition.", input$fmt),
      content  = function(file) { gg_device(file, input$fmt, 10, 6); print(comp_plot()); grDevices::dev.off() })

    # proportions table (samples x cell types)
    output$comp_tbl <- DT::renderDT({
      req(!is.null(dec))
      p <- dec$prop_wide[samp_order(), ct_order(), drop = FALSE]
      m <- dec$meta[match(samp_order(), dec$meta$sample), ]
      out <- cbind(data.frame(Sample = disp(m$sample), Condition = disp(m$condition), Sex = m$sex,
                              stringsAsFactors = FALSE),
                   round(as.data.frame(p), 3))
      DT::datatable(out, rownames = FALSE, extensions = "Buttons",
        options = list(pageLength = 14, scrollX = TRUE, dom = "Bfrtip",
          buttons = list(list(extend = "excel", text = "⤓ Excel", filename = "deconv_proportions"))))
    })

    # ============ Compare by contrast: box-plots + stats ============
    cname <- reactive(names(CONTRAST_LABELS)[match(input$contrast, CONTRAST_LABELS)])
    output$sni_note <- renderUI(sni_note(cname()))

    # side membership + palette for the two sides of the selected contrast
    sides <- reactive({
      tok <- strsplit(cname(), "_vs_")[[1]]
      m <- dec$meta
      sset <- function(t) m$sample[m$group == t | m$condition == t]
      condA <- sub("_.*", "", tok[1]); condB <- sub("_.*", "", tok[2])
      if (condA != condB) {                       # differ by condition -> condition palette
        col <- unname(COND_COLORS[c(condA, condB)])
      } else {                                     # differ by sex -> sex palette
        sx <- sub(".*_", "", tok); col <- unname(SEX_COLORS[sx])
      }
      labs <- gsub("_", " ", disp(tok))
      list(tok = tok, A = sset(tok[1]), B = sset(tok[2]), labs = labs, col = col)
    })

    cmp_plot <- function() {
      s <- sides()
      d <- dec$prop_long[dec$prop_long$sample %in% c(s$A, s$B), ]
      d$side <- factor(ifelse(d$sample %in% s$A, s$labs[1], s$labs[2]), levels = s$labs)
      d$cellType <- factor(d$cellType, levels = ct_order())
      pal <- stats::setNames(s$col, s$labs)
      ggplot2::ggplot(d, ggplot2::aes(cellType, proportion, fill = side)) +
        ggplot2::geom_boxplot(outlier.size = 0.6, position = ggplot2::position_dodge(0.8),
                              width = 0.7, alpha = 0.85) +
        ggplot2::scale_fill_manual(values = pal, name = NULL) +
        ggplot2::labs(title = sprintf("Cell-type proportions — %s", disp(cname())),
                      x = NULL, y = "Estimated proportion") +
        theme_pub(14) +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                       legend.position = "top")
    }
    output$cmp_wrap <- renderUI({
      if (is.null(dec)) return(processing())
      tagList(
        plotly::plotlyOutput(ns("cmp"), height = "520px"),
        tags$small(class = "text-muted",
          "With n=3 per group the exact Wilcoxon test cannot reach p<0.05 for a single-sex contrast; ",
          "the pooled MINP-vs-Sham contrast (6 vs 6) is the best-powered. Read differences alongside the effect size (diff / log2FC)."))
    })
    # interactive grouped box-plots (all points shown) — hover spells out the full name
    output$cmp <- plotly::renderPlotly({
      req(!is.null(dec)); s <- sides()
      d <- dec$prop_long[dec$prop_long$sample %in% c(s$A, s$B), ]
      d$side <- ifelse(d$sample %in% s$A, s$labs[1], s$labs[2])
      p <- plotly::plot_ly()
      for (i in seq_along(s$labs)) {
        dd <- d[d$side == s$labs[i], ]
        p <- plotly::add_trace(p, type = "box", x = dd$cellType, y = dd$proportion,
          name = s$labs[i], color = I(s$col[i]),
          boxpoints = "all", jitter = 0.5, pointpos = 0, marker = list(size = 4),
          text = sprintf("<b>%s</b><br>%s · %s<br>proportion: %.3f",
                         ct_full(dd$cellType), disp(dd$sample), s$labs[i], dd$proportion),
          hoverinfo = "text")
      }
      plotly::layout(p, boxmode = "group",
        title = sprintf("Cell-type proportions — %s", disp(cname())),
        xaxis = list(title = "", categoryorder = "array", categoryarray = ct_order(), tickangle = -45),
        yaxis = list(title = "Estimated proportion"),
        legend = list(orientation = "h", x = 0, y = 1.08)) |> ply_pub()
    })
    output$dl_cmp <- downloadHandler(
      filename = function() paste0("deconv_compare_", disp(cname()), ".", input$fmt),
      content  = function(file) { gg_device(file, input$fmt, 10, 6); print(cmp_plot()); grDevices::dev.off() })

    # per-cell-type stats for the selected contrast
    output$stat_tbl <- DT::renderDT({
      req(!is.null(dec)); s <- sides()
      d <- dec$stats[dec$stats$contrast == cname(), ]
      d <- d[order(d$p_wilcox), ]
      out <- data.frame(
        `Cell type` = d$cellType,
        `Mean (A)`  = round(d$mean_A, 3),
        `Mean (B)`  = round(d$mean_B, 3),
        Difference  = round(d$diff, 3),
        log2FC      = round(d$log2FC, 2),
        `p (Wilcoxon)` = signif(d$p_wilcox, 3),
        `adj.p (Wilcoxon)` = signif(d$padj_wilcox, 3),
        `p (t-test)` = signif(d$p_t, 3),
        check.names = FALSE, stringsAsFactors = FALSE)
      DT::datatable(out, rownames = FALSE, extensions = "Buttons",
        caption = sprintf("A = %s (n=%d), B = %s (n=%d)", s$labs[1], length(s$A), s$labs[2], length(s$B)),
        options = list(pageLength = 20, scrollX = TRUE, dom = "Bfrtip",
          buttons = list(list(extend = "excel", text = "⤓ Excel",
                              filename = paste0("deconv_stats_", disp(cname()))))))
    })
  })
}
