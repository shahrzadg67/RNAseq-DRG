# Heatmap tab — pick a comparison's DEGs at a cutoff; rows = genes, cols =
# samples, colour = per-gene z-score of corrected VST, with clustering. Save.
tab_heatmap_ui <- function(id) {
  ns <- NS(id)
  nav_panel(
    title = "Heatmap", icon = icon("table-cells"),
    div(class = "container-fluid py-3",
      h3("DEG expression heatmap"),
      div(class = "banner",
          "Corrected analysis (batch + sex chromosomes), DESeq2. Choose a comparison and cutoff; ",
          "the heatmap shows those DEGs (rows) across samples (columns), coloured by per-gene ",
          "z-score of batch-corrected expression, with hierarchical clustering."),
      uiOutput(ns("sni_note")),
      layout_sidebar(
        sidebar = sidebar(width = 340,
          radioButtons(ns("level"), "Level", c("Gene"="gene"), inline = TRUE),  # gene-only (transcript skipped)
          selectInput(ns("contrast"), "Comparison (its DEGs)", choices = unname(CONTRAST_LABELS)),
          radioButtons(ns("samples"), "Samples to show",
                       c("In this comparison"        = "comparison",
                         "All samples (corrected)"   = "all",
                         "All samples (uncorrected)" = "sni")),
          sliderInput(ns("padj"), "Adj. p cutoff", min = 0.01, max = 0.25, value = 0.05, step = 0.01),
          sliderInput(ns("lfc"), "|log2FC| cutoff", min = 0, max = 3, value = 1, step = 0.1),
          sliderInput(ns("ngenes"), "Max genes (top by adj.p)", min = 10, max = 100, value = 50, step = 5),
          checkboxInput(ns("crow"), "Cluster genes (rows)", TRUE),
          checkboxInput(ns("ccol"), "Cluster samples (columns)", TRUE),
          radioButtons(ns("hfmt"), "Save as", c("PNG"="png","SVG"="svg","PDF"="pdf"), inline = TRUE),
          downloadButton(ns("dl_hm"), "Download heatmap", class = "btn-sm btn-outline-secondary")
        ),
        navset_tab(
          nav_panel("Heatmap", uiOutput(ns("hm_wrap"))),
          nav_panel("Gene table", uiOutput(ns("tbl_wrap")))
        )
      )
    )
  )
}

tab_heatmap_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    de_long <- load_rds("de_long_noxy")
    expr    <- load_rds("expr_noxy")   # corrected VST, 12 samples (batch + sex-chr removed)
    expr_all<- load_rds("expr")        # uncorrected VST, 14 samples (incl SNI) — for the SNI option
    ann     <- ANNOT
    cname <- reactive(names(CONTRAST_LABELS)[match(input$contrast, CONTRAST_LABELS)])
    output$sni_note <- renderUI(sni_note(cname()))

    # expression source: corrected unless SNI requested (SNI has no corrected values)
    esrc <- reactive({ if (input$samples == "sni") expr_all[[input$level]] else expr[[input$level]] })
    # which sample columns to show
    sample_cols <- reactive({
      E <- esrc(); req(E); m <- E$meta
      if (input$samples == "comparison") rownames(m)[m$group %in% CONTRAST_GROUPS[[cname()]]]
      else rownames(m)
    })

    # top DEGs (by padj) present in the active expression matrix
    genes <- reactive({
      req(de_long); E <- esrc()
      d <- de_long[de_long$engine=="DESeq2" & de_long$level==input$level & de_long$contrast==cname(), ]
      d <- d[!is.na(d$padj) & d$padj < input$padj & abs(d$log2FC) >= input$lfc, ]
      d <- d[order(d$padj), ]
      d <- d[d$feature_id %in% rownames(E$vst), ]
      utils::head(d$feature_id, input$ngenes)
    })

    hm_inputs <- reactive({
      g <- genes(); req(length(g) >= 2)
      E <- esrc(); cols <- sample_cols(); req(length(cols) >= 2)
      mat <- E$vst[g, cols, drop = FALSE]; meta <- E$meta[cols, , drop = FALSE]
      sym <- ann$symbol[match(g, ann$feature_id)]; rownames(mat) <- make.unique(ifelse(is.na(sym), g, sym))
      ac <- data.frame(Sex = meta$sex, Condition = disp(meta$condition), row.names = rownames(meta))
      list(mat = mat, ann_col = ac,
           ann_colors = list(Sex = c(F="#941200", M="#021893"),
                             Condition = c(Sham="#000000", MINP="#E67E22", SNI="#7F8C8D")))
    })

    draw_hm <- function() {
      h <- hm_inputs()
      pheatmap::pheatmap(h$mat, scale = "row",
        cluster_rows = input$crow, cluster_cols = input$ccol,
        annotation_col = h$ann_col, annotation_colors = h$ann_colors,
        color = grDevices::colorRampPalette(c("#3B4CC0","#F7F7F7","#B40426"))(101),
        border_color = NA, fontsize = 11,
        show_rownames = nrow(h$mat) <= 80, show_colnames = TRUE, labels_col = disp(colnames(h$mat)),
        main = sprintf("%s — %s · %d DEGs × %d samples%s (z-score)",
                       disp(cname()), input$level, nrow(h$mat), ncol(h$mat),
                       if (input$samples == "sni") " · SNI shown on uncorrected VST" else ""))
    }

    output$hm_wrap <- renderUI({
      if (is.null(de_long) || is.null(expr)) return(coming_soon("Corrected DE / expression pending."))
      if (length(genes()) < 2)
        return(div(class = "banner",
                   "Fewer than 2 DEGs at this cutoff — loosen the adj.p / |log2FC| cutoff."))
      plotOutput(ns("hm"), height = "640px")
    })
    output$hm <- renderPlot({ draw_hm() })
    output$dl_hm <- downloadHandler(
      filename = function() sprintf("heatmap_%s_%s.%s", disp(cname()), input$level, input$hfmt),
      content  = function(file) { gg_device(file, input$hfmt, 8, 9); draw_hm(); grDevices::dev.off() })

    # ---- gene table with per-sample z-scores ----
    output$tbl_wrap <- renderUI({
      if (is.null(de_long)) return(coming_soon("Pending.")); DT::DTOutput(ns("tbl"))
    })
    output$tbl <- DT::renderDT({
      g <- genes(); req(length(g) >= 1)
      E <- esrc(); cols <- sample_cols(); z <- t(scale(t(E$vst[g, cols, drop = FALSE])))
      base <- data.frame(feature_id = g,
                         symbol = ann$symbol[match(g, ann$feature_id)],
                         chr = ann$chr[match(g, ann$feature_id)], stringsAsFactors = FALSE)
      colnames(z) <- disp(colnames(z))
      out <- cbind(base, round(as.data.frame(z), 2))
      DT::datatable(out, rownames = FALSE, filter = "top", extensions = "Buttons",
        options = list(pageLength = 15, scrollX = TRUE, dom = "Bfrtip",
          buttons = list(list(extend = "excel", text = "⤓ Excel",
                              filename = paste0("heatmap_genes_", disp(cname()), "_", input$level)))))
    })
  })
}
