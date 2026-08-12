# Cross-species gene lookup ("gene card") — one gene across mouse, rat & monkey.
# Type a mouse gene -> ortholog header, an all-species/all-contrast fold-change
# overview, and per-species expression-by-group boxplots + DE tables.
# Join key = mouse Ensembl gene id (feature_id) via ortholog_*_mouse.rds$mouse_gene_id.
# Data: de_long, expr, search_index, gene_annot, rat_de_long{,_noxy}, monkey_de_long{,_noxy},
#       rat_expr, monkey_expr, ortholog_rat_mouse, ortholog_monkey_mouse.
SPECIES_COL <- c(Mouse = "#E67E22", Rat = "#16A085", Macaque = "#8E44AD")

tab_gene_card_ui <- function(id, nested = FALSE) {
  ns <- NS(id)
  nav_panel(
    title = "Gene lookup", icon = icon("magnifying-glass-chart"),
    div(class = if (nested) "py-2" else "container-fluid py-3",
      if (!nested) h3("Cross-species gene lookup"),
      div(class = "banner",
          "Type one gene (mouse symbol or Ensembl id) to see its expression by group and its ",
          "differential-expression fold-change / adj.p across ", tags$b("every contrast in all three species"),
          " — mouse (MINP/SNI/Sham), rat (SNI/skin-lesion/naïve × DRG level) and macaque ",
          "(ipsi/contra/control). Species are joined through one-to-one orthologs."),
      layout_sidebar(
        sidebar = sidebar(width = 340,
          selectizeInput(ns("gene"), "Gene (mouse)", choices = NULL,
                         options = list(maxOptions = 50, placeholder = "type a symbol or ID…")),
          radioButtons(ns("xy"), "Sex chromosomes",
                       c("Include chrX/Y" = "all", "Exclude chrX/Y" = "noxy"), selected = "all", inline = TRUE),
          tags$small(class = "text-muted",
            "Sex-chromosome toggle swaps the rat & monkey DE between the all-genes and chrX/Y-excluded variants. Mouse DE shown is the uncorrected DESeq2 gene-level result.")
        ),
        div(
          uiOutput(ns("header")),
          tags$hr(),
          h5("Fold-change across all species & contrasts"),
          plotly::plotlyOutput(ns("overview"), height = "560px"),
          tags$small(class = "text-muted",
            "Each row = one contrast; x = log2 fold-change (dashed line = no change). Filled = significant (adj.p<0.05). Mouse SNI contrasts are exploratory (n=1/sex)."),
          tags$hr(),
          (if (nested) navset_pill else navset_tab)(
            nav_panel("Mouse",
              plotly::plotlyOutput(ns("mbox"), height = "420px"),
              tags$hr(), DT::DTOutput(ns("mtab"))),
            nav_panel("Rat",
              uiOutput(ns("rwrap"))),
            nav_panel("Macaque",
              uiOutput(ns("kwrap")))
          )
        )
      )
    )
  )
}

tab_gene_card_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    ns  <- session$ns
    mou <- load_rds("de_long"); if (!is.null(mou)) mou <- mou[mou$engine=="DESeq2" & mou$level=="gene", ]
    expr <- load_rds("expr")
    idx <- load_rds("search_index"); if (!is.null(idx)) idx <- idx[idx$level=="gene", ]
    rat <- load_rds("rat_de_long"); rat_noxy <- load_rds("rat_de_long_noxy")
    mon <- load_rds("monkey_de_long"); mon_noxy <- load_rds("monkey_de_long_noxy")
    re  <- load_rds("rat_expr"); ke <- load_rds("monkey_expr")
    ro  <- load_rds("ortholog_rat_mouse"); mo <- load_rds("ortholog_monkey_mouse")

    # server-side typeahead over mouse gene-level features (value = ENSMUSG id)
    observe({
      req(idx)
      updateSelectizeInput(session, "gene", choices = setNames(idx$feature_id, idx$label),
                           server = TRUE, selected = idx$feature_id[match("Atf3", idx$symbol)])
    })

    fid  <- reactive({ req(input$gene); input$gene })                    # mouse ENSMUSG id
    msym <- reactive({ s <- idx$symbol[match(fid(), idx$feature_id)]; if (is.na(s)) fid() else s })
    rid  <- reactive({ v <- as.character(ro$rat_gene_id[match(fid(), ro$mouse_gene_id)]); if (length(v)) v else NA_character_ })
    kid  <- reactive({ v <- as.character(mo$monkey_gene_id[match(fid(), mo$mouse_gene_id)]); if (length(v)) v else NA_character_ })
    ratd <- reactive({ if (identical(input$xy,"noxy") && !is.null(rat_noxy)) rat_noxy else rat })
    mond <- reactive({ if (identical(input$xy,"noxy") && !is.null(mon_noxy)) mon_noxy else mon })

    # ---------- header ----------
    output$header <- renderUI({
      req(fid())
      an <- annot_of(fid(), "gene")
      rsym <- ro$rat_symbol[match(fid(), ro$mouse_gene_id)];  rconf <- ro$confidence[match(fid(), ro$mouse_gene_id)]
      ksym <- mo$monkey_symbol[match(fid(), mo$mouse_gene_id)]; kconf <- mo$confidence[match(fid(), mo$mouse_gene_id)]
      orth_line <- function(lab, sym, conf) {
        if (is.na(sym)) span(tags$b(lab, ": "), tags$span(class="text-muted","no 1:1 ortholog"))
        else span(tags$b(lab, ": "), sym, if (!is.na(conf)) tags$small(class="text-muted", sprintf(" (conf %s)", conf)))
      }
      onXY <- !is.na(an$chr) && an$chr %in% c("X","Y")
      div(class = "banner",
        tags$h4(tags$b(msym()), tags$small(class="text-muted", paste0("  ", fid()))),
        if (!is.na(an$chr)) tags$span(sprintf("chromosome %s", an$chr),
          if (onXY) tags$span(class="badge bg-warning text-dark ms-1", "sex chromosome")),
        if (!is.na(an$description)) div(tags$small(em(an$description))),
        tags$hr(class="my-1"),
        div(class="d-flex flex-wrap gap-3",
            orth_line("Rat ortholog", rsym, rconf),
            orth_line("Macaque ortholog", ksym, kconf)))
    })

    # ---------- cross-species fold-change overview ----------
    overview_df <- reactive({
      req(fid())
      grab <- function(de, key_col, key, labels, species) {
        if (is.null(de) || is.na(key)) return(NULL)
        r <- de[de[[key_col]] == key, ]
        if (!nrow(r)) return(NULL)
        data.frame(species = species, contrast = as.character(r$contrast),
                   label = unname(labels[as.character(r$contrast)]),
                   log2FC = r$log2FC, padj = r$padj, stringsAsFactors = FALSE)
      }
      d <- rbind(
        grab(mou,    "feature_id", fid(),  CONTRAST_LABELS,        "Mouse"),
        grab(ratd(), "feature_id", rid(),  RAT_CONTRAST_LABELS,    "Rat"),
        grab(mond(), "feature_id", kid(),  MONKEY_CONTRAST_LABELS, "Macaque"))
      req(!is.null(d), nrow(d) > 0)
      d$species <- factor(d$species, levels = c("Mouse","Rat","Macaque"))
      d <- d[order(d$species, !is.na(d$padj) & d$padj < 0.05, abs(d$log2FC)), ]
      # spell the species out — abbreviating made Mouse and Monkey both read "Mo"
      d$row <- factor(paste0(as.character(d$species), " · ", d$label),
                      levels = paste0(as.character(d$species), " · ", d$label))
      d$sig <- !is.na(d$padj) & d$padj < 0.05
      d
    })
    output$overview <- plotly::renderPlotly({
      d <- overview_df()
      d$hover <- sprintf("<b>%s</b> — %s<br>log2FC %.2f · adj.p %s%s",
                         as.character(d$species), d$label, d$log2FC,
                         ifelse(is.na(d$padj), "NA", signif(d$padj,3)),
                         ifelse(d$species=="Mouse" & vapply(d$contrast, is_sni, logical(1)), "<br><i>SNI: exploratory (n=1/sex)</i>", ""))
      p <- plotly::plot_ly()
      for (sp in levels(d$species)) { dd <- d[d$species==sp, ]; if (!nrow(dd)) next
        p <- plotly::add_trace(p, data=dd, x=~log2FC, y=~row, type="scatter", mode="markers", name=sp,
              marker=list(size=11, color=SPECIES_COL[[sp]],
                          opacity=ifelse(dd$sig, 1, 0.35),
                          line=list(width=ifelse(dd$sig,1.1,0), color="#333")),
              text=~hover, hoverinfo="text") }
      plotly::layout(p, title=sprintf("%s — fold-change by contrast", msym()),
        xaxis=list(title="log2 fold-change", zeroline=FALSE),
        yaxis=list(title="", categoryorder="array", categoryarray=levels(d$row), automargin=TRUE),
        shapes=list(list(type="line", x0=0, x1=0, yref="paper", y0=0, y1=1, line=list(color="grey60", dash="dash"))),
        legend=list(orientation="h")) |> ply_pub()
    })

    # ---------- reusable expression boxplot + DE table ----------
    expr_box <- function(E, gid, title) {
      req(!is.null(E), !is.na(gid), gid %in% rownames(E$vst))   # blank if gene absent in this species
      df <- data.frame(expr = E$vst[gid, ], E$meta, stringsAsFactors = FALSE)
      df$.x <- factor(disp(df$group), levels = disp(sort(unique(as.character(df$group)))))
      sc <- color_spec(df$sex, "sex"); df$.sexf <- factor(df$sex, levels = sc$levels)
      plotly::plot_ly(df, x=~.x, y=~expr, color=~.sexf, colors=sc$colors, type="box",
                      boxpoints="all", jitter=0.5, pointpos=0,
                      text=~paste0(disp(sample_id)," (",disp(condition),", ",sex,")"), hoverinfo="text+y") |>
        plotly::layout(title=title, xaxis=list(title="", tickangle=-45),
                       yaxis=list(title="VST expression"), legend=list(title=list(text="Sex"))) |> ply_pub()
    }
    de_table <- function(de, key, labels, mouse=FALSE) {
      req(!is.null(de), !is.na(key))                            # blank if no ortholog in this species
      r <- de[de$feature_id == key, ]; req(nrow(r) > 0)
      r <- r[order(r$padj), ]
      out <- data.frame(Contrast = unname(labels[as.character(r$contrast)]),
                        log2FC = round(r$log2FC, 2), `adj.p` = signif(r$padj, 3),
                        check.names = FALSE, stringsAsFactors = FALSE)
      if (mouse) out$Note <- ifelse(vapply(as.character(r$contrast), is_sni, logical(1)), "exploratory (SNI n=1/sex)", "")
      DT::datatable(out, rownames = FALSE, extensions = "Buttons",
        options = list(pageLength = 12, scrollX = TRUE, dom = "Bfrtip",
          buttons = list(list(extend="excel", text="⤓ Excel", filename=paste0("gene_DE"))))) |>
        DT::formatStyle("adj.p", fontWeight = DT::styleInterval(0.05, c("bold","normal")))
    }

    # ---------- Mouse ----------
    output$mbox <- plotly::renderPlotly({ req(expr, fid()); expr_box(expr$gene, fid(), paste0("Mouse — ", msym())) })
    output$mtab <- DT::renderDT({ req(mou, fid()); de_table(mou, fid(), CONTRAST_LABELS, mouse=TRUE) })

    # ---------- Rat ----------
    output$rwrap <- renderUI({
      if (is.na(rid())) return(coming_soon("This mouse gene has no 1:1 rat ortholog."))
      tagList(plotly::plotlyOutput(ns("rbox"), height="420px"), tags$hr(), DT::DTOutput(ns("rtab")))
    })
    output$rbox <- plotly::renderPlotly({ req(re, rid()); expr_box(re, rid(), paste0("Rat — ", rid())) })
    output$rtab <- DT::renderDT({ req(rid()); de_table(ratd(), rid(), RAT_CONTRAST_LABELS) })

    # ---------- Monkey ----------
    output$kwrap <- renderUI({
      if (is.na(kid())) return(coming_soon("This mouse gene has no 1:1 macaque ortholog."))
      tagList(plotly::plotlyOutput(ns("kbox"), height="420px"), tags$hr(), DT::DTOutput(ns("ktab")))
    })
    output$kbox <- plotly::renderPlotly({ req(ke, kid()); expr_box(ke, kid(), paste0("Macaque — ", kid())) })
    output$ktab <- DT::renderDT({ req(kid()); de_table(mond(), kid(), MONKEY_CONTRAST_LABELS) })
  })
}
