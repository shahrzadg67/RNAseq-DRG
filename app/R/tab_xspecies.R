# Cross-species tab — mouse vs rat DRG. Sub-tabs: concordance scatter, ortholog-
# filtered Venn, rat volcano, rat deconvolution, rat DE table.
# Data: rat_de_long.rds, ortholog_rat_mouse.rds, de_long.rds (mouse), rat_deconv.rds.
tab_xspecies_ui <- function(id) {
  ns <- NS(id)
  nav_panel(
    title = "Cross-species (rat)", icon = icon("shuffle"),
    div(class = "container-fluid py-3",
      h3("Cross-species — mouse vs rat DRG"),
      div(class = "banner",
          "Rat DRG bulk RNA-seq (Ensembl GRCr8; SNI / skin-lesion / naïve × sex × DRG level ",
          tags$b("U/L"), ") mapped to mouse via ", tags$b("18,711 one-to-one orthologs"), ". ",
          "Note the injury response is concentrated in the ", tags$b("lower (L) DRG"),
          " — the L4/L5 level of the SNI model."),
      layout_sidebar(
        sidebar = sidebar(width = 340,
          selectInput(ns("mc"), "Mouse contrast", choices = unname(CONTRAST_LABELS),
                      selected = unname(CONTRAST_LABELS["SNI_vs_Sham"])),
          selectInput(ns("rc"), "Rat contrast", choices = unname(RAT_CONTRAST_LABELS),
                      selected = unname(RAT_CONTRAST_LABELS["SNI_vs_Naive_L"])),
          sliderInput(ns("padj"), "DEG adj.p cutoff", 0.01, 0.25, 0.05, 0.01),
          sliderInput(ns("lfc"), "|log2FC| cutoff", 0, 3, 1, 0.1),
          tags$hr(),
          radioButtons(ns("rxy"), "Rat sex chromosomes",
                       c("Include chrX/Y"="all", "Exclude chrX/Y"="noxy"), selected = "all", inline = TRUE),
          tags$small(class = "text-muted",
            "Excluding chrX/Y recomputes the rat DE (DESeq2) and PCA without sex-chromosome genes, so male/female contrasts and the composition axes are not driven by sex. Affects concordance, Venn, volcano, PCA and the DE table (not deconvolution)."),
          tags$hr(),
          radioButtons(ns("fmt"), "Save figures as", c("PNG"="png","SVG"="svg","PDF"="pdf"), inline = TRUE)
        ),
        navset_tab(
          nav_panel("Rat PCA",
            div(class = "row g-2 mt-1",
              div(class="col-md-4", radioButtons(ns("pcol"), "Colour by",
                   c("Condition"="condition","Sex"="sex","DRG level"="level"), selected="condition", inline=TRUE)),
              div(class="col-md-3", selectInput(ns("pcx"), "x-axis", paste0("PC",1:6), "PC1")),
              div(class="col-md-3", selectInput(ns("pcy"), "y-axis", paste0("PC",1:6), "PC2"))),
            plotly::plotlyOutput(ns("pca"), height = "500px"),
            tags$small(class = "text-muted",
              "VST-transformed rat counts, top-1000 most-variable genes. Toggle chrX/Y in the sidebar to see whether the dominant axes are sex-chromosome driven. ",
              "Symbol = DRG level (● upper / ▲ lower).")),
          nav_panel("Concordance",
            uiOutput(ns("stats")),
            plotly::plotlyOutput(ns("scatter"), height = "500px"),
            tags$small(class = "text-muted", "Each point = one ortholog. Dotted = identity, dashed = fit. ",
              "Orange = DEG in both species; positive correlation + shared same-direction DEGs = conserved.")),
          nav_panel("Ortholog Venn",
            uiOutput(ns("venn_wrap")),
            downloadButton(ns("dl_venn"), "Download Venn", class = "btn-sm btn-outline-secondary mt-1"),
            tags$hr(), h6("Shared / species-specific orthologous DEGs"), DT::DTOutput(ns("venn_tbl"))),
          nav_panel("Rat volcano",
            div(class = "row g-2 mt-1",
              div(class="col-md-4", sliderInput(ns("nlab"), "Label top-N", 0, 60, 20, 5)),
              div(class="col-md-4", sliderInput(ns("vx"), "x-limit |log2FC|", 1, 12, 8, 0.5)),
              div(class="col-md-4", sliderInput(ns("vy"), "y-limit -log10 p", 2, 60, 30, 2))),
            plotOutput(ns("volcano"), height = "500px"),
            downloadButton(ns("dl_volcano"), "Download volcano", class = "btn-sm btn-outline-secondary")),
          nav_panel("Rat deconvolution",
            uiOutput(ns("rdec_wrap"))),
          nav_panel("Rat DE table", DT::DTOutput(ns("rat_tbl")))
        )
      )
    )
  )
}

tab_xspecies_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    ns   <- session$ns
    rat  <- load_rds("rat_de_long")                 # with chrX/Y (default)
    rat_noxy <- load_rds("rat_de_long_noxy")         # chrX/Y-removed variant
    rpca <- load_rds("rat_pca")                       # list(withXY, noXY)
    orth <- load_rds("ortholog_rat_mouse")
    mou  <- load_rds("de_long"); if (!is.null(mou)) mou <- mou[mou$engine=="DESeq2" & mou$level=="gene", ]
    rdec <- load_rds("rat_deconv")
    mck  <- reactive(names(CONTRAST_LABELS)[match(input$mc, CONTRAST_LABELS)])
    rck  <- reactive(names(RAT_CONTRAST_LABELS)[match(input$rc, RAT_CONTRAST_LABELS)])
    r2m  <- function(ids) orth$mouse_symbol[match(ids, orth$rat_gene_id)]   # rat gene id -> mouse symbol
    ratd <- reactive({ if (identical(input$rxy,"noxy") && !is.null(rat_noxy)) rat_noxy else rat })  # active rat DE

    # ---------- joined ortholog table for the selected pair ----------
    xs <- reactive({
      req(rat, orth, mou, input$mc, input$rc)
      rd <- ratd()
      r <- rd[rd$contrast == rck(), c("rat_gene_id","symbol","log2FC","padj")]
      r$mouse_symbol <- r2m(r$rat_gene_id); r <- r[!is.na(r$mouse_symbol), ]
      m <- mou[mou$contrast == mck(), c("symbol","log2FC","padj")]
      j <- merge(r, m, by.x="mouse_symbol", by.y="symbol", suffixes=c("_rat","_mou"))
      j[!is.na(j$log2FC_rat) & !is.na(j$log2FC_mou), ]
    })

    # ================= 0. Rat PCA =================
    output$pca <- plotly::renderPlotly({
      req(rpca)
      pv <- if (identical(input$rxy,"noxy") && !is.null(rpca$noXY)) rpca$noXY else rpca$withXY
      d <- pv$coords; ve <- pv$var_explained
      xc <- input$pcx %||% "PC1"; yc <- input$pcy %||% "PC2"; req(xc %in% names(d), yc %in% names(d))
      cby <- input$pcol %||% "condition"
      d$x <- d[[xc]]; d$y <- d[[yc]]; d$grp <- as.character(d[[cby]])
      d$hover <- sprintf("<b>%s</b><br>%s / %s / %s<br>%s: %.1f · %s: %.1f",
                         d$sample, d$condition, d$sex, d$level, xc, d$x, yc, d$y)
      xi <- as.integer(sub("PC","",xc)); yi <- as.integer(sub("PC","",yc))
      # distinct colour + shape per group of the chosen variable, with a legend
      pca_scatter(d, "x", "y", "grp",
        title=sprintf("Rat DRG PCA — %s (%s genes%s)", cby, format(pv$ntop, big.mark=","),
                      if (identical(input$rxy,"noxy")) ", chrX/Y excluded" else ""),
        xlab=sprintf("%s (%.1f%%)", xc, ve[xi]), ylab=sprintf("%s (%.1f%%)", yc, ve[yi]),
        hovertext=d$hover, legend_title=cby)
    })

    # ================= 1. Concordance =================
    output$stats <- renderUI({
      if (is.null(rat) || is.null(mou)) return(coming_soon("Rat data pending."))
      j <- xs(); req(nrow(j) >= 3)
      pe <- suppressWarnings(cor.test(j$log2FC_mou, j$log2FC_rat))
      sp <- suppressWarnings(cor.test(j$log2FC_mou, j$log2FC_rat, method="spearman"))
      sm <- !is.na(j$padj_mou)&j$padj_mou<input$padj; sr <- !is.na(j$padj_rat)&j$padj_rat<input$padj
      both <- j[sm&sr,]; same <- sum(sign(both$log2FC_mou)==sign(both$log2FC_rat))
      div(class="banner", sprintf("n = %d orthologs · Pearson r = %.3f (p = %.2g) · Spearman ρ = %.3f · shared DEGs = %d (%d same-direction)",
        nrow(j), unname(pe$estimate), pe$p.value, unname(sp$estimate), nrow(both), same))
    })
    output$scatter <- plotly::renderPlotly({
      j <- xs(); req(nrow(j) >= 3)
      sm <- !is.na(j$padj_mou)&j$padj_mou<input$padj; sr <- !is.na(j$padj_rat)&j$padj_rat<input$padj
      cols <- c("DEG in both"="#E67E22","DEG in one"="#3477eb","n.s."="#B0B8BF")
      j$class <- factor(ifelse(sm&sr,"DEG in both",ifelse(sm|sr,"DEG in one","n.s.")), levels=names(cols))
      j$hover <- paste0("<b>",j$mouse_symbol,"</b> (rat ",j$symbol,")<br>mouse: ",round(j$log2FC_mou,2),"<br>rat: ",round(j$log2FC_rat,2))
      rng <- range(c(j$log2FC_mou,j$log2FC_rat),na.rm=TRUE); fit<-stats::lm(log2FC_rat~log2FC_mou,j); fy<-as.numeric(predict(fit,data.frame(log2FC_mou=rng)))
      plotly::plot_ly(j,x=~log2FC_mou,y=~log2FC_rat,type="scatter",mode="markers",color=~class,colors=cols,text=~hover,hoverinfo="text",marker=list(size=6,line=list(width=.3,color="#fff"))) |>
        plotly::layout(title=sprintf("mouse %s vs rat %s", input$mc, input$rc),
          xaxis=list(title=paste0("mouse log2FC · ",input$mc)), yaxis=list(title=paste0("rat log2FC · ",input$rc)),
          shapes=list(list(type="line",x0=rng[1],x1=rng[2],y0=rng[1],y1=rng[2],line=list(color="grey70",dash="dot")),
                      list(type="line",x0=rng[1],x1=rng[2],y0=fy[1],y1=fy[2],line=list(color="black",dash="dash")))) |> ply_pub()
    })

    # ================= 2. Ortholog Venn =================
    venn_sets <- reactive({
      req(rat, orth, mou)
      md <- mou[mou$contrast==mck(),]; md <- md[!is.na(md$padj)&md$padj<input$padj&abs(md$log2FC)>=input$lfc,]
      md <- md[md$symbol %in% orth$mouse_symbol, ]                 # mouse DEGs that HAVE a rat ortholog
      rd <- ratd(); rd <- rd[rd$contrast==rck(),]; rd <- rd[!is.na(rd$padj)&rd$padj<input$padj&abs(rd$log2FC)>=input$lfc,]
      rd_ms <- unique(r2m(rd$rat_gene_id)); rd_ms <- rd_ms[!is.na(rd_ms)]
      list(Mouse = unique(md$symbol), Rat = rd_ms)
    })
    venn_plot <- reactive({
      s <- venn_sets(); req(length(s$Mouse)+length(s$Rat) > 0)
      ggVennDiagram::ggVennDiagram(list(Mouse=s$Mouse, Rat=s$Rat), label="count", label_alpha=0, set_size=5) +
        ggplot2::scale_fill_gradient(low="#f7fbff", high="#E67E22", name="orthologous DEGs") +
        ggplot2::scale_x_continuous(expand=ggplot2::expansion(mult=.25)) +
        ggplot2::labs(title=sprintf("Orthologous DEGs — mouse %s vs rat %s\n(adj.p<%.3g, |log2FC|>=%g)", input$mc, input$rc, input$padj, input$lfc))
    })
    output$venn_wrap <- renderUI({ if (is.null(rat)||is.null(mou)) return(coming_soon("Rat data pending.")); plotOutput(ns("venn"), height="440px") })
    output$venn <- renderPlot({ print(venn_plot()) })
    output$dl_venn <- downloadHandler(filename=function() sprintf("xspecies_venn_%s_vs_%s.%s", mck(), rck(), input$fmt),
      content=function(file){ gg_device(file,input$fmt,7,6); print(venn_plot()); grDevices::dev.off() })
    output$venn_tbl <- DT::renderDT({
      s <- venn_sets(); u <- union(s$Mouse, s$Rat); req(length(u))
      region <- ifelse(u %in% s$Mouse & u %in% s$Rat, "Both", ifelse(u %in% s$Mouse, "Mouse only", "Rat only"))
      out <- data.frame(`Mouse symbol`=u, Region=region, check.names=FALSE, stringsAsFactors=FALSE)
      out <- out[order(out$Region!="Both", out$`Mouse symbol`), ]
      DT::datatable(out, rownames=FALSE, filter="top", extensions="Buttons",
        options=list(pageLength=15, scrollX=TRUE, dom="Bfrtip", buttons=list(list(extend="excel", text="⤓ Excel", filename="xspecies_venn_genes"))))
    })

    # ================= 3. Rat volcano =================
    rat_volc <- reactive({
      req(rat); rd <- ratd(); d <- rd[rd$contrast==rck(), ]
      volcano_gg_app(d, input$padj, input$lfc, input$nlab, input$vx, input$vy,
                     sprintf("Rat — %s%s", input$rc, if (identical(input$rxy,"noxy")) " (chrX/Y excluded)" else ""))
    })
    output$volcano <- renderPlot({ req(rat); print(rat_volc()) })
    output$dl_volcano <- downloadHandler(filename=function() sprintf("rat_volcano_%s.%s", rck(), input$fmt),
      content=function(file){ gg_device(file,input$fmt,8,6); print(rat_volc()); grDevices::dev.off() })

    # ================= 4. Rat deconvolution =================
    rct <- reactive({ p <- rdec$prop_wide; names(sort(colMeans(p), decreasing=TRUE)) })  # cell types by mean
    rgrp_order <- c("Naive","SL","SNI")
    output$rdec_wrap <- renderUI({
      if (is.null(rdec)) return(coming_soon("Rat deconvolution is being computed; it will appear here shortly."))
      tagList(
        div(class="text-muted small mb-2",
            sprintf("Reference: %s. %s cells, %d cell types, %s genes shared. %s.",
                    rdec$ref_info$source, format(rdec$ref_info$n_cells,big.mark=","), rdec$ref_info$n_celltypes,
                    format(rdec$ref_info$n_genes_shared,big.mark=","), rdec$ref_info$method)),
        radioButtons(ns("rmode"), "Composition — group by",
                     c("Sample"="sample","Group (cond×sex×level)"="group","Condition"="condition"), inline=TRUE),
        navset_tab(
          nav_panel("Composition", plotly::plotlyOutput(ns("rcomp"), height="480px")),
          nav_panel("Compare (rat contrast)",
            plotly::plotlyOutput(ns("rcmp"), height="460px"),
            tags$hr(), h6("Per-cell-type stats (selected rat contrast)"), DT::DTOutput(ns("rstat")))
        ))
    })
    rcomp_mat <- reactive({
      p <- rdec$prop_wide[, rct(), drop=FALSE]; meta <- rdec$meta[match(rownames(p), rdec$meta$sample), ]
      if (identical(input$rmode,"group")) key <- meta$group
      else if (identical(input$rmode,"condition")) key <- meta$condition else key <- NULL
      if (!is.null(key)) { ord <- if (identical(input$rmode,"condition")) intersect(rgrp_order, unique(key)) else sort(unique(key))
        agg <- t(vapply(ord, function(k) colMeans(p[key==k,,drop=FALSE]), numeric(ncol(p)))); rownames(agg)<-ord
        list(mat=agg, x=ord) } else { so <- meta$sample[order(meta$condition, meta$sex, meta$level)]
        list(mat=p[so,,drop=FALSE], x=so) }
    })
    output$rcomp <- plotly::renderPlotly({ req(rdec); cm<-rcomp_mat(); cts<-rct(); pw<-cm$mat[,cts,drop=FALSE]
      pal<-stats::setNames(scales::hue_pal()(length(cts)),cts); p<-plotly::plot_ly()
      for(ct in cts){ y<-as.numeric(pw[,ct]); p<-plotly::add_bars(p,x=cm$x,y=y,name=ct,marker=list(color=pal[[ct]],line=list(color="white",width=.4)),
        text=sprintf("<b>%s</b><br>%s<br>%.3f",ct_full(ct),cm$x,y),hoverinfo="text") }
      plotly::layout(p,barmode="stack",title=sprintf("Rat cell-type composition per %s", input$rmode),
        xaxis=list(title="",categoryorder="array",categoryarray=cm$x,tickangle=-45),
        yaxis=list(title="Estimated proportion",range=c(0,1.001)),legend=list(title=list(text="Cell type"))) |> ply_pub() })
    output$rcmp <- plotly::renderPlotly({ req(rdec); s<-rat_sides(rck(), rdec$meta); req(length(s$A)>0,length(s$B)>0)
      d<-rdec$prop_long[rdec$prop_long$sample %in% c(s$A,s$B),]; labs<-strsplit(input$rc," vs ")[[1]]
      if(length(labs)<2) labs<-c("A","B"); d$side<-ifelse(d$sample %in% s$A, labs[1], labs[2])
      p<-plotly::plot_ly()
      for(i in seq_along(labs)){ dd<-d[d$side==labs[i],]; col<-c("#E67E22","#3477eb")[i]
        p<-plotly::add_trace(p,type="box",x=dd$cellType,y=dd$proportion,name=labs[i],color=I(col),boxpoints="all",jitter=.4,pointpos=0,marker=list(size=4),
          text=sprintf("<b>%s</b><br>%s<br>%.3f",ct_full(dd$cellType),dd$sample,dd$proportion),hoverinfo="text") }
      plotly::layout(p,boxmode="group",title=sprintf("Rat cell-type proportions — %s", input$rc),
        xaxis=list(title="",categoryorder="array",categoryarray=rct(),tickangle=-45),yaxis=list(title="Estimated proportion"),legend=list(orientation="h")) |> ply_pub() })
    output$rstat <- DT::renderDT({ req(rdec); d<-rdec$stats[rdec$stats$contrast==rck(),]; d<-d[order(d$p_wilcox),]
      out<-data.frame(`Cell type`=d$cellType,`Mean A`=round(d$mean_A,3),`Mean B`=round(d$mean_B,3),Diff=round(d$diff,3),
        log2FC=round(d$log2FC,2),`p (Wilcoxon)`=signif(d$p_wilcox,3),`adj.p`=signif(d$padj_wilcox,3),check.names=FALSE,stringsAsFactors=FALSE)
      DT::datatable(out,rownames=FALSE,extensions="Buttons",options=list(pageLength=20,scrollX=TRUE,dom="Bfrtip",
        buttons=list(list(extend="excel",text="⤓ Excel",filename=paste0("rat_deconv_",rck()))))) })

    # ================= 5. Rat DE table =================
    output$rat_tbl <- DT::renderDT({
      req(rat); rd <- ratd(); r <- rd[rd$contrast==rck(), ]; r$mouse_ortholog <- r2m(r$rat_gene_id); r <- r[order(r$padj), ]
      out <- data.frame(`Rat symbol`=r$symbol, `Mouse ortholog`=r$mouse_ortholog, `Rat gene`=r$rat_gene_id,
                        log2FC=round(r$log2FC,2), `adj.p`=signif(r$padj,3), check.names=FALSE, stringsAsFactors=FALSE)
      DT::datatable(out, rownames=FALSE, filter="top", extensions="Buttons",
        options=list(pageLength=15, scrollX=TRUE, dom="Bfrtip", buttons=list(list(extend="excel", text="⤓ Excel", filename=paste0("rat_DE_", rck())))))
    })
  })
}
