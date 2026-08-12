# Cross-species (monkey) tab — mouse vs macaque DRG. Sub-tabs: PCA, concordance
# scatter, ortholog-filtered Venn, monkey volcano, monkey deconvolution, DE table.
# Data: monkey_de_long{,_noxy}.rds, monkey_pca.rds, ortholog_monkey_mouse.rds,
#       de_long.rds (mouse), monkey_deconv.rds.  Transection axotomy model:
#       ipsi = injured L4-6 side, contra = contralateral, control = naïve.
tab_xspecies_monkey_ui <- function(id) {
  ns <- NS(id)
  nav_panel(
    title = "Cross-species (monkey)", icon = icon("shuffle"),
    div(class = "container-fluid py-3",
      h3("Cross-species — mouse vs macaque DRG"),
      div(class = "banner",
          "Macaque (", tags$i("Macaca fascicularis"), ") DRG bulk RNA-seq (Ensembl Macaca_fascicularis_6.0; ",
          tags$b("ipsi / contra / control"), " × sex, L4-6 nerve-transection axotomy) mapped to mouse via ",
          tags$b("one-to-one orthologs"), ". ipsi = injured side, contra = contralateral, control = naïve."),
      layout_sidebar(
        sidebar = sidebar(width = 340,
          selectInput(ns("mc"), "Mouse contrast", choices = unname(CONTRAST_LABELS),
                      selected = unname(CONTRAST_LABELS["SNI_vs_Sham"])),
          selectInput(ns("kc"), "Monkey contrast", choices = unname(MONKEY_CONTRAST_LABELS),
                      selected = unname(MONKEY_CONTRAST_LABELS["ipsi_vs_control"])),
          sliderInput(ns("padj"), "DEG adj.p cutoff", 0.01, 0.25, 0.05, 0.01),
          sliderInput(ns("lfc"), "|log2FC| cutoff", 0, 3, 1, 0.1),
          tags$hr(),
          radioButtons(ns("kxy"), "Monkey sex chromosomes",
                       c("Include chrX/Y"="all", "Exclude chrX/Y"="noxy"), selected = "all", inline = TRUE),
          tags$small(class = "text-muted",
            "Excluding chrX/Y recomputes the monkey DE (DESeq2) and PCA without sex-chromosome genes (the macaque assembly is female-derived: X only). Affects concordance, Venn, volcano, PCA and the DE table (not deconvolution)."),
          tags$hr(),
          radioButtons(ns("fmt"), "Save figures as", c("PNG"="png","SVG"="svg","PDF"="pdf"), inline = TRUE)
        ),
        navset_tab(
          nav_panel("Monkey PCA",
            div(class = "row g-2 mt-1",
              div(class="col-md-4", radioButtons(ns("pcol"), "Colour by",
                   c("Condition"="condition","Sex"="sex"), selected="condition", inline=TRUE)),
              div(class="col-md-3", selectInput(ns("pcx"), "x-axis", paste0("PC",1:6), "PC1")),
              div(class="col-md-3", selectInput(ns("pcy"), "y-axis", paste0("PC",1:6), "PC2"))),
            plotly::plotlyOutput(ns("pca"), height = "500px"),
            tags$small(class = "text-muted",
              "VST-transformed macaque counts, top-1000 most-variable genes. Toggle chrX/Y in the sidebar to see whether the dominant axes are sex-chromosome driven.")),
          nav_panel("Concordance",
            uiOutput(ns("stats")),
            plotly::plotlyOutput(ns("scatter"), height = "500px"),
            tags$small(class = "text-muted", "Each point = one ortholog. Dotted = identity, dashed = fit. ",
              "Orange = DEG in both species; positive correlation + shared same-direction DEGs = conserved.")),
          nav_panel("Ortholog Venn",
            uiOutput(ns("venn_wrap")),
            downloadButton(ns("dl_venn"), "Download Venn", class = "btn-sm btn-outline-secondary mt-1"),
            tags$hr(), h6("Shared / species-specific orthologous DEGs"), DT::DTOutput(ns("venn_tbl"))),
          nav_panel("Monkey volcano",
            div(class = "row g-2 mt-1",
              div(class="col-md-4", sliderInput(ns("nlab"), "Label top-N", 0, 60, 20, 5)),
              div(class="col-md-4", sliderInput(ns("vx"), "x-limit |log2FC|", 1, 12, 8, 0.5)),
              div(class="col-md-4", sliderInput(ns("vy"), "y-limit -log10 p", 2, 60, 30, 2))),
            plotOutput(ns("volcano"), height = "500px"),
            downloadButton(ns("dl_volcano"), "Download volcano", class = "btn-sm btn-outline-secondary")),
          nav_panel("Monkey deconvolution",
            uiOutput(ns("kdec_wrap"))),
          nav_panel("Monkey DE table", DT::DTOutput(ns("mon_tbl")))
        )
      )
    )
  )
}

tab_xspecies_monkey_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    ns   <- session$ns
    mon  <- load_rds("monkey_de_long")                 # with chrX/Y (default)
    mon_noxy <- load_rds("monkey_de_long_noxy")         # chrX/Y-removed variant
    kpca <- load_rds("monkey_pca")                       # list(withXY, noXY)
    orth <- load_rds("ortholog_monkey_mouse")
    mou  <- load_rds("de_long"); if (!is.null(mou)) mou <- mou[mou$engine=="DESeq2" & mou$level=="gene", ]
    kdec <- load_rds("monkey_deconv")
    mck  <- reactive(names(CONTRAST_LABELS)[match(input$mc, CONTRAST_LABELS)])
    kck  <- reactive(names(MONKEY_CONTRAST_LABELS)[match(input$kc, MONKEY_CONTRAST_LABELS)])
    k2m  <- function(ids) orth$mouse_symbol[match(ids, orth$monkey_gene_id)]   # monkey gene id -> mouse symbol
    mond <- reactive({ if (identical(input$kxy,"noxy") && !is.null(mon_noxy)) mon_noxy else mon })  # active monkey DE

    # ---------- joined ortholog table for the selected pair ----------
    xs <- reactive({
      req(mon, orth, mou, input$mc, input$kc)
      kd <- mond()
      k <- kd[kd$contrast == kck(), c("monkey_gene_id","symbol","log2FC","padj")]
      k$mouse_symbol <- k2m(k$monkey_gene_id); k <- k[!is.na(k$mouse_symbol), ]
      m <- mou[mou$contrast == mck(), c("symbol","log2FC","padj")]
      j <- merge(k, m, by.x="mouse_symbol", by.y="symbol", suffixes=c("_mon","_mou"))
      j[!is.na(j$log2FC_mon) & !is.na(j$log2FC_mou), ]
    })

    # ================= 0. Monkey PCA =================
    output$pca <- plotly::renderPlotly({
      req(kpca)
      pv <- if (identical(input$kxy,"noxy") && !is.null(kpca$noXY)) kpca$noXY else kpca$withXY
      d <- pv$coords; ve <- pv$var_explained
      xc <- input$pcx %||% "PC1"; yc <- input$pcy %||% "PC2"; req(xc %in% names(d), yc %in% names(d))
      cby <- input$pcol %||% "condition"
      d$x <- d[[xc]]; d$y <- d[[yc]]; d$grp <- as.character(d[[cby]])
      d$hover <- sprintf("<b>%s</b><br>%s / %s<br>%s: %.1f · %s: %.1f",
                         d$sample, d$condition, d$sex, xc, d$x, yc, d$y)
      xi <- as.integer(sub("PC","",xc)); yi <- as.integer(sub("PC","",yc))
      pca_scatter(d, "x", "y", "grp",
        title=sprintf("Macaque DRG PCA — %s (%s genes%s)", cby, format(pv$ntop, big.mark=","),
                      if (identical(input$kxy,"noxy")) ", chrX/Y excluded" else ""),
        xlab=sprintf("%s (%.1f%%)", xc, ve[xi]), ylab=sprintf("%s (%.1f%%)", yc, ve[yi]),
        hovertext=d$hover, legend_title=cby)
    })

    # ================= 1. Concordance =================
    output$stats <- renderUI({
      if (is.null(mon) || is.null(mou)) return(coming_soon("Monkey data pending."))
      j <- xs(); req(nrow(j) >= 3)
      pe <- suppressWarnings(cor.test(j$log2FC_mou, j$log2FC_mon))
      sp <- suppressWarnings(cor.test(j$log2FC_mou, j$log2FC_mon, method="spearman"))
      sm <- !is.na(j$padj_mou)&j$padj_mou<input$padj; sk <- !is.na(j$padj_mon)&j$padj_mon<input$padj
      both <- j[sm&sk,]; same <- sum(sign(both$log2FC_mou)==sign(both$log2FC_mon))
      div(class="banner", sprintf("n = %d orthologs · Pearson r = %.3f (p = %.2g) · Spearman ρ = %.3f · shared DEGs = %d (%d same-direction)",
        nrow(j), unname(pe$estimate), pe$p.value, unname(sp$estimate), nrow(both), same))
    })
    output$scatter <- plotly::renderPlotly({
      j <- xs(); req(nrow(j) >= 3)
      sm <- !is.na(j$padj_mou)&j$padj_mou<input$padj; sk <- !is.na(j$padj_mon)&j$padj_mon<input$padj
      cols <- c("DEG in both"="#E67E22","DEG in one"="#3477eb","n.s."="#B0B8BF")
      j$class <- factor(ifelse(sm&sk,"DEG in both",ifelse(sm|sk,"DEG in one","n.s.")), levels=names(cols))
      j$hover <- paste0("<b>",j$mouse_symbol,"</b> (monkey ",j$symbol,")<br>mouse: ",round(j$log2FC_mou,2),"<br>monkey: ",round(j$log2FC_mon,2))
      rng <- range(c(j$log2FC_mou,j$log2FC_mon),na.rm=TRUE); fit<-stats::lm(log2FC_mon~log2FC_mou,j); fy<-as.numeric(predict(fit,data.frame(log2FC_mou=rng)))
      plotly::plot_ly(j,x=~log2FC_mou,y=~log2FC_mon,type="scatter",mode="markers",color=~class,colors=cols,text=~hover,hoverinfo="text",marker=list(size=6,line=list(width=.3,color="#fff"))) |>
        plotly::layout(title=sprintf("mouse %s vs monkey %s", input$mc, input$kc),
          xaxis=list(title=paste0("mouse log2FC · ",input$mc)), yaxis=list(title=paste0("monkey log2FC · ",input$kc)),
          shapes=list(list(type="line",x0=rng[1],x1=rng[2],y0=rng[1],y1=rng[2],line=list(color="grey70",dash="dot")),
                      list(type="line",x0=rng[1],x1=rng[2],y0=fy[1],y1=fy[2],line=list(color="black",dash="dash")))) |> ply_pub()
    })

    # ================= 2. Ortholog Venn =================
    venn_sets <- reactive({
      req(mon, orth, mou)
      md <- mou[mou$contrast==mck(),]; md <- md[!is.na(md$padj)&md$padj<input$padj&abs(md$log2FC)>=input$lfc,]
      md <- md[md$symbol %in% orth$mouse_symbol, ]                 # mouse DEGs that HAVE a monkey ortholog
      kd <- mond(); kd <- kd[kd$contrast==kck(),]; kd <- kd[!is.na(kd$padj)&kd$padj<input$padj&abs(kd$log2FC)>=input$lfc,]
      kd_ms <- unique(k2m(kd$monkey_gene_id)); kd_ms <- kd_ms[!is.na(kd_ms)]
      list(Mouse = unique(md$symbol), Monkey = kd_ms)
    })
    venn_plot <- reactive({
      s <- venn_sets(); req(length(s$Mouse)+length(s$Monkey) > 0)
      ggVennDiagram::ggVennDiagram(list(Mouse=s$Mouse, Monkey=s$Monkey), label="count", label_alpha=0, set_size=5) +
        ggplot2::scale_fill_gradient(low="#f7fbff", high="#E67E22", name="orthologous DEGs") +
        ggplot2::scale_x_continuous(expand=ggplot2::expansion(mult=.25)) +
        ggplot2::labs(title=sprintf("Orthologous DEGs — mouse %s vs monkey %s\n(adj.p<%.3g, |log2FC|>=%g)", input$mc, input$kc, input$padj, input$lfc))
    })
    output$venn_wrap <- renderUI({ if (is.null(mon)||is.null(mou)) return(coming_soon("Monkey data pending.")); plotOutput(ns("venn"), height="440px") })
    output$venn <- renderPlot({ print(venn_plot()) })
    output$dl_venn <- downloadHandler(filename=function() sprintf("xspecies_monkey_venn_%s_vs_%s.%s", mck(), kck(), input$fmt),
      content=function(file){ gg_device(file,input$fmt,7,6); print(venn_plot()); grDevices::dev.off() })
    output$venn_tbl <- DT::renderDT({
      s <- venn_sets(); u <- union(s$Mouse, s$Monkey); req(length(u))
      region <- ifelse(u %in% s$Mouse & u %in% s$Monkey, "Both", ifelse(u %in% s$Mouse, "Mouse only", "Monkey only"))
      out <- data.frame(`Mouse symbol`=u, Region=region, check.names=FALSE, stringsAsFactors=FALSE)
      out <- out[order(out$Region!="Both", out$`Mouse symbol`), ]
      DT::datatable(out, rownames=FALSE, filter="top", extensions="Buttons",
        options=list(pageLength=15, scrollX=TRUE, dom="Bfrtip", buttons=list(list(extend="excel", text="⤓ Excel", filename="xspecies_monkey_venn_genes"))))
    })

    # ================= 3. Monkey volcano =================
    mon_volc <- reactive({
      req(mon); kd <- mond(); d <- kd[kd$contrast==kck(), ]   # has feature_id + symbol for volcano_gg_app
      volcano_gg_app(d, input$padj, input$lfc, input$nlab, input$vx, input$vy,
                     sprintf("Macaque — %s%s", input$kc, if (identical(input$kxy,"noxy")) " (chrX/Y excluded)" else ""))
    })
    output$volcano <- renderPlot({ req(mon); print(mon_volc()) })
    output$dl_volcano <- downloadHandler(filename=function() sprintf("monkey_volcano_%s.%s", kck(), input$fmt),
      content=function(file){ gg_device(file,input$fmt,8,6); print(mon_volc()); grDevices::dev.off() })

    # ================= 4. Monkey deconvolution =================
    kct <- reactive({ p <- kdec$prop_wide; names(sort(colMeans(p), decreasing=TRUE)) })
    kgrp_order <- c("control","contra","ipsi")
    output$kdec_wrap <- renderUI({
      if (is.null(kdec)) return(coming_soon("Monkey deconvolution is being computed; it will appear here shortly."))
      tagList(
        div(class="text-muted small mb-2",
            sprintf("Reference: %s. %s cells, %d cell types, %s genes shared. %s.",
                    kdec$ref_info$source, format(kdec$ref_info$n_cells,big.mark=","), kdec$ref_info$n_celltypes,
                    format(kdec$ref_info$n_genes_shared,big.mark=","), kdec$ref_info$method)),
        radioButtons(ns("kmode"), "Composition — group by",
                     c("Sample"="sample","Group (cond×sex)"="group","Condition"="condition"), inline=TRUE),
        navset_tab(
          nav_panel("Composition", plotly::plotlyOutput(ns("kcomp"), height="480px")),
          nav_panel("Compare (monkey contrast)",
            plotly::plotlyOutput(ns("kcmp"), height="460px"),
            tags$hr(), h6("Per-cell-type stats (selected monkey contrast)"), DT::DTOutput(ns("kstat")))
        ))
    })
    kcomp_mat <- reactive({
      p <- kdec$prop_wide[, kct(), drop=FALSE]; meta <- kdec$meta[match(rownames(p), kdec$meta$sample), ]
      if (identical(input$kmode,"group")) key <- meta$group
      else if (identical(input$kmode,"condition")) key <- meta$condition else key <- NULL
      if (!is.null(key)) { ord <- if (identical(input$kmode,"condition")) intersect(kgrp_order, unique(key)) else sort(unique(key))
        agg <- t(vapply(ord, function(k) colMeans(p[key==k,,drop=FALSE]), numeric(ncol(p)))); rownames(agg)<-ord
        list(mat=agg, x=ord) } else { so <- meta$sample[order(meta$condition, meta$sex)]
        list(mat=p[so,,drop=FALSE], x=so) }
    })
    output$kcomp <- plotly::renderPlotly({ req(kdec); cm<-kcomp_mat(); cts<-kct(); pw<-cm$mat[,cts,drop=FALSE]
      pal<-stats::setNames(scales::hue_pal()(length(cts)),cts); p<-plotly::plot_ly()
      for(ct in cts){ y<-as.numeric(pw[,ct]); p<-plotly::add_bars(p,x=cm$x,y=y,name=ct,marker=list(color=pal[[ct]],line=list(color="white",width=.4)),
        text=sprintf("<b>%s</b><br>%s<br>%.3f",ct_full(ct),cm$x,y),hoverinfo="text") }
      plotly::layout(p,barmode="stack",title=sprintf("Macaque cell-type composition per %s", input$kmode),
        xaxis=list(title="",categoryorder="array",categoryarray=cm$x,tickangle=-45),
        yaxis=list(title="Estimated proportion",range=c(0,1.001)),legend=list(title=list(text="Cell type"))) |> ply_pub() })
    output$kcmp <- plotly::renderPlotly({ req(kdec); s<-monkey_sides(kck(), kdec$meta); req(length(s$A)>0,length(s$B)>0)
      d<-kdec$prop_long[kdec$prop_long$sample %in% c(s$A,s$B),]; labs<-strsplit(input$kc," vs ")[[1]]
      if(length(labs)<2) labs<-c("A","B"); d$side<-ifelse(d$sample %in% s$A, labs[1], labs[2])
      p<-plotly::plot_ly()
      for(i in seq_along(labs)){ dd<-d[d$side==labs[i],]; col<-c("#E67E22","#3477eb")[i]
        p<-plotly::add_trace(p,type="box",x=dd$cellType,y=dd$proportion,name=labs[i],color=I(col),boxpoints="all",jitter=.4,pointpos=0,marker=list(size=4),
          text=sprintf("<b>%s</b><br>%s<br>%.3f",ct_full(dd$cellType),dd$sample,dd$proportion),hoverinfo="text") }
      plotly::layout(p,boxmode="group",title=sprintf("Macaque cell-type proportions — %s", input$kc),
        xaxis=list(title="",categoryorder="array",categoryarray=kct(),tickangle=-45),yaxis=list(title="Estimated proportion"),legend=list(orientation="h")) |> ply_pub() })
    output$kstat <- DT::renderDT({ req(kdec); d<-kdec$stats[kdec$stats$contrast==kck(),]; d<-d[order(d$p_wilcox),]
      out<-data.frame(`Cell type`=d$cellType,`Mean A`=round(d$mean_A,3),`Mean B`=round(d$mean_B,3),Diff=round(d$diff,3),
        log2FC=round(d$log2FC,2),`p (Wilcoxon)`=signif(d$p_wilcox,3),`adj.p`=signif(d$padj_wilcox,3),check.names=FALSE,stringsAsFactors=FALSE)
      DT::datatable(out,rownames=FALSE,extensions="Buttons",options=list(pageLength=20,scrollX=TRUE,dom="Bfrtip",
        buttons=list(list(extend="excel",text="⤓ Excel",filename=paste0("monkey_deconv_",kck()))))) })

    # ================= 5. Monkey DE table =================
    output$mon_tbl <- DT::renderDT({
      req(mon); kd <- mond(); k <- kd[kd$contrast==kck(), ]; k$mouse_ortholog <- k2m(k$monkey_gene_id); k <- k[order(k$padj), ]
      out <- data.frame(`Monkey symbol`=k$symbol, `Mouse ortholog`=k$mouse_ortholog, `Monkey gene`=k$monkey_gene_id,
                        log2FC=round(k$log2FC,2), `adj.p`=signif(k$padj,3), check.names=FALSE, stringsAsFactors=FALSE)
      DT::datatable(out, rownames=FALSE, filter="top", extensions="Buttons",
        options=list(pageLength=15, scrollX=TRUE, dom="Bfrtip", buttons=list(list(extend="excel", text="⤓ Excel", filename=paste0("monkey_DE_", kck())))))
    })
  })
}
