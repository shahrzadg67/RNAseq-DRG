# Cross-species meta-analysis — mouse + rat + macaque DRG in shared ortholog space.
# Sub-tabs:
#   PCA             : cross_pca.rds (raw / aligned); colour = variable, shape = species
#   Cell-type comp. : deconv.rds / rat_deconv.rds / monkey_deconv.rds — same 18 Renthal cell
#                     types. Every injury is named explicitly (which model, and for rat which
#                     DRG level); pooled panel on top, male/female panel below.
#   Injury markers  : the 47 conserved SNI+/MINP- genes as a within-species z-score heatmap
#                     across all 93 samples, plus the downloadable table.
#   Gene lookup     : the tab_gene_card module, nested here because it is a cross-species view.
CPCA_SP_LAB  <- c(Mouse="Mouse", Rat="Rat", Monkey="Macaque")   # "Monkey" is the internal key
# Marker shapes for the PCA "Shape by" selector — 7 distinct glyphs, enough for the
# widest encoding (7 conditions). Species keeps the app-wide ● mouse / ▲ rat / ■ macaque
# ordering because Mouse/Rat/Macaque sort into the first three slots.
CPCA_SHAPES  <- c("circle","triangle-up","square","diamond","star","cross","hexagon")
CPCA_DIVERGE <- list(list(0,"#2166AC"), list(0.5,"#F7F7F7"), list(1,"#B2182B"))

# ---- Injury contrasts, both sexes pooled -----------------------------------------------
# grp/ctl are `meta$group` values; key is the deconv stats contrast used for the ✱ p-value.
INJ_POOLED <- list(
  list(sp="Mouse",  inj=c("NP_M","NP_F"),           ctl=c("Sham_M","Sham_F"),
       key="NP_vs_Sham",      lab="Mouse MINP<br>vs Sham"),
  list(sp="Mouse",  inj=c("SNI_M","SNI_F"),         ctl=c("Sham_M","Sham_F"),
       key="SNI_vs_Sham",     lab="Mouse SNI<br>vs Sham"),
  list(sp="Rat",    inj=c("SL_M_L","SL_F_L"),       ctl=c("Naive_M_L","Naive_F_L"),
       key="SL_vs_Naive_L",   lab="Rat SL vs Naïve<br>lower DRG (L4/5)"),
  list(sp="Rat",    inj=c("SL_M_U","SL_F_U"),       ctl=c("Naive_M_U","Naive_F_U"),
       key="SL_vs_Naive_U",   lab="Rat SL vs Naïve<br>upper DRG"),
  list(sp="Rat",    inj=c("SNI_M_L","SNI_F_L"),     ctl=c("Naive_M_L","Naive_F_L"),
       key="SNI_vs_Naive_L",  lab="Rat SNI vs Naïve<br>lower DRG (L4/5)"),
  list(sp="Rat",    inj=c("SNI_M_U","SNI_F_U"),     ctl=c("Naive_M_U","Naive_F_U"),
       key="SNI_vs_Naive_U",  lab="Rat SNI vs Naïve<br>upper DRG"),
  list(sp="Monkey", inj=c("ipsi_M","ipsi_F"),       ctl=c("control_M","control_F"),
       key="ipsi_vs_control", lab="Macaque SN transection<br>ipsi vs control"))

# ---- Same injuries, male and female adjacent -------------------------------------------
# Mouse SNI has no sex-specific deconv test (SNI is n=1/sex, so `SNI_vs_Sham_M` pools both
# SNI animals against Sham males). Δ is still computed from the group means, but key=NA so
# no ✱ is drawn — the n= in the hover makes the n=1 explicit.
INJ_SEX <- list(
  list(sp="Mouse",  inj="NP_M",     ctl="Sham_M",    key="NP_M_vs_Sham_M",     lab="Mouse MINP<br>M"),
  list(sp="Mouse",  inj="NP_F",     ctl="Sham_F",    key="NP_F_vs_Sham_F",     lab="Mouse MINP<br>F"),
  list(sp="Mouse",  inj="SNI_M",    ctl="Sham_M",    key=NA,                   lab="Mouse SNI<br>M"),
  list(sp="Mouse",  inj="SNI_F",    ctl="Sham_F",    key=NA,                   lab="Mouse SNI<br>F"),
  list(sp="Rat",    inj="SL_M_L",   ctl="Naive_M_L", key="SL_M_vs_Naive_M_L",  lab="Rat SL lower<br>M"),
  list(sp="Rat",    inj="SL_F_L",   ctl="Naive_F_L", key="SL_F_vs_Naive_F_L",  lab="Rat SL lower<br>F"),
  list(sp="Rat",    inj="SL_M_U",   ctl="Naive_M_U", key="SL_M_vs_Naive_M_U",  lab="Rat SL upper<br>M"),
  list(sp="Rat",    inj="SL_F_U",   ctl="Naive_F_U", key="SL_F_vs_Naive_F_U",  lab="Rat SL upper<br>F"),
  list(sp="Rat",    inj="SNI_M_L",  ctl="Naive_M_L", key="SNI_M_vs_Naive_M_L", lab="Rat SNI lower<br>M"),
  list(sp="Rat",    inj="SNI_F_L",  ctl="Naive_F_L", key="SNI_F_vs_Naive_F_L", lab="Rat SNI lower<br>F"),
  list(sp="Rat",    inj="SNI_M_U",  ctl="Naive_M_U", key="SNI_M_vs_Naive_M_U", lab="Rat SNI upper<br>M"),
  list(sp="Rat",    inj="SNI_F_U",  ctl="Naive_F_U", key="SNI_F_vs_Naive_F_U", lab="Rat SNI upper<br>F"),
  list(sp="Monkey", inj="ipsi_M",   ctl="control_M", key="ipsi_M_vs_control_M", lab="Macaque SN<br>M"),
  list(sp="Monkey", inj="ipsi_F",   ctl="control_F", key="ipsi_F_vs_control_F", lab="Macaque SN<br>F"))

# ---- Named conditions for the stacked-composition views --------------------------------
COMP_GROUPS <- list(
  list(sp="Mouse",  grp=c("Sham_M","Sham_F"),       lab="Mouse<br>Sham",       ctl=TRUE),
  list(sp="Mouse",  grp=c("NP_M","NP_F"),           lab="Mouse<br>MINP",       ctl=FALSE),
  list(sp="Mouse",  grp=c("SNI_M","SNI_F"),         lab="Mouse<br>SNI",        ctl=FALSE),
  list(sp="Rat",    grp=c("Naive_M_L","Naive_F_L"), lab="Rat Naïve<br>lower",  ctl=TRUE),
  list(sp="Rat",    grp=c("SL_M_L","SL_F_L"),       lab="Rat SL<br>lower",     ctl=FALSE),
  list(sp="Rat",    grp=c("SNI_M_L","SNI_F_L"),     lab="Rat SNI<br>lower",    ctl=FALSE),
  list(sp="Rat",    grp=c("Naive_M_U","Naive_F_U"), lab="Rat Naïve<br>upper",  ctl=TRUE),
  list(sp="Rat",    grp=c("SL_M_U","SL_F_U"),       lab="Rat SL<br>upper",     ctl=FALSE),
  list(sp="Rat",    grp=c("SNI_M_U","SNI_F_U"),     lab="Rat SNI<br>upper",    ctl=FALSE),
  list(sp="Monkey", grp=c("control_M","control_F"), lab="Macaque<br>control",  ctl=TRUE),
  list(sp="Monkey", grp=c("contra_M","contra_F"),   lab="Macaque<br>contra",   ctl=FALSE),
  list(sp="Monkey", grp=c("ipsi_M","ipsi_F"),       lab="Macaque<br>ipsi",     ctl=FALSE))

# ---- Sample ordering for the 47-marker heatmap -----------------------------------------
MK_BLOCKS <- list(
  Mouse = list(rds="expr", nest=TRUE, orth=NULL, grp=list(
    list(f=function(m) m$condition=="Sham", lab="Mouse<br>Sham"),
    list(f=function(m) m$condition=="NP",   lab="Mouse<br>MINP"),
    list(f=function(m) m$condition=="SNI",  lab="Mouse<br>SNI"))),
  Rat = list(rds="rat_expr", nest=FALSE, orth="ortholog_rat_mouse", grp=list(
    list(f=function(m) m$condition=="Naive" & m$level=="DRGL", lab="Rat Naïve<br>lower"),
    list(f=function(m) m$condition=="SL"    & m$level=="DRGL", lab="Rat SL<br>lower"),
    list(f=function(m) m$condition=="SNI"   & m$level=="DRGL", lab="Rat SNI<br>lower"),
    list(f=function(m) m$condition=="Naive" & m$level=="DRGU", lab="Rat Naïve<br>upper"),
    list(f=function(m) m$condition=="SL"    & m$level=="DRGU", lab="Rat SL<br>upper"),
    list(f=function(m) m$condition=="SNI"   & m$level=="DRGU", lab="Rat SNI<br>upper"))),
  Monkey = list(rds="monkey_expr", nest=FALSE, orth="ortholog_monkey_mouse", grp=list(
    list(f=function(m) m$condition=="control", lab="Macaque<br>control"),
    list(f=function(m) m$condition=="contra",  lab="Macaque<br>contra"),
    list(f=function(m) m$condition=="ipsi",    lab="Macaque<br>ipsi"))))

tab_cross_pca_ui <- function(id, show_gene_card = TRUE) {
  ns <- NS(id)
  panels <- list(
    # ---------------- PCA ----------------
    nav_panel("PCA",
      div(class = "banner mt-2",
          "All ", tags$b("93 DRG samples"), " (14 mouse, 67 rat, 12 macaque) in one PCA over ",
          tags$b("11,296 one-to-one orthologs"), ". ",
          tags$b("Aligned"), " z-scores each gene within its species to remove the between-species baseline, ",
          "so injury/condition structure can line up across species; ", tags$b("Raw"),
          " keeps it, so PC1 mostly separates species."),
      layout_sidebar(
        sidebar = sidebar(width = 320,
          radioButtons(ns("mode"), "View",
            c("Aligned (remove species offset)"="aligned", "Raw (species separate)"="raw"), selected="aligned"),
          radioButtons(ns("color"), "Colour by",
            c("Condition"="condition", "Unified condition"="unified", "Species"="species",
              "Sex"="sex", "DRG level"="level", "Group (cond×sex×level)"="group"),
            selected="condition"),
          radioButtons(ns("shape"), "Shape by",
            c("Species"="species", "Sex"="sex", "Condition"="condition", "DRG level"="level"),
            selected="species"),
          selectInput(ns("pcx"), "x-axis", paste0("PC",1:6), "PC1"),
          selectInput(ns("pcy"), "y-axis", paste0("PC",1:6), "PC2"),
          tags$hr(),
          tags$small(class="text-muted",
            "Colour and shape are chosen independently, each with its own legend — pick two different ",
            "variables to separate every group at once (e.g. colour = Condition, shape = Sex). ",
            "Sex is always red = female / blue = male, matching the rest of the app. ",
            "Unified: nerve injury / pain model = mouse SNI + MINP, rat SNI, macaque ipsi; ",
            "control = Sham, Naïve/control, macaque contra; other = rat skin-lesion.")
        ),
        div(
          plotly::plotlyOutput(ns("pca"), height = "560px"),
          tags$small(class="text-muted",
            "In the Aligned view, nerve-injury samples from different species drifting together = a conserved injury axis.")
        )
      )
    ),
    # ---------------- Cell-type composition ----------------
    nav_panel("Cell-type composition",
      div(class = "banner mt-2",
          "The same ", tags$b("18 Renthal DRG cell types"), " were deconvolved (MuSiC) in all three species ",
          "(rat & macaque mapped to the mouse reference via orthologs), so proportions are directly comparable. ",
          "Every column names its own injury model and control — mouse ", tags$b("MINP"), " and ", tags$b("SNI"),
          " (both vs Sham), rat ", tags$b("skin-lesion (SL)"), " and ", tags$b("SNI"), " (both vs Naïve, and ",
          tags$b("lower L4/5 and upper DRG shown separately"), "), macaque ", tags$b("SN transection"),
          " (ipsi vs control)."),
      radioButtons(ns("dview"), "Show",
        c("Injury shift (Δ heatmap)"="delta", "Composition by condition"="comp", "Baseline (control) composition"="base"),
        selected="delta", inline=TRUE),
      checkboxGroupInput(ns("dclust"), "Cluster Δ heatmap", inline=TRUE,
        c("Cell types (rows)"="rows", "Injuries (cols)"="cols")),
      uiOutput(ns("dec_wrap")),
      tags$small(class="text-muted",
        "Δ = injury − control mean proportion for that exact pair of groups. A row the same colour across ",
        "columns = a compositional shift shared by those injuries. ",
        tags$b("✱"), " = adj.p<0.05 and ", tags$b("·"), " = nominal p<0.05 (Wilcoxon) — with 18 cell types ",
        "and these group sizes, FDR significance is only reachable in the well-powered rat, so the nominal ",
        "tier is shown to indicate where a trend exists. Blank/grey = that cell type is not detected in that ",
        "species (rather than detected and unchanged). Hover shows Δ, the group sizes and the exact adj.p.")
    ),
    # ---------------- Conserved injury markers ----------------
    nav_panel("Injury markers (SNI⁺/MINP⁻)",
      div(class = "banner mt-2",
          "Genes that are ", tags$b("up-regulated by nerve injury in ≥2 species"), " (mouse SNI, rat SNI, macaque ipsi; ",
          "adj.p<0.05, log2FC>1) yet ", tags$b("flat in mouse MINP"), " (|log2FC|<0.5, n.s.) — a conserved ",
          "nerve-injury signature that cleanly separates true nerve injury from the nucleus-pulposus (MINP) model. ",
          tags$b("Sprr1a is the archetype; Atf3 is the best pan-species marker"), " (Sprr1a has no macaque ortholog, ",
          "so it cannot enter a panel that requires all three species and is absent from the 47 below)."),
      checkboxGroupInput(ns("mclust"), "Cluster heatmap", inline=TRUE,
        c("Genes (rows)"="rows", "Samples within group (cols)"="cols")),
      plotly::plotlyOutput(ns("mk_heat"), height="720px"),
      tags$small(class="text-muted",
        "Each cell is one sample. Expression is VST, z-scored ", tags$b("per gene within each species"),
        " — mouse, rat and macaque differ in library size, depth and ortholog set, so raw values are not ",
        "comparable across species; z-scoring within a species centres every gene on its own baseline, which ",
        "is what makes the colour mean the same thing in all three blocks. Colour clipped at z = ±2. ",
        "Solid lines separate species, dotted lines separate conditions."),
      tags$hr(),
      DT::DTOutput(ns("markers")),
      downloadButton(ns("mk_dl"), "Download table (CSV)", class="btn-sm btn-outline-secondary mt-2"),
      tags$small(class="text-muted d-block mt-1",
        "log2FC orientation: injury vs its control (SNI vs Sham / SNI vs Naïve-lower / ipsi vs control), and MINP vs Sham.")
    ),
    # ---------------- MINP-only genes ----------------
    nav_panel("MINP-only genes",
      div(class="banner mt-2",
          "Genes associated with ", tags$b("MINP (nucleus pulposus)"), " — significant in MINP-vs-Sham ",
          "(pooled or within a sex) — flagged for whether each is ", tags$b("well-annotated"), ", ",
          tags$b("injury-negative"), " (not up in SNI), on a ", tags$b("sex chromosome"),
          ", differs in the baseline ", tags$b("Sham male-vs-female"), " comparison, or is ",
          tags$b("female-biased"), " within MINP. Take-home: the MINP response is small, ",
          tags$b("overwhelmingly female-biased"), ", mostly in poorly-annotated genes, and overlaps the ",
          "sex signature more than the nerve-injury program — only ", tags$b("Myh7"), " and ", tags$b("Slc15a2"),
          " are well-annotated, injury-negative and MINP-up (and both only reach significance in females)."),
      DT::DTOutput(ns("minp_tbl")),
      downloadButton(ns("minp_dl"), "Download (CSV)", class="btn-sm btn-outline-secondary mt-2"),
      tags$small(class="text-muted d-block mt-1",
        "Flags: well_annotated (named gene + description) · injury_negative (SNI |log2FC|<1 or n.s.) · ",
        "sex_chr (chrX/Y) · sex_DE (differs in Sham M vs F) · female_biased (sig in MINP females, not males).")
    )
  )
  # Gene lookup is a cross-species view, so it lives here rather than in the top navbar.
  if (isTRUE(show_gene_card)) panels <- c(panels, list(tab_gene_card_ui("gene_card", nested = TRUE)))

  nav_panel(
    title = "Cross-species meta-analysis", icon = icon("braille"),
    div(class = "container-fluid py-3",
      h3("Cross-species meta-analysis — mouse · rat · macaque"),
      do.call(navset_tab, panels)
    )
  )
}

tab_cross_pca_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    # ===================== PCA =====================
    cp <- load_rds("cross_pca")
    LABS <- c(unified="Unified condition", species="Species", condition="Condition",
              sex="Sex", level="DRG level", group="Group")
    # Explicit palette per colour-by option so every level is visually distinct and stable.
    # Anything without an entry (group, 24 levels) falls back to hue_pal.
    PCA_PAL <- list(
      unified   = c("Control / baseline"="#3477eb", "Nerve injury / pain model"="#E74C3C",
                    "Other (skin lesion)"="#F1C40F"),
      species   = c(Mouse="#E67E22", Rat="#16A085", Macaque="#8E44AD"),
      sex       = c(Female=unname(SEX_COLORS[["F"]]), Male=unname(SEX_COLORS[["M"]])),
      condition = c(Sham="#000000", MINP="#E67E22", SNI="#C0392B", "Naïve/control"="#3477eb",
                    SL="#F1C40F", contra="#16A085", ipsi="#8E44AD"),
      level     = c("Lower DRG (L4/5)"="#B2182B", "Upper DRG"="#2166AC",
                    "n.a. (mouse/monkey)"="#95A5A6"))
    # display values for a given encoding key
    pca_field <- function(d, key) switch(key %||% "condition",
      unified   = as.character(d$unified),
      species   = unname(CPCA_SP_LAB[as.character(d$species)]),
      sex       = unname(c(F="Female", M="Male")[as.character(d$sex)]),
      condition = as.character(d$condition_disp),
      level     = as.character(d$level_disp),
      group     = as.character(d$group_disp),
      as.character(d$condition_disp))
    # palette order first (so Female/Male, Lower/Upper read naturally), then any extras
    pca_levels <- function(v, key) {
      u <- unique(v); p <- names(PCA_PAL[[key]])
      c(p[p %in% u], sort(setdiff(u, p)))
    }

    output$pca <- plotly::renderPlotly({
      req(cp)
      pv <- if (identical(input$mode,"raw")) cp$raw else cp$aligned
      d <- pv$coords; ve <- pv$var_explained
      xc <- input$pcx %||% "PC1"; yc <- input$pcy %||% "PC2"; req(xc %in% names(d), yc %in% names(d))
      cby <- input$color %||% "condition"; sby <- input$shape %||% "species"
      d$x <- d[[xc]]; d$y <- d[[yc]]

      cv  <- pca_field(d, cby); lv <- pca_levels(cv, cby)
      pal <- PCA_PAL[[cby]]
      if (is.null(pal) || !all(lv %in% names(pal)))
        pal <- stats::setNames(scales::hue_pal()(length(lv)), lv)
      d$pcol <- unname(pal[cv])

      sv  <- pca_field(d, sby); slv <- pca_levels(sv, sby)
      sym <- stats::setNames(CPCA_SHAPES[(seq_along(slv) - 1) %% length(CPCA_SHAPES) + 1], slv)

      d$hover <- sprintf("<b>%s</b><br>%s · %s · %s%s", d$sample, CPCA_SP_LAB[as.character(d$species)],
                         d$condition_disp, ifelse(d$sex=="F", "female", "male"),
                         ifelse(d$level=="n.a.", "", paste0(" · ", d$level_disp)))
      p <- plotly::plot_ly()
      for (s in slv) { dd <- d[sv == s, ]; if (!nrow(dd)) next
        p <- plotly::add_trace(p, data=dd, x=~x, y=~y, type="scatter", mode="markers",
              marker=list(size=11, symbol=unname(sym[s]), color=dd$pcol, line=list(width=.5, color="#fff")),
              text=~hover, hoverinfo="text", showlegend=FALSE) }
      for (l in lv)
        p <- plotly::add_trace(p, x=c(NA), y=c(NA), type="scatter", mode="markers",
              name=l, legendgroup="colour", legendgrouptitle=list(text=paste0("Colour — ", LABS[[cby]] %||% cby)),
              marker=list(size=11, symbol="circle", color=unname(pal[l])), showlegend=TRUE, hoverinfo="skip")
      if (!identical(cby, sby)) for (s in slv)
        p <- plotly::add_trace(p, x=c(NA), y=c(NA), type="scatter", mode="markers",
              name=s, legendgroup="shape", legendgrouptitle=list(text=paste0("Shape — ", LABS[[sby]] %||% sby)),
              marker=list(size=11, symbol=unname(sym[s]), color="#666"), showlegend=TRUE, hoverinfo="skip")
      xi <- as.integer(sub("PC","",xc)); yi <- as.integer(sub("PC","",yc))
      out <- plotly::layout(p, title=sprintf("Cross-species PCA — %s view (%s genes)",
                       input$mode, format(pv$ngenes, big.mark=",")),
        xaxis=list(title=sprintf("%s (%.1f%%)", xc, ve[xi])),
        yaxis=list(title=sprintf("%s (%.1f%%)", yc, ve[yi])),
        legend=list(itemsizing="constant")) |> ply_pub()
      suppressWarnings(plotly::plotly_build(out))   # legend-proxy NA points -> silence "Ignoring N observations"
    })

    # ===================== Cell-type composition =====================
    md <- load_rds("deconv"); rd <- load_rds("rat_deconv"); kd <- load_rds("monkey_deconv")
    DEC <- list(Mouse=md, Rat=rd, Monkey=kd)

    # meta$group aligned to the rows of prop_wide, for one species
    grp_of <- function(sp) {
      d <- DEC[[sp]]; if (is.null(d)) return(NULL)
      mta <- d$meta
      sid <- if (!is.null(mta$sample)) mta$sample else if (!is.null(mta$sample_id)) mta$sample_id else rownames(mta)
      as.character(mta$group)[match(rownames(d$prop_wide), sid)]
    }
    # mean proportion vector over a set of meta$group values
    grp_mean <- function(sp, groups) {
      d <- DEC[[sp]]; g <- grp_of(sp); if (is.null(d) || is.null(g)) return(NULL)
      i <- which(g %in% groups); if (!length(i)) return(NULL)
      colMeans(d$prop_wide[i, , drop=FALSE])
    }
    # Δ (injury − control) + padj for one INJ_* entry
    inj_block <- function(s) {
      d <- DEC[[s$sp]]; g <- grp_of(s$sp); if (is.null(d) || is.null(g)) return(NULL)
      ii <- which(g %in% s$inj); ci <- which(g %in% s$ctl)
      if (!length(ii) || !length(ci)) return(NULL)
      pw <- d$prop_wide
      mi <- colMeans(pw[ii, , drop=FALSE]); mc <- colMeans(pw[ci, , drop=FALSE])
      padj <- pnom <- rep(NA_real_, ncol(pw))
      if (!is.na(s$key)) {
        st <- d$stats[d$stats$contrast == s$key, ]
        if (nrow(st)) { m <- match(colnames(pw), st$cellType)
          padj <- st$padj_wilcox[m]; pnom <- st$p_wilcox[m] }
      }
      list(cellType = colnames(pw), delta = mi - mc, padj = padj, pnom = pnom,
           detected = (mi + mc) > 1e-6, n_inj = length(ii), n_ctl = length(ci))
    }
    inj_mat <- function(spec, cts) {
      cols <- vapply(spec, function(s) s$lab, character(1))
      Z <- matrix(NA_real_, length(cts), length(cols), dimnames=list(cts, cols))
      P <- Q <- Z; N <- D <- matrix("", length(cts), length(cols))
      for (j in seq_along(spec)) {
        b <- inj_block(spec[[j]]); if (is.null(b)) next
        m <- match(cts, b$cellType)
        z <- b$delta[m]
        z[!b$detected[m]] <- NA_real_   # not detected in this species -> blank, not "no change"
        Z[, j] <- z; P[, j] <- b$padj[m]; Q[, j] <- b$pnom[m]
        N[, j] <- sprintf("n = %d injury vs %d control", b$n_inj, b$n_ctl)
        D[, j] <- ifelse(b$detected[m], "", "not detected in this species")
      }
      list(Z=Z, P=P, Q=Q, N=N, D=D, cols=cols)
    }

    # cell types ordered by mouse Sham baseline; cLTMR2 is all-zero in every species
    ct_order <- reactive({
      b <- grp_mean("Mouse", c("Sham_M","Sham_F")); req(b)
      b <- b[setdiff(names(b), "cLTMR2")]
      names(b)[order(-b)]
    })
    dmats <- reactive({
      cts <- ct_order(); req(length(cts) > 0)
      pooled <- inj_mat(INJ_POOLED, cts); sex <- inj_mat(INJ_SEX, cts)
      if ("rows" %in% (input$dclust %||% character(0))) {           # cluster cell types on the pooled Δ
        ro <- cluster_order(pooled$Z, do_row=TRUE)$row; cts <- cts[ro]
        pooled <- inj_mat(INJ_POOLED, cts); sex <- inj_mat(INJ_SEX, cts)
      }
      list(cts=cts, pooled=pooled, sex=sex)
    })

    delta_heat <- function(mm, cts, title, zm, showscale = TRUE, clust_cols = FALSE) {
      Z <- mm$Z; P <- mm$P; Q <- mm$Q; N <- mm$N; D <- mm$D; cols <- mm$cols
      if (isTRUE(clust_cols) && ncol(Z) > 2) {                     # cluster injury columns
        co <- cluster_order(Z, do_col=TRUE)$col
        Z<-Z[,co,drop=FALSE]; P<-P[,co,drop=FALSE]; Q<-Q[,co,drop=FALSE]
        N<-N[,co,drop=FALSE]; D<-D[,co,drop=FALSE]; cols<-cols[co]
      }
      sig <- function(v) ifelse(!is.na(v) & v < 0.05, TRUE, FALSE)
      txt <- matrix(sprintf("<b>%s</b><br>%s<br>%s<br>%s%s",
                     rep(gsub("<br>", " ", cols), each=length(cts)),
                     rep(ct_full(cts), times=length(cols)),
                     ifelse(nzchar(as.character(D)), as.character(D),
                            sprintf("Δ %+.3f · %s", as.numeric(Z), as.character(N))),
                     ifelse(is.na(as.numeric(P)), "no test",
                            sprintf("adj.p %.3g", as.numeric(P))),
                     ifelse(sig(as.numeric(P)), " ✱", ifelse(sig(as.numeric(Q)), " · (nominal p<0.05)", ""))),
                   nrow=length(cts))
      p <- plotly::plot_ly(x=cols, y=cts, z=Z, type="heatmap", zmid=0, zmin=-zm, zmax=zm,
             colorscale=CPCA_DIVERGE, text=txt, hoverinfo="text", showscale=showscale,
             colorbar=list(title="Δ prop\n(injury−control)"))
      ann <- list()
      for (i in seq_along(cts)) for (j in seq_along(cols)) {
        if (sig(P[i,j]))
          ann[[length(ann)+1]] <- list(x=cols[j], y=cts[i], text="✱", showarrow=FALSE,
                                       font=list(color="black", size=15))
        else if (sig(Q[i,j]))
          ann[[length(ann)+1]] <- list(x=cols[j], y=cts[i], text="·", showarrow=FALSE,
                                       font=list(color="#444", size=22))
      }
      shp <- list()
      if (!isTRUE(clust_cols)) {                                   # species separators only in fixed order
        sp_of <- sub("^(Mouse|Rat|Macaque).*", "\\1", gsub("<br>", " ", cols))
        brk <- utils::head(cumsum(rle(sp_of)$lengths), -1)
        shp <- lapply(brk, function(b) list(type="line", x0=b-0.5, x1=b-0.5, yref="paper",
                                            y0=0, y1=1, line=list(color="black", width=2)))
      }
      plotly::layout(p, title=title,
        xaxis=list(title="", categoryorder="array", categoryarray=cols, tickangle=-30),  # labels at bottom
        yaxis=list(title="", autorange="reversed"),
        annotations=ann, shapes=shp) |> ply_pub()
    }

    output$dec_wrap <- renderUI({
      if (is.null(md) || is.null(rd) || is.null(kd)) return(coming_soon("Deconvolution pending."))
      if (identical(input$dview, "delta"))
        tagList(
          h6(class="mt-2", "Both sexes pooled"),
          plotly::plotlyOutput(ns("dec_pooled"), height="430px"),
          tags$hr(),
          h6("Split by sex — male and female adjacent within each injury"),
          plotly::plotlyOutput(ns("dec_sex"), height="430px"),
          sni_note("SNI_vs_Sham"))
      else plotly::plotlyOutput(ns("dec_one"), height="560px")
    })

    zmax <- reactive({ x <- dmats(); max(abs(c(x$pooled$Z, x$sex$Z)), na.rm=TRUE) })
    output$dec_pooled <- plotly::renderPlotly({
      x <- dmats(); cc <- "cols" %in% (input$dclust %||% character(0))
      delta_heat(x$pooled, x$cts, "Injury − control cell-type proportion (both sexes)", zmax(), clust_cols=cc)
    })
    output$dec_sex <- plotly::renderPlotly({
      x <- dmats(); cc <- "cols" %in% (input$dclust %||% character(0))
      delta_heat(x$sex, x$cts, "Same injuries, split by sex (same colour scale as above)", zmax(), clust_cols=cc)
    })

    output$dec_one <- plotly::renderPlotly({
      cts <- ct_order(); req(length(cts) > 0)
      spec <- if (identical(input$dview, "base")) Filter(function(g) isTRUE(g$ctl), COMP_GROUPS) else COMP_GROUPS
      cols <- vapply(spec, function(g) g$lab, character(1))
      M <- vapply(spec, function(g) { v <- grp_mean(g$sp, g$grp)
                    if (is.null(v)) rep(NA_real_, length(cts)) else unname(v[cts]) }, numeric(length(cts)))
      pal <- stats::setNames(scales::hue_pal()(length(cts)), cts)
      p <- plotly::plot_ly()
      for (i in seq_along(cts))
        p <- plotly::add_bars(p, x=cols, y=M[i, ], name=cts[i],
               marker=list(color=pal[[cts[i]]], line=list(color="white", width=.4)),
               text=sprintf("<b>%s</b><br>%s<br>%.3f", ct_full(cts[i]), gsub("<br>", " ", cols), M[i, ]),
               hoverinfo="text")
      plotly::layout(p, barmode="stack",
        title=if (identical(input$dview, "base")) "Baseline (control) cell-type composition"
              else "Cell-type composition by condition",
        xaxis=list(title="", tickangle=-20, categoryorder="array", categoryarray=cols),
        yaxis=list(title="Estimated proportion", range=c(0, 1.001)),
        legend=list(title=list(text="Cell type"))) |> ply_pub()
    })

    # ===================== Conserved injury markers =====================
    # 47 genes x 93 samples of VST expression, z-scored per gene WITHIN each species so the
    # three blocks are on a common scale (see the caption under the plot).
    mk_z <- reactive({
      mk <- load_rds("conserved_injury_markers"); req(mk)
      gsym <- as.character(mk$Gene)
      ag  <- ANNOT[ANNOT$level == "gene", ]
      mid <- as.character(ag$feature_id[match(gsym, as.character(ag$symbol))])

      blocks <- lapply(names(MK_BLOCKS), function(sp) {
        o <- MK_BLOCKS[[sp]]; E <- load_rds(o$rds); if (is.null(E)) return(NULL)
        vst <- if (o$nest) E$gene$vst else E$vst
        mta <- if (o$nest) E$gene$meta else E$meta
        ids <- if (is.null(o$orth)) mid else {
          orth <- load_rds(o$orth); if (is.null(orth)) return(NULL)
          idc <- setdiff(grep("_gene_id$", names(orth), value=TRUE), "mouse_gene_id")[1]
          as.character(orth[[idc]][match(mid, orth$mouse_gene_id)])
        }
        keep <- !is.na(ids) & ids %in% rownames(vst)
        sel  <- Filter(function(x) length(x$i) > 0,
                       lapply(o$grp, function(g) list(i=which(g$f(mta)), lab=g$lab)))
        if (!length(sel)) return(NULL)
        ord <- unlist(lapply(sel, `[[`, "i"))
        M <- matrix(NA_real_, length(gsym), length(ord))
        if (any(keep)) M[keep, ] <- vst[ids[keep], ord, drop=FALSE]
        Z <- t(scale(t(M)))                       # per gene, WITHIN this species
        Z[!is.finite(Z)] <- NA_real_
        list(sp=sp, Z=Z, raw=M,
             grp=rep(vapply(sel, `[[`, character(1), "lab"),
                     vapply(sel, function(s) length(s$i), integer(1))),
             sample=as.character(mta$sample_id)[ord],
             sex=as.character(mta$sex)[ord])
      })
      blocks <- Filter(Negate(is.null), blocks); req(length(blocks) > 0)
      list(genes = gsym,
           Z    = do.call(cbind, lapply(blocks, `[[`, "Z")),
           raw  = do.call(cbind, lapply(blocks, `[[`, "raw")),
           grp  = unlist(lapply(blocks, `[[`, "grp")),
           sample = unlist(lapply(blocks, `[[`, "sample")),
           sex  = unlist(lapply(blocks, `[[`, "sex")),
           species = rep(vapply(blocks, `[[`, character(1), "sp"),
                         vapply(blocks, function(b) ncol(b$Z), integer(1))))
    })

    output$mk_heat <- plotly::renderPlotly({
      z <- mk_z(); Z <- z$Z; raw <- z$raw; g <- z$genes
      grp <- z$grp; sample <- z$sample; sex <- z$sex; species <- z$species
      cl <- input$mclust %||% character(0)
      if ("rows" %in% cl) {                                         # cluster the 47 genes
        ro <- cluster_order(Z, do_row=TRUE)$row
        g <- g[ro]; Z <- Z[ro, , drop=FALSE]; raw <- raw[ro, , drop=FALSE]
      }
      if ("cols" %in% cl) {                                         # cluster samples WITHIN each group block
        ord <- integer(0); rgg <- rle(grp); st <- 1L
        for (k in seq_along(rgg$lengths)) { len <- rgg$lengths[k]; idx <- st:(st+len-1L)
          oo <- if (len > 2) cluster_order(t(Z[, idx, drop=FALSE]), do_row=TRUE)$row else seq_len(len)
          ord <- c(ord, idx[oo]); st <- st + len }
        Z<-Z[,ord,drop=FALSE]; raw<-raw[,ord,drop=FALSE]; grp<-grp[ord]; sample<-sample[ord]; sex<-sex[ord]; species<-species[ord]
      }
      n <- ncol(Z); ng <- length(g)
      hov <- matrix(sprintf("<b>%s</b><br>%s<br>%s · %s<br>z %+.2f · VST %.2f",
                     rep(g, times=n), rep(sample, each=ng),
                     rep(gsub("<br>", " ", grp), each=ng), rep(sex, each=ng),
                     as.numeric(Z), as.numeric(raw)), nrow=ng)
      rg  <- rle(grp); brk <- cumsum(rg$lengths)
      mid <- brk - rg$lengths / 2 + 0.5
      spb <- utils::head(cumsum(rle(species)$lengths), -1)
      shp <- c(
        lapply(spb, function(b) list(type="line", x0=b+0.5, x1=b+0.5, yref="paper", y0=0, y1=1,
                                     line=list(color="black", width=2.5))),
        lapply(setdiff(utils::head(brk, -1), spb), function(b)
          list(type="line", x0=b+0.5, x1=b+0.5, yref="paper", y0=0, y1=1,
               line=list(color="grey45", width=1, dash="dot"))))
      p <- plotly::plot_ly(x=seq_len(n), y=g, z=Z, type="heatmap",
             zmid=0, zmin=-2, zmax=2, colorscale=CPCA_DIVERGE,
             text=hov, hoverinfo="text", colorbar=list(title="z-score\n(within\nspecies)"))
      plotly::layout(p, title="47 conserved injury markers — VST z-scored within each species",
        xaxis=list(title="", tickmode="array", tickvals=mid,
                   ticktext=gsub("<br>", " ", rg$values), tickangle=-35),
        yaxis=list(title="", autorange="reversed", categoryorder="array", categoryarray=rev(g),
                   tickfont=list(size=9)),
        shapes=shp) |> ply_pub()
    })

    output$markers <- DT::renderDT({
      mk <- load_rds("conserved_injury_markers"); req(mk)
      DT::datatable(mk, rownames=FALSE, filter="top", extensions="Buttons",
        options=list(pageLength=15, scrollX=TRUE, dom="Bfrtip", order=list(list(9,"desc")),
          buttons=list(list(extend="excel", text="⤓ Excel", filename="conserved_injury_markers_SNIpos_MINPneg")))) |>
        DT::formatStyle(c("Mouse SNI log2FC","Rat SNI log2FC","Macaque ipsi log2FC","Mean SNI log2FC"),
                        color=DT::styleInterval(0, c("#B0B8BF","#c0392b"))) |>
        DT::formatStyle("MINP log2FC", fontWeight="bold", color="#333")
    })
    output$mk_dl <- downloadHandler(
      filename=function() "conserved_injury_markers_SNIpos_MINPneg.csv",
      content=function(f){ mk <- load_rds("conserved_injury_markers"); utils::write.csv(mk, f, row.names=FALSE) })

    # ===================== MINP-only genes =====================
    output$minp_tbl <- DT::renderDT({
      d <- load_rds("minp_specific"); req(d)
      yn <- function(x) ifelse(isTRUE(x) | x %in% TRUE, "✓", "")
      out <- data.frame(
        Gene=d$Gene, chr=d$chr, Description=d$description,
        `MINP log2FC`=round(d$MINP_l2fc,2), `MINP adj.p`=signif(d$MINP_padj,2),
        `MINP-F log2FC`=round(d$MINP_F_l2fc,2), `MINP-F adj.p`=signif(d$MINP_F_padj,2),
        `MINP-M log2FC`=round(d$MINP_M_l2fc,2), `MINP-M adj.p`=signif(d$MINP_M_padj,2),
        `SNI log2FC`=round(d$SNI_l2fc,2), `SNI adj.p`=signif(d$SNI_padj,2),
        `Sex(M-F) log2FC`=round(d$Sex_l2fc,2), `Sex adj.p`=signif(d$Sex_padj,2),
        `well-annot`=yn(d$well_annotated), `injury-neg`=yn(d$injury_negative),
        `female-biased`=yn(d$female_biased), `sex-DE`=yn(d$sex_DE), `sex-chr`=yn(d$sex_chr),
        check.names=FALSE, stringsAsFactors=FALSE)
      DT::datatable(out, rownames=FALSE, filter="top", extensions="Buttons",
        options=list(pageLength=20, scrollX=TRUE, dom="Bfrtip",
          buttons=list(list(extend="excel", text="⤓ Excel", filename="MINP_specific_genes")))) |>
        DT::formatStyle("MINP-F log2FC", fontWeight="bold")
    })
    output$minp_dl <- downloadHandler(
      filename=function() "MINP_specific_genes.csv",
      content=function(f){ d <- load_rds("minp_specific"); utils::write.csv(d, f, row.names=FALSE) })
  })
}
