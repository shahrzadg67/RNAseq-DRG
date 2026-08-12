# =====================================================================
# TUY35595 RNA-seq report app  (shinylive / WebR static build)
# Reads precomputed artifacts from data/ and images from www/.
# Tab visibility is gated by data/build_config.json (written by publish.sh):
#   target = "staging" -> show every tab
#   target = "prod"    -> show only stages with released = true
# The Methods tab ALWAYS shows the full per-stage status, transparently.
# =====================================================================
library(shiny); library(bslib)
library(jsonlite)
have <- function(p) requireNamespace(p, quietly = TRUE)
if (have("plotly")) library(plotly)
if (have("DT")) library(DT)
if (have("ggplot2")) library(ggplot2)          # static volcano / venn / theme_pub
if (have("ggrepel")) library(ggrepel)           # volcano gene labels
# pheatmap, ggVennDiagram, UpSetR are referenced with :: in their tabs

source("R/helpers.R", local = TRUE)
source("R/tab_methods.R", local = TRUE)
source("R/tab_qc.R",      local = TRUE)
source("R/tab_pca.R",     local = TRUE)
source("R/tab_de.R",      local = TRUE)
source("R/tab_gsea.R",    local = TRUE)
source("R/tab_pca_noxy.R",  local = TRUE)
source("R/tab_de_noxy.R",   local = TRUE)
source("R/tab_gsea_noxy.R", local = TRUE)
source("R/tab_venn.R",      local = TRUE)
source("R/tab_heatmap.R",   local = TRUE)
source("R/tab_summary.R",   local = TRUE)
source("R/tab_deconv.R",    local = TRUE)
source("R/tab_xspecies.R",  local = TRUE)
source("R/tab_xspecies_monkey.R", local = TRUE)
source("R/tab_gene_card.R", local = TRUE)
source("R/tab_cross_pca.R", local = TRUE)
source("R/tab_sc_injury.R", local = TRUE)

CFG    <- read_build_config()          # target + stage status
STAGES <- CFG$stages

# decide which stage keys are visible in THIS build
visible <- function(key) {
  st <- STAGES[[key]]
  if (is.null(st)) return(FALSE)
  if (identical(CFG$target, "staging")) return(TRUE)
  isTRUE(st$released)
}

# assemble tabs (methods always present)
tabs <- list(tab_methods_ui("methods", CFG))
if (visible("summary")) tabs <- c(tabs, list(tab_summary_ui("summary")))
if (visible("qc"))   tabs <- c(tabs, list(tab_qc_ui("qc")))
# group the analysis tabs into two dropdown menus to declutter the navbar
unc <- list()
if (visible("pca"))  unc <- c(unc, list(tab_pca_ui("pca")))
if (visible("de"))   unc <- c(unc, list(tab_de_ui("de")))
if (visible("gsea")) unc <- c(unc, list(tab_gsea_ui("gsea")))
if (length(unc)) tabs <- c(tabs, list(do.call(nav_menu, c(list("Uncorrected"), unc))))
# Corrected group: PCA/DE/GSEA (batch + sex-chr) plus the downstream views that
# run on the corrected data (Venn, Heatmap) and the deconvolution.
cor <- list()
if (visible("pca_noxy"))  cor <- c(cor, list(tab_pca_noxy_ui("pca_noxy")))
if (visible("de_noxy"))   cor <- c(cor, list(tab_de_noxy_ui("de_noxy")))
if (visible("gsea_noxy")) cor <- c(cor, list(tab_gsea_noxy_ui("gsea_noxy")))
if (visible("venn"))    cor <- c(cor, list(tab_venn_ui("venn")))
if (visible("heatmap")) cor <- c(cor, list(tab_heatmap_ui("heatmap")))
if (visible("deconvolution")) cor <- c(cor, list(tab_deconv_ui("deconvolution")))
if (length(cor)) tabs <- c(tabs, list(do.call(nav_menu, c(list("Corrected (batch + sex-chr)"), cor))))
if (visible("xspecies")) tabs <- c(tabs, list(tab_xspecies_ui("xspecies")))
if (visible("xspecies_monkey")) tabs <- c(tabs, list(tab_xspecies_monkey_ui("xspecies_monkey")))
# Gene lookup is a cross-species view, so it is nested inside the meta-analysis tab
# (tab_cross_pca_ui) rather than sitting in the top navbar. Its server stays registered below.
if (visible("cross_pca"))
  tabs <- c(tabs, list(tab_cross_pca_ui("cross_pca", show_gene_card = visible("gene_card"))))
if (visible("sc_injury")) tabs <- c(tabs, list(tab_sc_injury_ui("sc_injury")))

ui <- page_navbar(
  title = CFG$project$client_visible_title %||% "RNA-seq Analysis",
  theme = bs_theme(version = 5, preset = "flatly"),
  header = tags$head(tags$style(HTML(
    ".banner{background:#f4f6f8;border-left:4px solid #18bc9c;padding:.6rem 1rem;margin:.4rem 0;border-radius:4px}
     .status-chip{font-size:.75rem;padding:.1rem .5rem;border-radius:10px;color:#fff}
     .s-planned{background:#95a5a6}.s-in_progress{background:#f39c12}.s-complete{background:#18bc9c}.s-finalized{background:#27ae60}
     .app-id{display:flex;flex-direction:column;align-items:flex-end;justify-content:center;line-height:1.15;padding:.1rem .9rem;text-align:right}
     .app-id .who{font-weight:600;font-size:.9rem;white-space:nowrap}
     .app-id .meta{font-size:.7rem;opacity:.7;white-space:nowrap}"))),
  !!!tabs,
  nav_spacer(),
  nav_item(div(class = "app-id",
    tags$span(class = "who", "Shahrzad Ghazisaeidi"),
    tags$span(class = "meta", sprintf("nf-core/rnaseq v3.19.0 · %s", CFG$project$last_updated %||% ""))))
)

server <- function(input, output, session) {
  tab_methods_server("methods", CFG)
  if (visible("summary")) tab_summary_server("summary")
  if (visible("qc"))   tab_qc_server("qc")
  if (visible("pca"))  tab_pca_server("pca")
  if (visible("de"))   tab_de_server("de")
  if (visible("gsea")) tab_gsea_server("gsea")
  if (visible("pca_noxy"))  tab_pca_noxy_server("pca_noxy")
  if (visible("de_noxy"))   tab_de_noxy_server("de_noxy")
  if (visible("gsea_noxy")) tab_gsea_noxy_server("gsea_noxy")
  if (visible("venn"))    tab_venn_server("venn")
  if (visible("heatmap")) tab_heatmap_server("heatmap")
  if (visible("deconvolution")) tab_deconv_server("deconvolution")
  if (visible("xspecies")) tab_xspecies_server("xspecies")
  if (visible("xspecies_monkey")) tab_xspecies_monkey_server("xspecies_monkey")
  if (visible("gene_card")) tab_gene_card_server("gene_card")
  if (visible("cross_pca")) tab_cross_pca_server("cross_pca")
  if (visible("sc_injury")) tab_sc_injury_server("sc_injury")
}

shinyApp(ui, server)
