# build_cross_pca.R — one PCA over mouse + rat + monkey DRG samples in shared 1:1
# ortholog space (mouse gene-id). Two variants:
#   raw     : combined VST (PC1 will separate species)
#   aligned : per-species per-gene z-score (removes species baseline -> lets condition
#             structure align across species)
# Output app/data/cross_pca.rds = list(raw=..., aligned=...), each list(coords, var_explained, ntop, ngenes).
setwd("/hpf/projects/msalter/sghazis/rnaseq_TUY35595/app")
suppressMessages(library(matrixStats))
me <- readRDS("data/expr.rds")$gene          # mouse: vst (ENSMUSG) + meta
re <- readRDS("data/rat_expr.rds")           # rat:   vst (ENSRNOG) + meta
ke <- readRDS("data/monkey_expr.rds")        # monkey:vst (ENSMFAG) + meta
ro <- readRDS("data/ortholog_rat_mouse.rds")
mo <- readRDS("data/ortholog_monkey_mouse.rds")

# map a species VST to mouse-gene-id rows (1:1); collapse dups by rowMeans
to_mouse <- function(vst, ids2mouse) {
  mm <- ids2mouse(rownames(vst)); keep <- !is.na(mm)
  v <- vst[keep, , drop=FALSE]; mm <- mm[keep]
  if (any(duplicated(mm))) v <- rowsum(v, mm) / as.numeric(table(mm)[match(rownames(rowsum(v,mm)), names(table(mm)))])
  else rownames(v) <- mm
  v
}
rat_m <- to_mouse(re$vst, function(x) as.character(ro$mouse_gene_id[match(x, ro$rat_gene_id)]))
mky_m <- to_mouse(ke$vst, function(x) as.character(mo$mouse_gene_id[match(x, mo$monkey_gene_id)]))
mou_m <- me$vst
common <- Reduce(intersect, list(rownames(mou_m), rownames(rat_m), rownames(mky_m)))
cat("shared 1:1 orthologs across 3 species:", length(common), "\n")

blocks <- list(Mouse=mou_m[common,,drop=FALSE], Rat=rat_m[common,,drop=FALSE], Monkey=mky_m[common,,drop=FALSE])
comb <- do.call(cbind, blocks)
species <- rep(names(blocks), vapply(blocks, ncol, integer(1)))

# unified meta across species (mapping confirmed with the user)
#   Nerve injury / pain model : mouse SNI + MINP(NP), rat SNI, monkey ipsi
#   Control / baseline        : mouse Sham, rat Naive, monkey control + contra
#   Other (skin lesion)       : rat SL
uni_cond <- function(sp, cond) {
  inj  <- (sp=="Mouse" & cond %in% c("SNI","NP")) | (sp=="Rat" & cond=="SNI") | (sp=="Monkey" & cond=="ipsi")
  base <- (sp=="Mouse" & cond=="Sham") | (sp=="Rat" & cond=="Naive") | (sp=="Monkey" & cond %in% c("control","contra"))
  out <- rep("Other (skin lesion)", length(cond)); out[base] <- "Control / baseline"; out[inj] <- "Nerve injury / pain model"; out
}
dispf  <- function(x) gsub("NP", "MINP", as.character(x))
lv_map <- c(DRGU="Upper DRG", DRGL="Lower DRG (L4/5)", "n.a."="n.a. (mouse/monkey)")
meta <- rbind(
  data.frame(sample=colnames(mou_m), species="Mouse",  condition=me$meta$condition, sex=me$meta$sex, group=me$meta$group, level="n.a.", stringsAsFactors=FALSE),
  data.frame(sample=colnames(re$vst), species="Rat",    condition=re$meta$condition, sex=re$meta$sex, group=re$meta$group, level=re$meta$level, stringsAsFactors=FALSE),
  data.frame(sample=colnames(ke$vst), species="Monkey", condition=ke$meta$condition, sex=ke$meta$sex, group=ke$meta$group, level="n.a.", stringsAsFactors=FALSE))
meta$unified <- uni_cond(meta$species, meta$condition)
meta$condition_disp <- dispf(meta$condition)
meta$condition_disp[meta$condition %in% c("Naive","control")] <- "Naïve/control"
meta$group_disp <- dispf(meta$group)
meta$level_disp <- unname(lv_map[meta$level])
stopifnot(nrow(meta)==ncol(comb))

pca_of <- function(m, ntop=2000){
  v <- rowVars(m); sel <- order(v, decreasing=TRUE)[seq_len(min(ntop, sum(v>0)))]
  pc <- prcomp(t(m[sel,]), center=TRUE, scale.=FALSE)
  ve <- round(100*pc$sdev^2/sum(pc$sdev^2),1)
  cd <- as.data.frame(pc$x[,1:min(10,ncol(pc$x)),drop=FALSE]); cd <- cbind(cd, meta)
  list(coords=cd, var_explained=ve, ntop=length(sel), ngenes=nrow(m))
}
# aligned: z-score each gene WITHIN each species block, then recombine
zblock <- function(mat){ t(scale(t(mat))) }
comb_z <- do.call(cbind, lapply(blocks, zblock))
comb_z[is.na(comb_z)] <- 0

cross_pca <- list(raw = pca_of(comb), aligned = pca_of(comb_z))
saveRDS(cross_pca, "data/cross_pca.rds")
cat("cross_pca.rds: samples", ncol(comb), "| raw PC1-2", paste0(cross_pca$raw$var_explained[1:2],"%",collapse="/"),
    "| aligned PC1-2", paste0(cross_pca$aligned$var_explained[1:2],"%",collapse="/"), "\n")
cat("size MB:", round(as.numeric(object.size(cross_pca))/1e6,2), "\n")
