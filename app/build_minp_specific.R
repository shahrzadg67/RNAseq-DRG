# MINP-specific genes: DE in MINP but NOT in nerve injury (SNI). Flag annotation
# quality and whether they are sex-linked / sex-different. -> minp_specific.rds
setwd("/hpf/projects/msalter/sghazis/rnaseq_TUY35595/app")
mou<-readRDS("data/de_long.rds"); mou<-mou[mou$engine=="DESeq2"&mou$level=="gene",]
ann<-readRDS("data/gene_annot.rds"); ann<-ann[ann$level=="gene",]
gM<-function(c){d<-mou[mou$contrast==c,];data.frame(fid=as.character(d$feature_id),sym=as.character(d$symbol),l=d$log2FC,p=d$padj,stringsAsFactors=FALSE)}
minp<-gM("NP_vs_Sham"); minp_f<-gM("NP_F_vs_Sham_F"); minp_m<-gM("NP_M_vs_Sham_M")
sni <-gM("SNI_vs_Sham"); sex<-gM("Sham_M_vs_Sham_F")
key<-function(d) setNames(seq_len(nrow(d)),d$fid)
# MINP-associated = sig in pooled OR either sex MINP contrast
sigset<-function(d) d$fid[!is.na(d$p)&d$p<0.05]
minp_ids<-unique(c(sigset(minp),sigset(minp_f),sigset(minp_m)))
cat("MINP-associated features (any MINP contrast padj<.05):",length(minp_ids),"\n")
# build table
row<-function(fid){
  gm<-minp[match(fid,minp$fid),]; gf<-minp_f[match(fid,minp_f$fid),]; gmm<-minp_m[match(fid,minp_m$fid),]
  gs<-sni[match(fid,sni$fid),]; gx<-sex[match(fid,sex$fid),]; a<-ann[match(fid,ann$feature_id),]
  sym<-gm$sym; if(is.na(sym)||sym=="") sym<-fid
  data.frame(Gene=sym, feature_id=fid, chr=a$chr, description=a$description,
    MINP_l2fc=gm$l, MINP_padj=gm$p, MINP_F_l2fc=gf$l, MINP_F_padj=gf$p, MINP_M_l2fc=gmm$l, MINP_M_padj=gmm$p,
    SNI_l2fc=gs$l, SNI_padj=gs$p, Sex_l2fc=gx$l, Sex_padj=gx$p, stringsAsFactors=FALSE)
}
tab<-do.call(rbind,lapply(minp_ids,row))
# annotation quality
predicted<-function(s) is.na(s)|grepl("Rik$|^Gm[0-9]|^LOC|^A[0-9]{3}[0-9]+|^BC[0-9]+|^RP[0-9]|orf",s)
tab$well_annotated<- !predicted(tab$Gene) & !is.na(tab$description)
# injury-negative = SNI not up (not sig-up); flat/absent in injury
tab$injury_negative<- is.na(tab$SNI_padj)|tab$SNI_padj>=0.05|abs(tab$SNI_l2fc)<1
# direction in MINP
tab$MINP_up<- !is.na(tab$MINP_l2fc)&tab$MINP_l2fc>0
# sex-linked / sex-different
tab$sex_chr<- !is.na(tab$chr)&tab$chr %in% c("X","Y")
tab$sex_DE <- !is.na(tab$Sex_padj)&tab$Sex_padj<0.05
tab$female_biased<- !is.na(tab$MINP_F_padj)&tab$MINP_F_padj<0.05 & (is.na(tab$MINP_M_padj)|tab$MINP_M_padj>=0.05)
tab<-tab[order(tab$MINP_padj),]
saveRDS(tab,"data/minp_specific.rds")

cat("\n=== MINP-associated genes: annotation & sex ===\n")
cat(sprintf("total=%d | well-annotated=%d | injury-negative=%d | well-annot & injury-neg & MINP-up=%d\n",
  nrow(tab),sum(tab$well_annotated),sum(tab$injury_negative),sum(tab$well_annotated&tab$injury_negative&tab$MINP_up)))
cat(sprintf("on sex chromosome=%d | also sex-DE(Sham M vs F)=%d | female-biased MINP=%d\n",
  sum(tab$sex_chr),sum(tab$sex_DE),sum(tab$female_biased)))
cat("\n=== the well-annotated, injury-negative, MINP-up set ===\n")
core<-tab[tab$well_annotated&tab$injury_negative&tab$MINP_up,]
for(i in seq_len(nrow(core))) cat(sprintf("  %-12s chr%-3s MINP %+5.2f(p=%.1g)  SNI %+5.2f(p=%.1g)  Sex %+5.2f(p=%.1g)  %s%s\n",
  core$Gene[i],core$chr[i],core$MINP_l2fc[i],core$MINP_padj[i],core$SNI_l2fc[i],core$SNI_padj[i],
  core$Sex_l2fc[i],core$Sex_padj[i], ifelse(core$sex_chr[i],"[sexchr]",""), ifelse(core$female_biased[i],"[F-biased]","")))
cat("\n(", nrow(core), "core genes)\n")
