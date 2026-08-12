# build_injury_markers.R — the conserved "SNI-positive / MINP-negative" panel:
# genes UP & significant in nerve injury in >=2 species and FLAT in mouse MINP.
# Output app/data/conserved_injury_markers.rds (data.frame for a downloadable table).
setwd("/hpf/projects/msalter/sghazis/rnaseq_TUY35595/app")
mou<-readRDS("data/de_long.rds"); mou<-mou[mou$engine=="DESeq2"&mou$level=="gene",]
rat<-readRDS("data/rat_de_long.rds"); mon<-readRDS("data/monkey_de_long.rds")
ro<-readRDS("data/ortholog_rat_mouse.rds"); mo<-readRDS("data/ortholog_monkey_mouse.rds")
gM<-function(c){d<-mou[mou$contrast==c,];data.frame(sym=as.character(d$symbol),l=d$log2FC,p=d$padj)}
gR<-function(c){d<-rat[rat$contrast==c,];s<-as.character(ro$mouse_symbol[match(as.character(d$feature_id),ro$rat_gene_id)]);k<-!is.na(s);data.frame(sym=s[k],l=d$log2FC[k],p=d$padj[k])}
gK<-function(c){d<-mon[mon$contrast==c,];s<-as.character(mo$mouse_symbol[match(as.character(d$feature_id),mo$monkey_gene_id)]);k<-!is.na(s);data.frame(sym=s[k],l=d$log2FC[k],p=d$padj[k])}
m<-Reduce(function(a,b)merge(a,b,by="sym"),list(
  setNames(gM("SNI_vs_Sham"),c("sym","sni_m_l","sni_m_p")),
  setNames(gM("NP_vs_Sham"), c("sym","minp_l","minp_p")),
  setNames(gR("SNI_vs_Naive_L"),c("sym","sni_r_l","sni_r_p")),
  setNames(gK("ipsi_vs_control"),c("sym","sni_k_l","sni_k_p"))))
sigup<-function(l,p) !is.na(p)&p<0.05&l>1
n_inj<-sigup(m$sni_m_l,m$sni_m_p)+sigup(m$sni_r_l,m$sni_r_p)+sigup(m$sni_k_l,m$sni_k_p)
flat_minp<-(is.na(m$minp_p)|m$minp_p>0.5)&abs(m$minp_l)<0.5
cand<-m[n_inj>=2 & flat_minp & !is.na(m$sym) & m$sym!="",]
cand$n_sig<-(sigup(cand$sni_m_l,cand$sni_m_p)+sigup(cand$sni_r_l,cand$sni_r_p)+sigup(cand$sni_k_l,cand$sni_k_p))
cand$mean_sni<-rowMeans(cand[,c("sni_m_l","sni_r_l","sni_k_l")],na.rm=TRUE)
cand<-cand[order(-cand$mean_sni),]
out<-data.frame(
  Gene=cand$sym,
  `Mouse SNI log2FC`=round(cand$sni_m_l,2), `Mouse SNI adj.p`=signif(cand$sni_m_p,2),
  `Rat SNI log2FC`=round(cand$sni_r_l,2),   `Rat SNI adj.p`=signif(cand$sni_r_p,2),
  `Macaque ipsi log2FC`=round(cand$sni_k_l,2), `Macaque ipsi adj.p`=signif(cand$sni_k_p,2),
  `MINP log2FC`=round(cand$minp_l,2), `MINP adj.p`=signif(cand$minp_p,2),
  `Mean SNI log2FC`=round(cand$mean_sni,2), `n species sig`=cand$n_sig,
  check.names=FALSE, stringsAsFactors=FALSE)
saveRDS(out,"data/conserved_injury_markers.rds")
cat("conserved_injury_markers.rds:",nrow(out),"genes\n"); print(utils::head(out,8))
