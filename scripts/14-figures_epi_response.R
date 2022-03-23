source("scripts/utils/new_utils.R")
out<-"outputs/figures_epi_response"
dir.create(out)

#FIGURE 1 : Epigenetic Memory ####

#1A : volcanoplots ####
#by cohort
res_coh<-fread("outputs/01-lga_vs_ctrl_limma_DMCs_analysis/res_limma_cohorts.tsv.gz")
res_coh<-res_coh[compa=="C.L"]
table(res_coh[adj.P.Val<=0.1,.(batch)])
ggplot(res_coh)+
  geom_point(aes(x=logFC,y=-log10(P.Value),col=P.Value<0.001&abs(logFC)>25),size=0.1)+
  facet_wrap("batch")+
  scale_color_manual(values = c("grey","red"))+
  theme_minimal()
ggsave(fp(out,"1A-volcanos_limma_cohorts_C.L.pdf"),width = 12,height = 5)

#all samples
res<-fread("outputs/01-lga_vs_ctrl_limma_DMCs_analysis/res_limma.tsv.gz")
res[,meth.change:=logFC][,p_val_adj:=adj.P.Val][,p_val:=P.Value]

res[p_val_adj<0.05&abs(meth.change)>25] #43 DMCs
res[p_val_adj<0.1&abs(meth.change)>25] #1255 DMCs
res[p_val_adj<0.1&meth.change>25] #1250 DMCs

res[p_val<0.001&abs(meth.change)>25] #4815 DMCs
res[p_val<0.001&meth.change>25] #4787 
4815-4787 
#28 hypo
cpgs<-fread("ref/all_HELP_tagging_cpgs_hg19_pos.csv")
res<-merge(res,cpgs)
fwrite(res[p_val<0.001&abs(meth.change)>25],fp(out,"supptabl2-res_limma_pval0.001_meth.change25.csv"))

ggplot(res)+
  geom_point(aes(x=meth.change,y=-log10(p_val),col=p_val<0.001&abs(meth.change)>25),size=0.1)+
  scale_color_manual(values = c("grey","red"))+scale_x_continuous(limits = c(-100,100))+theme_minimal()
ggsave(fp(out,"1A-volcano_plot_hypermet_LGA.pdf"))

#table : genescore, gene
resg<-fread("outputs/02-gene_score_calculation_and_validation/res_genes.csv.gz")
fwrite(resg[,.(gene,gene_score_add,n.cpg.gene,n.cpg.sig.gene,chr,pos,pval,padj)],fp(out,"table-genes_methylation_genescore.csv"))
summary(resg[n.cpg.sig.gene>3]$gene_score_add)

#% of the total number of analyzed CpGs fall in putative regulatory regions
cpgs_annot<-fread("outputs/02-gene_score_calculation_and_validation/cpgs_genes_annot_and_weight.csv.gz")
res_meth<-fread("outputs/01-lga_vs_ctrl_limma_DMCs_analysis/res_limma.tsv.gz")
cpgs_annotf<-cpgs_annot[cpg_id%in%res_meth$cpg_id]
unique(cpgs_annotf[ensembl_reg_score>0|chromatin_feature%in%c(4,5,6)],by="cpg_id")
381075/756470#50%
cpgs_annotf[in_eQTR==T]
253927/756470

res_methg<-fread("outputs/02-gene_score_calculation_and_validation/res_genes.csv.gz")

#1C- validation GS based on DEGs ####
res_anno<-fread("outputs/02-gene_score_calculation_and_validation/res_anno.csv.gz")
resg<-unique(res_anno[order(gene,pval)],by="gene")
resg[,mlog10pval:=-log10(pval)]
meth_scores<-melt(resg[,.SD,.SDcols=c("gene","gene_score_add","min.pval","avg.pval","mlog10pval","meth.change","avg.mlog10.pval","avg.meth.change","max.dmc_score", "avg.dmc_score")],id.vars = "gene",variable.name = "meth_metric",value.name = "score")
meth_scores[score==Inf,score:=NA]
meth_scores[score==-Inf,score:=NA]
meth_scores<-meth_scores[!is.na(score)]
res_de_cl<-fread("outputs/09-LGA_vs_Ctrl_Activated/res_pseudobulkDESeq2_all_cbps.csv.gz")

res_de_cl<-res_de_cl[!is.na(padj)]
meth_scores_de<-merge(meth_scores,res_de_cl,by=c("gene"))
unique(meth_scores_de[padj<0.05]$gene)#362 DEGS


meth_scores_de[,score_scaled:=scale(score),by="meth_metric"]
meth_scores_de$meth_metric<-factor(meth_scores_de$meth_metric,levels = c("min.pval","mlog10pval","meth.change","avg.pval","avg.mlog10.pval","avg.meth.change","max.dmc_score", "avg.dmc_score","gene_score_add"))
ggplot(meth_scores_de)+
  geom_boxplot(aes(fill=padj<0.05,y=score_scaled,x=meth_metric),outlier.shape = NA)+
  coord_cartesian(ylim = c(-3,3))+
  theme_classic()+
  scale_fill_manual(values = c("white","grey"))
ggsave(fp(out,"1C-boxplot_valid_genescore_with_degs_association.pdf"))

meth_scores_de[,pval:=wilcox.test(score[padj<0.05],score[padj>=0.05])$p.value,by="meth_metric"]
meth_scores_de[,score_diff:=mean(score_scaled[padj<0.05])-mean(score_scaled[padj>=0.05]),by="meth_metric"]
meth_scores_de[,FC:=mean(score[padj<0.05])/mean(score[padj>=0.05]),by="meth_metric"]

unique(meth_scores_de[,.(meth_metric,pval,score_diff,FC)],by="meth_metric")

#        meth_metric         pval
# 1:  gene_score_add 0.0006818586
# 2:        min.pval 0.0249488972
# 3:        avg.pval 0.7968159120
# 4:      mlog10pval 0.0757731002
# 5:     meth.change 0.0180098893
# 6: avg.mlog10.pval 0.3534897816
# 7: avg.meth.change 0.0368265523
# 8:   max.dmc_score 0.0192803029
# 9:   avg.dmc_score 0.4128018069

ggplot(unique(meth_scores_de[!meth_metric%in%c('min.pval','avg.pval')],by="meth_metric"))+
  geom_col(aes(y=-log10(pval),x=meth_metric,fill=FC))+
  scale_fill_gradient(low = "grey",high = "red")+
  theme(axis.text.x = element_text(angle = 45,hjust = 1))
ggsave(fp(out,"1C-barplot_valid_genescore_with_degs_association.pdf"))

#1D : pathway GSEA 
library(enrichplot)
res_go<-readRDS("outputs/03-pathway_analysis/res_gsea_go.rds")
max(data.table(as.data.frame(res_go))$p.adjust) #699 go term with padj < 0,001

pdf(fp(out,"1D-emapplot_gsea_go.pdf"),width = 14,height = 8)
emapplot(pairwise_termsim(res_go,showCategory = 40),showCategory = 40)
dev.off()

gsea_go<-fread("outputs/03-pathway_analysis/res_gsea_go.csv")



#supp : PC plot et cohort####
#PC plot
library(pheatmap)
source("scripts/utils/methyl_utils.R")
pval_cov<-fread("outputs/01-lga_vs_ctrl_limma_DMCs_analysis/pval_correl_covar_pcs.csv")

pval_covf<-pval_cov[!covar%in%c("groupbatch_complexity_fac","seq.depth","library_complexity","group_complexity","groupbatch_complexity")]
pval_covf[covar=="group_complexity_fac",covar:="library_complexity"]
pca<-readRDS("outputs/01-lga_vs_ctrl_limma_DMCs_analysis/pca_lgactrl.rds")

pval_mat<-as.matrix(data.frame(pval_covf,row.names = 1))

pct.varPCs<-pctPC(pca,1:10)*100
break_cols<-c(1,0.1,0.05,0.01,0.001,0.0001,0.00001)

pdf(fp(out,"supp10-pcs_covars_plot.pdf"))
pheatmap(-log10(pval_mat),cluster_cols = F,
           labels_col= paste0("PC",1:10,"(",round(pct.varPCs[as.character(1:10)]),"%)"),
           display_numbers = T,
           color = colorRampPalette(c("white", "red"))(length(break_cols)+1),
         breaks =sort(-log10(break_cols)))
dev.off()



####FIGURE 2 LGA vs CTRL scRNA-seq####
#3A-hematomap ####
hmap<-readRDS("../singlecell/outputs/02-hematopo_datasets_integration/hematomap_ctrls_sans_stress/hematomap_ctrls_sans_stress.rds")
hmap<-subset(hmap,cell_type!="18")

#reorder lineage levels
Idents(hmap)<-"lineage"

levels(hmap)<-c("LT-HSC",
                "HSC",
                "MPP/LMPP",
                "Lymphoid",
                "B cell",
                "T cell",
                "Erythro-Mas",
                "Mk/Er",
                "Myeloid",
                "DC")
#3A-UMAP 
hmap[["lineage"]]<-Idents(hmap)
DimPlot(hmap,group.by = c("lineage"),label = T)
ggsave(fp(out,"3A-hematomap.pdf"))


#+% OF EACH SUBPOPULATION
round(table(Idents(hmap))/ncol(hmap)*100)

#3B-markers ####
#dotplot
markers<-fread("../singlecell/ref/hematopo_markers.csv")
fwrite(markers[order(cell_type)],fp(out,"markers_hematopo_lineage.csv"))

markers[gene%in%c("CTSG","AZU1")]
markers[gene%in%c("CTSG","AZU1")]


markers[cell_type=="ErP"]

markers<-fread("outputs/05-make_hematomap/hematomap_ctrls_sans_stressSCT_Leiden_res0.6_markers.csv.gz")
markers[gene%in%c("CDK6","MLLT3")]
m_of_int<-markers[score>=3][MarqueurPops!=""][order(lineage)][cell_type!="18"]


m_anno<-fread("outputs/05-make_hematomap/markers_subpop")
setdiff(m_anno$marker,m_of_int$gene) #"ID2"   "DUSP2"  "FOS" "CDK6"
markers[gene%in%setdiff(m_anno$marker,m_of_int$gene)]
m_of_int<-rbind(m_of_int,markers)

m_of_int

fwrite(,fp(out,"markers_clusters_subpop.csv"))

m_lin<-fread("outputs/05-make_hematomap/markers_lineage_annotated.csv.gz")
m_lin


m_of_int<-rbind(m_lin[cluster!="18"][order(-score,p_val_adj)][,.SD[1:5],by="cluster"][,topspe:=T],m_lin[cluster!="18"][order(-score,p_val_adj)][MarqueurPops!=""][,.SD[1:5],by="cluster"][,topanno:=T],fill=T)[order(cluster)]
m_of_int<-unique(m_of_int[,topspe:=gene%in%gene[topspe==T],by="cluster"][,topanno:=gene%in%gene[topanno==T],by="cluster"])
m_of_int<-m_of_int[!is.na(gene)]
m_of_int[,MarqueurPops:=str_replace_all(MarqueurPops,";","|")]
fwrite(m_of_int,fp(out,"3B-markers_lineage.csv"))


DefaultAssay(hmap)<-"SCT"
key_genes_lin<-c("ID1","DUSP2", #LT-HSC
           "AVP","FOS", #HSC
           "MLLT3","CDK6", #MPP
           "CD99", #LMPP
           "LTB", #CLP
           "IGHM","CD37", #B cell
           "TNFAIP3","CD7", #T cell
           "GATA2", #MPP-Ery
           "GATA1", #EMP
           "PLEK","HBD", #Mk/Er
           "MPO","CEBPA","CTSG","AZU1", #GMP
         "CST3","CD83") #DC
DotPlot(hmap,features = key_genes_lin,
        group.by = "lineage",cols = c("white","black"))
ggsave(fp(out,"3B-dotplot_markers.pdf"))
DimPlot(hmap,group.by = "cell_type",label=T)

FeaturePlot(hmap,c("IGLL1"))


#featureplot
ps<-FeaturePlot(hmap,
                features = c("ID1","AVP","CDK6","LTB","GATA1","MPO"),
                max.cutoff = "q99",
                combine = F)
ps[[3]]<-FeaturePlot(hmap,
                features = c("CDK6"),
                min.cutoff = "q60")
ps<-lapply(ps, function(x)x+NoAxes()+NoLegend())
wrap_plots(ps)
ggsave(fp(out,"3B-featureplot_markers.pdf"))


#2C-distr LGA / Ctrl
cbps<-readRDS("outputs/06-integr_singlecell_cbps/cbps_light.rds")
cbps<-subset(cbps,lineage_hmap!="18")

#dimplot
DimPlot(cbps,group.by = c("group"))
ggsave(fp(out,"2C-dimplot_group.pdf"))

#barplot
ggplot(cbps@meta.data)+geom_bar(aes(x=group,fill=lineage_hmap),position = "fill")
ggsave(fp(out,"2C-barplot_group.pdf"))


#2C-volcano HSC ####
res_lin<-fread("outputs/09-LGA_vs_Ctrl_Activated/res_pseudobulkDESeq2_by_lineage.csv.gz")

res_hsc<-res_lin[lineage=="HSC"]
ggplot(res_hsc,aes(x=log2FoldChange,y=-log10(padj),col=padj<0.05&abs(log2FoldChange)>0.5))+
  geom_point()+
  scale_color_manual(values = c("grey","red")) +
  scale_x_continuous(limits=c(-4.5,4.5))+
  theme_minimal() +
  theme(legend.position = "bottom")
ggsave(fp(out,"4A-volcano_HSC.pdf"),height = 6.6)

#MA plot
genes_of_interest<-c("SESN2","SOCS3","HES1","JUN","FOS","JUNB","ZFP36","EGR1",
                      "DUSP2","DUSP1","FOSB","SOCS1","KLF2","KLF4",
                       "PLK2","PLK3","ID1","MYC","","ID2","IDS","RGCC")

ggplot(res_hsc,aes(y=log2FoldChange,x=baseMean,col=padj<0.05&abs(log2FoldChange)>0.5))+
  geom_point()+
  geom_label_repel(aes(label=ifelse(padj<0.05&abs(log2FoldChange)>0.5&gene%in%genes_of_interest,gene,"")),

                   max.overlaps = 3000,
                   size=3)+
  scale_color_manual(values = c("grey","red")) +scale_x_log10()+
  theme_minimal() +
  theme(legend.position = "bottom")
ggsave(fp(out,"2C-ma_plot_HSC.pdf"))

#all progens
progens<-c("LT-HSC","HSC","MPP/LMPP","Lymphoid","Myeloid","Erythro-Mas")
res_lin<-res_lin[lineage%in%progens][,lineage:=factor(lineage,levels = progens)]
ggplot(res_lin,aes(y=log2FoldChange,x=baseMean,alpha=padj<0.05&abs(log2FoldChange)>0.5,
                   col=padj<0.05&abs(log2FoldChange)>0.5))+
  geom_point()+
  facet_wrap("lineage")+
  scale_color_manual(values = c("grey","red")) +scale_x_log10()+
    scale_alpha_manual(values = c(0.6,1)) +
  theme_minimal() +
  theme(legend.position = "bottom")

ggsave(fp(out,"2Csupp-ma_plot_all_lin.pdf"))


#2D-pathway_down ####

#gsea go bp
library(clusterProfiler)
library(enrichplot)
res_gsea_go<-readRDS("outputs/09-LGA_vs_Ctrl_Activated/res_gsea_go_bp.rds")
pdf(fp(out,"4B-emapplot_gsea_go_bp_dn_top50.pdf"),width = 10,height = 6)
emapplot(pairwise_termsim(res_gsea_go),showCategory = 50,cex_label_category=0.66)
dev.off()


library(rrvgo)
go_dn<-readRDS("outputs/09-LGA_vs_Ctrl_Activated/res")
res_go_dn<-fread("outputs/09-LGA_vs_Ctrl_Activated/res_go_bp_dn.csv.gz")
res_go_dn[p.adjust<0.05]#47
res_go_dn[,n.overlap:=as.numeric(Count)]
res_go_dn[,n.query:=as.numeric(str_extract(GeneRatio,"[0-9]+$"))]
res_go_dn[,n.gene_set:=as.numeric(str_extract(BgRatio,"[0-9]+"))]
res_go_dn[,n.background:=as.numeric(str_extract(BgRatio,"[0-9]+$"))]

res_go_dn[,pct.overlap.query:=n.overlap/n.gene_set]
res_go_dn[,pct.overlap.background:=n.gene_set/n.background]
res_go_dn[,fold.enrichment:=pct.overlap.query/pct.overlap.background]
res_go_dn[,log2FE:=log2(fold.enrichment)]
fwrite(res_go_dn[p.adjust<0.05],fp(out,"4B-res_go_bp_down_padj0.05.csv"))

pdf(fp(out,"4B-emapplot_go_bp_dn_padj0.05.pdf"),width = 10,height = 6)
emapplot(pairwise_termsim(go_dn),showCategory = 47,cex_label_category=0.66)
dev.off()

####FIGURE 3 : correl Methylation / Expression####
#1A-MAP plot meth expr ####
quants<-quantile(res_mg$gene_score_add,1:100/100)
res_mg[,gs_percentile:=sum(gene_score_add>=quants),by="gene"]
res_hsc_dn<-fread("outputs/09-LGA_vs_Ctrl_Activated/res_pseudobulkDESeq2_by_lineage.csv.gz")[lineage=="HSC"]
res_meth_expr<-merge(res_mg,res_hsc_dn,by="gene")

ggplot(res_meth_expr,aes(y=log2FoldChange*-log10(padj.y),x=gs_percentile,col=padj.y<0.05&abs(log2FoldChange)>0.5))+
  geom_point()+
  geom_label_repel(aes(label=ifelse(padj.y<0.05&abs(log2FoldChange)>0.5&gene%in%genes_of_interest,gene,"")),
                    size=3,
                   max.overlaps = 3000)+
  scale_color_manual(values = c("grey","red")) +
  theme_minimal() +
  theme(legend.position = "bottom")
ggsave(fp(out,"4D-plot_expr_and_meth_change_.pdf"))


#1B-correl meth expr (boxplot)####
res_meth<-fread("outputs/02-gene_score_calculation_and_validation/res_genes.csv.gz")
res_merge<-merge(res_hsc,res_meth[,.(gene,gene_score_add)])
res_merge[,degs:=ifelse(padj<0.05&abs(log2FoldChange)>0.5,ifelse(log2FoldChange>0,"upregulated","downregulated"),"non-regulated")]
res_merge[,genes:=factor(degs,levels = c("non-regulated","upregulated","downregulated"))]
res_merge[,genes.:=paste0(genes,"\n(n=",.N,")"),by="genes"]

ggplot(res_merge)+
  geom_boxplot(aes(x=genes.,y=gene_score_add,fill=genes),outlier.shape = NA,col="bisque4")+
  scale_x_discrete(limits=unique(res_merge[order(genes)]$genes.))+
  scale_fill_manual(values = c("white","grey","black"))+
  coord_cartesian(ylim = c(0,1800))+theme_minimal()
ggsave(fp(out,"4C-boxplot_correl_meth_expr.pdf"))

res_merge[,p.up:=wilcox.test(gene_score_add[genes=='upregulated'],
                                                      gene_score_add[genes=='non-regulated'])$p.value]
res_merge[,p.down:=wilcox.test(gene_score_add[genes=='downregulated'],
                                                      gene_score_add[genes=='non-regulated'])$p.value]

unique(res_merge,by="p.down")
#  p.up       p.down
# 0.398549 2.373248e-05


#3C-overlap####
#pathways overlap

#emaplot meth expr

meth<-fread("outputs/03-pathway_analysis/res_gsea_go_bp_all.csv")[p.adjust<0.05]

meth[,mlog10padj.meth:=-log10(p.adjust)]
meth[mlog10padj.meth>10,mlog10padj.meth:=10]

hsc_dn<-fread("outputs/09-LGA_vs_Ctrl_Activated/res_gsea_go_bp.csv.gz")
hsc_dn<-hsc_dn[order(pvalue)][1:50]
meth_dn<-merge(hsc_dn,meth,all.x=T,by=c("ID","Description"),suffixes = c("",".meth"))

res_hsc_dn<-readRDS("outputs/09-LGA_vs_Ctrl_Activated/res_gsea_go_bp.rds")

res_hsc_dnm<-res_hsc_dn
res_hsc_dnm@result<-data.frame(meth_dn,row.names = "ID")
res_hsc_dnm@result$ID<-rownames(res_hsc_dnm@result)

pdf(fp(out,"3C-emapplot_overlap_methyl_on_expr_go_bp_top50.pdf"),width = 10,height = 6)
emapplot(pairwise_termsim(res_hsc_dnm),
         color="mlog10padj.meth",
         showCategory = 50,
         cex_label_category=0.66)+
  scale_color_gradient2(low = "grey",mid = "red",high = "darkred",midpoint = 5,limits=c(-log10(0.05),10))
dev.off()


####FIGURE 4 : ATAC ####
#4A : atac umap####
library(Seurat)
library(Signac)
inp<-"outputs/14-DMCs_atac_integr/"
atacs<-readRDS("outputs/14-DMCs_atac_integr/cbps_atacs.rds")
atacs[["lin_peaks"]]<-readRDS("outputs/14-DMCs_atac_integr/cbps_lin_spe_peaks_assay.rds")

DimPlot(atacs, group.by = "predicted.id", label = TRUE,reduction = "humap")
ggsave(fp(out,"4A-cbps_atac_umaps.pdf"))

#4B : TF motif enrichment in lineage spe peak ####
res_motif_lin<-fread(fp(inp,"enriched_motifs_by_lineage_specific_peaks.csv"))
resf<-res_motif_lin[pvalue<10^-50&lineage!="MPP/LMPP"]

resf[,motif.name:=factor(motif.name,levels=unique(motif.name[order(pvalue)]))]

res_motif_lin[,]
ggplot(resf)+
  geom_point(aes(x=motif.name,col=fold.enrichment,y=-log10(pvalue),size=observed))+
  facet_grid(~lineage,scales = "free",space="free_x")+
  scale_color_gradient2(low = "grey",mid = "red",high = "darkred")+
  scale_y_continuous()+scale_x_discrete(guide = guide_axis(angle = 60))+
  theme(axis.text=element_text(size=8))

ps<-lapply(unique(resf$lineage),function(lin)ggplot(resf[lineage==lin])+
  geom_point(aes(x=motif.name,col=fold.enrichment,y=-log10(pvalue),size=observed))+
  scale_color_gradient2(low = "white",mid = "darkgrey",high = "black",midpoint = 2,limits=c(1,8))+
  scale_y_continuous(expand = c(0.2,0.2))+
    scale_x_discrete(guide = guide_axis(angle = 60))+
    scale_size_continuous(limits = c(100,1200))+
  theme(axis.text=element_text(size=8))+ggtitle(lin))

wrap_plots(ps)+plot_layout(guides = 'collect')
ggsave(fp(out,"4B-lineage_spe_tf_mlog10pval50.pdf"),height = 8)


#4C : DMCs enrichment in HSC peak ####
res_dmcs_peaks_enr<-fread(fp(inp,"res_dmcs_peaks_enrichment_for_lineage_specific_peaks.csv"))
res_dmcs_peaks_enr[,percent.observed:=precision*100]
ggplot(res_dmcs_peaks_enr)+
  geom_point(aes(x=query,size=n.overlap,col=-log10(padj),y=percent.observed))+
  scale_color_gradient(low = "grey",high = "red")+
    scale_y_continuous(expand = c(0,1))+scale_x_discrete(guide = guide_axis(angle = 0))+
    theme(axis.title.x =element_blank())

ggsave(fp(out,"6C-dmcs_enrichment_in_lineage_spe_peaks.pdf"))


#4D : TF motif enrichment in peak down ####
hsc_lga_tf_all<-fread(fp(inp,"motif_enrichment_in_peaks_up_and_down_lga_vs_ctrl_hsc.csv.gz"))

hsc_lga_tf_all[,top20:=rank(pvalue)<=20,by="accessibility"]
hsc_lga_tf_all_20<-hsc_lga_tf_all[top20==T]
hsc_lga_tf_all_20[,motif.name:=factor(motif.name,levels =hsc_lga_tf_all_20[order(pvalue)]$motif.name,ordered = T )]
ggplot(hsc_lga_tf_all_20)+
  geom_point(aes(x=motif.name,size=fold.enrichment,y=-log10(pvalue),col=percent.observed))+
  facet_wrap("accessibility",scales = "free_x",ncol = 2)+
  scale_color_gradient2(low = "grey",mid = "red",high = "darkred",limits=c(0,100),midpoint = 50)+
  scale_y_continuous(limits=c(0,55),expand = c(0,0))+
  scale_x_discrete(guide = guide_axis(angle = 60))+
  theme(axis.text=element_text(size=8))
ggsave(fp(out,"4D-dotplot_TF_motif_enrichment_in_up_or_down_peaks_LGA_vs_Ctrl_HSC_top20.pdf"))
fwrite(hsc_lga_tf_all[pvalue<10^-6],fp(out,"4D-res_TF_motif_enrichment_in_up_or_down_peaks_LGA_vs_Ctrl_HSC_pval10m6.csv"))


#4E : DMCs and DEGs enrichment in DA peaks####
res_or_da_peaks_dmcs_degs<-fread(fp(inp,"res_da_peaks_enrichment_in_dmcs_degs_hsc_peaks.csv"))

res_or_da_peaks_dmcs_degs[,padj_b:=ifelse(padj<0.001,"***",ifelse(padj<0.01,"**",ifelse(padj<0.05,'*','ns')))]
res_or_da_peaks_dmcs_degs[,padj_b:=factor(padj_b,levels=c("ns","*","**","***"))]
res_or_da_peaks_dmcs_degs[,candidat_peaks:=paste0(term,"\n(n=",term.size,")")]
res_or_da_peaks_dmcs_degs[,da_peaks:=paste0(query,"\n(n=",n.query,")")]
res_or_da_peaks_dmcs_degs[,candidat_peaks:=factor(candidat_peaks,levels =res_or_da_peaks_dmcs_degs[query=="down_peaks"][order(fold.enrichment)]$candidat_peaks )]

ggplot(res_or_da_peaks_dmcs_degs)+geom_col(aes(y=log2(fold.enrichment),x=candidat_peaks,fill=padj_b),position = "dodge")+
  facet_wrap("da_peaks")+
  scale_fill_manual(values=c("grey","orange","red","darkred"))+
  scale_x_discrete(guide =guide_axis(angle = 66))+
  theme(axis.text.x = element_text(size = 8))
ggsave(fp(out,"Ah-barplot_dmcs_degs_peaks_enrichment_in_hsc_da_peaks.pdf"))


####FIGURE 5 : Regulons ####
#5A : Regulon activity LGA vs Ctrl####
cbps<-readRDS("outputs/10-SCENIC/cbps_with_regulons_activity.rds")
res_tf_diff<-fread("outputs/10-SCENIC/regulon_activity_lga_vs_ctrl_HTO_by_lineage.csv.gz")
res_tf_sig<-res_tf_diff[p_val_adj<0.001&lineage=="HSC"&hto==T&!str_detect(regulon,"e$")]
ggplot(res_tf_sig)+geom_density(aes(x=abs(avg_diff)))


res_tf_sig<-res_tf_diff[p_val_adj<0.001&abs(avg_diff)>0.03&lineage=="HSC"&hto==T&!str_detect(regulon,"e$")]

fwrite(res_tf_sig,fp(out,"5A-res_tf_activity_alteration_lga_hto_hsc_padj0.001_avg_diff0.3.csv"))

ggplot(res_tf_sig)+
  geom_col(aes(x=regulon,y=avg_diff,fill=-log10(p_val_adj)))+
 scale_x_discrete(limits = res_tf_sig[order(p_val_adj)]$regulon)+
  scale_fill_gradient(low="white",high="black",limits=c(10,60))

ggsave(fp(out,"5A-barplot_tf_change_hsc_lga_vs_ctrl_hto.pdf"))

mtd<-data.table(cbps@meta.data,keep.rownames = "cell")
tfs_alt<-res_hsc_htof$regulon
tf_alt_act<-data.table(t(as.matrix(cbps@assays$TF_AUC@data[tfs_alt,])),keep.rownames = "cell")
tf_alt_act<-melt(tf_alt_act,id.vars = "cell",variable.name ="regulon",value.name = "activity" )
tf_alt_act<-merge(tf_alt_act,mtd)
ggplot(tf_alt_act[lineage_hmap=="HSC"&hto==T&regulon!="KLF2e"])+
  geom_boxplot(aes(x=regulon,y=activity,fill=group))+theme_bw()
ggsave(fp(out,"5A-boxplot_tf_activity_change_padj0.001_avg_diff0.3_lga_vs_ctrl_hsc_hto.pdf"))

#5B : Regulons annotation####
library(ComplexHeatmap)
inp<-"outputs/08-HTO_signature/regulon_level/"
reg_enrich<-readRDS(fp(inp,"reg_enrich_gprofil.rds"))

res_bp<-fread(fp(inp,"res_regulons_enrich_for_BP_of_interest.csv"),select = c(1,11,3),col.names = c("regulon","term","padj"))
res_hsc<-fread(fp(inp,"res_regulons_enrich_hsc_state_signatures_geneset.csv"),select = c(9,1,8),col.names = c("regulon","term","padj"))
res_hscf<-res_hsc[term!='5FU_treatment']
terms_of_interest<-c("cellular response to stress",
                    "regulation of cellular response to stress" ,
                     "cellular response to growth factor stimulus",
                     "regulation of cellular response to growth factor stimulus" ,
                      "regulation of cell differentiation",
                     "regulation of hemopoiesis" ,
                     "hematopoietic progenitor cell differentiation"   ,
                    "cell population proliferation",
                     "regulation of cell population proliferation",
                     "regulation of cell cycle" )

enrich_mat<-as.matrix(data.frame(dcast(Reduce(rbind,list(res_bp[term%in%terms_of_interest],res_hscf)),formula = term~regulon),row.names = "term"))
enrich_mat<-(-log10(enrich_mat))
dim(enrich_mat)

col_fun4 = colorRamp2(c(0,3, 8), c("white", "grey","black"))
pdf(fp(out,"5B-regulon_of_int_annotation.pdf"),height = 12)
Heatmap(t(enrich_mat[,res_tf_sig$regulon]), 
                col=col_fun4,
                row_names_gp = grid::gpar(fontsize = 8),
                column_names_gp = grid::gpar(fontsize = 8),
                cluster_rows = T,
                cluster_columns=T,
                name="-log10(adjusted p-value)",
                row_title = "Regulon",
                column_title = "Gene sets")
dev.off()

#5C : Correlation Methylation Expression in Regulons####
res_gsea_gs<-fread("outputs/11-regulons_enrichment_genescore/res_gsea_genescore_regulons_hiconf.csv.gz",
                   select = c(1:3,6,9),
                     col.names = c("regulon","pval.gs","padj.gs","NES.gs","size.regulon"))
res_gsea_degs<-fread("outputs/11-regulons_enrichment_genescore/res_gsea_regulons_degs.csv",
                     select = c(1:3,6,9),
                     col.names = c("regulon","pval.degs","padj.degs","NES.degs","size.regulon"))

res_merge<-merge(res_gsea_gs,res_gsea_degs,all.x=T)

res_merge[padj.gs<0.01&NES.gs>1.6] #33


ggplot(res_merge,aes(y=NES.gs,x=NES.degs,col=padj.gs<0.01&NES.gs>1.6&padj.degs<0.01&abs(NES.degs)>1.6))+
  geom_point(aes(size=size.regulon))+
    geom_label_repel(aes(label=ifelse(padj.gs<0.01&NES.gs>1.6&padj.degs<0.01&abs(NES.degs)>1.6,regulon,"")),
                   max.overlaps = 3000)+
  scale_color_manual(values = c("grey","red"))+
  theme_minimal()+theme(legend.position = "bottom")

ggsave(fp(out,"5C-dotplot_NES_gsea_regulons_genescore_degs.pdf"))


#5D : TF motif enrichment in DMCs regions####
res<-fread("outputs/03B-motif_analysis/knownResults.txt",
           select = c(1,2,3,5,6,7,8,9),
           col.names = c("motif","consensus","pval","padj","n_dmc_with_motif","pct_dmc_with_motif","n_background_with_motif","pct_background_with_motif"))
res[padj<0.05]
res[,pct_dmc_with_motif:=as.numeric(str_remove(pct_dmc_with_motif,"%"))]
res[,pct_background_with_motif:=as.numeric(str_remove(pct_background_with_motif,"%"))]
res[,motif:=str_remove(motif,"/Homer")]
res[,fold.enrichment:=pct_dmc_with_motif/pct_background_with_motif]

res_perm<-fread("outputs/03B-motif_analysis/res_known_motif_all_perm_bgrandom.csv")
res_perm[,motif:=str_remove(motif,"/Homer")]

res[,permut:=0]
res_perm_merge<-merge(res,res_perm,all=T)
res_perm_merge[,p.perm:=sum(pct_dmc_with_motif[permut==0]<=pct_dmc_with_motif[permut!=0])/sum(permut!=0,na.rm = T),by=.(motif)]

res_perm_merge[padj<=0.05&n_dmc_with_motif>30&p.perm<0.05&permut==0]
res_perm_<-res_perm_merge[permut==0]

res_perm_[padj<=0.05] #26
res_perm_[padj<=0.05&n_dmc_with_motif>30] #26
res_perm_[padj<=0.05&n_dmc_with_motif>30&p.perm<0.01] #26

res_perm_[,motif.name:=sapply(motif,function(x)strsplit(x,"/")[[1]][1])]
ggplot(res_perm_[padj<=0.05&p.perm<0.01&fold.enrichment>1.30][order(pval)])+
  geom_point(aes(x=motif.name,y=-log10(pval),col=fold.enrichment,size=n_dmc_with_motif))+
  scale_x_discrete(limits=res_perm_[padj<=0.05&p.perm<0.01&fold.enrichment>1.30][order(pval)]$motif.name)+
  scale_color_gradient2(low = "white",mid = "grey",high = "black",midpoint = 1.25,limits=c(1,2.1))+
  scale_y_continuous(limits = c(0,25))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 80, hjust=1))

ggsave(fp(out,"5D-motif_enrichment_homer_own_background_no_norm_padj0.05_p.perm0.01_fold.enrichment1.3.pdf"))


#5E : TF motif enrichment in DMCs containing HSC peaks####

res_motif_hsc_dmcs<-fread(fp(inp,"res_tf_motif_enrichment_in_DMCs_containing_vs_CpGs_containing_hsc_peaks.csv"))
res_motif_hsc_dmcs[,pvalue:=pvalue+10^-316]
ggplot(res_motif_hsc_dmcs[1:15])+
  geom_point(aes(x=motif.name,size=fold.enrichment,col=-log10(pvalue),y=percent.observed))+
  scale_color_gradient2(low = "grey",mid = "red",high = "darkred",limits=c(1,317),midpoint = 150)+
  scale_y_continuous(limits=c(0,100),expand = c(0,0))+scale_x_discrete(guide = guide_axis(angle = 60),
                                                                       limits=res_motif_hsc_dmcs[1:15][order(pvalue,-fold.enrichment)]$motif.name)+
  scale_size(limits = c(1,2.7))+
  theme(axis.text=element_text(size=9))

ggsave(fp(out,"5E-dotplot_tf_motif_dmcs_vs_cpg_hsc_peaks_top15.pdf"))


#### FIGURE 6 : Regulatory Network ####
#6A : Coregulatory network SCENIC####
reg_egr1r<-fread(fp(out,"egr1_network_plus_tf_regulators_tf_target_interactions.csv"))
egr1_modul<-c("KLF2","FOSB","JUN","EGR1","KLF4","ARID5A","KLF10","JUNB")
reg_egr1_o<-reg_egr1r[tf%in%egr1_modul&target%in%egr1_modul] 

net_egr1<-as.network(reg_egr1_o[,.(tf,target)],loops = T,directed = T)
#add a vertex attributes wich indicates if the gene is a tf or not
net_egr1 %v% "type" = ifelse(network.vertex.names(net_egr1) %in% regf$tf, "tf", "gene")

#add methyl / expression info
res_m<-fread("outputs/02-gene_score_calculation_and_validation/res_genes.csv.gz")
res_e<-fread("outputs/09-LGA_vs_Ctrl_Activated/res_pseudobulkDESeq2_by_lineage.csv.gz")[lineage=="HSC"]

res_e[padj>0.05,deg:="black"]
res_e[padj<=0.05&log2FoldChange>0.5,deg:="red"]
res_e[padj<=0.05&log2FoldChange<(-0.5),deg:="blue"]

net_egr1 %v% "deg" = res_e[network.vertex.names(net_egr1),on="gene"]$deg

res_m[gene_score_add>500,meth:="cadetblue"]
res_m[gene_score_add<=500,meth:="cornsilk3"]

net_egr1 %v% "meth" = sapply(res_m[network.vertex.names(net_egr1),on="gene"]$meth,function(x)ifelse(is.na(x),"cornsilk3",x))

ggnet2(net_egr1,
       color = "meth",
       label = T,label.color = "deg",label.size = 3,
       size = "degree",size.min = 2,size.cut = 4,
       shape = "type" ,
       edge.alpha = 0.7,
       arrow.size = 5, arrow.gap =0.02) +
  theme(panel.background = element_rect(fill = "white"))

ggsave(fp(out,"6A-network_altered_regulons_and_genes_mincon2.pdf"))


#6B : genome tracks altered genes####
#SOCS3
gene<-"SOCS3"
peaks_hsc_anno<-fread("outputs/14-DMCs_atac_integr/peaks_hsc_genes_anno.csv.gz")
peaks_hsc_anno[gene_name==gene]
region<-"chr17-78352543-78361026"

cov_plot<-CoveragePlot(
  object = atacs_hsc,
  region = region,
  group.by = "group",
  annotation = F,
  peaks = F
)

gene_plot <- AnnotationPlot(
  object = atacs_hsc,
  region = region
)

peak_plot <- PeakPlot(
  object = atacs_hsc,assay = "lin_peaks",
  region = region
)

expr_plot <- ExpressionPlot(
  object = subset(cbps_hsc,hto==T),
  features = gene,
  assay = "SCT"
)
res_meth_reg<-MethChangeReg(res_meth,region)
res_meth_reg[-log10(P.Value)>5]
res_meth_reg[order(P.Value)]
meth_plot<-MethChangePlot(res_meth,region = region,limits = c(0,4))

tfs_plot<-TFsMotifPlot(atacs_hsc,
                      region = region,
                      motif.names= c("EGR1","KLF2"),
                      assay = "lin_peaks",
                      size=4,alpha = 0.6)

plots<-CombineTracks(
  plotlist = list(cov_plot, meth_plot, peak_plot,tfs_plot, gene_plot),
  expression.plot = expr_plot,
  heights = c(10, 6, 1,2,3),
  widths = c(9, 2)
)
plots
ggsave(fp(out,ps("6B-genome_track_",gene,"_v2.pdf")),plot = plots,height = 7)

#KLF2
gene<-"KLF2"
peaks_da_anno[gene_name==gene]
peaks_hsc_anno[gene_name==gene]
region<-"chr19-16319352-16330749"

cov_plot<-CoveragePlot(
  object = atacs_hsc,
  region = region,
  group.by = "group",
  annotation = F,
  peaks = F
)

gene_plot <- AnnotationPlot(
  object = atacs_hsc,
  region = region
)

peak_plot <- PeakPlot(
  object = atacs_hsc,assay = "lin_peaks",
  region = region
)

expr_plot <- ExpressionPlot(
  object = subset(cbps_hsc,hto==T),
  features = gene,
  assay = "SCT"
)
res_meth_reg<-MethChangeReg(res_meth,region)
res_meth_reg[-log10(P.Value)>5]
res_meth_reg[order(P.Value)]
meth_plot<-MethChangePlot(res_meth,region = region,limits = c(0,4))

tfs_plot<-TFsMotifPlot(atacs_hsc,
                      region = region,
                      motif.names= c("KLF4","KLF2","JUNB"),
                      assay = "lin_peaks",pad=20,
                      size=4,alpha = 0.6)
plots<-CombineTracks(
  plotlist = list(cov_plot, meth_plot, peak_plot,tfs_plot, gene_plot),
  expression.plot = expr_plot,
  heights = c(10, 6, 1,2,3),
  widths = c(9, 2)
)
plots
ggsave(fp(out,ps("6F-genome_track_",gene,"_v2.pdf")),plot = plots,height = 7)

#JUNB
res_meth_anno[gene=="JUNB"][order(pval)]
gene<-"JUNB"
peaks_hsc_anno<-fread("outputs/14-DMCs_atac_integr/peaks_hsc_genes_anno.csv.gz")
peaks_hsc_anno[gene_name==gene]
region<-"chr19-12790128-12793822"

cov_plot<-CoveragePlot(
  object = atacs_hsc,
  region = region,
  group.by = "group",
  annotation = F,
  peaks = F
)

gene_plot <- AnnotationPlot(
  object = atacs_hsc,
  region = region
)

peak_plot <- PeakPlot(
  object = atacs_hsc,assay = "lin_peaks",
  region = region
)

expr_plot <- ExpressionPlot(
  object = subset(cbps_hsc,hto==T),
  features = gene,
  assay = "SCT"
)

res_meth_reg<-MethChangeReg(res_meth,region)
res_meth_reg[,pval:=sapply(P.Value,function(x)ifelse(x<0.001,"***",ifelse(x<0.01,"**",ifelse(x<0.05,"*","ns"))))]
res_meth_reg[,pval:=factor(pval,levels = c("ns","*","**","***"))]
meth_plot<-ggplot(data = res_meth_reg) + geom_segment(aes(x = start, y = 0, 
        xend = end, yend = logFC,col=pval), size = 2, data = res_meth_reg)+
    scale_color_manual(values = c("grey","yellow","orange","red"))
  
meth_plot<-meth_plot+ theme_classic() + ylab(label = "Methylation change") + 
    xlab(label = paste0(seqid(region), " position (bp)")) + 
    xlim(c(start(region), end(region)))


tfs_plot<-TFsMotifPlot(atacs_hsc,
                      region = region,
                      motif.names= c("EGR1","KLF2"),
                      assay = "lin_peaks",
                      size=4,alpha = 0.6)

plots<-CombineTracks(
  plotlist = list(cov_plot, meth_plot, peak_plot,tfs_plot, gene_plot),
  expression.plot = expr_plot,
  heights = c(10, 6, 1,2,3),
  widths = c(9, 2)
)
plots
ggsave(fp(out,ps("6B-genome_track_",gene,"_v2.pdf")),plot = plots,height = 7)


#6C : multimodal GRN####

regulonsf<-fread(fp(out,"tf_target_interactions.csv"))[!(extended)]

#start build network only with tf> interact with high conf
regf<-regulons[(!extended)]
regf<-regf[!is.na(target)]
#only for tf of interest
egr1_modul<-c("KLF2","EGR1","KLF4")
reg_egr1<-regf[tf%in%c(egr1_modul)] #add only targets of the tfs altered
fwrite(reg_egr1,fp(out,"egr1_KLF2_KLF4_network_tf_target_interactions.csv"))
tfs<-unique(reg_egr1r1$tf)
net_egr1<-as.network(reg_egr1r1[,.(tf,target)],loops = T,directed = T)
net_egr1
#add a vertex attributes wich indicates if the gene is a tf or not
net_egr1 %v% "type" = ifelse(network.vertex.names(net_egr1) %in% regf$tf, "tf", "gene")

#add methyl info
res_m<-fread("outputs/02-gene_score_calculation_and_validation/res_genes.csv.gz")
res_m[gene_score_add>500,meth:="darkred"]
res_m[gene_score_add<=500,meth:="black"]

net_egr1 %v% "meth" = sapply(res_m[network.vertex.names(net_egr1),on="gene"]$meth,function(x)ifelse(is.na(x),"cornsilk3",x))


#add expr info 
res_e<-fread("outputs/09-LGA_vs_Ctrl_Activated/res_pseudobulkDESeq2_by_lineage.csv.gz")[lineage=="HSC"]

res_e[padj>0.05,deg:="cornsilk3"]
res_e[padj<=0.05&log2FoldChange>0,deg:="coral2"]
res_e[padj<=0.05&log2FoldChange>0.5,deg:="coral3"]
res_e[padj<=0.05&log2FoldChange<(0),deg:="cadetblue3"]
res_e[padj<=0.05&log2FoldChange<(-0.5),deg:="cadetblue4"]
res_e[padj<=0.05&log2FoldChange<(-0.25)]
net_egr1 %v% "deg" = res_e[network.vertex.names(net_egr1),on="gene"]$deg


#add atac info 
#on vertice
#need add target info
res_a<-fread("outputs/15-chromatin_change_LGA_vs_Ctrl/differential_peaks_accessibility_lga_vs_ctrl_hsc_logFC0.csv.gz")
res_a<-res_a[!str_detect(peak,"chr[XY]")]
peaks_hsc_genes[,peak:=query_region]
res_at<-merge(res_a,peaks_hsc_genes,by="peak")
res_at[,target:=gene_name]

res_at[p_val_adj<0.001&avg_log2FC>0.25,da:="red"]
res_at[p_val_adj<0.001&avg_log2FC<(-0.25),da:="blue"]
res_at[is.na(da),da:="grey75"]

net_egr1 %v% "da" = res_at[network.vertex.names(net_egr1),on="target"]$da

#on edge : df tf-target link with if peak with motif found, FC / pval of the change
#need merge network df with res_atac df
#add TF info on res_atac => for each peak, merge with tf(of the network)-peak dt 
tfs<-unique(reg_egr1r1[,.(tf,target)]$tf) #"KLF4" "EGR1" "KLF2" "ATF4" "ATF3" "JUN"  "FOS"  "JUNB"
peaks<-unique(res_at[target%in%reg_egr1r1$target]$peak)
length(peaks)#690

motif.all <- GetMotifData(
    object = atacs, assay = "lin_peaks", slot = "data"
  )

motifs_peaks_tfs <- motif.all[peaks,GetMotifIDs(atacs,tfs) , drop = FALSE]
tf_peak_dt<-melt(data.table(data.frame(as.matrix(motifs_peaks_tfs==1)),keep.rownames = "peak"),id.vars = "peak",variable.name ="motif.id" ,value.name = "is_present")
tf_peak_dt<-merge(tf_peak_dt,GetMotifIDs(atacs,tfs,return_dt=TRUE))
tf_peak_dt<-tf_peak_dt[is_present==TRUE]

res_at_tf<-merge(res_at,tf_peak_dt,by="peak")
res_at_tf[,tf:=motif.name]

#merge with network df
reg_egr1r1_peaks<-merge(reg_egr1r1,
                        res_at_tf[,.(tf,target,peak,p_val,p_val_adj,avg_log2FC,pct.1,pct.2,type,da)],
                        by = c("tf","target"),
                        all.x = T)

reg_egr1r1_peaks[,n.tf.target.peaks:=.N,by=.(tf,target)]
reg_egr1r1_peaks[,biggest_change:=p_val==min(p_val),.(tf,target)]
reg_egr1r1_peaks[(biggest_change)|is.na(biggest_change)]
unique(reg_egr1r1_peaks[(biggest_change)|is.na(biggest_change)],by=c("tf","target"))#ok


reg_egr1r1_peaks[,da.peak:=p_val_adj<0.001&abs(avg_log2FC)>0.25]

reg_egr1r1_peaks[,n.da.tf.target.peaks:=sum(da.peak),.(tf,target)]

reg_egr1r1_peaks[da.peak==T] #9

reg_egr1r1_peak1<-reg_egr1r1_peaks[(biggest_change)|is.na(biggest_change)]

#add DMCs infos on edge
#need merge peaks DMCs df with reg_egr1 df
peaks_cpgs<-fread("outputs/14-DMCs_atac_integr/cpgs_in_lin_OCRs.csv.gz")
peaks_meth<-merge(peaks_cpgs,fread("outputs/01-lga_vs_ctrl_limma_DMCs_analysis/res_limma.tsv.gz"))
peaks_meth[,peak:=peaks]
peaks_meth_hsc<-peaks_meth[peak%in%peaks_hsc_genes$peak]

reg_egr1r1_peak1_meth<-merge(reg_egr1r1_peak1,
                             peaks_meth_hsc[,.(peak,cpg_id,logFC,P.Value,adj.P.Val)],
                             by="peak",
                             all.x=T)
reg_egr1r1_peak1_meth[,n.cpg.peak:=.N,by=.(peak,tf)]
reg_egr1r1_peak1_meth[,biggest_meth_change:=P.Value==min(P.Value),.(peak,tf)]
reg_egr1r1_peak1_meth[(biggest_meth_change)|is.na(biggest_meth_change)] #ok


reg_egr1r1_peak1_meth[,dmcs:=P.Value<0.001&abs(logFC)>25]

reg_egr1r1_peak1_meth[,n.dmcs.peak:=sum(dmcs),.(peak,tf)]

reg_egr1r1_peak1_meth[dmcs==T] #24

reg_egr1r1_peak1_meth1<-reg_egr1r1_peak1_meth[(biggest_meth_change)|is.na(biggest_meth_change)]


#ADD edge atttibute (tf> target) : color depend of if atac based tf> target interact is altered by chromatin change
#if tf-gene peak dn : blue if tf-gene peak up :  red

net_egr1_a<-network(reg_egr1r1_peak1_meth1[,-c("extended","biggest_change","biggest_meth_change","peak")],loops = T,directed = T)
list.edge.attributes(net_egr1_a)
as.matrix(net_egr1_a,attrname='da')
net_egr1_a %e% "da"
net_egr1_a %e% "da"=sapply(net_egr1_a %e% "da",function(x)ifelse(x=="grey75","darkgrey",x))

net_egr1_a %e% "dmc_line"=net_egr1_a %e% "n.dmcs.peak"+1
net_egr1_a %e% "dmc_line"=sapply(net_egr1_a %e% "dmc_line",function(x)ifelse(is.na(x),1,x))

#set.edge.attribute(net_egr1, "color", ifelse(net_egr1 %e% "dap" > 1, "black", "grey75"))
#add vertices attributes
net_egr1_a %v% "type" = ifelse(network.vertex.names(net_egr1_a) %in% regf$tf, "tf", "gene")
net_egr1_a %v% "deg" = res_e[network.vertex.names(net_egr1_a),on="gene"]$deg
net_egr1_a %v% "deg" = sapply(net_egr1_a %v% "deg",function(x)ifelse(is.na(x),"cornsilk3",x))

net_egr1_a %v% "meth" = sapply(res_m[network.vertex.names(net_egr1_a),on="gene"]$meth,function(x)ifelse(is.na(x),"black",x))

#genes_of_interest<-union(res_e[padj<=0.05&abs(log2FoldChange)>0.5]$gene,union(res_m[gene_score_add>500]$gene,unique(reg_egr1$tf)))


#GRN sans selection,label all genes
ggnet2(net_egr1_a,
       color = "deg",
       label = T,label.color = "meth",label.size = 2,
       size = "type" ,size.palette = c("tf"=3,"gene"=1),
       shape = "type",
       edge.alpha = 0.8,
       edge.size=0.5,
       edge.color = "da",
       arrow.size = 5,
       edge.lty = "dmc_line",
       arrow.gap =0.02) +
  theme(panel.background = element_rect(fill = "white"))

ggsave(fp(out,"6C-final_network_EGR1_KLF2_KLF4_tf_targets_2.pdf"))


#### FIGURE 7 : Functional alteration ####
library(Seurat)
cbps<-readRDS("outputs/06-integr_singlecell_cbps/cbps_filtered.rds")

#7C : lineage bias
mtd<-data.table(cbps@meta.data,keep.rownames = "bc")
mtd$lineage_hmap<-factor(mtd$lineage_hmap,
                         levels = c("LT-HSC",
                                    "HSC",
                                    "MPP/LMPP",
                                    "Lymphoid",
                                    "B cell",
                                    "T cell",
                                    "Erythro-Mas",
                                    "Mk/Er",
                                    "Myeloid",
                                    "DC"))

mtd[,n.cells:=.N,"sample_hto"]
mtd[,n.cells.lin:=.N,c("sample_hto","lineage_hmap")]
mtd[,pct.lin:=.N/n.cells,c("sample_hto","lineage_hmap")]
mtsl<-unique(mtd[lineage_hmap%in%c("LT-HSC","HSC","MPP/LMPP","Erythro-Mas","Myeloid","Lymphoid")],by=c("sample_hto","lineage_hmap"))
mtsl$lineage_hmap<-factor(mtsl$lineage_hmap,levels = c("LT-HSC","HSC","MPP/LMPP","Erythro-Mas","Myeloid","Lymphoid"))

#rm outlyers
ggplot(mtsl)+geom_boxplot(aes(x=hto,y=pct.lin,fill=group))+facet_wrap("lineage_hmap")
  # scale_x_discrete(guide = guide_axis(angle = 45))+
  # scale_y_continuous(expand = c(0,0))

mtsl[,is.outlier:=pct.lin%in%boxplot.stats(pct.lin)$out,by=.(group,hto,lineage_hmap)]
mtsl[(is.outlier)]$sample_hto 
mtsl[,is.outlier:=pct.lin%in%boxplot.stats(pct.lin)$out,by=.(group,hto,lineage_hmap)]
mtsl[,outlier:=sample_hto%in%sample_hto[(is.outlier)]]
unique(mtsl[(outlier)]$sample_hto) #"lgaF552FALSE"  "ctrlF523TRUE"  "ctrlM537FALSE" "ctrlM530TRUE"
mtslf<-mtsl[(!outlier)]
ggplot(mtslf)+
  geom_boxplot(aes(x=lineage_hmap,y=pct.lin,fill=group))+facet_wrap("hto")
ggsave(fp(out,"7C-boxplot_lineage_proport_in_samples_by_group.pdf"))


mtslf[,pval:=wilcox.test(pct.lin[group=="lga"],pct.lin[group=="ctrl"])$p.value,by=c("hto","lineage_hmap")]

mthl<-unique(mtslf[,.(hto,lineage_hmap,pval)])
mthl[,padj:=p.adjust(pval),by="hto"]
mthl
#       hto lineage_hmap       pval       padj
#  1: FALSE     MPP/LMPP 0.66233766 1.00000000
#  2: FALSE          HSC 1.00000000 1.00000000
#  3: FALSE  Erythro-Mas 0.42857143 1.00000000
#  4: FALSE      Myeloid 0.79220779 1.00000000
#  5: FALSE     Lymphoid 0.79220779 1.00000000
#  6: FALSE       LT-HSC 0.28571429 1.00000000
#  7:  TRUE     MPP/LMPP 0.13203463 0.66017316
#  8:  TRUE          HSC 0.01515152 0.09090909
#  9:  TRUE  Erythro-Mas 1.00000000 1.00000000
# 10:  TRUE     Lymphoid 0.69913420 1.00000000
# 11:  TRUE      Myeloid 0.24025974 0.96103896
# 12:  TRUE       LT-HSC 0.73015873 1.00000000

#7D : UMAP pseudotime
pseudo_mtd<-fread("outputs/12A-Pseudotime_integrated/metadata_pseudotime_ComputRoot.csv")
cbps<-AddMetaData(cbps,
                  metadata = data.frame(pseudo_mtd[,.(bc,pseudotime)],row.names = "bc"),
                  col.name = "pseudotime")

Idents(cbps)<-"lineage_hmap"
FeaturePlot(cbps,"pseudotime",label=T,reduction = "ref.umap",cols = c("aquamarine1","darkorange2"))
ggsave(fp(out,"umap_pseudotime.pdf"))

mtd<-data.table(cbps@meta.data,keep.rownames = "bc")
mtd$lineage_hmap<-factor(mtd$lineage,levels = c("LT-HSC","HSC","MPP/LMPP","Erythro-Mas","Myeloid","Lymphoid"))

ggplot(mtd)+geom_boxplot(aes(y=pseudotime,x=lineage_hmap),alpha=0.7)+theme_minimal()
ggsave(fp(out,"boxplot_pseudotime_by_lineage.pdf"))

#7E : cumulative proportion of cells / pseudotime
#rm outliers
mtdf<-mtd[sample_hto%in%unique(mtslf$sample_hto)]
mtdf
n_bins<-30
pseudotime_thr<-1:n_bins*(max(mtdf$pseudotime,na.rm = T)/n_bins)
pseudo_bins<-data.table(bin=0:(n_bins-1),
                        pseudotime_thr=pseudotime_thr)
mtdf[,bin:=sum(pseudotime>=pseudotime_thr),by="bc"]
summary(mtdf$bin)

pseudo_bins_mtd<-merge(pseudo_bins,mtdf)

pseudo_bins_mtd[,n.cells:=.N,by="sample_hto"]
pseudo_bins_mtd[,n.cells.bin:=.N,by=.(sample_hto,bin)]
pseudo_bins_<-unique(pseudo_bins_mtd,by=c("sample_hto","bin"))
for(gr in c("ctrl","lga")){
  for(h in c(FALSE,TRUE)){
    samples<-unique(pseudo_bins_[group==gr&hto==h]$sample)
    for(i in 1:max(pseudo_bins_$bin)){
      for(s in setdiff(samples,pseudo_bins_[bin==i&group==gr&hto==h]$sample)){
        pseudo_bins_<-rbind(pseudo_bins_,pseudo_bins_[bin==i-1&group==gr&hto==h&sample==s][,bin:=bin+1][,n.cells.bin:=0])
        
      }
      
    }
    
  }
  
}

pseudo_bins_[,pct.cells.bin:=n.cells.bin/n.cells,by=.(sample_hto,bin)]

pseudo_bins_[,n.cells.below.thr:=NA]
pseudo_bins_[,n.cells.below.thr:=as.numeric(n.cells.below.thr)]

for(binn in 0:n_bins){
  pseudo_bins_[bin<=binn,n.cells.below.thr:=ifelse(is.na(n.cells.below.thr),sum(n.cells.bin),n.cells.below.thr),by="sample_hto"]

  }

pseudo_bins_[,pct.cells.below.thr:=n.cells.below.thr/n.cells]
ggplot(pseudo_bins_,aes(x=factor(bin),y=pct.cells.below.thr,fill=group))+
  geom_boxplot()+
facet_wrap("hto")
ggsave(fp(out,"7E-boxplot_cumulated_pseudotime_by_group.pdf"))

#visuellement le pique de l'influence semble être à bin = 7 wich corresponf to HSC cells
#so, is there less cells with pseudotim <7 in stimulated LGA compared to stimulated control ?

wilcox.test(pseudo_bins_[hto==T&bin==7&group=="lga"]$pct.cells.below.thr,pseudo_bins_[hto==T&bin==7&group=="ctrl"]$pct.cells.below.thr)
#p-value = 0.04113


#+ density plot lineage by pseudotime
ggplot(mtdf)+geom_density(aes(x=pseudotime,fill=lineage),alpha=0.5)+theme_minimal()
ggsave(fp(out,"7E-density_pseudotime_by_lineage.pdf"))



