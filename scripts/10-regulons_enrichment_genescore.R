#TF signature enriched in DMC  ?
source("../methyl/scripts/utils/new_utils.R")
out<-here("outputs/11-regulons_enrichment_genescore")
dir.create(out)
regulons_list<-readRDS("../singlecell/outputs/05-SCENIC/cbps0-8_clean/regulons_list.rds")

#in GeneScore
res_methg<-fread("outputs/02-gene_score_calculation_and_validation/res_genes.csv.gz")

#GSEA
library(fgsea)

res_methg[,gs_rank:=rank(gene_score_add)]

genes_rank<-res_methg$gs_rank
names(genes_rank)<-res_methg$gene
genes_rank<-sort(genes_rank,decreasing = T)
head(genes_rank)
regulons_listf<-regulons_list[!str_detect(names(regulons_list),"e$")]
length(regulons_listf) #106
res_gsea<-fgsea(pathways=regulons_listf,
      stats=genes_rank,eps=0,scoreType = "pos",minSize = 10,maxSize = 500)

res_gsea[,size.regulon:=length(regulons_listf[[pathway]]),by="pathway"]
head(res_gsea[order(padj)],30)
res_gsea[pathway=="EGR1"]
fwrite(res_gsea,fp(out,"res_gsea_genescore_regulons_hiconf.csv.gz"))
res_gsea[padj<0.002]
res_gsea[padj<0.002]$regulon
res_gsea[order(padj)]$pathway[1:10]
res_gsea[pathway=="STAT3e"]
res_gsea[pathway=="STAT3"]
res_gsea[pathway=="EGR1"]


#with clusterprofiler
library(clusterProfiler)
library(enrichplot)
?GSEA

regulons_df<-Reduce(rbind,lapply(names(regulons_list), function(tf)data.table(term=tf,gene=regulons_list[[tf]])))
regulons_df<-regulons_df[,tf:=term][,.(tf,gene)]
fwrite(regulons_df,"outputs/10-SCENIC/regulons_extended.csv")
fwrite(regulons_df[!str_detect(tf,"e$")],"outputs/10-SCENIC/regulons.csv")

res_gsea_cl<-GSEA(geneList = genes_rank,
                  TERM2GENE =data.frame(regulons_df) ,scoreType="pos",
                  minGSSize = 10,maxGSSize = 500,eps=0,pvalueCutoff = 1)
res_gsea_cl_dt<-data.table(as.data.frame(res_gsea_cl))

saveRDS(res_gsea_cl,fp(out,"res_gsea_genescore_regulon_clusterprofiler.rds"))
fwrite(res_gsea_cl_dt,fp(out,"res_gsea_genescore_regulon_clusterprofiler.csv"))

res_gsea_cl_dt[p.adjust<0.01]
res_gsea_cl_dt[ID=="EGR1"]
head(res_gsea_cl_dt,100)
head(res_gsea_cl_dt[order(p.adjust)]$ID,100)
emapplot(pairwise_termsim(res_gsea_cl,showCategory = 30),showCategory = 30)

pdf("outputs/07-regulons_enrichment_genescore/emapplot_gsea_genescore_regulon_clusterprofiler.pdf",width=14,height=8)
emapplot(pairwise_termsim(res_gsea_cl,showCategory = 100),showCategory = 100)
dev.off()

res_gsea_cl2<-GSEA(geneList = genes_rank,TERM2GENE =data.frame(regulons_df[!str_detect(term,"e$")]) ,
                   scoreType="pos",
                   minGSSize = 10,maxGSSize = 500,eps=0)
res_gsea_cl2_dt<-data.table(as.data.frame(res_gsea_cl2))

saveRDS(res_gsea_cl2,fp(out,"res_gsea_genescore_regulon_non_extended_clusterprofiler.rds"))
fwrite(res_gsea_cl2_dt,fp(out,"res_gsea_genescore_regulon_non_extended_clusterprofiler.csv"))

res_gsea_cl2_dt[p.adjust<0.1]
res_gsea_cl2_dt[ID=="EGR1"]
head(res_gsea_cl2_dt[order(p.adjust)]$ID,100) #
res_gsea_cl2<-pairwise_termsim(res_gsea_cl2,showCategory = 69)

emapplot(res_gsea_cl2,showCategory = 50,min_edge=0.15)

pdf("outputs/07-regulons_enrichment_genescore/emapplot_gsea_genescore_regulon_non_extended_clusterprofiler.pdf",width=14,height=8)
emapplot(res_gsea_cl2,showCategory = 69)
dev.off()

res_gsea_cl2_sub = res_gsea_cl2[res_gsea_cl2$ID%in%c("ATF3","FOS","JUNB","KLF2","FOSB","JUN","EGR1","JUND","KLF10","KLF4","ARID5A"), asis=T]
class(res_gsea_cl2_sub)

res_degs<-fread("outputs/09-LGA_vs_Ctrl_Activated/res_scEdgeR_by_lineage.csv.gz")[lineage_hmap=="HSC"]
fold_changes<-res_degs[p_val_adj<0.001&abs(avg_logFC)>0.3]$avg_logFC
names(fold_changes)<-res_degs[p_val_adj<0.001&abs(avg_logFC)>0.3]$gene
fold_changes_bin<-fold_changes
fold_changes_bin[fold_changes<0]<--1
fold_changes_bin[fold_changes>0]<-+1

cnetplot(res_gsea_cl2_sub,showCategory = 11,
         foldChange = fold_changes_bin,cex_gene=0.5,cex_label_gene=0.5,
         )

pdf("outputs/07-regulons_enrichment_genescore/cnetplot_gsea_genescore_EGR1JUNB_regulons_network_dn_upreg_genes_highlight.pdf",width=14,height=8)

cnetplot(res_gsea_cl2_sub,showCategory = 11,
         foldChange = fold_changes_bin,cex_category=0.8,cex_label_category=0.8,
         node_label = "category"
         )
dev.off()


#GSEA - degs
library(fgsea)
res_degs<-fread("outputs/09-LGA_vs_Ctrl_Activated/res_pseudobulkDESeq2_by_lineage.csv.gz")[lineage=="HSC"]


res_degs[log2FoldChange>0,fc_rank:=rank(log2FoldChange*-log10(pvalue))]
res_degs[log2FoldChange<0,fc_rank:=-rank(abs(log2FoldChange)*-log10(pvalue))]

genes_rank<-res_degs$fc_rank
names(genes_rank)<-res_degs$gene
genes_rank<-sort(genes_rank,decreasing = T)
head(genes_rank)
res_gsea<-fgsea(pathways=regulons_listf,
      stats=genes_rank,eps=0,scoreType = "std",minSize = 10,maxSize = 500)

res_gsea[,size.regulon:=length(regulons_listf[[pathway]]),by="pathway"]



fwrite(res_gsea,fp(out,"res_gsea_regulons_degs.csv"))


#old : 
res_gsea_cl2_sub = res_gsea_cl2[res_gsea_cl2$ID%in%c("ATF3","FOS","JUNB","KLF2","FOSB","JUN","EGR1","JUND","KLF10","KLF4","ARID5A"), asis=T]
class(res_gsea_cl2_sub)

res_degs<-fread("outputs/09-LGA_vs_Ctrl_Activated/res_scEdgeR_by_lineage.csv.gz")[lineage_hmap=="HSC"]
fold_changes<-res_degs[p_val_adj<0.001&abs(avg_logFC)>0.3]$avg_logFC
names(fold_changes)<-res_degs[p_val_adj<0.001&abs(avg_logFC)>0.3]$gene
fold_changes_bin<-fold_changes
fold_changes_bin[fold_changes<0]<--1
fold_changes_bin[fold_changes>0]<-+1

cnetplot(res_gsea_cl2_sub,showCategory = 11,
         foldChange = fold_changes_bin,cex_gene=0.5,cex_label_gene=0.5,
         )

pdf("outputs/07-regulons_enrichment_genescore/cnetplot_gsea_genescore_EGR1JUNB_regulons_network_dn_upreg_genes_highlight.pdf",width=14,height=8)

cnetplot(res_gsea_cl2_sub,showCategory = 11,
         foldChange = fold_changes_bin,cex_category=0.8,cex_label_category=0.8,
         node_label = "category"
         )
dev.off()



?bitr
renv::install("ggnewscale")

#over representation test 
res_methg[,gene_altered:=abs(gene_score_add)>150]

res_or<-data.table(regulon=names(regulons_list),regulon.size=sapply(regulons_list,length))
genes_altered<-res_methg[gene_altered==T]$gene
res_or[,n.genes.altered:=length(genes_altered)]
res_or[,n.enriched:=sum(genes_altered%in%regulons_list[[regulon]]),by="regulon"]
res_or[,genes.enriched:=paste(genes_altered[genes_altered%in%regulons_list[[regulon]]],collapse="|"),by="regulon"]
res_or[,pct.enriched:=n.enriched/regulon.size]

size_universe<-length(res_methg$gene)

res_or[,pval:=phyper(q=n.enriched-1, 
     m=n.genes.altered, 
     n=size_universe-n.genes.altered, 
     k=regulon.size, 
     lower.tail=FALSE),
     by="regulon"]

res_or[,padj:=p.adjust(pval,method = 'BH')]
  
fwrite(res_or,fp(out,"res_or_genescore150_regulons.csv.gz"))


#over representation test with universe = genes of regulons and padj
res_methg[,gene_altered:=abs(gene_score_add)>300&pval_gs<0.01]

universe<-Reduce(union,regulons_list)
size_universe<-length(universe)
size_universe

res_or<-data.table(regulon=names(regulons_list),regulon.size=sapply(regulons_list,length))
genes_altered_uni<-intersect(res_methg[gene_altered==T]$gene,universe)
res_or[,n.genes.altered:=length(genes_altered_uni)]
res_or[,n.enriched:=sum(genes_altered_uni%in%regulons_list[[regulon]]),by="regulon"]
res_or[,genes.enriched:=paste(genes_altered_uni[genes_altered_uni%in%regulons_list[[regulon]]],collapse="|"),by="regulon"]
res_or[,pct.enriched:=n.enriched/regulon.size]

res_or[,pval:=phyper(q=n.enriched-1, 
     m=n.genes.altered, 
     n=size_universe-n.genes.altered, 
     k=regulon.size, 
     lower.tail=FALSE),
     by="regulon"]

res_or[,padj:=p.adjust(pval,method = 'BH')]
res_or[padj<0.1]



#pca on the cpgs determining TF #[to update]

genes_determining_tf<-Reduce(union,regulons_list[unique(res_or_all[padj<0.001&!str_detect(regulon,"e$")]$regulon)])

res_meth[,tf_determining_genes:=gene%in%genes_determining_tf]

res_meth[,tf_determining_cpgs:=tf_determining_genes&pval<0.001&abs(meth.change)>25]

cpgs_determining_tf<-unique(res_meth[tf_determining_cpgs==T]$locisID)
length(cpgs_determining_tf)
meth<-fread("../methyl/datasets/cd34/2020-05-25_methyl_data_before_limma.csv")
meth_deter<-meth[locisID%in%cpgs_determining_tf] #1413 cpgs

pca<-prcomp(t(as.matrix(meth_deter[,.SD,.SDcols=mtd[Group_name%in%c("C","L")]$sample])))

pca$x
vars_pcs<-GetVarPCs(pca)
round(vars_pcs*100)

res_pca_deter<-as.data.table(pca$x,keep.rownames = "sample")

mtd<-fread("../methyl/datasets/cd34/cleaned_batch_CD34_250121.csv")

res_pca_deter<-merge(res_pca_deter,mtd,by="sample")

ggplot(res_pca_deter)+geom_point(aes(x=PC1,y=PC2,col=Group_name))
ggplot(res_pca_deter)+geom_point(aes(x=PC1,y=PC2,col=Group_Sex))
ggplot(res_pca_deter)+geom_point(aes(x=PC1,y=PC3,col=Group_Sex))
ggplot(res_pca_deter)+geom_point(aes(x=PC1,y=PC3,col=Group_Sex))

ggplot(res_pca_deter)+geom_boxplot(aes(x=Group_Sex,y=PC1,fill=Group_Sex))
