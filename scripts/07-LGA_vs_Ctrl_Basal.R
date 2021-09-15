
out<-"outputs/07-LGA_vs_Ctrl_Basal"
dir.create(out)
source("scripts/utils/new_utils.R")
library(Seurat)
library(Matrix)
library(Matrix.utils)
library(DESeq2)

cbps<-readRDS("outputs/06-integr_singlecell_cbps/cbps_filtered.rds")

cbps_b<-subset(cbps,hto==F)

#pseudobulk all cbps
#get mtd of interest
mtd<-data.table(cbps_b@meta.data,keep.rownames = "bc")
mts<-unique(mtd,by=c("sample"))
#get counts and filter genes lowly express
counts<-as.matrix(cbps_b@assays$RNA@counts)
dim(counts) 

counts <- counts[rowSums(counts > 0) >= 100|rowSums(counts > 0)>=ncol(counts)*0.1, ] 
message(nrow(counts)," genes kept after filtering") 

  # Aggregate across cluster-sample groups
sample_counts <- t(aggregate.Matrix(t(counts[,mtd$bc]), 
                     groupings = mtd$sample, fun = "sum"))
#DEseq2_analysis
dds <- DESeqDataSetFromMatrix(sample_counts, 
                               colData = data.frame(mts,row.names="sample")[colnames(sample_counts),], 
                               design = ~ group+orig.ident+sex)

dds <- DESeq(dds)

res <- results(dds,contrast = c("group","lga","ctrl"),alpha = 0.05)

res<-data.table(as.data.frame(res),keep.rownames="gene")

res[padj<0.05]$gene
fwrite(res,fp(out,"res_pseudobulkDESeq2_all_cbps.csv.gz"))


#Distribution pareil
cbps_b$lineage_hmap<-factor(cbps_b$lineage_hmap,levels = c("LT-HSC","HSC","MPP/LMPP","Lymphoid","B cell","T cell","Erythro-Mas","Mk/Er","Myeloid","DC"))


mtd<-data.table(cbps_b@meta.data,keep.rownames="bc")
mtd[,n.sample:=.N,by="sample_hto"]
mtd[,pct.lin:=.N/n.sample,by=c("sample_hto","lineage_hmap")]

fwrite(mtd,fp(out,"metadata_basal.csv"))
ggplot(unique(mtd,by=c("sample","lineage_hmap")))+
  geom_boxplot(aes(x=lineage_hmap,y=pct.lin,fill=group))
ggsave("outputs/figures_epi_response/figure2/2supp-distribution_lineage_control_lga_basal.pdf")

#No expr change in lineage
#pseudobulk
Idents(cbps_b)<-"lineage_hmap"
res_lin<-Reduce(rbind,lapply(levels(cbps_b),function(lin){
  print(lin)
  cbps_sub<-subset(cbps_b,lineage_hmap==lin)
  #get mtd of interest
  mtd<-data.table(cbps_sub@meta.data,keep.rownames = "bc")
  mts<-unique(mtd,by=c("sample"))
  #get counts and filter genes lowly express
  counts<-as.matrix(cbps_sub@assays$RNA@counts)
  dim(counts) 

  counts <- counts[rowSums(counts > 0) >= 100|rowSums(counts > 0)>=ncol(counts)*0.1, ] 
  message(nrow(counts)," genes kept after filtering") 
  if(nrow(counts)>0){
    # Aggregate across cluster-sample groups
  sample_counts <- t(aggregate.Matrix(t(counts[,mtd$bc]), 
                       groupings = mtd$sample, fun = "sum"))
  #DEseq2_analysis
  dds <- DESeqDataSetFromMatrix(sample_counts, 
                                 colData = data.frame(mts,row.names="sample")[colnames(sample_counts),], 
                                 design = ~ group+orig.ident+sex)
  
  dds <- DESeq(dds)
  
  res <- results(dds,contrast = c("group","lga","ctrl"),alpha = 0.05)
  
  return(data.table(as.data.frame(res),keep.rownames="gene")[,lineage:=lin])
  }else{
    return(data.table())
      }
  

  }))

fwrite(res_lin,fp(out,"res_pseudobulkDESeq2_by_lineage.csv.gz"),sep=";")


#sc 
#run 07A #[to redo]


