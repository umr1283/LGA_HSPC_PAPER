#02A1a-meta_eQTL
library(data.table)
library(stringr)
library(biomartr)
out<-"outputs/02A1a-meta_eQTL"
dir.create(out)

eqtl<-fread("ref/GTEx_Analysis_v8.metasoft.txt.gz")
#I) recover col gene, chr, pos 
#1) in hg38 :
eqtl[,chr.hg38:=str_extract(RSID,"chr[0-9XY]+")]
eqtl[,pos.hg38:=str_extract(RSID,"_[0-9]+")]
eqtl[,pos.hg38:=as.numeric(str_extract(pos.hg38,"[0-9]+"))]
eqtl[,gene_id:=str_extract(RSID,"ENSG[0-9]+")]
eqtl

#trans in symbol :  with hsapiens_gene_ensembl when biomartr server will work
genes_translator <- biomart( genes      = unique(eqtl$gene_id), # genes that we wanted info
                                mart       = "ENSEMBL_MART_ENSEMBL", # marts were selected with biomartr::getMarts()
                                dataset    = "hsapiens_gene_ensembl", # datasets were selected with biomartr::getDatasets()
                                attributes = "hgnc_symbol", # attributes were selected with biomartr::getAttributes()
                                filters    = "ensembl_gene_id")

genes_translator<-data.table(genes_translator)
genes_translator<-genes_translator[,gene_id:=ensembl_gene_id][,-"ensembl_gene_id"]
genes_translator<-genes_translator[,gene:=hgnc_symbol][,-"hgnc_symbol"]
eqtl<-merge(eqtl,genes_translator,all.x=T,by="gene_id")


#2) trans in hg19 :
translator<-fread("ref/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.lookup_table.txt.gz",
                  select = c(1,8))
translator

eqtl[,variant_id:=str_remove(RSID,",ENSG[0-9]+\\.[0-9]+")]


eqtl<-merge(eqtl,translator,all.x=T,by="variant_id")


eqtl[,chr.hg19:=paste0("chr",str_extract(variant_id_b37,"^[0-9XY]{1,2}"))]

eqtl[,pos.hg19:=as.numeric(str_sub(str_extract(variant_id_b37,"_[0-9]+"),2))]
eqtl

#anno for tissue_wide eQTL : 
# - eQTL with low Heterogeneity : Qvalue <0.1 ;
# - eQTL in >30/49 tissue =>  mval>0.8 & pval <10^-4 in >30 tissue
# - PVALUE_FE < 10^-30
# -at least 20 tissues with mval > 0.8
eqtl[,n.tissue_sig:=rowSums(eqtl[,.SD,.SDcols=str_detect(colnames(eqtl),"mval")]>0.8 & eqtl[,.SD,.SDcols=str_detect(colnames(eqtl),"pval")]>0.0001,na.rm = T)]
eqtl[,tissue_wide:=PVALUE_Q>0.1&n.tissue_sig>30&PVALUE_FE<10e-30]
nrow(eqtl[tissue_wide==T]) #552921
unique(eqtl[tissue_wide==T&!is.na(gene)],by="gene") #5641
# #anno for liver
# eqtl[,liver:=mval_Liver>0.8&pval_Liver>10^-4]
# nrow(eqtl[liver==T]) #3578315

fwrite(eqtl,fp(out,"meta_eQTL_annotated.csv.gz"))

#formatage tissue wide eqtls :
eqtl[,chr:=chr.hg19]
eqtl[,pos:=pos.hg19]
eqtl[gene=="",gene:=NA]
eqtl_meta<-eqtl[tissue_wide==T&!is.na(gene),.(chr,pos,gene,PVALUE_FE)]

#add tss distance
tss<-fread("ref/hg19/tss_genes.bed")
tss[,tss_pos:=start]

eqtl_meta<-merge(eqtl_meta,tss[,.(gene,chr,tss_pos,strand)],allow.cartesian=TRUE)
eqtl_meta[strand=="+",tss_dist:=pos-tss_pos][strand=="-",tss_dist:=tss_pos-pos]
eqtl_meta[,n.gene.eqtl:=length(unique(gene)),by=c("chr","pos")]
summary(eqtl_meta$n.gene.eqtl)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#   1.000   1.000   1.000   1.416   2.000   9.000 
summary(abs(eqtl_meta$tss_dist))
   # Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
   #     0   188748   430644   640813   806979 28986909 
#keep only most closest eqtl-tss link
eqtl_meta[,closest:=abs(tss_dist)==min(abs(tss_dist)),by=c("chr","pos","gene")]
eqtl_meta<-eqtl_meta[closest==T]
#and save
fwrite(eqtl_meta,fp(out,"tissue_wide_signif_variant_gene_pairs_hg19.csv.gz"))
