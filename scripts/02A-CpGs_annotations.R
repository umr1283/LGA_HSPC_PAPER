#ANNOT CPGs

source("scripts/utils/new_utils.R")
cpgs<-fread("datasets/cd34/CD34_angle_119_noEmptyLocis_withConfScore_withoutChrXY.txt",
            select = c("id","chr","start"),
            col.names = c("cpg_id","chr","pos"))
fwrite(cpgs,"ref/all_HELP_tagging_cpgs_hg19_pos.csv")
out<-"outputs/02A-CpGs_annotations"
dir.create(out)

#1)cpgs_genes_links
#by distance to tss
#first gene within +/-200kb

genes_in_200kb<-bed_inter(cpgs[,start:=pos-200e3][start<0,start:=0][,end:=pos+200e3][,.(chr,start,end,cpg_id)][order(chr,start)],
          "ref/hg19/tss_genes.bed",
          select = c(5,6,8,9,4),col.names = c("chr","tss_pos","gene","strand","cpg_id"))

genes_in_200kb
cpgs_genes_tss<-merge(cpgs,genes_in_200kb,all.x=T)[,.(chr,pos,cpg_id,gene,tss_pos,strand)]
cpgs_genes_tss[strand=="+",tss_dist:=pos-tss_pos][strand=="-",tss_dist:=tss_pos-pos]
cpgs_genes_tss[,n.gene.cpg:=length(unique(gene)),by="cpg_id"]
summary(cpgs_genes_tss$n.gene.cpg)
  # Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  #  1.00    6.00   10.00   12.36   17.00  101.00 
summary(abs(cpgs_genes_tss$tss_dist))
   # Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA
   #    0   45620   96440   96946  147918  200000  131093 

#why a lot of NA ?
cpgs_genes_tss[is.na(tss_dist)] #131093 CpGs with no gene at +/-200kb

cpgs_genes_tss[!is.na(tss_dist),closest.rank:=rank(abs(tss_dist)),by="cpg_id"]

summary(abs(cpgs_genes_tss[closest.rank==1]$tss_dist))
  # Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
  #     0    4455   18193   35849   50208  200000  131093 

summary(abs(cpgs_genes_tss[closest.rank==2]$tss_dist))
   # Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
   #    0   10139   31030   48263   72131  199996 

fwrite(cpgs_genes_tss,fp(out,"cpg_tss_200kb.csv.gz"))


cpgs_genes_1<-cpgs_genes_tss[!is.na(tss_dist)&closest.rank==1]

#for >1 cpg_gene link (because multiple TSS), keep most closest cpg_gene link
cpgs_genes_1<-unique(cpgs_genes_1[order(cpg_id,gene,abs(tss_dist))],by=c("cpg_id","gene")) 

cpgs_genes_1
rm(cpgs_genes_tss)

fwrite(cpgs_genes_1,fp(out,"cpgs_closest_gene_tss_linked_within_200kb_around.tsv"),sep="\t")
cpgs_genes_1<-fread(fp(out,"cpgs_closest_gene_tss_linked_within_200kb_around.tsv"),sep="\t")

#by presence in eQTL region
#need first create eQTR
#see 02A1-create_eQTR
eqtrs<-fread("outputs/02A1-create_eQTR/all_eQTR_symbol.bed")
eqtrs[,eqtr_id:=paste(eqtr_id,tissue,sep="-")]
eqtrs<-unique(eqtrs[!is.na(start.eQTR)],by="eqtr_id")
table(eqtrs$tissue)
# tissue_wide whole_blood 
#       37112      118529 

table(unique(eqtrs,by=c("gene","tissue"))$tissue)
# tissue_wide whole_blood 
#        5254       12342 


cpgs_eQTR<-bed_inter(a=cpgs[,start:=pos][,end:=pos+1][,.(chr,start,end,cpg_id)][order(chr,start)],
          b=eqtrs[,start:=start.eQTR][,end:=end.eQTR][,.(chr,start,end,eqtr_id)][order(chr,start)],
          select = c(4,8),col.names = c("cpg_id","eqtr_id"))

cpgs_eQTR #549k  cpgs - eQTR match
cpgs_eQTR<-merge(cpgs_eQTR,unique(eqtrs[,.(eqtr_id,gene,tissue,avg.mlog10pval)]),all.x=T)
unique(cpgs_eQTR,by="cpg_id")#322155 cpgs linked to a gene
unique(cpgs_eQTR,by="gene")#12k genes linked to a cpg
cpgs_eQTR[,in_eQTR:=T]
cpgs_eQTR[,in_both_eQTR:=all(c("whole_blood","tissue_wide")%in%tissue),by=.(cpg_id,gene)]
unique(cpgs_eQTR[in_both_eQTR==T],by="cpg_id") #(only) 4k cpg in both tissue
unique(cpgs_eQTR[in_both_eQTR==T],by="gene") #in (only) 406 genes
cpgs_eQTR<-unique(cpgs_eQTR)

#merge with coord
cpgs_eQTR<-merge(cpgs_eQTR,cpgs[,.(cpg_id,chr,pos)],by="cpg_id")
#add tss dist
  #need strand
gene_strands<-fread("ref/eQTL/GTEx_Analysis_v8_eQTL/Whole_Blood.v8.egenes.txt.gz",select = c(2,4:6),col.names = c("gene","start","end","strand"))
gene_strands[strand=="+",tss_pos:=start]
gene_strands[strand=="-",tss_pos:=end]
cpgs_eQTR<-merge(cpgs_eQTR,gene_strands[,.(gene,tss_pos,strand)],,by="gene",all.x = T)
cpgs_eQTR[strand=="+",tss_dist:=pos-tss_pos][strand=="-",tss_dist:=tss_pos-pos]


#merge the 2 links
cpgs_genes_1[,in_eQTR:=F]
cpgs_genes_1[,in_both_eQTR:=F]
cpgs_genes<-rbind(cpgs_genes_1[,-c("n.gene.cpg","closest.rank")],cpgs_eQTR,fill=T)
cpgs_genes[in_eQTR==T]

fwrite(cpgs_genes,fp(out,"all_cpgs_gene_links.csv.gz"))

#2)cpgs_regulatory region
#ensembl regulatory reg matching
cpgs_ensembl<-bed_inter(a=unique(cpgs,by="cpg_id")[,chr:=str_remove(chr,'chr')][,start:=pos][,end:=pos+1][,.(chr,start,end,cpg_id)][order(chr,start)],
          b="ref/ensembl_regulatory/ensembl_regulatory_hg19.bed",
          select = c(4,8),col.names = c("cpg_id","ensembl_regulatory_domain"))



cpgs_ensembl[,ensembl_regulatory_domain:=paste(ensembl_regulatory_domain,collapse = "/"),by="cpg_id"] 
cpgs_ensembl<-unique(cpgs_ensembl)
cpgs_ensembl

chrine_feat<-fread("ref/Chromatin_Annot/CBP/CD34_all_chromatin_feature.csv")

cpgs_chrine<-bed_inter(a=unique(cpgs,by="cpg_id")[,start:=pos][,end:=pos+1][,.(chr,start,end,cpg_id)][order(chr,start)],
          b=chrine_feat[,.(chr,start,end,type)][order(chr,start)],
          select = c(4,8),col.names = c("cpg_id","chromatin_feature"))

cpgs_reg<-merge(cpgs_chrine,cpgs_ensembl,all.x=T,by="cpg_id")

fwrite(cpgs_reg,fp(out,"cpgs_annot_ensembl_regulatory_domain_and_chromHMM_chromatin_features.csv.gz"))
cpgs_reg<-fread(fp(out,"cpgs_annot_ensembl_regulatory_domain_and_chromHMM_chromatin_features.csv.gz"))

cpgs_anno<-merge(cpgs_reg,cpgs_genes,all=T,by="cpg_id")

fwrite(unique(cpgs_anno[order(cpg_id,gene)][,.(cpg_id,chr,pos,gene,tss_pos,tss_dist,in_eQTR,eqtr_id,tissue,avg.mlog10pval,in_both_eQTR,strand,chromatin_feature,ensembl_regulatory_domain)]),
       fp(out,"cpgs_annot.csv.gz"))


