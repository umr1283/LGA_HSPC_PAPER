
#Gene score calculation and validation
source("scripts/utils/new_utils.R")
source("scripts/utils/methyl_utils.R")
out<-"outputs/02-gene_score_calculation_and_validation"
dir.create(out)


res<-fread("outputs/01-lga_vs_ctrl_limma_DMCs_analysis/res_limma.tsv.gz",sep="\t",
           select = c("cpg_id","P.Value","adj.P.Val","AveExpr","logFC"),
           col.names = c("cpg_id","pval","padj","avg.meth","meth.change"))
mtd<-fread("datasets/cd34/metadata_pcs_cl_190421.csv",sep=";")

#link cpg to gene and weitgh link confidence
# 1)annots cpg ====
#see 02A-

# 2) calculate cpg weight====
cpgs_anno<-fread("outputs/02A-CpGs_annotations/cpgs_annot.csv.gz")
cpgs_genes<-cpgs_anno[!is.na(gene)&gene!=""]
rm(cpgs_anno)
unique(cpgs_genes,by="cpg_id") #1,2M/~1.7M cpgs link to a gene

unique(cpgs_genes[in_eQTR==T],by="cpg_id") #dont 320k linked thx to eQTR
cpgs_genes[,double.linked:=any(in_eQTR==T)&any(in_eQTR==F),by=c("cpg_id","gene")]
unique(cpgs_genes[double.linked==T],by="cpg_id")#dont 67k also gene linked based on tss
#=> + ~250k cpgs gene link thx to eQTR

# a) linksWeight
cpgs_genes[in_eQTR==F,links_score:=sapply(abs(tss_dist),function(x){
  if(x<1000)return(1)
  else if(x<20000)return(0.5+0.5*sqrt(1000/x))
  else return(0.5*sqrt(20000/x))
    })]

summary(cpgs_genes[in_eQTR==F]$links_score)
 #   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
 # 0.1581  0.3156  0.6172  0.5612  0.7369  1.0000 

summary(unique(cpgs_genes[in_eQTR==T&tissue=="whole_blood"],by="eqtr_id")$avg.mlog10pval)
  # Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  # 3.255   7.455  12.600  20.560  23.741 277.112 

summary(cpgs_genes[tissue=="tissue_wide"]$avg.mlog10pval) #not unique because crash idk why
  #  Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  # 29.02   53.02   73.91   79.29  103.04  191.19

cpgs_genes[in_eQTR==T,links_score:=ifelse(avg.mlog10pval>quantile(avg.mlog10pval,0.9),
                                          1,
                                          avg.mlog10pval/quantile(avg.mlog10pval,0.9))
           ,by="tissue"]

summary(cpgs_genes[in_eQTR==T&tissue=="whole_blood"]$links_score)
summary(cpgs_genes[in_eQTR==T&tissue=="tissue_wide"]$links_score)
cpgs_genes[,links_weight:=max(links_score),by=c("cpg_id","gene")]


#b) Regulatory Weights
cpgs_reg<-unique(cpgs_genes,by=c("cpg_id"))
unique(cpgs_reg$chromatin_feature)
cpgs_reg[is.na(chromatin_feature)]
cpgs_reg[,chromatin_score:=sapply(chromatin_feature,function(x){
      if(is.na(x))return(0)
      else if(x%in%c(6,4))return(1)
      else if(x==5)return(0.75)
      else if(x%in%c(1,2,3))return(0.5)
      else return(0)
    })]


unique(cpgs_reg$chromatin_feature)
cpgs_reg[,ensembl_reg_score:=sapply(ensembl_regulatory_domain,function(x){
    vecX<-strsplit(as.character(x),"/")[[1]]
    if(any(c("CTCF_binding_site","promoter","enhancer")%in%vecX)){
      score<-0.5
    }else if(any(c("open_chromatin_region","promoter_flanking_region")%in%vecX)){
      score<-0.25
    }else{
      score<-0
    }
    if("TF_binding_site"%in%vecX){
      score<-score+0.5
    }
    return(score)
  })]

cpgs_reg[,regul_weight:=(0.5+1.5*((chromatin_score+ensembl_reg_score)/2))]

cpgs_score<-merge(cpgs_genes,cpgs_reg[,.(cpg_id,chromatin_score,ensembl_reg_score,regul_weight)],by="cpg_id")
fwrite(cpgs_score,fp(out,"cpgs_genes_annot_and_weight.csv.gz"))


res_anno<-merge(res,cpgs_score,by="cpg_id")

#GeneScore calculation====
# pval
res_anno[,min.pval:=min(pval[which(in_eQTR==F)],na.rm = T),by=c("gene")]
res_anno[,avg.pval:=mean(-log10(pval[in_eQTR==F]),na.rm = T),by=c("gene")]
res_anno[,avg.m.log10.pval:=mean(-log10(pval[in_eQTR==F]),na.rm = T),by=c("gene")]
meth_metrics<-c("min.pval","avg.pval","avg.m.log10.pval")

# dmcscore
res_anno[,max.dmc_score:=max(-log10(pval[which(in_eQTR==F)])*abs(meth.change[which(in_eQTR==F)]),na.rm = T),by=c("gene")]
res_anno[,avg.dmc_score:=mean(-log10(pval[in_eQTR==F])*abs(meth.change[in_eQTR==F]),na.rm = T),by=c("gene")]
meth_metrics<-c(meth_metrics,"max.dmc_score","avg.dmc_score")

# +LinksScore
res_anno[,max.dmc_score.tss_pond:=max(-log10(pval[which(in_eQTR==F)])*abs(meth.change[which(in_eQTR==F)])*links_score[which(in_eQTR==F)],na.rm = T),by=c("gene")]
res_anno[,avg.dmc_score.tss_pond:=sum(-log10(pval[in_eQTR==F])*abs(meth.change[in_eQTR==F])*links_score[in_eQTR==F],na.rm = T)/sum(links_score[in_eQTR==F]),by=c("gene")]
meth_metrics<-c(meth_metrics,"max.dmc_score.tss_pond","avg.dmc_score.tss_pond")

# +eQTL
res_anno[,max.dmc_score.tss_eqtl_pond:=max(-log10(pval)*abs(meth.change)*links_weight,na.rm = T),by=c("gene")]
res_anno[,avg.dmc_score.tss_eqtl_pond:=sum(-log10(pval)*abs(meth.change)*links_weight,na.rm = T)/sum(links_weight),by=c("gene")]
meth_metrics<-c(meth_metrics,"max.dmc_score.tss_eqtl_pond","avg.dmc_score.tss_eqtl_pond")


# +annot_ensembl
res_anno[,max.dmc_score.tss_ens_pond:=max(-log10(pval)*abs(meth.change)*links_weight*ensembl_reg_score,na.rm = T),by=c("gene")]
res_anno[,avg.dmc_score.tss_ens_pond:=sum(-log10(pval)*abs(meth.change)*links_weight*ensembl_reg_score,na.rm = T)/sum(links_weight*ensembl_reg_score),by=c("gene")]
meth_metrics<-c(meth_metrics,"max.dmc_score.tss_ens_pond","avg.dmc_score.tss_ens_pond")


# +chrine
res_anno[,max.dmc_score.tss_ens_chrine_pond:=max(-log10(pval)*abs(meth.change)*links_weight*regul_weight),by=c("gene")]
res_anno[,avg.dmc_score.tss_ens_chrine_pond:=sum(-log10(pval)*abs(meth.change)*links_weight*regul_weight,na.rm = T)/sum(links_weight*regul_weight),by=c("gene")]
meth_metrics<-c(meth_metrics,"max.dmc_score.tss_ens_chrine_pond","avg.dmc_score.tss_ens_chrine_pond")


# +nCpGweight
res_anno[,cpg_score:=-log10(pval)*meth.change*links_weight*regul_weight] #divided by n_sample ro normalized gene score ~ n_sampleres_anno[,n.cpg_weight:=(1/sum(1/(abs(cpg_score)+1)))^(1/5),by="gene"] #n.cpg_weight to reduce the influence of the n_cpg by gene to the GeneScore
res_anno[,n.cpg_weight:=(1/sum(1/(abs(cpg_score)+1)))^(1/5),by="gene"] #n.cpg_weight to reduce the influence of the n_cpg by gene to the GeneScore
res_anno[,gene_score:=sum(cpg_score)*n.cpg_weight,by="gene"] 
meth_metrics<-c(meth_metrics,"gene_score")


# +sep prom/enh
res_anno[,region_type:=ifelse(abs(tss_dist)<=2000,"promoter","other")]
res_anno[is.na(tss_dist),region_type:="other"]
res_anno[is.na(region_type)]

res_anno[,n_cpg_weight_region:=(1/sum(1/(abs(cpg_score)+1)))^(1/3.8),by=c('region_type',"gene")]
res_anno[region_type=="promoter",gene_score_region:=sum(cpg_score)*n_cpg_weight_region,by=c("gene")]
res_anno[region_type=="other",gene_score_region:=sum(abs(cpg_score))*n_cpg_weight_region,by=c("gene")]

res_anno[,gene_score_add:=sum(unique(abs(gene_score_region)),na.rm = T),by="gene"]
meth_metrics<-c(meth_metrics,"gene_score_add")

fwrite(res_anno,"outputs/02-gene_score_calculation_and_validation/res_anno.csv.gz")
res_anno<-fread("outputs/02-gene_score_calculation_and_validation/res_anno.csv.gz")

#Rq : to see how we define the weight see 02B :

#valid that not correlated to ncpg.gene, but crrelated to n cpg sig :
res_anno[,n.cpg.gene:=.N,by=.(gene)]
res_anno[,n.cpg.sig.gene:=sum(pval<0.01),by=.(gene)]
resg<-unique(res_anno[order(gene,pval)],by=c("gene"))

ggplot(resg)+
      geom_boxplot(aes(x = as.factor(n.cpg.gene),y =gene_score_add )) 


#VALIDATION  Gene Score====
# 1)valid wieight
#see correl covariates with genescore
summary(lm(gene_score_add~n.cpg.gene+n.cpg.sig.gene+pval+meth.change+chromatin_feature+ensembl_reg_score+in_eQTR+abs(tss_dist),data = resg)) #best gene_score_add than gene_score


# 2)valid expression change prediction
meth_scores<-melt(resg[,.SD,.SDcols=c("gene",meth_metrics)],id.vars = "gene",variable.name = "meth_metric",value.name = "score")
meth_scores[score==Inf,score:=NA]
meth_scores[score==-Inf,score:=NA]
meth_scores<-meth_scores[!is.na(score)]
res_de_cl<-fread("../singlecell/outputs/04-DEG_in_LGA/2020-09-01_pseudo_bulk_DEseq2_LgaVsCtrl_CBP1andcbp558_559_samples_excluded_regr_on_batch_and_sex_all_genes.csv")

res_de_cl<-res_de_cl[!is.na(padj)]
meth_scores_de<-merge(meth_scores,res_de_cl,by=c("gene"))
unique(meth_scores_de[padj<0.1]$gene)#120 DEGS

meth_scores_de[,score_scaled:=scale(score),by="meth_metric"]
ggplot(meth_scores_de)+
  geom_boxplot(aes(fill=padj<0.1,y=score_scaled,x=meth_metric),outlier.shape = NA)+
  coord_cartesian(ylim = c(-2.5,3))+
  theme(axis.text.x = element_text(angle = 45,hjust = 1))


meth_scores_de[,p_diff:=wilcox.test(score[padj<0.1],score[padj>=0.1])$p.value,by="meth_metric"]

ggplot(unique(meth_scores_de,by="meth_metric"))+geom_col(aes(y=-log10(p_diff),x=meth_metric))+
  theme(axis.text.x = element_text(angle = 45,hjust = 1))

meth_scores_de[,degs:=padj<0.1]
meth_scores_de[,q75_dt:=quantile(score_scaled[padj<0.1],0.75,na.rm=T)-quantile(score_scaled[padj>=0.1],0.75,na.rm=T),by=.(meth_metric)]

ggplot(unique(meth_scores_de,by=c("meth_metric")))+
  geom_col(aes(y=q75_dt,x=meth_metric))+
    theme(axis.text.x = element_text(angle = 45,hjust = 1))


meth_metrics<-c("mlog10.pval.max","avg.mlog10.pval","meth.change.max","avg.meth.change","gene_score","gene_score_add")

res_anno[,meth.change.max:=max(meth.change),by=c("gene")]
res_anno[,avg.meth.change:=mean(meth.change),by=c("gene")]
res_anno[,mlog10.pval.max:=max(-log10(pval)),by=c("gene")]
res_anno[,avg.mlog10.pval:=mean(-log10(pval)),by=c("gene")]
resg<-unique(res_anno[order(gene,pval)],by=c("gene"))
meth_scores2<-melt(resg[,.SD,.SDcols=c("gene",meth_metrics)],id.vars = "gene",variable.name = "meth_metric",value.name = "score")

meth_scores_de2<-merge(meth_scores2,res_de_cl,by=c("gene"))
meth_scores_de2[,score_scaled:=scale(score),by="meth_metric"]

ggplot(meth_scores_de2)+
  geom_boxplot(aes(fill=padj<0.1,y=score_scaled,x=meth_metric),outlier.shape = NA)+
  coord_cartesian(ylim = c(-2.5,3))+
  theme(axis.text.x = element_text(angle = 45,hjust = 1))
  
meth_scores_de2[,p_diff:=wilcox.test(score[padj<0.1],score[padj>=0.1])$p.value,by="meth_metric"]
unique(meth_scores_de2,by="meth_metric")
ggplot(unique(meth_scores_de2,by="meth_metric"))+geom_col(aes(y=-log10(p_diff),x=meth_metric))


meth_metrics<-c("avg.mlog10.pval","avg.meth.change","gene_score_add")

meth_scores3<-melt(resg[,.SD,.SDcols=c("gene",meth_metrics)],id.vars = "gene",variable.name = "meth_metric",value.name = "score")
meth_scores3[meth_metric=="gene_score_add",meth_metric:="gene_score"]
meth_scores_de3<-merge(meth_scores3,res_de_cl,by=c("gene"))
meth_scores_de3[,score_scaled:=scale(score),by="meth_metric"]

ggplot(meth_scores_de3)+
  geom_boxplot(aes(fill=padj<0.1,y=score_scaled,x=meth_metric),outlier.shape = NA)+
  coord_cartesian(ylim = c(-3,3))+
  theme_classic()+
  scale_fill_manual(values = c("white","grey"))
  
meth_scores_de3[,p_diff:=wilcox.test(score[padj<0.1],score[padj>=0.1])$p.value,by="meth_metric"]
unique(meth_scores_de3,by="meth_metric")
ggplot(unique(meth_scores_de,by="meth_metric"))+geom_col(aes(y=-log10(p_diff),x=meth_metric))

#=> on voit que les degs ont un gene_score + élevé que les autres genes alors quils ont les meme pvaleurs et methchange min
#par contre max(-log10(padj))  semble suffire pour predire gene expr change. Really ?
#apres permut est ce tjr signif ?

meth_metrics<-c("max.dmc_score","gene_score_add")


meth_scores4<-melt(resg[,.SD,.SDcols=c("gene",meth_metrics)],id.vars = "gene",variable.name = "meth_metric",value.name = "score")
meth_scores4[meth_metric=="max.dmc_score",meth_metric:="max(-log10(pval)*meth.change)"]
meth_scores4[meth_metric=="gene_score_add",meth_metric:="gene_score"]

meth_scores4[score==Inf,score:=NA]
meth_scores4[score==-Inf,score:=NA]

meth_scores4<-meth_scores4[!is.na(score)]
meth_scores_de4<-merge(meth_scores4,res_de_cl,by=c("gene"))
meth_scores_de4[,score_scaled:=scale(score),by="meth_metric"]

ggplot(meth_scores_de4)+
  geom_boxplot(aes(fill=padj<0.1,y=score_scaled,x=meth_metric),outlier.shape = NA)+
  coord_cartesian(ylim = c(-3,3))+
  theme_classic()+
  scale_fill_manual(values = c("white","grey"))
  
meth_scores_de4[,p_diff:=wilcox.test(score[padj<0.1],score[padj>=0.1])$p.value,by="meth_metric"]
    
meth_scores_de4[,q75_dt:=quantile(score_scaled[padj<0.1],0.75,na.rm=T)-quantile(score_scaled[padj>=0.1],0.75,na.rm=T),by=.(meth_metric)]

res_pred<-unique(meth_scores_de4,by="meth_metric")


res_pred<-res_pred[,.(meth_metric,p_diff,)]

pvals_perms<-sapply(1:100, function(i){
  print(i)

  mtd_p<-copy(mtd_f)
  mtd_p[,group:=sample(group)]
  print(paste(sum(mtd_p$group==mtd_f$group),"matchs/",nrow(mtd_f)))
  mtd_p[,group_sex:=paste(group,sex,sep="_")]
  design<-model.matrix(formule,data = data.frame(mtd_p,row.names = "sample"))
  fit <- lmFit(data.frame(methf,row.names = "cpg_id")[,mtd_p$sample], design)
  cont.matrix <- makeContrasts(C.L = "(group_sexCTRL_F+group_sexCTRL_M)-(group_sexLGA_F+group_sexLGA_M)",
                             levels=design)

  fit2  <- contrasts.fit(fit, cont.matrix)
  fit2  <- eBayes(fit2)
  resp<-data.table(topTable(fit2,coef = "C.L",n = Inf),keep.rownames = "cpg_id")
  resp[,cpg_id:=as.numeric(cpg_id)]
  setnames(resp,c("P.Value","adj.P.Val","AveExpr","logFC"),c("pval","padj","avg.meth","meth.change"))
  print(table(resp[padj<0.01]$meth.change>0))
  resp<-merge(resp,cpgs_score,by="cpg_id")

  resp[,max.dmc_score:=max(-log10(pval[which(in_eQTR==F)])*abs(meth.change[which(in_eQTR==F)]),na.rm = T),by=c("gene")]
  resp[,cpg_score:=-log10(pval)*meth.change*links_weight*regul_weight]
  resp[,region_type:=ifelse(abs(tss_dist)<=2000,"promoter","other")]
  resp[is.na(tss_dist),region_type:="other"]

  resp[,n_cpg_weight_region:=(1/sum(1/(abs(cpg_score)+1)))^(1/3.9),by=c('region_type',"gene")]
  resp[region_type=="promoter",gene_score_region:=sum(cpg_score)*n_cpg_weight_region,by=c("gene")]
  resp[region_type=="other",gene_score_region:=sum(abs(cpg_score))*n_cpg_weight_region,by=c("gene")]
  resp[,gene_score:=sum(unique(abs(gene_score_region)),na.rm = T),by="gene"]

  meth_metrics<-c("max.dmc_score","gene_score")
  resgp<-unique(resp[order(gene,pval)],by=c("gene"))

  meth_scores4p<-melt(resgp[,.SD,.SDcols=c("gene",meth_metrics)],id.vars = "gene",variable.name = "meth_metric",value.name = "score")
  meth_scores4p[score==Inf,score:=NA]
  meth_scores4p[score==-Inf,score:=NA]
  meth_scores4p<-meth_scores4p[!is.na(score)]
  meth_scores_de4p<-merge(meth_scores4p,res_de_cl,by=c("gene"))
  meth_scores_de4p[,score_scaled:=scale(score),by="meth_metric"]
  meth_scores_de4p[,p_diff:=wilcox.test(score[padj<0.1],score[padj>=0.1])$p.value,by="meth_metric"]
  meth_scores_de4p[,q75_dt:=quantile(score_scaled[padj<0.1],0.75,na.rm=T)-quantile(score_scaled[padj>=0.1],0.75,na.rm=T),by=.(meth_metric)]
  res_predp<-unique(meth_scores_de4p,by="meth_metric")

  return(list(p_diff=res_predp$p_diff,
              q75_dt=res_predp$q75_dt))

  })
pvals_perms
p_perm<-Reduce(rbind,pvals_perms[1,])
q_perm<-Reduce(rbind,pvals_perms[2,])

# max(dmc_score)
sum(p_perm[,1]<0.00192)/100 #p_perm p = 0.48
sum(q_perm[,1]>0.448)/100 #p_perm q75 = 0.1


# gene_score
sum(p_perm[,2]<0.0007759)/100 #p_perm p = 0.19
sum(q_perm[,2]>0.5135)/100 #p_perm q75 = 0.03

#save res_by_gene
res_anno[,n.cpg.in.eQTR:=sum(in_eQTR),by="gene"]
res_anno[,gene_score_prom:=gene_score_region[region_type=="promoter"][1],by="gene"]
res_anno[,gene_score_enh:=gene_score_region[region_type=="other"][1],by="gene"]

resg<-unique(res_anno[order(gene,pval)],by=c("gene"))
fwrite(resg[,.(gene,chr,tss_pos,gene_score_add,gene_score_prom,gene_score_enh,n.cpg.gene,n.cpg.sig.gene,n.cpg.in.eQTR,cpg_id,pos,tss_dist,pval,padj,meth.change)],"outputs/02-gene_score_calculation_and_validation/res_genes.csv.gz")
resg[,gene_score_add_scaled:=scale(gene_score_add)]
#plot dmgs plot
ggplot(resg,aes(x=gene_score_add,y=-log10(pval)))+
  geom_point()+
  theme_bw()

ggplot(resg,aes(x=gene_score_add,y=-log10(pval),col=padj<0.1&gene_score_add>300))+
  geom_point()+
  scale_color_manual(values = c("grey","red"))+theme_minimal()+
  geom_label_repel(aes(label = ifelse(gene_score_add>1500&padj<0.05,gene,"")),
                   max.overlaps = 3000,
                   box.padding   = 0.35, 
                   point.padding = 0.5,
                   segment.color = 'grey50')

###PERM GENESCORE====
library(limma)
library(parallel)
source("scripts/utils/new_utils.R")
set.seed(1234)
methf<-fread("datasets/cd34/meth_data_filtered.csv.gz")

mtd_f<-fread("datasets/cd34/metadata_cl_pcs_040521.csv")
cpgs_score<-fread("outputs/02-gene_score_calculation_and_validation/cpgs_genes_annot_and_weight.csv.gz")

res_perm<-Reduce(rbind,mclapply( 1:1000,function(i){

  mtd_f$group_sex<-sample(mtd_f$group_sex)
  
  #limma
  design<-model.matrix(~0 + group_sex   + batch+ group_complexity_fac +mat.age  + latino + PC2,
                       data = data.frame(mtd_f,row.names = "sample"))
  
  fit <- lmFit(data.frame(methf,row.names = "cpg_id")[,mtd_f$sample], design)
  
  cont.matrix <- makeContrasts(C.L = "(group_sexCTRL_F+group_sexCTRL_M)-(group_sexLGA_F+group_sexLGA_M)",
                               levels=design)
  fit2  <- contrasts.fit(fit, cont.matrix)
  fit2  <- eBayes(fit2)
  
  res<-data.table(topTable(fit2,coef = "C.L",n = Inf),keep.rownames = "cpg_id")
  res[,cpg_id:=as.numeric(cpg_id)]
  cols1<-c("cpg_id","P.Value","adj.P.Val","AveExpr","logFC")
  cols2<-c("cpg_id","pval","padj","avg.meth","meth.change")
  res<-res[,(cols2):=.SD,.SDcols=cols1][,.SD,.SDcols=cols2]
  
  res_anno<-merge(res,cpgs_score,by="cpg_id")

  #genescore
  res_anno[,cpg_score:=-log10(pval)*meth.change*links_weight*regul_weight] #divided by n_sample ro normalized gene score ~ n_sampleres_anno[,n.cpg_weight:=(1/sum(1/(abs(cpg_score)+1)))^(1/5),by="gene"] #n.cpg_weight to reduce the influence of the n_cpg by gene to the GeneScore
  
  res_anno[,region_type:=ifelse(abs(tss_dist)<=2000,"promoter","other")]
  res_anno[is.na(tss_dist),region_type:="other"]
  res_anno[,n_cpg_weight_region:=(1/sum(1/(abs(cpg_score)+1)))^(1/3.8),by=c('region_type',"gene")]
  
  res_anno[region_type=="promoter",gene_score_region:=sum(cpg_score)*n_cpg_weight_region,by=c("gene")]
  res_anno[region_type=="other",gene_score_region:=sum(abs(cpg_score))*n_cpg_weight_region,by=c("gene")]
  res_anno[,gene_score_add:=sum(unique(abs(gene_score_region)),na.rm = T),by="gene"]
  
  return(unique(res_anno,by="gene")[,perm:=i][order(gene)][,.(gene,gene_score_add,perm)])
  
  
},mc.cores = 20))

fwrite(res_perm,"outputs/02-gene_score_calculation_and_validation/res_1000perm_genescore_add.csv.gz")

res_perm<-fread("outputs/02-gene_score_calculation_and_validation/res_1000perm_genescore_add.csv.gz")
resg<-fread("outputs/02-gene_score_calculation_and_validation/res_genes.csv.gz")
resg[,perm:=0]
resg_perm<-rbind(resg[,.(gene,gene_score_add,perm)],res_perm)
resg_perm[,pval_gs:=sum(gene_score_add[perm==0]<=gene_score_add[perm!=0])/1000,by="gene"]
resgp<-resg_perm[perm==0][,-"perm"]
resg2<-merge(resg[,-"perm"],resgp[,.(gene,pval_gs)])
fwrite(resg2[,.SD,.SDcols=c(1,4,ncol(resg2),2:3,5:(ncol(resg2)-1))],"outputs/02-gene_score_calculation_and_validation/res_genes.csv.gz")

resg2[gene_score_add>1200&pval_gs<0.01]$gene
plot(density(resg$gene_score_add))
abline(v=200)
abline(v=quantile(resg2$gene_score_add,0.90))
resg2[gene%in%c("SOCS3","HES1")]
#NEXT, PAthway analysis, see 03-
