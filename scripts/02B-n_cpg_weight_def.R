
#deter if ^1/4 for n_cpg_weight were the opti one to reduce influence of n_cpg_gene (werease kept influences of n_cpg_sig)
source("scripts/utils/new_utils.R")
out<-"outputs/02B-n_cpg_weight_def"
dir.create(out)
res_anno<-fread("outputs/02-gene_score_calculation_and_validation/res_anno.csv.gz")

res_test<-copy(res_anno)

#for gene_score
rs_ncpg<-rep(0,10)
rs_ncpg_sig<-rep(0,10)
ps<-list()
for(i in 1:10){
  print(i)
  res_test[,n.cpg_weight:=(1/sum(1/(abs(cpg_score)+1)))^(1/i),by="gene"] #n.cpg_weight to reduce the influence of the n_cpg by gene to the GeneScore
  res_test[,gene_score:=sum(cpg_score)*n.cpg_weight,by="gene"] 
  resg_test<-unique(res_test[order(gene,pval)],by="gene")
  print(summary(lm(gene_score~n.cpg.gene+n.cpg.sig.gene+pval+meth.change+chromatin_feature+ensembl_reg_score+in_eQTR+abs(tss_dist),data = resg_test))) 

  rs_ncpg[i]<-summary(lm(gene_score~n.cpg.gene,data = resg_test))$r.squared
  rs_ncpg_sig[i]<-summary(lm(gene_score~n.cpg.sig.gene,data = resg_test))$r.squared
  ps[[i]]<-ggplot(resg_test)+geom_boxplot(aes(x = as.factor(n.cpg.gene),y =gene_score_add )) +ggtitle(paste("x =",i))


}

plot(1:10,rs_ncpg)
plot(1:10,rs_ncpg_sig)
wrap_plots(ps)

rs_ncpg2<-rep(0,10)
rs_ncpg_sig2<-rep(0,10)
xs<-4+(1:10/10)
for(i in 1:length(xs)){
  x<-xs[i]
  print(x)
  res_test[,n.cpg_weight:=(1/sum(1/(abs(cpg_score)+1)))^(1/x),by="gene"] #n.cpg_weight to reduce the influence of the n_cpg by gene to the GeneScore
  res_test[,gene_score:=sum(cpg_score)*n.cpg_weight,by="gene"] 
  resg_test<-unique(res_test[order(gene,pval)],by="gene")
  print(summary(lm(gene_score~n.cpg.gene+n.cpg.sig.gene+pval+meth.change+chromatin_feature+ensembl_reg_score+in_eQTR+abs(tss_dist),data = resg_test))) 

  rs_ncpg2[i]<-summary(lm(gene_score~n.cpg.gene,data = resg_test))$r.squared
  rs_ncpg_sig2[i]<-summary(lm(gene_score~n.cpg.sig.gene,data = resg_test))$r.squared


}
plot(xs,rs_ncpg2)
plot(xs,rs_ncpg_sig2)

#gene_score_add
rs_ncpg<-rep(0,10)
rs_ncpg_sig<-rep(0,10)
ps<-list()
for(i in 1:10){
  print(i)
  res_test[,n_cpg_weight_region:=(1/sum(1/(abs(cpg_score)+1)))^(1/i),by=c('region_type',"gene")]
  res_test[region_type=="promoter",gene_score_region:=sum(cpg_score)*n_cpg_weight_region,by=c("gene")]
  res_test[region_type=="other",gene_score_region:=sum(abs(cpg_score))*n_cpg_weight_region,by=c("gene")]
  
  res_test[,gene_score_add:=sum(unique(abs(gene_score_region)),na.rm = T),by="gene"]
  resg_test<-unique(res_test[order(gene,region_type,pval)],by="gene")
  n_cpg_max<-10^max(boxplot.stats(log10(resg_test$n.cpg.gene))$stats)
  resg_test<-resg_test[n.cpg.gene<=n_cpg_max]
  print(summary(lm(gene_score_add~n.cpg.gene+n.cpg.sig.gene+pval+meth.change+chromatin_feature+ensembl_reg_score+in_eQTR+abs(tss_dist),data = resg_test))) 

  rs_ncpg[i]<-summary(lm(gene_score_add~n.cpg.gene,data = resg_test))$r.squared
  rs_ncpg_sig[i]<-summary(lm(gene_score_add~n.cpg.sig.gene,data = resg_test))$r.squared
  ps[[i]]<-ggplot(resg_test)+geom_boxplot(aes(x = as.factor(n.cpg.gene),y =gene_score_add )) +ggtitle(paste("x =",i))


}

plot(1:10,rs_ncpg)
plot(1:10,rs_ncpg_sig)
wrap_plots(ps)

rs_ncpg2<-rep(0,10)
rs_ncpg_sig2<-rep(0,10)
xs<-3+(1:10/10)
for(i in 1:length(xs)){
  x<-xs[i]
  print(x)
  res_test[,n_cpg_weight_region:=(1/sum(1/(abs(cpg_score)+1)))^(1/x),by=c('region_type',"gene")]
  res_test[region_type=="promoter",gene_score_region:=sum(cpg_score)*n_cpg_weight_region,by=c("gene")]
  res_test[region_type=="other",gene_score_region:=sum(abs(cpg_score))*n_cpg_weight_region,by=c("gene")]
  
  res_test[,gene_score_add:=sum(unique(abs(gene_score_region)),na.rm = T),by="gene"]
  resg_test<-unique(res_test[order(gene,region_type,pval)],by="gene")
  print(summary(lm(gene_score_add~n.cpg.gene+n.cpg.sig.gene+pval+meth.change+chromatin_feature+ensembl_reg_score+in_eQTR+abs(tss_dist),data = resg_test))) 

  rs_ncpg2[i]<-summary(lm(gene_score_add~n.cpg.gene,data = resg_test))$r.squared
  rs_ncpg_sig2[i]<-summary(lm(gene_score_add~n.cpg.sig.gene,data = resg_test))$r.squared


}
plot(xs,rs_ncpg2)
plot(xs,rs_ncpg_sig2)
#3.4 is best : 
# [1] 3.4
# 
# Call:
# lm(formula = gene_score_add ~ n.cpg.gene + n.cpg.sig.gene + pval + 
#     meth.change + chromatin_feature + ensembl_reg_score + in_eQTR + abs(tss_dist), 
#     data = resg_test)
# 
# Residuals:
#     Min      1Q  Median      3Q     Max 
# -421.86  -14.97   -6.11   11.93  202.40 
# 
# Coefficients:
#                 Estimate Std. Error t value Pr(>|t|)    
# (Intercept)    1.108e+01  6.097e-01  18.176  < 2e-16 ***
# n.cpg.gene    -9.621e-04  8.930e-03  -0.108   0.9142    
# n.cpg.sig.gene  2.334e+01  1.250e-01 186.726  < 2e-16 ***
# pval          -8.026e+00  1.084e+00  -7.405 1.37e-13 ***
# meth.change    6.203e-02  1.142e-02   5.434 5.58e-08 ***
# chromatin_feature           7.337e-01  9.591e-02   7.650 2.09e-14 ***
# ensembl_reg_score    1.452e+01  1.074e+00  13.520  < 2e-16 ***
# in_eQTRTRUE   -1.304e+00  5.737e-01  -2.273   0.0231 *  
# abs(tss_dist) -4.914e-06  1.160e-06  -4.237 2.28e-05 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 27.83 on 21189 degrees of freedom
#   (54 observations deleted due to missingness)
# Multiple R-squared:  0.728,	Adjusted R-squared:  0.7279 
# F-statistic:  7090 on 8 and 21189 DF,  p-value: < 2.2e-16


xs<-1:50/10
ps_cov<-sapply(xs, function(x){
  print(x)
  res_test[,n_cpg_weight_region:=(1/sum(1/(abs(cpg_score)+1)))^(1/x),by=c('region_type',"gene")]
  res_test[region_type=="promoter",gene_score_region:=sum(cpg_score)*n_cpg_weight_region,by=c("gene")]
  res_test[region_type=="other",gene_score_region:=sum(abs(cpg_score))*n_cpg_weight_region,by=c("gene")]
  res_test[,gene_score_add:=sum(unique(abs(gene_score_region)),na.rm = T),by="gene"]
  resg_test<-unique(res_test[order(gene,pval)],by="gene")
  res<-summary(lm(gene_score_add~n.cpg.gene+n.cpg.sig.gene+pval+meth.change+chromatin_feature+ensembl_reg_score+in_eQTR+abs(tss_dist),data = resg_test)) 
  return(res$coefficients[,4])

})
colnames(ps_cov)<-paste("x",xs,sep="=")
ps_cov<-data.table(ps_cov,keep.rownames = "covariates")
ps_cov<-melt(ps_cov,variable.name = "x",value.name = "p_val")
ps_cov[,x:=as.numeric(str_remove(x,"x="))]
fwrite(ps_cov,fp(out,"signif_covariates_effect_on_gene_score.csv"))
ggplot(ps_cov[x==3.9&covariates!="(Intercept)"])+geom_col(aes(x=covariates,y=-log10(p_val)))+scale_y_log10()

ggplot(ps_cov[covariates=="n.cpg.gene"])+geom_point(aes(x=x,y=-log10(p_val)))+scale_y_log10()

plot(xs,rs_ncpg_sig2)
