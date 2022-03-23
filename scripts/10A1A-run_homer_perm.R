#random background 

library(data.table)
library(stringr)

fp<-function(...)file.path(...)

set.seed(1234)
old_path <- Sys.getenv("PATH")
Sys.setenv(PATH = paste(old_path, "/home/apelletier/HSPC_EpiStress/Data/tools/homer/bin/", sep = ":"))

out0<-"outputs/03B-motif_analysis"
dir.create(out0)

#bed_path<-fp(out,"40bp_win_dmc_pvalnom0.001_meth.change30")

res<-fread("outputs/01-lga_vs_ctrl_limma_DMCs_analysis/res_limma.tsv.gz")
res[P.Value<0.001&abs(logFC)>25] 
cpgs<-fread("ref/all_HELP_tagging_cpgs_hg19_pos.csv")

res<-merge(res,cpgs)
res[,start:=pos-20][start<1,start:=1][,end:=pos+20]
res<-res[!is.na(chr)&!is.na(pos)]

n_cpgs<-nrow(res[P.Value<0.001&abs(logFC)>25])

n_perm=1000

for(i in 1:n_perm){
  message("perm ",i,"/",n_perm)
  out<-fp(out0,paste0("PERM",i))
  dir.create(out)
  
  bed_path<-fp(out,paste0("random",i,".bed"))
  
  random_cpgs<-res[cpg_id%in%sample(cpg_id,n_cpgs),.(chr,start,end)][order(chr,start)]

  fwrite(random_cpgs[,.(chr,start,end)][order(chr,start)],
       bed_path,sep="\t",
       col.names=F)
  
  cmd<-paste("findMotifsGenome.pl",bed_path,"hg19",out,"-size 41 -nomotif -bg outputs/03B-motif_analysis/random1_HOMER.bed -noweight -nlen 0 -p 4")
  system(cmd)
  
  
  
  try({
    
    known<-fread(fp(out,"knownResults.txt"),
           select = c(1,3,5,6,7,8,9),
           col.names = c("motif","pval","padj","n_dmc_with_motif","pct_dmc_with_motif","n_background_with_motif","pct_background_with_motif"))
    known[,pct_dmc_with_motif:=as.numeric(str_remove(pct_dmc_with_motif,"%"))]
    known[,pct_background_with_motif:=as.numeric(str_remove(pct_background_with_motif,"%"))]
    known[,permut:=i]
     if(i==1){
    res_all<-copy(known)
  }else{
    res_all<-rbind(res_all,known)
  }
  fwrite(res_all,fp(out0,"res_known_motif_all_perm_bgrandom.csv"))
    })
  
  if(i>10){
    message("removing res")
    unlink(out,recursive = T)
  }
  
  
}

fwrite(res_all,fp(out0,"res_known_motif_all_perm_bgrandom.csv"))
message("Success !")
