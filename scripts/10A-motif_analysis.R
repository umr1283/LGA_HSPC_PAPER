#Tf motif enrichmnent around DMC ?
source("scripts/utils/new_utils.R")
out<-"outputs/03B-motif_analysis"
dir.create(out)

#simple motif analysis on methyl : HOMER puis FIMO
# Finding Enriched Motifs in Genomic Regions 
#use: findMotifsGenome.pl <peak/BED file> <genome> <output directory> -size # [options]
#By default this will perform de novo motif discovery as well as check the enrichment of known motifs
#all parameter in http://homer.ucsd.edu/homer/ngs/peakMotifs.html

#LGA vs Ctrl DMC pval<0.001 and abs(meth.change)>30 with background
#run 03B1

res<-fread("outputs/03B-motif_analysis/knownResults.txt",
           select = c(1,2,3,5,6,7,8,9),
           col.names = c("motif","consensus","pval","padj","n_dmc_with_motif","pct_dmc_with_motif","n_background_with_motif","pct_background_with_motif"))
res[,pct_dmc_with_motif:=as.numeric(str_remove(pct_dmc_with_motif,"%"))]
res[,pct_background_with_motif:=as.numeric(str_remove(pct_background_with_motif,"%"))]
res[,motif:=str_remove(motif,"/Homer")]
res[,fold.enrichment:=pct_dmc_with_motif/pct_background_with_motif]
ggplot(res[padj<=0.2&  n_dmc_with_motif>30][order(pval)])+geom_point(aes(x=motif,col=-log10(pval),size=fold.enrichment,y=n_dmc_with_motif))+
  scale_x_discrete(limits=res[padj<=0.2&  n_dmc_with_motif>30][order(pval)]$motif)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 70, hjust=1))

#Some TF motif are similar to CCGG (our HELP-tagging enzymes cleaving site) 
#so perform also autobckgroun to see if enrichment in such TF
#run 03B2 



