library(Seurat)
source("scripts/utils/new_utils.R")
out<-"12-Pseudotime"
dir.create(out)

#Pseudotime analysis on ref umap (hematomap) root= LT-HSC
#run 12A

cbps.cds<-readRDS(fp(out,"cbps.cds.rds"))
plot_cells(
  cds = cbps.cds,
  color_cells_by = "pseudotime",
  show_trajectory_graph = TRUE,
  label_branch_points= TRUE,
  label_leaves=FALSE
)

cbps<-readRDS(file="outputs/12-Pseudotime/cbps_RNA_Pseudotime_ComputRoot.rds")

cbps_inf<-subset(cbps,pseudo_time!=Inf)
Idents(cbps_inf)<-"cell_type_hmap"
FeaturePlot(cbps_inf,"pseudotime",reduction = 'umap',label=T)

ggplot(mtd)+geom_boxplot(aes(y=pseudotime,x=lineage_hmap),alpha=0.7)+theme_minimal()


ggplot(mtd)+geom_density(aes(x=pseudotime,fill=group,col=group),alpha=0.7)+facet_wrap("hto")+theme_minimal()



#increase lo pseudotime in lga, due to lo gene cells (allegement du seuil de filtration) ?

Idents(cbps_hm)<-"lineage_hmap"
p1<-FeaturePlot(subset(cbps_hm,Pseudotime_ALL_LTHSC_CTRLonly!=Inf), "Pseudotime_ALL_LTHSC_CTRLonly",label=T,reduction = "ref.umap")
p2<-DimPlot(subset(cbps_hm,Pseudotime_ALL_LTHSC_CTRLonly!=Inf), group.by="lineage_hmap",label=T,reduction = "ref.umap")
p1+p2

mtd<-data.table(cbps_hm@meta.data,keep.rownames = "bc")



unique(mtd$batch)
ggplot(mtd)+geom_density(aes(x=Pseudotime_ALL_LTHSC_CTRLonly,fill=group,col=group),alpha=0.7)+
  facet_wrap("hto")+theme_minimal()

mtd[,n.sample:=.N,"sample_hto"]

mtd[,pct.lin:=.N/n.sample,c("sample_hto","lineage_hmap")]

ggplot(unique(mtd,by=c("sample_hto","lineage_hmap")))+geom_boxplot(aes(x=hto,y=pct.lin,fill=group))+facet_wrap("lineage_hmap")+theme_minimal()

ggplot(mtd[Pseudotime_ALL_LTHSC_CTRLonly<10])+geom_boxplot(aes(x=hto,y=nFeature_RNA,fill=group))

ggplot(mtd[Pseudotime_ALL_LTHSC_CTRLonly<10])+geom_boxplot(aes(x=hto,y=percent.mt,fill=group))


ggplot(mtd[Pseudotime_ALL_LTHSC_CTRLonly<10])+geom_bar(aes(x=hto,fill=group),position = "dodge")

ggplot(mtd[Pseudotime_ALL_LTHSC_CTRLonly<10])+geom_bar(aes(x=hto,fill=group),position = "dodge")+facet_wrap("lineage_hmap")

ggplot(mtd[Pseudotime_ALL_LTHSC_CTRLonly>10])+geom_bar(aes(x=hto,fill=group),position = "dodge")+facet_wrap("lineage_hmap")


ggplot(mtd[Pseudotime_ALL_LTHSC_CTRLonly<10])+geom_boxplot(aes(x=hto,y=predicted.lineage.score,fill=group))

