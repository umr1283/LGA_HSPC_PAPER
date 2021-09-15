### Project Setup ==================================================================================
library(here)
out<- here("outputs", "12A-Pseudotime_integrated")
dir.create(out, recursive = TRUE, showWarnings = FALSE, mode = "0775")


### Load Packages ==================================================================================
# renv::install("bioc::batchelor")
# renv::install('cole-trapnell-lab/leidenbase')
# renv::install("cole-trapnell-lab/monocle3")
# renv::install("satijalab/seurat-wrappers")

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratWrappers)
  library(monocle3)
  library(Matrix)
  library(ggplot2)
  library(patchwork)
  library(data.table)
  set.seed(1234)

  
})


### Tables and Figures Theme =======================================================================
# theme_set(theme_light())


### Functions ======================================================================================
fp<-function(...)file.path(...)


# a helper function to identify the root principal points:
get_earliest_principal_node <- function(cbp.cds, time_bin="LT-HSC"){
  cell_ids <- which(colData(cbp.cds)[, "cell_type_hmap"] == time_bin )
  
  closest_vertex <-
  cbp.cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cbp.cds), ])
  root_pr_nodes <-
  igraph::V(principal_graph(cbp.cds)[["UMAP"]])$name[as.numeric(names
  (which.max(table(closest_vertex[cell_ids,]))))]
  
  root_pr_nodes
}


### Analysis =======================================================================================

cbps<-readRDS("outputs/10A-classical_integr/cbps0-8_clean.rds")
DefaultAssay(cbps) <- "integrated"
Idents(cbps) <- "cell_type_hmap"
cbps <- RunPCA(cbps,reduction.name = "integrated.pca",features = VariableFeatures(cbps))

cbps <- RunUMAP(cbps, reduction = "integrated.pca",n.components = 3,reduction.name = "integrated.umap", dims = c(1:50))

DimPlot(cbps,reduction="integrated.umap",label=T)


cbps[["pca"]]<-cbps@reductions$integrated.pca
cbps[["umap"]]<-cbps@reductions$integrated.umap

cbps.cds <- as.cell_data_set(cbps)

cbps.cds <- cluster_cells(cds = cbps.cds, reduction_method = "UMAP")
cbps.cds <- learn_graph(cbps.cds, use_partition = TRUE)

cbps.cds <- order_cells(cbps.cds, root_pr_nodes=get_earliest_principal_node(cbps.cds))


plot_cells(
  cds = cbps.cds,
  color_cells_by = "pseudotime",
  show_trajectory_graph = TRUE,
  label_branch_points= TRUE,
  label_leaves=FALSE
)

cbps <- AddMetaData(
  object = cbps,
  metadata = cbps.cds@principal_graph_aux@listData$UMAP$pseudotime,
  col.name = "pseudotime"
)

cbps$pseudotime[cbps$pseudotime==Inf]<-NA
Idents(cbps)<-"cell_type_hmap"
FeaturePlot(cbps,"pseudotime",reduction = 'umap',label=T)
FeaturePlot(cbps,"pseudotime",dims=c(1,3),reduction = 'umap',label=T)


cbps_f<-subset(cbps,lineage_hmap!="18"&ambigous==F&group!="iugr")

mts<-unique(data.table(cbps_f@meta.data,keep.rownames = "bc"),by="sample_hto")
print(table(mts$group_hto))

mtd<-data.table(cbps_f@meta.data,keep.rownames = "bc")

fwrite(mtd,fp(out,"metadata_pseudotime_ComputRoot.csv"))

ggplot(mtd)+geom_boxplot(aes(y=pseudotime,x=lineage_hmap),alpha=0.7)+theme_minimal()

ggplot(mtd)+geom_density(aes(x=pseudotime,col=group))+facet_wrap("hto")

lga_b<-ggplot(mtd[group=="lga"&hto==F])+geom_density(aes(x=pseudotime,fill=group,col=group),alpha=0.7)+
  facet_wrap("sample")+theme_minimal()

ctrl_b<-ggplot(mtd[group=="ctrl"&hto==F])+geom_density(aes(x=pseudotime,fill=group,col=group),alpha=0.7)+
  facet_wrap("sample")+theme_minimal()
lga_b/ctrl_b

saveRDS(cbps_f,file=fp(out,"cbps_filtered_pseudotime_computroot.rds"))

cbps_f<-readRDS(file=fp(out,"cbps_filtered_pseudotime_computroot.rds"))


### Complete =======================================================================================
message("Success!", appendLF = TRUE)


