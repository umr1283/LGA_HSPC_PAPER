library("BiocManager")
library(R.utils)
renv::hydrate()
renv::install(c("YuLab-SMU/clusterProfiler"))
renv::install("missMDA")
renv::snapshot() #git commit
renv::install("bioc::batchelor")
renv::install('cole-trapnell-lab/leidenbase')
renv::install("cole-trapnell-lab/monocle3")
renv::install("satijalab/seurat-wrappers")

renv::install(c("bioc::scater",
                "bioc::GENIE3",
                "bioc::RcisTarget")) #for SCENIC
                
renv::install("doRNG") #for SCENIC
renv::install("doMC") #for SCENIC
renv::install("aertslab/SCENIC")


renv::snapshot() #git commit
