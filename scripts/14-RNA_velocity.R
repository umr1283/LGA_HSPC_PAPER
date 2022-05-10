out<-"outputs/20-RNA_velocity"
dir.create(out)

source("scripts/utils/new_utils.R")

#infer RNA velocity using scVelo to 1) find if LGA HSC have a differentiation bias ?
#2) study which genes (IEGs ?) are freshly transcribed (++ premRNA/ mRNA) 
#and 3) if there is transcription differences between LGA and control

####need first install python####
renv::snapshot()
#rm python files

renv::use_python()

library(reticulate)

reticulate::py_install(packages = c("numpy" ,"scipy", "cython" ,"numba" ,"matplotlib", "scikit-learn", "h5py", "click"))
reticulate::py_install(packages ="umap-learn")
reticulate::py_install(packages ="pysam",pip=T)
reticulate::py_install(packages = c("velocyto"),pip = T)
reticulate::py_install(packages ="scvelo",pip=T)

#need downgrade numba to 0.52
reticulate::py_install(packages ="numba==0.52")
reticulate::py_run_string("
import scvelo as scv
scv.logging.print_version()
                          ") ##Running scvelo 0.2.4 (python 3.7.3) on 2022-04-07 13:23

renv::install("bioc::LoomExperiment")
renv::snapshot()

#need also samtools
system("../singlecell/tools/samtools-1.12/samtools")


#run velocyto counting pipeline
  #first, need dl genome repeat sequence to mask here : https://genome.ucsc.edu/cgi-bin/hgTables?hgsid=611454127_NtvlaW6xBSIRYJEBI0iRDEWisITa&clade=mammal&org=Human&db=0&hgta_group=allTracks&hgta_track=rmsk&hgta_table=rmsk&hgta_regionType=genome&position=&hgta_outputType=gff&hgta_outFileName=mm10_rmsk.gtf
dir.create("outputs/20-RNA_velocity/velocyto_counts")
#run 14A files


#then use scVelo (https://scvelo.readthedocs.io/getting_started.html) in python to generate the RNA velocity matrix
#need first merge all velocito object #see https://github.com/basilkhuder/Seurat-to-RNA-Velocity#integrating-loom-file-and-meta-data
#with cells kept in seurat analysis
mtd<-fread("outputs/06-integr_singlecell_cbps/metadata_cbps_filtered.csv.gz")

mtd[,cell_id:=paste(str_extract(bc,"[ATCG]+"),orig.ident,sep="_")]
fwrite(mtd,"outputs/06-integr_singlecell_cbps/metadata_cbps_filtered.csv.gz")
#update reticulat because of bug
renv::install("rstudio/reticulate")

#need get umap coord
mtd<-fread("outputs/06-integr_singlecell_cbps/metadata_cbps_filtered.csv.gz")

cbps<-readRDS("outputs/06-integr_singlecell_cbps/cbps_filtered.rds")
umap_coord<-data.table(cbps@reductions$ref.umap@cell.embeddings,keep.rownames = "bc")
umap_coord<-merge(umap_coord,mtd[,.(bc,cell_id)])
fwrite(umap_coord,"outputs/06-integr_singlecell_cbps/umap_cbps.csv")
table(mtd[hto==T]$batch)

#finally run scvelo, 14B- file
# see https://scvelo.readthedocs.io/DynamicalModeling/

#add velocity matrix on seurat object
velo<-fread("outputs/20-RNA_velocity/cbps_hto_dynamical_velocity_matrix.csv")
velo[1:10,.(cell_id)]
velo<-as.matrix(t(data.frame(velo,row.names = "cell_id")))

mtdf<-mtd[cell_id%in%colnames(velo)]
cbps_h<-cbps[,mtdf$bc]

cbps_h[["velocity"]]<-CreateAssayObject(data =velo[2:nrow(velo),] )

mtdvelo<-fread("outputs/14-RNA_velocity/cbps_hto_dynamical_velocity_metadata.csv")
mtdvelo[,cell_id:=V1]
cbps_h<-AddMetaData(cbps_h,metadata = data.frame(mtdvelo[,-c("V1","batch")],row.names = "cell_id"))
saveRDS(cbps_h@assays$velocity,fp(out,"cbps_hto_dynamical_velocity_assay.rds"))

#Differentiation Analysis with pseudotime####
trans<-fread("outputs/14-RNA_velocity/cbps_hto_dynamical_transition_matrix.csv")
trans<-as.matrix(t(data.frame(trans,row.names = "cell_id")))
dim(trans) #12684 12683
trans[1:10,1:10]
trans<-trans[2:nrow(trans),]

mtd<-fread("outputs/06-integr_singlecell_cbps/metadata_cbps_filtered.csv.gz")
mtdf<-mtd[cell_id%in%colnames(trans)]

#add pseudotime
pseudo<-fread("outputs/12-Pseudotime/metadata_pseudotime.csv")
pseudo[,cell_id:=ps(str_extract(bc,"[ATCG]+"),orig.ident,sep='_')]
mtdf<-merge(mtdf,pseudo[,.(cell_id,pseudotime)])

mtdff<-mtdf[pseudotime!=Inf]
summary(mtdff$pseudotime)
mtdff[,pseudo_pred:=sapply(1:.N,function(i){
  cell<-cell_id[i]
  return(sum(pseudotime*trans[cell,cell_id]))
  })]

mtdff[,pseudo_bias:=pseudo_pred-pseudotime]


fwrite(mtdff,fp(out,"pseudo_bias_rna_velo_based_lineages.csv.gz"))

lins<-c("LT-HSC","HSC","MPP/LMPP","Myeloid","Lymphoid","Erythro-Mas")


mtdff[,lineage_hmap:=factor(lineage_hmap,levels=c("LT-HSC","HSC","MPP/LMPP","Myeloid","Lymphoid","Erythro-Mas"))]
ggplot(mtdff[lineage_hmap%in%lins])+
  geom_boxplot(aes(x=group,y=pseudo_bias,fill=group,group=sample),outlier.shape = NA)+
  facet_wrap("lineage_hmap")

ggplot(mtdff[lineage_hmap%in%lins])+
  geom_boxplot(aes(x=lineage_hmap,y=pseudo_bias,fill=group))

ggplot(mtdff[lineage_hmap%in%lins])+
  geom_boxplot(aes(x=group,y=pseudo_bias,fill=group),outlier.shape = NA)+
  facet_wrap("lineage_hmap")

lins<-c("LT-HSC","HSC","MPP/LMPP","Myeloid","Lymphoid","Erythro-Mas")

mtdffl<-mtdff[lineage_hmap%in%lins]
mtdffl[,lineage_hmap:=factor(lineage_hmap,levels = lins)]
mtdffl[,avg_pseudobias:=mean(pseudo_bias),by=.(lineage_hmap,sample)]

ggplot(unique(mtdffl,by=c("sample","lineage_hmap")))+geom_boxplot(aes(x=lineage_hmap,y=avg_pseudobias,fill=group))
unique(mtdffl[lineage_hmap=="MPP/LMPP"],by=c("sample","lineage_hmap"))

#latent time (~pseudotime) analysis ####
#what lineage/ct are the bigger latent time
mtdvelo<-fread("outputs/14-RNA_velocity/cbps_hto_dynamical_velocity_metadata.csv")
mtdvelo[,cell_id:=V1]

mtdcbps<-fread("outputs/14-integr_singlecell_cbps/metadata_cbps_filtered.csv.gz")
mtd<-merge(mtdvelo[,-"batch"],mtdcbps)
ggplot(mtd)+geom_boxplot(aes(x=cell_type_hmap,y=latent_time))
fwrite(mtd,fp(out,"metadata_cbps_hto_dynamical_merged.csv"))

#keys Genes following main dynamics ####
#Top-likelihood genes
#=> assume Driver genes because display pronounced dynamic behavior

genes_mtd<-fread("outputs/14-RNA_velocity/cbps_hto_dynamical_genes_metadata.csv")
genes_mtd<-genes_mtd[!is.na(fit_likelihood)]
genes_mtd#1742
head(genes_mtd[order(-fit_likelihood)],20)
head(genes_mtd[order(-fit_likelihood),.(Gene,fit_likelihood)],50)
#             Gene fit_likelihood
#  1:        RPL41      0.7024196
#  2:         FTH1      0.4409182
#  3:        RPS21      0.4193488
#  4:        RPL26      0.4011444
#  5:   AC002454.1      0.3567726
#  6:        SERF2      0.3304627
#  7:         TMA7      0.3147169
#  8:        MYADM      0.3128698
#  9:        SNHG5      0.3084447
# 10:         CPA3      0.3069082
# 11:        SSBP1      0.3048434
# 12:       SH2D4B      0.3040029
# 13:         MYL6      0.2988126
# 14:         AREG      0.2983997
# 15:       RPL27A      0.2977637
# 16:          FOS      0.2969011
# 17:         CD99      0.2942007
# 18:        PPM1G      0.2932404
# 19:      SEPTIN6      0.2924571
# 20:       TIPARP      0.2923209
# 21:         FOSB      0.2871259 *
# 22: ADAMTSL4-AS1      0.2866284
# 23:       ATP5MG      0.2844941
# 24:      POU2AF1      0.2822485
# 25:       IGFBP7      0.2806364
# 26:       ADGRG6      0.2798706
# 27:        COPS9      0.2797676
# 28:    LINC01220      0.2766412
# 29:        DUSP1      0.2734219 *
# 30:        MLLT3      0.2731625
# 31:        RPS4X      0.2726536
# 32:   AC008440.1      0.2719394
# 33:         KLF2      0.2718106 *
# 34:       NDUFS5      0.2700264
# 35:         ACTB      0.2697666
# 36:       AKAP12      0.2694474
# 37:       COX6B1      0.2688626
# 38:        RBPMS      0.2686618
# 39:       SPINK2      0.2678411
# 40:         ATF3      0.2665551
# 41:        SNHG7      0.2664279
# 42:        MYBL2      0.2653777
# 43:        TSTD1      0.2645509
# 44:       RNF130      0.2644003
# 45:       PET100      0.2629900
# 46:        SNRPF      0.2615482
# 47:   SPACA6P-AS      0.2603058
# 48:        IFFO2      0.2587123
# 49:         NASP      0.2567542
# 50:         GYPC      0.2563279
#             Gene fit_likelihood

#likelihood by lineage (lineage specific driver genes)
dyna_genes<-fread("outputs/20-RNA_velocity/dynamical_genes_by_lineage.csv",header = T)
dyna_genes<-melt(dyna_genes[,-"V1"],value.name = "gene",variable.name = "lineage",measure.vars = colnames(dyna_genes)[-1])
dyna_genes[lineage=="HSC"]
#     lineage         gene
#   1:     HSC        RPL41
#   2:     HSC        RPS21
#   3:     HSC   AC002454.1
#   4:     HSC         CD99
#   5:     HSC      SEPTIN6
#   6:     HSC       ATP5MG
#   7:     HSC         TMA7
#   8:     HSC       SPINK2
#   9:     HSC        MLLT3
#  10:     HSC         NASP
#  11:     HSC        RBPMS
#  12:     HSC         RHOH
#  13:     HSC        SNHG5
#  14:     HSC        MYADM
#  15:     HSC          TKT
#  16:     HSC          FOS
#  17:     HSC        SNHG7
#  18:     HSC         FTH1
#  19:     HSC       COX6B1
#  20:     HSC        SNRPF
#  21:     HSC       TMSB10
#  22:     HSC       EPSTI1
#  23:     HSC       RPL27A
#  24:     HSC       RNF130
#  25:     HSC      PIP4K2A
#  26:     HSC         GBP2
#  27:     HSC        YPEL5
#  28:     HSC         PTMA
#  29:     HSC       SMIM24
#  30:     HSC         GYPC
#  31:     HSC         ATF3
#  32:     HSC         KLF4 *
#  33:     HSC         ZNF3
#  34:     HSC        ZFP36
#  35:     HSC         KLF2 *
#  36:     HSC ADAMTSL4-AS1
#  37:     HSC         FOSB
#  38:     HSC     SERPINB1
#  39:     HSC         DTD1
#  40:     HSC        STX11
#  41:     HSC         BEX1
#  42:     HSC         JUND
#  43:     HSC        H3F3B
#  44:     HSC         AREG
#  45:     HSC       TUBA1A
#  46:     HSC         YBX3
#  47:     HSC        SOCS2
#  48:     HSC         MCL1
#  49:     HSC        DUSP1
#  50:     HSC       DYNLL1
#  51:     HSC         FDX1
#  52:     HSC        RPL26
#  53:     HSC         EGR1 *
#  54:     HSC        PTGS2
#  55:     HSC         MAFF
#  56:     HSC          IDS
#  57:     HSC       TIPARP
#  58:     HSC        KLF10
#  59:     HSC         SRGN
#  60:     HSC         BTG2
#  61:     HSC        BUD31
#  62:     HSC      TSC22D2
#  63:     HSC         YRDC
#  64:     HSC       NFKBIA
#  65:     HSC       CHST11
#  66:     HSC      CHORDC1
#  67:     HSC         EMP3
#  68:     HSC       PTGER4
#  69:     HSC        NR4A1
#  70:     HSC     PPP1R15A
#  71:     HSC      SERTAD3
#  72:     HSC         LMNA
#  73:     HSC        HSPA8
#  74:     HSC         IER2
#  75:     HSC       DNAJB1
#  76:     HSC         NEU1
#  77:     HSC        RPS4X
#  78:     HSC        RUNX3
#  79:     HSC        IFFO2
#  80:     HSC        HSPH1
#  81:     HSC         EREG
#  82:     HSC          AVP
#  83:     HSC          ID2
#  84:     HSC        ARL4A
#  85:     HSC          PNP
#  86:     HSC          KIN
#  87:     HSC        CXCL8
#  88:     HSC        UBE2S
#  89:     HSC       MAGED2
#  90:     HSC         TOB2
#  91:     HSC      CCDC173
#  92:     HSC       ARRDC2
#  93:     HSC        PPM1N
#  94:     HSC        NR4A2
#  95:     HSC   AL021155.5
#  96:     HSC         MLF1
#  97:     HSC         SKIL
#  98:     HSC       NDUFS5
#  99:     HSC         ICA1
# 100:     HSC         CD74
#      lineage         gene


#gene level velocity analysis [todo, redo]####
#IEGs/ EGRns
#IEGs / EGRns LGAs bias


#Velocity on control cells (wihtout HTO)  ####
#est ce que la dynamique de velocity est inversée ? (va vers les cellules differenciées?)
#run 20D
# => non

#Top-likelihood genes compared to HTO cells





