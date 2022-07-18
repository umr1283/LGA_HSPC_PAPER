# LGA_HSPC_PAPER

[![DOI](https://zenodo.org/badge/404302229.svg)](https://zenodo.org/badge/latestdoi/404302229)

This repository contains all scripts used to perform the analyses described in the manuscript entitled: "_**Epigenetic and Transcriptomic Programming of HSC Quiescence Signaling in Large for Gestational Age Neonates**_" from Pelletier _et al._, 2022 (doi: 10.3390/ijms23137323).

Scripts were executed using the following Docker images stored on Docker Hub:

* [umr1283/stat:R403](https://hub.docker.com/layers/umr1283/stat/R403/images/sha256-2b0fb490ba31e186f8b8f53c8e11e79177b6b7eeb8e8267b8c90f5f4c4d6c4f1)
* [umr1283/stat:4.0.5](https://hub.docker.com/layers/stat/umr1283/stat/4.0.5/images/sha256-b4a4e46c4586841e0c3b8b3a4b0c8d0c69c19f85207bfe64ab60f06ef2ba40ec)

Data sharing statement: *see manuscript*.

Analyses:
- Optimized methylation gene set analysis reveals association between LGA DNA hypermethylation and stem cell differentiation pathways
  1. [01-lga_vs_ctrl_limma_DMCs_analysis.R](scripts/01-lga_vs_ctrl_limma_DMCs_analysis.R)
  2. [02-gene_score_calculation_and_validation.R](scripts/02-gene_score_calculation_and_validation.R)
  3. [03-pathway_analysis.R](scripts/03-pathway_analysis.R)
- Single-cell transcriptomic analysis confirms alteration of hyper-methylated genes in pathways regulating stem cell differentiation among LGA HSCs
  1. [04-make_hematomap.R](scripts/04-make_hematomap.R)
  2. [05-integr_singlecell_cbps.R](scripts/05-integr_singlecell_cbps.R)
  3. [06-LGA_vs_Ctrl_RNA.R](scripts/06-LGA_vs_Ctrl_RNA.R)
- DNA methylation changes occurs in HSCs and DEGs associated open chromatin regions
  1. [07-DMCs_ATAC_integr.R](scripts/07-DMCs_ATAC_integr.R)
  2. [08-chromatin_change_LGA_vs_Ctrl.R](scripts/08-chromatin_change_LGA_vs_Ctrl.R)
- EGR1, KLF2 and KLF4 are key upstream regulators influenced by early epigenetic programming in LGA
  1. [09-SCENIC.R](scripts/09-SCENIC.R)
  2. [10-regulons_enrichment_genescore.R](scripts/10-regulons_enrichment_genescore.R)
  3. [10A-motif_analysis.R](scripts/10A-motif_analysis.R)
- Multimodal co-regulatory network recapitulating TF-gene interactions influenced by early epigenetic programming in LGA
  1. [11-GRN_regulons.R](scripts/11-GRN_regulons.R)
  2. [12-GRN_final.R](scripts/12-GRN_final.R)
- In vitro analysis confirms the alteration of HSPCs differentiation capacities in LGA
  1. [13-Pseudotime.R](scripts/13-Pseudotime.R)
  2. [14-RNA_velocity.R](scripts/14-RNA_velocity.R)
  3. [15-figures_epi_response.R](scripts/15-figures_epi_response.R)

Figures and statistical tests: [15-figures_epi_response.R](scripts/15-figures_epi_response.R)
