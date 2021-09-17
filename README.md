# LGA_HSPC_PAPER
This repository contains all scripts used to perform the analyses described in the manuscript entitled: "Excessive fetal growth affects HSC homeostasis through epigenetic programming of EGR1 transcriptional network" from Pelletier et al., 2021.

Scripts were executed using the following Docker images stored on Docker Hub:

* [umr1283/stat:R403](https://hub.docker.com/layers/umr1283/stat/R403/images/sha256-2b0fb490ba31e186f8b8f53c8e11e79177b6b7eeb8e8267b8c90f5f4c4d6c4f1?context=repo)
* [umr1283/stat:4.0.4](https://hub.docker.com/layers/umr1283/stat/4.0.4/images/sha256-774404184836d1dafeada1b4635adb7fa20d8520d27c38751867d52c0c2c09bd?context=repo)

Data sharing statement: see manuscript.

## Extreme fetal growth is associated with hypermethylation of key genes regulating hematopoietic homeostasis 
See scripts: 
  1. [01-lga_vs_ctrl_limma_DMCs_analysis.R](scripts/01-lga_vs_ctrl_limma_DMCs_analysis.R)
  2. [02-gene_score_calculation_and_validation.R](scripts/02-gene_score_calculation_and_validation.R)
  3. [03-pathway_analysis.R](scripts/03-pathway_analysis.R)
  4. [04-motif_analysis.R](scripts/04-motif_analysis.R)

## Epigenetic programming impacts hematopoietic stem cell response to stimulation at transcriptomic level
See scripts:
  1. [05-lga_vs_ctrl_limma_DMCs_analysis.R](scripts/)
  2. [06-integr_singlecell_cbps.R](scripts/06-integr_singlecell_cbps.R)
  3. [07-LGA_vs_Ctrl_Basal.R](scripts/07-LGA_vs_Ctrl_Basal.R)
  4. [08-HTO_signature.R](scripts/08-HTO_signature.R)
  9. [09-LGA_vs_Ctrl_Activated.R](scripts/09-LGA_vs_Ctrl_Activated.R)

## Epigenetic programming affects key regulons essential to HSCs self-renewal and differentiation
See scripts:
  1. [10-SCENIC.R](scripts/10-SCENIC.R)
  2. [11-regulons_enrichment_genescore.R](scripts/11-regulons_enrichment_genescore.R)
  3. [13-GRN_integr.R](scripts/13-GRN_integr.R)

## Epigenetic programming is associated with a shift in differentiation process
See scripts:
  1. [12-Pseudotime_integrated.R](scripts/12-Pseudotime_integrated.R)

## Epigenetic programming is associated with alteration of the hematopoietic compartment function and integrity
See scripts:
  1. [14-figures_Paper_LGA_HSPC.R](scripts/14-figures_Paper_LGA_HSPC.R)

## Figures and statistical tests
See scripts:
  1. [14-figures_Paper_LGA_HSPC.R](scripts/14-figures_Paper_LGA_HSPC.R)
