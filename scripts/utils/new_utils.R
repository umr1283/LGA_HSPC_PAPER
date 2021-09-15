
library("data.table")
library("stringr")
library("ggplot2")
library("ggrepel")
library("patchwork")
library("here")
fp<-function(...)file.path(...)

ps<-function(...,sep="",collapse = NULL)paste(...,sep=sep,collapse = collapse)

bed_inter<- function(a, b, opt1="-wa", opt2="-wb",out_dir=".", select=NULL, col.names=NULL){
  require(data.table)
  l<-list(a,b)
  files_to_rm<-c(FALSE,FALSE)
  file_paths<-sapply(1:2, function(i){
    x<-l[[i]]
    if(is.data.frame(x)){
      files_to_rm[i]<-TRUE
      file_path<-fp(out_dir,paste0("temp",i,".bed"))
      fwrite(x,file_path,sep="\t")
      
    }else{
      file_path<-x
    }
    return(file_path)
  })
  
  out_file<-fp(out_dir,"temp_inter.bed")
  cmd<-paste("bedtools intersect -a",file_paths[1],"-b",file_paths[2], opt1, opt2,">",out_file)
  message("run in shell : ",cmd)
  system(cmd)
  message("done.")
  
  if(!is.null(col.names)){
    dt<-fread(out_file,select = select,col.names = col.names)
    file.remove(out_file)
    file.remove(file_paths[files_to_rm])
    return(dt)
  }else{
    dt<-fread(out_file,select = select)
    file.remove(out_file)
    file.remove(file_paths[files_to_rm])
    return(dt)
  }
}




GetMartGenes<-function()biomaRt::useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
GetMartReg<-function()biomaRt::useEnsembl(biomart = "regulation", dataset = "hsapiens_regulatory_feature")
GetMartMotif<-function()biomaRt::useEnsembl(biomart = "regulation", dataset = "hsapiens_motif_feature")



GetBiomartAttrs<-function(mart)data.table::data.table(biomaRt::listAttributes(mart))
GetBiomartFilter<-function()data.table::data.table(biomaRt::listFilter(mart))

GetBiomartAttrs_reg<-function()data.table::data.table(listAttributes(biomaRt::useEnsembl(biomart = "regulation",dataset = "hsapiens_regulatory_feature")))
GetBiomartFilter_reg<-function()data.table::data.table(biomaRt::listFilter(biomaRt::useEnsembl(biomart = "regulation", dataset = "hsapiens_regulatory_feature")))

TransNMtoSymbol<-function(refseq_ids){
  require(biomaRt)
  require(data.table)
      return(data.table::data.table(getBM(attributes = c('refseq_mrna', 'hgnc_symbol'),
      filters = 'refseq_mrna', 
      values = refseq_ids, 
      mart = GetMartGenes())))
  }

TransSymboltoNM<-function(hgnc_symbols){
  require(biomaRt)
  require(data.table)
      return(data.table::data.table(biomaRt::getBM(attributes = c('refseq_mrna', 'hgnc_symbol'),
      filters = 'hgnc_symbol', 
      values = hgnc_symbols, 
      mart = GetMartGenes())))
}

TransEnsembltoSymbol<-function(ensembl_ids){
  require(biomaRt)
  require(data.table)
      return(data.table::data.table(biomaRt::getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol'),
      filters = 'ensembl_gene_id', 
      values = ensembl_ids, 
      mart = GetMartGenes())))
  }

TransSymboltoEnsembl<-function(hgnc_symbols){
  require(biomaRt)
  require(data.table)
      return(data.table::data.table(biomaRt::getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol'),
      filters = 'hgnc_symbol', 
      values = hgnc_symbols, 
      mart = GetMartGenes())))
}

TransEnsemblVerstoSymbol<-function(ensembl_ids){
  require(biomaRt)
  require(data.table)
      return(data.table::data.table(biomaRt::getBM(attributes = c('ensembl_gene_id_version', 'hgnc_symbol'),
      filters = 'ensembl_gene_id_version', 
      values = ensembl_ids, 
      mart = GetMartGenes())))
  }

TransSymboltoEnsemblVers<-function(hgnc_symbols){
  require(biomaRt)
  require(data.table)
      return(data.table::data.table(biomaRt::getBM(attributes = c('ensembl_gene_id_version', 'hgnc_symbol'),
      filters = 'hgnc_symbol', 
      values = hgnc_symbols, 
      mart = GetMartGenes())))
}



TransEnsembltoNM<-function(ensembl_ids){
  require(biomaRt)
  require(data.table)
      return(data.table::data.table(biomaRt::getBM(attributes = c('ensembl_gene_id', 'refseq_mrna'),
      filters = 'ensembl_gene_id', 
      values = ensembl_ids, 
      mart = GetMartGenes())))
  }

TransNMtoEnsembl<-function(refseq_ids){
  require(biomaRt)
  require(data.table)
      return(data.table::data.table(biomaRt::getBM(attributes = c('ensembl_gene_id', 'refseq_mrna'),
      filters = 'refseq_mrna', 
      values = refseq_ids, 
      mart = GetMartGenes())))
}

TransEnsemblVerstoNM<-function(ensembl_ids){
  require(biomaRt)
  require(data.table)
      return(data.table::data.table(biomaRt::getBM(attributes = c('ensembl_gene_id_version', 'refseq_mrna'),
      filters = 'ensembl_gene_id_version', 
      values = ensembl_ids, 
      mart = GetMartGenes())))
  }

TransNMtoEnsemblVers<-function(refseq_ids){
  require(biomaRt)
  require(data.table)
      return(data.table::data.table(biomaRt::getBM(attributes = c('ensembl_gene_id_version', 'refseq_mrna'),
      filters = 'refseq_mrna', 
      values = refseq_ids, 
      mart = GetMartGenes())))
  }



tr<-function(ids_sepBySlash,retourne="all",sep="/",tradEntrezInSymbol=FALSE,uniqu=TRUE){
  IDs<-as.vector(strsplit(ids_sepBySlash,sep)[[1]])
  if(retourne=="all"){
    ret<-1:length(IDs)
  }else{
    ret<-retourne
  }
  if(tradEntrezInSymbol){
    require(clusterProfiler)
    library(org.Hs.eg.db)
    if(retourne=="all"){
      return(clusterProfiler::bitr(IDs, fromType = "ENTREZID",toType =  "SYMBOL",OrgDb = org.Hs.eg.db)$SYMBOL)
    }
    else{
      return(clusterProfiler::bitr(IDs, fromType = "ENTREZID",toType =  "SYMBOL",OrgDb = org.Hs.eg.db)$SYMBOL[ret])
    }
    
  }else{
    if(uniqu){
      return(unique(IDs[ret]))
    }else{
      return(IDs[ret])
    }
    
    
  }
  
}


GetVarPCs<-function(pca,rngPCs="all"){
  if(is.character(rngPCs)){
    rngPCs<-1:length(pca$sdev)
  }
  pct.varPCs<-pca$sdev[rngPCs]^2/sum(pca$sdev^2)
  names(pct.varPCs)<-rngPCs
  return( pct.varPCs)
}





# Basic function to convert human to mouse gene names
convertHumanGeneList <- function(x){
  
  require("biomaRt")
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  
  genesV2 = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = x , mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=T)
  
  mousex <- unique(genesV2[, 2])
  
  # Print the first 6 genes found to the screen
  print(head(mousex))
  return(mousex)
}

convertMouseGeneList <- function(x){
  
  require("biomaRt")
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  
  genesV2 = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = x , mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)
  humanx <- unique(genesV2[, 2])
  
  # Print the first 6 genes found to the screen
  print(head(humanx))
  return(humanx)
}


#Over repre / GSEA



over_repr_test_simple<-function(set1,set2,size_universe){
  if(any(duplicated(set1))){
    set1<-unique(set1)
  }
  if(any(duplicated(set2))){
    set2<-unique(set2)
    }
  phyper(q=sum(set1%in%set2)-1, 
       #number of white balls drawn without replacement from an urn which contains both black and white balls. 
       #<=> number of tcf4 target in DEGs. "-1" because normally give P[X > x] but we want P[X >= x])
       m=length(set2), #the number of white balls in the urn. <=> number of DEGs 
       n=size_universe-length(set2), #the number of black balls in the urn. <=> number of genes tested - number of DEGs
       k=length(set1), #the number of balls drawn from the urn <=> numbers of tcf4 targets 
       lower.tail=FALSE)  # if TRUE (default), probabilities are P[X ≤ x] (under-representation), otherwise, P[X > x] (over-representation).

}

over_repr_test_multi<-function(genes_of_interest,terms_list,size_universe,min.term.size=10,max.term.size=500){
  res_or<-data.table(term=names(terms_list),term.size=sapply(terms_list,length))
  res_or<-res_or[term.size<=max.term.size]
  n_terms<-nrow(res_or)
  message(length(terms_list)-n_terms, " terms were filtered due to term.size above the limit of ",max.term.size," genes")
  res_or<-res_or[term.size>=min.term.size]
  message(n_terms-nrow(res_or), " terms were filtered due to term.size below the limit of ",min.term.size," genes")
  
  res_or[,n.genes.altered:=length(genes_of_interest)]
  res_or[,n.enriched:=sum(genes_of_interest%in%terms_list[[term]]),by="term"]
  res_or[,genes.enriched:=paste(genes_of_interest[genes_of_interest%in%terms_list[[term]]],collapse="|"),by="term"]
  res_or[,pct.enriched:=n.enriched/term.size]
  res_or[,pval:=phyper(q=n.enriched-1, 
                     m=n.genes.altered, 
                     n=size_universe-n.genes.altered, 
                     k=term.size, 
                     lower.tail=FALSE),
       by="term"]
  res_or[,padj:=p.adjust(pval,method = 'BH')]
  message(nrow(res_or[padj<0.05])," terms enriched in your genes of interest with padj<0.05")
  return(res_or)
}





