# Use the CRISPRScreenAnalysis.R output files, put this script in the same folder, adjust the settings, and run the script in Rstudio.
# This script performs gene set enrichment analysis using the ClusterProfiler package
# Author: M.F.M. de Rooij PhD, Amsterdam UMC, Spaargaren Lab, 2023, info: m.f.derooij@amsterdamumc.nl
################################################################################
## Install package
#install.packages("BiocManager")
#BiocManager::install(c("clusterProfiler","msigdbr","enrichplot","org.Hs.eg.db","AnnotationDbi","ReactomePA","HGNChelper","pathview", "ggplot2"))
#install.packages("writexl")
################################################################################
# 
# Folder and files
setwd("H:/BioWin/CRISPRscreen/")
comparisons<- c("DESeq2 T0vsT1 Genes", "DESeq2 T0vsT2 Genes", "DESeq2 T1vsT2 Genes")
# You can also rename the files to cell line or drug names

# Minimal size of each geneSet for analyzing, and p.adjust cutoff
minGSSize<- 10
pvalueCutoff<- 1

###########################################################################################################################
library("clusterProfiler")
library("msigdbr")
library("enrichplot")
library("org.Hs.eg.db")
library("AnnotationDbi")
library("ReactomePA")
library("HGNChelper")
library("writexl")
library("pathview")
library("ggplot2")

for (com in comparisons) {
  data <- read.csv(paste0(com,".csv"))
  
  # Data
  Hugo <- checkGeneSymbols(data$GeneSymbol, unmapped.as.na=T)
  data$Hugo <- Hugo$Suggested.Symbol
  data$EntrezID <- as.character(mapIds(org.Hs.eg.db, keys = data$Hugo, column = "ENTREZID", keytype = "SYMBOL"))
  data<- data[!is.na(data$EntrezID),]
  data<- data[data$EntrezID!="NULL",]
  
  # Gene list with gene symbols
  geneList<- log2(data$MedianFoldChange)
  names(geneList)<- data$Hugo
  geneList<- geneList[order(geneList, decreasing=T)]
  # To prevent a bug
  geneList[geneList==0]<- 10^-10
  
  # Gene list with entrez ids
  geneListID<- log2(data$MedianFoldChange)
  names(geneListID)<- data$EntrezID
  geneListID<- geneListID[order(geneListID, decreasing=T)]
  # To prevent a bug
  geneList[geneList==0]<- 10^-10
  
  sheets<- list()
  
  set.seed(211081)
  
  # KEGG
  cat<- "KEGG"
  sheets[[paste0("GSEA_",cat)]]<- vector()
  try(sheets[[paste0("GSEA_",cat)]] <- as.data.frame(setReadable(gseKEGG(geneList=geneListID, organism='hsa', minGSSize=minGSSize, pvalueCutoff=pvalueCutoff, seed=T), 'org.Hs.eg.db', 'ENTREZID')@result))

  # REACTOME
  cat<- "REACTOME"
  sheets[[paste0("GSEA_",cat)]]<- vector()
  try(sheets[[paste0("GSEA_",cat)]] <- as.data.frame(setReadable(gsePathway(geneList=geneListID, organism = "human", minGSSize=minGSSize, pvalueCutoff=pvalueCutoff, seed=T), 'org.Hs.eg.db', 'ENTREZID')@result))

  # WIKI
  cat<- "WIKI"
  sheets[[paste0("GSEA_",cat)]]<- vector()
  try(sheets[[paste0("GSEA_",cat)]] <- as.data.frame(setReadable(gseWP(geneList=geneListID, organism="Homo sapiens", minGSSize=minGSSize, pvalueCutoff=pvalueCutoff, seed=T), 'org.Hs.eg.db', 'ENTREZID')@result))

  # GO
  for (cat in c("BP", "MF", "CC")) {
    sheets[[paste0("GSEA_GO_",cat)]]<- vector()
    try(sheets[[paste0("GSEA_GO_",cat)]] <- as.data.frame(gseGO(geneList=geneList, OrgDb=org.Hs.eg.db, ont=cat,minGSSize=minGSSize, pvalueCutoff=pvalueCutoff, seed=T, keyType="SYMBOL")@result))
  }

  # MSIGDB
  for (cat in c("H", "C2", "C1")) {
    sheets[[paste0("GSEA_MSIGDB_",cat)]]<- vector()
    sig <- msigdbr(species="Homo sapiens", category=cat)[,c("gs_name", "entrez_gene")]
    try(sheets[[paste0("GSEA_MSIGDB_",cat)]] <- as.data.frame(setReadable(GSEA(geneList=geneListID, TERM2GENE = sig, minGSSize=minGSSize, pvalueCutoff=pvalueCutoff, seed=T), 'org.Hs.eg.db', 'ENTREZID')@result))
  }

  # Write to Excel
  for (i in 1:length(sheets)){
    if (!is.data.frame(sheets[[i]])){
      sheets[[i]]<- data.frame(ID="No Matches")
    }
  }
  write_xlsx(sheets, paste0("GSEA_ClusterProfiler_", com,".xlsx"))

  GSEAplots <- setReadable(gseKEGG(geneList=geneListID, organism='hsa', minGSSize=minGSSize, pvalueCutoff=pvalueCutoff, seed=T), 'org.Hs.eg.db', 'ENTREZID')
  if (GSEAplots@setType=="KEGG"){
    dir.create(paste0("GSEA_KEGGpathview_",com))
    setwd(paste0("GSEA_KEGGpathview_",com))
    for (i in 1:20){
      try(pathview(gene.data  = geneListID,
                   pathway.id = GSEAplots$ID[i],
                   species    = "hsa"))
    }
    setwd("../")
  }
  pdf(paste0("GSEA_",GSEAplots@setType,"_", com,".pdf"), 13,10)
  
    print(dotplot(GSEAplots, showCategory=20) + ggtitle("Dotplot")+xlab("Gene ratio (# leading egde genes / # available genes)"))
    
    print(treeplot(pairwise_termsim(GSEAplots), cluster.params = list(method = "average"))+ggtitle("Treeplot"))
    
    print(emapplot(pairwise_termsim(GSEAplots))+ggtitle("Enrichment map"))
    
    print(cnetplot(GSEAplots, color.params=list(foldChange=geneList, edge=T))+ggtitle("CNET plot"))
    
    print(ridgeplot(GSEAplots, showCategory=20)+ggtitle("Ridgeplot")+xlab("Log2 fold change"))
    
    print(heatplot(GSEAplots, foldChange=geneList, showCategory=20)+ggtitle("Heatplot"))
    
    for (i in 1:20){
      print(gseaplot2(GSEAplots, geneSetID = i, title = paste0(com, ": ",GSEAplots$Description[i], " (",GSEAplots@setType,":", GSEAplots$ID[i],")", 
                                                              "\nES: ", round(GSEAplots$enrichmentScore[i],3), ",",
                                                              " NES: ", round(GSEAplots$NES[i],3), ",",
                                                              " pvalue: ", round(GSEAplots$pvalue[i],3), ",",
                                                              " p.adjust: ", round(GSEAplots$p.adjust[i],3)
                                                              
      )))
    }
  dev.off()
}