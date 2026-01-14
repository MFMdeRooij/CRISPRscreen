# Use the GeneExpressionRNAseqDESeq2.R output files, put this script in the same folder, adjust the settings, and run the script in Rstudio.
# This script performs STRING analysis with STRINGdb, and gene set enrichment analysis with ClusterProfiler
# Author: M.F.M. de Rooij PhD, Amsterdam UMC, Spaargaren Lab, 2023, info: m.f.derooij@amsterdamumc.nl
################################################################################
## Install package
#install.packages("BiocManager")
#BiocManager::install(c("clusterProfiler","msigdbr","enrichplot","org.Hs.eg.db","AnnotationDbi","ReactomePA","HGNChelper","pathview", "STRINGdb","ggtree"))
#install.packages("writexl","ggridges")
################################################################################

# Put this script in the folder where the files of  gene statistics are located
folder<- dirname(rstudioapi::getActiveDocumentContext()$path)

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
library("STRINGdb")
library("ggridges")
library("ggtree")

setwd(folder)
dir.create("GSEA_STRING")

comparisons <- list.files(pattern="DESeq2_RNAseq_.*.csv$")
comparisons <- comparisons[!grepl("kinome", comparisons)]
#comparisons<- c("DESeq2_RNAseq_AvsB.csv","DESeq2_RNAseq_AvsC.csv","DESeq2_RNAseq_BvsC.csv")
for (com in comparisons) {
  # Read data
  data <- read.csv(paste0(com))
  
  # Adjust gene symbols to recent hugo symbols
  Hugo <- checkGeneSymbols(data$hgnc_symbol, unmapped.as.na=T)
  data$Hugo <- Hugo$Suggested.Symbol
  data$EntrezID <- as.character(mapIds(org.Hs.eg.db, keys = data$Hugo, column = "ENTREZID", keytype = "SYMBOL"))
  data<- data[!is.na(data$EntrezID),]
  data<- data[data$EntrezID!="NULL",]
  
  # Remove duplicate Hugo symbols, take the one with highest average expression if the original is not there 
  duplicates<- data.frame(h=data$Hugo, n=1)
  duplicates<- aggregate(duplicates$n, by=list(duplicates$h), FUN=sum)
  duplicates<- duplicates[duplicates$x>1,]
  dataDupl<- data[data$Hugo %in% duplicates$Group.1,]
  
  data<- data[!data$Hugo %in% duplicates$Group.1,]
  for (g in unique(dataDupl$Hugo)){
    temp<- dataDupl[dataDupl$Hugo==g,]
    temp$counts<- temp$BaseMeanA+temp$BaseMeanB
    idmax<- temp$ensembl_gene_id[temp$count==max(temp$count)][1]
    
    if (sum(temp$hgnc_symbol==g)==1){
      data<- rbind(data, dataDupl[dataDupl$hgnc_symbol==g,])
    } else{
      data<- rbind(data, dataDupl[dataDupl$ensembl_gene_id==idmax,])
    }
  }
  
  # Gene list with gene symbols
  geneList<- log2(data$FoldChange)
  names(geneList)<- data$Hugo
  geneList<- geneList[order(geneList, decreasing=T)]
  # To prevent a bug
  geneList[geneList==0]<- 10^-10
  
  # Gene list with entrez ids
  geneListID<- log2(data$FoldChange)
  names(geneListID)<- data$EntrezID
  geneListID<- geneListID[order(geneListID, decreasing=T)]
  # To prevent a bug
  geneListID[geneListID==0]<- 10^-10
  
  # STRING
  data$logFC<- log2(data$FoldChange)
  data$fdr<- data$padj
  
  hitsDown<- data$Hugo[data$logFC < -1 & data$fdr < 0.1 & !is.na(data$fdr)]
  hitsUp<- data$Hugo[data$logFC > 1 & data$fdr < 0.1 & !is.na(data$fdr)]
  hits<- c(hitsDown,hitsUp)
  hitsDown50<- data$Hugo[order(data$FoldChange)][1:50]
  hitsUp50<- data$Hugo[order(data$FoldChange, decreasing = T)][1:50]

  string_db <- STRINGdb$new(version="12.0", species=9606, score_threshold=400, input_directory="")
  example1_mapped = string_db$map(data, "Hugo", removeUnmappedRows = TRUE)
  example1_mapped_fdr = string_db$add_diff_exp_color(subset(example1_mapped, fdr<0.1), logFcColStr="logFC")
  payload_id = string_db$post_payload(example1_mapped_fdr$STRING_id, colors=example1_mapped_fdr$color)
  
  pdf(paste0("GSEA_STRING/STRING_",com,".pdf"),10,10)
    if (length(hits)>0){
      string_db$plot_network(hits, payload_id=payload_id)
      title(main = "Hits All", adj = 0, cex.main = 1.5, col.main = "black")
    }
    if (length(hitsDown)>0){
      string_db$plot_network(hitsDown, payload_id=payload_id)
      title(main = "Hits Down", adj = 0, cex.main = 1.5, col.main = "black")
    }
    if (length(hitsUp)>0){
      string_db$plot_network(hitsUp, payload_id=payload_id)
      title(main = "Hits Up", adj = 0, cex.main = 1.5, col.main = "black")
    }
    string_db$plot_network(hitsDown50, payload_id=payload_id)
    title(main = "Down 50", adj = 0, cex.main = 1.5, col.main = "black")
    string_db$plot_network(hitsUp50, payload_id=payload_id)
    title(main = "Up 50", adj = 0, cex.main = 1.5, col.main = "black")
  dev.off()
  
  # GSEA
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
  for (cat in c("H", "C1", "C2", "C3")) {
    sheets[[paste0("GSEA_MSIGDB_",cat)]]<- vector()
    sig <- msigdbr(species="Homo sapiens", collection=cat)[,c("gs_name", "ncbi_gene")]
    try(sheets[[paste0("GSEA_MSIGDB_",cat)]] <- as.data.frame(setReadable(GSEA(geneList=geneListID, TERM2GENE = sig, minGSSize=minGSSize, pvalueCutoff=pvalueCutoff, seed=T), 'org.Hs.eg.db', 'ENTREZID')@result))
  }

  # Write to Excel
  for (i in 1:length(sheets)){
    if (!is.data.frame(sheets[[i]])){
      sheets[[i]]<- data.frame(ID="No Matches")
    }
  }
  write_xlsx(sheets, paste0("GSEA_STRING/GSEA_ClusterProfiler_", com,".xlsx"))

  sig <- msigdbr(species="Homo sapiens", collection="H")[,c("gs_name", "ncbi_gene")]
  GSEAplots <- setReadable(GSEA(geneList=geneListID, TERM2GENE = sig, minGSSize=minGSSize, pvalueCutoff=pvalueCutoff, seed=T), 'org.Hs.eg.db', 'ENTREZID')

  if (GSEAplots@setType=="KEGG"){
    dir.create(paste0("GSEA_STRING/GSEA_KEGGpathview_",com))
    setwd(paste0("GSEA_STRING/GSEA_KEGGpathview_",com))
    for (i in 1:20){
      try(pathview(gene.data  = geneListID,
                   pathway.id = GSEAplots$ID[i],
                   species    = "hsa"))
    }
    setwd("../../")
  }
  pdf(paste0("GSEA_STRING/GSEA_",GSEAplots@setType,"_", com,".pdf"), 13,10)

    print(dotplot(GSEAplots, x = "NES", color = "p.adjust", size = "GeneRatio", orderBy = "NES", showCategory = 10, split = ".sign") +
            ggtitle("Dotplot") + theme(aspect.ratio = 2)#+ facet_grid(.~.sign)
    )

    print(ridgeplot(GSEAplots, showCategory=20)+ggtitle("Ridgeplot")+xlab("Log2 fold change")+ theme(aspect.ratio = 2))

    print(treeplot(pairwise_termsim(GSEAplots), cladelab_offset=8, tiplab_offset=.3, fontsize_cladelab=5)+hexpand(.1)+ggtitle("Treeplot"))

    print(emapplot(pairwise_termsim(GSEAplots))+ggtitle("Enrichment map"))

    print(cnetplot(GSEAplots, foldChange=geneList)+ggtitle("CNET plot"))

    print(heatplot(GSEAplots, foldChange=geneList, showCategory=20)+ggtitle("Heatplot"))

    for (i in 1:20){
      print(gseaplot2(GSEAplots, geneSetID = i, title = paste0(com, ": ",GSEAplots$Description[i], " (",GSEAplots@setType,":", GSEAplots$ID[i],")",
                                                               "\nES: ", round(GSEAplots$enrichmentScore[i],3), ",",
                                                               " NES: ", round(GSEAplots$NES[i],3), ",",
                                                               " pvalue: ",if (GSEAplots$pvalue[i]>= 0.001) {round(GSEAplots$pvalue[i],3)}
                                                               else{formatC(GSEAplots$pvalue[i],format = "e", digits = 2)},",",
                                                               " p.adjust: ",if (GSEAplots$p.adjust[i]>= 0.001) {round(GSEAplots$p.adjust[i],3)}
                                                               else{formatC(GSEAplots$p.adjust[i],format = "e", digits = 2)},",",
                                                               " leading edge: (", GSEAplots$leading_edge[i],")" )))
    }
    dev.off()
}
