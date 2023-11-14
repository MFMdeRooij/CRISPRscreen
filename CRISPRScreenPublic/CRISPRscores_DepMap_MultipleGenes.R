# Plot CRISPR scores of the public CRISPR screens from DepMap.org (Broad Institute)
# Download CRISPR data (23q2) "CRISPRGeneEffect.csv, and model.csv from https://depmap.org/portal/download/all/
# All cell lines are plotted on page 1, and the lymphoid ones are plotted including labels on page 2
# Author: M.F.M. de Rooij PhD, Amsterdam UMC, Spaargaren Lab, 2023, info: m.f.derooij@amsterdamumc.nl
##################################################################################
# Workdirectory (folder where "CRISPR_(DepMap_22Q2_Public+Score,_Chronos).csv" is located)
setwd("G:/divg/Pathologie-ADHESIE/ALL/TheScreeningTeam/DepMap/")

# Genes of interest (add to line 15)
q <- function(...) {
  sapply(match.call()[-1], deparse)
}
interestingGenes <- toupper(q(
  # Which genes?  
  MCL1,BCL2L1,BCL2,BCL2L2,BCL2A1
))

disease<- "Plasma Cell Myeloma"
##################################################################################
# # R packages
#install.packages(c("ggplot2", "factoextra", "ggrepel"))
library("ggplot2")
library("factoextra")
library("ggrepel")
#install.packages("devtools")
#devtools::install_github("JosephCrispell/basicPlotteR")
library("basicPlotteR")
if (!exists("dataDepMap")) {
  dataDepMap<-read.csv("CRISPRGeneEffect.csv",stringsAsFactors = F)
  chronos<-as.data.frame(t(dataDepMap[,2:ncol(dataDepMap)]))
  colnames(chronos)<- dataDepMap$ModelID
  rownames(chronos)<- unlist(lapply(strsplit(as.character(rownames(chronos)),"\\.\\."), "[", 1))
  cells<- read.csv("Model.csv",stringsAsFactors = F)
  UsedIDs<- cells[match(colnames(chronos), cells$ModelID),]
  colnames(chronos)<- UsedIDs$StrippedCellLineName
  
  # Sample metadata
  df_designAll<-cells[cells$StrippedCellLineName %in% UsedIDs$StrippedCellLineName,c("StrippedCellLineName","OncotreeLineage", "OncotreeSubtype")]
  colnames(df_designAll)<- c("Samples", "Group", "Rep")
  df_designAll$Rep[df_designAll$Rep==""]<-"Unspecified"
  df_designAll<-df_designAll[order(df_designAll$Rep, decreasing = T),]
  df_designAll<-df_designAll[order(df_designAll$Group, decreasing = T),]
  
  # Disease selection
  cellsDisease <- cells$StrippedCellLineName[cells$OncotreeSubtype %in% disease]
  df_designDisease<- df_designAll[df_designAll$Samples %in% cellsDisease,]
}

# Selected Data table
df_gene<- chronos[toupper(interestingGenes), df_designDisease$Samples]

df_geneStack<- stack(df_gene)
colnames(df_geneStack)<-c("Chronos", "CellLine")
df_geneStack$Gene<- rep(rownames(df_gene), ncol(df_gene))
df_geneStack$color<- "black"
df_geneStack$Gene<- as.factor(df_geneStack$Gene)
for (g in interestingGenes){
  df_geneStack$Gene<- relevel(df_geneStack$Gene, g)
}

# Boxplot
pdf(paste0("CRISPRscores_",disease,paste(substr(interestingGenes,1,1),collapse=""),".pdf"), width=12, height=2*length(interestingGenes))  
  g<-ggplot(df_geneStack, aes(x=Chronos, y=Gene, fill=Gene)) + geom_boxplot()+ geom_jitter(color=df_geneStack$color, size=2, alpha=0.9, height = 0) + 
    theme(axis.text.x = element_text(size=20, angle = 90, vjust = 0.5, hjust=1), axis.text.y = element_text(size = 20), axis.title.x = element_text(size = 20)) + labs(y = "", x="CRISPR score (Chronos)") +
    ggtitle(paste0(disease, " (DepMap)")) + theme(legend.position="none") + theme(plot.title = element_text(size=20, face="bold", hjust = 0.5)) +
    geom_vline(xintercept = 0, color = "blue", linewidth=0.5) + geom_vline(xintercept = -1, color = "red", linewidth=0.5)
  print(g+geom_text_repel(label=df_geneStack$CellLine, angle=90, max.overlaps = Inf))
dev.off()
