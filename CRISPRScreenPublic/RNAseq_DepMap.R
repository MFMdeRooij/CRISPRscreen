# Plot RNAseq data from DepMap.org (Broad Institute)
# Download RNAseq data (22q4) OmicsExpressionProteinCodingGenesTPMLogp1.csv, and model.csv from https://depmap.org/portal/download/all/
# All cell lines are plotted on page 1, and the lymphoid ones are plotted including labels on page 2
# Author: M.F.M. de Rooij PhD, Amsterdam UMC, Spaargaren Lab, 2023, info: m.f.derooij@amsterdamumc.nl
##################################################################################
# Workdirectory (folder where "Expression_22Q2_Public.csv" is located)
setwd("H:/BioWin/DepMap/")

# Genes of interest (add to line 15)
q <- function(...) {
  sapply(match.call()[-1], deparse)
}
interestingGenes <- toupper(q(
  # Which genes?  
  BTK, PIK3CD, CD79A, BCL6, PRDM1
))
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
  dataDepMap<-read.csv("OmicsExpressionProteinCodingGenesTPMLogp1.csv",stringsAsFactors = F)
  RNAseq<-as.data.frame(t(dataDepMap[,2:ncol(dataDepMap)]))
  colnames(RNAseq)<- dataDepMap$X
  rownames(RNAseq)<- unlist(lapply(strsplit(as.character(rownames(RNAseq)),"\\.\\."), "[", 1))
  cells<- read.csv("Model.csv",stringsAsFactors = F)
  UsedIDs<- cells[match(colnames(RNAseq), cells$ModelID),]
  colnames(RNAseq)<- UsedIDs$StrippedCellLineName
  
  # Sample metadata
  df_designAll<-cells[cells$StrippedCellLineName %in% UsedIDs$StrippedCellLineName,c("StrippedCellLineName","OncotreeLineage", "OncotreeSubtype")]
  colnames(df_designAll)<- c("Samples", "Group", "Rep")
  df_designAll$Rep[df_designAll$Rep==""]<-"Unspecified"
  df_designAll<-df_designAll[order(df_designAll$Rep, decreasing = T),]
  df_designAll<-df_designAll[order(df_designAll$Group, decreasing = T),]
  
  # Lymphoid selection
  cellsHemato <- cells$StrippedCellLineName[cells$OncotreeLineage %in% c("Lymphoid")]
  df_designHemato<- df_designAll[df_designAll$Samples %in% cellsHemato,]
}

# Generate plots
interestingGenes<-toupper(interestingGenes)
for (gene in interestingGenes){
  # data table
  if (gene %in% rownames(RNAseq)){
    pdf(paste0("RNAseq_",gene,".pdf"), width=10, height=20)
    i<-1
    for (df_design in list(df_designAll, df_designHemato)){
      df_gene<- RNAseq[rownames(RNAseq)==gene,df_design$Sample]
      df_gene<-as.data.frame(t(df_gene))
      colnames(df_gene)<-"Gene"
      df_gene$Group<-factor(df_design$Group, levels=unique(df_design$Group))
      df_gene$color<- "black"
      df_gene$Rep<-factor(df_design$Rep, levels=unique(df_design$Rep))
      df_gene$Group<-interaction(df_gene$Group,df_gene$Rep, lex.order = T)
      
      # Barplot
      g<-ggplot(df_gene, aes(x=Gene, y=Group, fill=Group, label=rownames(df_gene))) + geom_boxplot() + geom_jitter(color=df_gene$color, size=2, alpha=0.9, height = 0) + 
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + labs(y = "", x="Log2 TPM+1") +
        ggtitle(paste0(gene, " expression (DepMap)")) + theme(legend.position="none") + theme(plot.title = element_text(hjust = 0.5))
      if (i==1){  
        print(g)
      } else {
        print(g+geom_text_repel(angle=90, max.overlaps = Inf))
      }
      i<-i+1
    }
    dev.off()
  }else{
    print(paste0(gene, " is not available"))
  }
}
