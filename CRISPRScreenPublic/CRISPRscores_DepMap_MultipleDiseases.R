# Plot CRISPR scores of the public CRISPR screens from DepMap.org (Broad Institute)
# Download CRISPR data (24q4) "CRISPRGeneEffect.csv, and model.csv from https://depmap.org/portal/download/all/
# All cell lines are plotted on page 1, and the lymphoid ones are plotted including labels on page 2
# Author: M.F.M. de Rooij PhD, Amsterdam UMC, Spaargaren Lab, 2023, info: m.f.derooij@amsterdamumc.nl
##################################################################################
# Workdirectory (folder where "CRISPR_(CRISPRGeneEffect.csv" and "model.csv" are located)
setwd("G:/divg/Pathologie-ADHESIE/ALL/G_ScreeningTeam/DepMap/")

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
#install.packages(c("ggplot2", "factoextra", "ggrepel", "RColorBrewer"))
library("ggplot2")
library("factoextra")
library("ggrepel")
library("RColorBrewer")
#install.packages("devtools")
#devtools::install_github("JosephCrispell/basicPlotteR")
library("basicPlotteR")
if (!exists("dataDepMap")) {
  dataDepMap<-read.csv("CRISPRGeneEffect.csv",stringsAsFactors = F)
  chronos<-as.data.frame(t(dataDepMap[,2:ncol(dataDepMap)]))
  colnames(chronos)<- dataDepMap$X
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
  
  # Lymphoid selection
  cellsHemato <- cells$StrippedCellLineName[cells$OncotreeLineage %in% c("Lymphoid")]
  df_designHemato<- df_designAll[df_designAll$Samples %in% cellsHemato,]
}

# Generate plots
interestingGenes<-toupper(interestingGenes)
for (gene in interestingGenes){
  # data table
  if (gene %in% rownames(chronos)){
    pdf(paste0("CRISPRscores_",gene,".pdf"), width=12, height=20)
    i<-1
    for (df_design in list(df_designAll, df_designHemato)){
      df_gene<- chronos[rownames(chronos)==gene,df_design$Sample]
      df_gene<- as.data.frame(t(df_gene))
      colnames(df_gene)<- "Gene"
      df_gene$Group<- factor(df_design$Group, levels=unique(df_design$Group))
      df_gene$Rep<- factor(df_design$Rep, levels=unique(df_design$Rep))
      df_gene$Group<- interaction(df_gene$Group,df_gene$Rep, lex.order = T)
      
      df_gene$Color<- rep(brewer.pal(8,"Dark2"),100)[df_gene$Rep]
      color<- unique(df_gene[,c("Rep", "Color")])
      colorVec<- as.character(color$Color)
      names(colorVec) <- color$Rep
      
      # Boxplot
      g<-ggplot(df_gene, aes(x=Gene, y=Group, fill=Rep)) + theme_classic() + scale_fill_manual(values = colorVec) + 
        geom_boxplot(color="snow3", alpha=0.05, show.legend = FALSE) + geom_jitter(color=df_gene$Color, size=2, alpha=1, height = 0, show.legend = FALSE) + 
        theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5), axis.text.y = element_text(size=8,color=colorVec)) + labs(y = "", x="CRISPR score (Chronos)") +
        ggtitle(paste0(gene, " knockout (DepMap)")) + theme(plot.title = element_text(hjust = 0.5)) +
        geom_vline(xintercept = 0, color = "blue", linewidth=0.5) + geom_vline(xintercept = -1, color = "red", linewidth=0.5)
      if (i==1){  
        print(g)
      } else {
        print(g + theme(axis.text.y = element_text(size = 10)) + geom_text_repel(label=rownames(df_gene), angle=90, box.padding = 0.3, max.overlaps = Inf, size=3))
      }
      i<-i+1
    }
    dev.off()
  }else{
    print(paste0(gene, " is not available"))
  }
}