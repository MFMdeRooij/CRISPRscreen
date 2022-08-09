# Plot RNAseq data from DepMap.org (Broad Institute)
# Download RNAseq data (Expression_22Q2_Public)" from https://depmap.org/portal/download/custom/, and include "Add cell line metadata to download"
# All cell lines are plotted on page 1, and the hematological ones are plotted including labels on page 2
# Author: M.F.M. de Rooij PhD, Amsterdam UMC, Spaargaren Lab, 2022, info: m.f.derooij@amsterdamumc.nl
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
  dataDepMap<-read.csv("Expression_22Q2_Public.csv",stringsAsFactors = F)
  RNAseq<-as.data.frame(t(dataDepMap[,c(7:ncol(dataDepMap))]))
  colnames(RNAseq)<- dataDepMap$cell_line_display_name
  cells<- dataDepMap[,c(2,1,3:6)]
  cellsHemato <- cells$cell_line_display_name[cells$lineage_1 %in% c("Blood", "Lymphocyte", "Plasma Cell")]
}

# Sample metadata
df_designAll<-cells[,c(1,3,4)]
colnames(df_designAll)<- c("Samples", "Group", "Rep")
df_designAll$Rep[df_designAll$Rep==""]<-"Unspecified"
df_designAll<-df_designAll[order(df_designAll$Rep, decreasing = T),]
df_designAll<-df_designAll[order(df_designAll$Group, decreasing = T),]
df_designHemato<- df_designAll[df_designAll$Group %in% c("Plasma Cell", "Lymphocyte", "Blood"),]

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
