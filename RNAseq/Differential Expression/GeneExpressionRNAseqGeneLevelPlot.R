## RNAseq gene level plots:
# - Make count tables with the GeneExpressionRNAseqHisat2.R script
# - Perform DESeq2 analysis with the GeneExpressionRNAseqDESeq2.R script 
# - Use the NormalizedCounts_rlog.csv table derived from the DESeq2 script
# - Fill in the Design table in RNAseqDesign.csv
# - Adjust the settings and run the script
# Author: M.F.M. de Rooij PhD, Amsterdam UMC, Spaargaren Lab, 2022, info: m.f.derooij@amsterdamumc.nl
######################################################################################
# # R packages
#install.packages(c("ggplot2", "factoextra","RColorBrewer", "matrixStats"))
library("ggplot2")
library("factoextra")
library("RColorBrewer")
library('matrixStats')
#install.packages("devtools")
#devtools::install_github("JosephCrispell/basicPlotteR")
library("basicPlotteR")
#install.packages("BiocManager")
#BiocManager::install("DESeq2")
######################################################################################
#                                     SETTINGS

# Put this script in the folder where the count tables are located
Workdirectory<- dirname(rstudioapi::getActiveDocumentContext()$path)
## Fill in workdirectory (folder in which the count tables are located, use always slash (/) instead of backslash)
#Workdirectory<- "C:/BioWin/RNAseq/"

# Which count table?
Filename<- "NormalizedCounts_rlog.csv"

# Title of dataset
Title<- "B cell lines"

# Start y axis of boxplot (gene expression) at 0: 0 = yes, 1 = no
zero<- 0

# Use Group and Rep as combined factors (Like different cells (-> Group) and different drugs (-> Rep)): 0 = Yes, 1 = No
GRC<- 1

# Genes of interest
q <- function(...) {
  sapply(match.call()[-1], deparse)
}
interestingGenes <- toupper(q(
  # Which genes?  
  BTK, SYK, PIK3CA, PIK3CD, PIK3R1, LYN
))
##################################################################################
setwd(Workdirectory)
CountTableNor<-read.csv(file=Filename, sep=",", header=TRUE, stringsAsFactors = FALSE)

# Sample metadata
df_design<-read.csv("RNAseqDesign.csv", stringsAsFactors = F)

# Select samples from design table
CountTableNor<- CountTableNor[,c("ensembl_gene_id", "hgnc_symbol", colnames(CountTableNor)[colnames(CountTableNor) %in% df_design$Sample])]

# PCA analysis
df_pr<-CountTableNor[3:ncol(CountTableNor)]
rownames(df_pr)<- CountTableNor$ensembl_gene_id

variances <- rowVars(as.matrix(df_pr))
names(variances) <- rownames(df_pr)
top_500_indices <- order(variances, decreasing = TRUE)[1:500]
df_prSel <- df_pr[top_500_indices, ]

# # Optioneel: Bekijk de namen van de top 10 meest variabele genen
head(sort(variances, decreasing = TRUE), 10)

# Generate plots
interestingGenes<- interestingGenes[toupper(interestingGenes) %in% toupper(CountTableNor$hgnc_symbol)]

for (gene in interestingGenes){
  # data table
  res.pca<- prcomp(t(df_prSel))
  df_pca<- data.frame(PC1=res.pca$x[,1], PC2=res.pca$x[,2])
  df_gene<- CountTableNor[toupper(CountTableNor$hgnc_symbol)==gene,df_design$Sample]
  if (nrow(df_gene)>1){
    df_gene<-df_gene[which.max(rowMeans(df_gene)),]
  }
  df_gene<- as.data.frame(t(df_gene))
  colnames(df_gene)<- "Gene"
  df_pca<- cbind(df_pca, df_gene)
  df_pca$Group<- factor(df_design$Group, levels=unique(df_design$Group))
  
  if (GRC==0){
    df_pca$Rep<- factor(df_design$Rep, levels=unique(df_design$Rep))
    df_pca$Group<- interaction(df_pca$Group,df_pca$Rep, lex.order = T)
  }
  
  df_pca$Color<- rep(brewer.pal(8,"Dark2"),100)[df_pca$Group]
  color<- unique(df_pca[,c("Group", "Color")])
  colorVec<- as.character(color$Color)
  names(colorVec) <- color$Group
  
  pdf(paste0("GE_",gene,".pdf"), width=7, height=7)
  # Boxplot
  g<- ggplot(df_pca, aes(x=Group, y=Gene, fill=Group)) + theme_classic() + scale_fill_manual(values = colorVec) + geom_boxplot(color="snow3", alpha=0.05, show.legend = FALSE) + 
    geom_jitter(color=df_pca$Color, size=3, alpha=1, height = 0, show.legend = FALSE) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + labs(x = "", y=paste0(gene," expression (rlog)"))
  if(zero==0) {
    g<- g + expand_limits(y = 0)
  }
  set.seed(101)
  print(g + ggtitle(Title) + theme(
    plot.title = element_text(size = 20, hjust = 0.5, face = "bold"),
    axis.title.x = element_text(size = 20),
    axis.text.x = element_text(angle=45, , vjust = 1, hjust=1, size = 20, color = colorVec),
    axis.title.y = element_text(size = 20),
    axis.text.y = element_text(size = 20)
  )
  )
  
  var1<-round(get_eig(res.pca)[1,2],1)
  var2<-round(get_eig(res.pca)[2,2],1)
  
  # PCA plot
  par(mar=c(5,5,5,1))
  par(bty="l", lwd=2)
  par(fig=c(0,0.83,0,0.83))
  palette <- colorRampPalette(c("dodgerblue","red2"))
  plot(df_pca$PC1,df_pca$PC2, pch=16, cex=1, col=palette((max(df_pca$Gene)-min(df_pca$Gene))*100+1)[(df_pca$Gene-min(df_pca$Gene))*100+1],
       main=paste0(gene, " expression"), xlab=paste0("PC1 (", var1,")"), ylab=paste0("PC2 (", var2,")"), cex.main=1.6, cex.axis=1.5, cex.lab=1.5)
  addTextLabels(df_pca$PC1,df_pca$PC2,paste0(rownames(df_pca),"\n(", as.character(df_pca$Group),")"), avoidPoints = TRUE,
                keepLabelsInside = TRUE, col.label=df_pca$Color, cex.label=0.5)
  x=1
  y=seq((min(df_pca$Gene)-min(df_pca$Gene))*100+1,(max(df_pca$Gene)-min(df_pca$Gene))*100+1,len=100)
  z=matrix(1:100,nrow=1)
  par(fig=c(0.8,1,0,0.83),new=TRUE)
  image(x,((y-1)/100)+min(df_pca$Gene),z,col=palette(max(y))[y], xaxt='n', xlab="",ylab=paste0(gene, " expression (rlog)"))
  dev.off()
}