## RNAseq gene level plots:
# - Make count tables with the GeneExpressionRNAseqHisat2.R script
# - Use the RNAseqCountTableNorTPX.csv count table
# - Fill in the Design table in RNAseqDesign.csv
# - Adjust the settings and run the script
# Author: M.F.M. de Rooij PhD, Amsterdam UMC, Spaargaren Lab, 2022, info: m.f.derooij@amsterdamumc.nl
######################################################################################
# # R packages
#install.packages(c("ggplot2", "factoextra"))
library("ggplot2")
library("factoextra")
#install.packages("devtools")
#devtools::install_github("JosephCrispell/basicPlotteR")
library("basicPlotteR")
#install.packages("BiocManager")
#BiocManager::install("DESeq2")
######################################################################################
#                                     SETTINGS

# Workdirectory (folder in which the count tables are located, use always slash (/) instead of backslash)
Workdirectory<- "H:/BioWin/RNAseq/"

# Which count table?
Filename<- "RNAseqCountTableNorTPX.csv"

# Is the count table sample normalized: 0 = Yes, 1 = No
Norm<- 0

# Title of dataset
Title<- "B cell lines"

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

# Sample normalize count table
if (Norm==1){
  library("DESeq2")
  numericColumns <- unlist(lapply(CountTableNor, is.numeric))  
  sf<- estimateSizeFactorsForMatrix(CountTableNor[,numericColumns])
  CountTableNor[,numericColumns]<- as.data.frame(round(t(t(CountTableNor[,numericColumns])/sf),1))
}

# PCA analysis
df_pr<-CountTableNor[3:ncol(CountTableNor)]
rownames(df_pr)<- CountTableNor$ensembl_gene_id
df_pr<- log2(df_pr+1)
df_pr$meanExpr<- apply(df_pr, 1, FUN=mean)
df_prSel <- df_pr[df_pr$meanExpr>1,]
df_prSel$meanExpr<- NULL

# Generate plots
interestingGenes<-toupper(interestingGenes)
for (gene in interestingGenes){
  # data table
  res.pca<- prcomp(t(df_prSel))
  df_pca<- data.frame(PC1=res.pca$x[,1], PC2=res.pca$x[,2])
  df_gene<- CountTableNor[toupper(CountTableNor$hgnc_symbol)==gene,df_design$Sample]
  if (nrow(df_gene)>1){
    df_gene<-df_gene[which.max(rowMeans(df_gene)),]
  }
  df_gene<-as.data.frame(t(df_gene))
  colnames(df_gene)<-"Gene"
  df_pca<- cbind(df_pca, df_gene)
  df_pca$Group<-factor(df_design$Group, levels=unique(df_design$Group))
  df_pca$color<- "white"
  #df_pca$color<- as.numeric(df_pca$Group)
  if (GRC==0){
    df_pca$Rep<-factor(df_design$Rep, levels=unique(df_design$Rep))
    df_pca$Group<-interaction(df_pca$Group,df_pca$Rep, lex.order = T)
  }
  
  pdf(paste0("GE_",gene,".pdf"), width=10, height=10)
    # Barplot
    g<-ggplot(df_pca, aes(x=Group, y=Gene, fill=Group)) + geom_boxplot() + geom_jitter(color=df_pca$color, size=2, alpha=0.9, height = 0) 
    g<-g + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + labs(x = "", y=paste0(gene," expression (nTPX)"))
    set.seed(101)
    print(g + ggtitle(Title) + theme(
      plot.title = element_text(size = 20, hjust = 0.5),
      axis.title.x = element_text(size = 20),
      axis.text.x = element_text(size = 20),
      axis.title.y = element_text(size = 20),
      axis.text.y = element_text(size = 20),
      legend.title=element_text(size=20),
      legend.text=element_text(size=20))
    )
    
    var1<-round(get_eig(res.pca)[1,2],1)
    var2<-round(get_eig(res.pca)[2,2],1)
    
    # PCA plot
    par(mar=c(4,5,2,1))
    par(fig=c(0.1,0.83,0.1,0.83))
    palette <- colorRampPalette(c("blue","red"))
    plot(df_pca$PC1,df_pca$PC2, pch=16, cex=2, col=palette((max(df_pca$Gene)-min(df_pca$Gene))*100+1)[(df_pca$Gene-min(df_pca$Gene))*100+1],
         main=paste0(gene, " expression"), xlab=paste0("PC1 (", var1,")"), ylab=paste0("PC2 (", var2,")"))
    addTextLabels(df_pca$PC1,df_pca$PC2,paste0(rownames(df_pca),"\n(", as.character(df_pca$Group),")"), avoidPoints = TRUE,
                  keepLabelsInside = TRUE, col.label=as.numeric(df_pca$Group), cex.label=0.5)
    x=1
    y=seq((min(df_pca$Gene)-min(df_pca$Gene))*100+1,(max(df_pca$Gene)-min(df_pca$Gene))*100+1,len=100)
    z=matrix(1:100,nrow=1)
    par(fig=c(0.8,0.95,0.1,0.83),new=TRUE)
    image(x,((y-1)/100)+min(df_pca$Gene),z,col=palette(max(y))[y], xaxt='n', xlab="",ylab=paste0(gene, " (nTPX)"))
  dev.off()
}