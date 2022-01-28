## RNAseq differential expression:
# - Make count tables with the GeneExpressionRNAseqHisat2.R script
# - Use the RNAseqCountTableRawProteinCoding.csv count table
# - This script also uses the CRISPRScreenAnalysisLibraries.csv file to highlight the essentials/nonessential genes, 
#   you can change these genes by for example certain pathway genes; this script also extract the kinome in the end
# - Fill in the Design table in RNAseqDesign.csv
# - Adjust the settings and run the script
# Author: M.F.M. de Rooij PhD, Amsterdam UMC, Spaargaren Lab, 2022, info: m.f.derooij@amsterdamumc.nl
######################################################################################
## Install DESeq2
#install.packages("BiocManager")
#install.packages("gtools")
#BiocManager::install(c("DESeq2", "pheatmap", "RColorBrewer"), dependencies=TRUE)
library("DESeq2")
library("pheatmap")
library("RColorBrewer")
library("gtools")
######################################################################################
#                                     SETTINGS

# Workdirectory (folder in which the count tables are located, use always slash (/) instead of backslash)
# Windows
Workdirectory<- "H:/BioWin/RNAseq/"
# Linux
#Workdirectory <- "~/BioLin/RNAseq/"

# Which count table?
Filename <- "RNAseqCountTableRawProteinCoding.csv"

# Fill in the table in RNAseqDesign.csv (Group = tumor-subtypes, Rep = replicates (when paired, this should be matched))

# Round numbers in output tables: 0 = Yes, 1 = No
RoundNumbers <- 0

# Paired replicates: 0 = Paired, 1 = Unpaired
Paired <- 0

# Minimal fold change of guides to be a hit: 1 = No minimal fold change,  >1: The minimal fold change (linear scale, 2^abs(l2fc))
minimalFoldChange <- 1

# MA plots of all genes: 0 = Yes, 1 = No, 2 = Top 10, 3 = Genes of interest
MA_all_genes <- 3
# Genes of interest
if (MA_all_genes==3){
  q <- function(...) {
    sapply(match.call()[-1], deparse)
  }
  interestingGenes <-q(
    # Which genes?  
    BTK, SYK, PIK3CA, PIK3CD, PIK3R1, LYN
  )
}
######################################################################################
setwd(Workdirectory)

# Make a data folder
dirname2<- paste0(Filename,Sys.time())
dirname1<-gsub("[[:punct:]]", "", dirname2) 
dirname<- gsub("\\s", "", dirname1) 
dir.create(dirname)

# Read count table
df_raw <- read.csv(file=Filename, sep=",", header=TRUE, stringsAsFactors = FALSE)

# Gene annotation
df_Gene_ID<- df_raw[,1:2]  

# Sample metadata
df_design<-read.csv("RNAseqDesign.csv")
df_colData <- data.frame(Group=df_design$Group, Rep=df_design$Rep) 
df_colData$Group <- as.factor(df_colData$Group)
df_colData$Rep <- as.factor(df_colData$Rep)

# Count table for DESEq2
counts <- df_raw[,-1:-2]
rownames(counts) <- df_raw[,1]
counts2<-counts[,df_design$Sample]

# DESeq2 pipeline
if (Paired==0) {
  dds <- DESeqDataSetFromMatrix(countData = counts, colData = df_colData, design = ~ Rep + Group)
}
if (Paired==1) {
  dds <- DESeqDataSetFromMatrix(countData = counts, colData = df_colData, design = ~ Group)
}  
dds <- DESeq(dds, betaPrior = TRUE)

# Define comparisons
Groups<-unique(df_design$Group)
combi<-as.data.frame(combinations(length(Groups),2,Groups))

# Analyse comparisons
for (i in 1:nrow(combi)){
  con1 <- combi[i,1]
  con2 <- combi[i,2]
  res <- results(dds, contrast=c("Group",con1,con2), addMLE=TRUE)
  
  # DESeq2 data table
  df_res <- as.data.frame(res)
  
  df_baseMeanPerLvl <- as.data.frame(sapply(levels(dds$Group), function(lvl) rowMeans(counts(dds, normalized=TRUE)[,dds$Group==lvl])))
  df_res$BaseMeanA <- df_baseMeanPerLvl[[con1]]
  df_res$logBaseMeanA <- log(df_baseMeanPerLvl[[con1]]+1)/log(10)
  df_res$BaseMeanB <- df_baseMeanPerLvl[[con2]]

  # GeneIDs
  df_res$ensembl_gene_id <- rownames(df_res)
  df_res <- merge(df_Gene_ID, df_res, by='ensembl_gene_id', all.y=T)
  df_res$log2FoldChange[is.na(df_res$log2FoldChange)] < -0
  df_res$FoldChange <- 2^df_res$log2FoldChange
    
  if (minimalFoldChange > 1) {
    df_res$padj[!is.na(df_res$padj) & df_res$FoldChange >= (1/minimalFoldChange) & df_res$FoldChange <= minimalFoldChange] <- 1
  }
  
  # CRISPR screen controls
  df_Gene_IDX<- read.csv("CRISPRScreenAnalysisLibraries.csv", sep=',', header=TRUE, stringsAsFactors = FALSE)
  df_Control<- df_Gene_IDX[1:927,1:3]
  PC<- "Essential"
  NC<- "NonEssential"
  df_res$Type<- "x"
  df_res[df_res$hgnc_symbol %in% df_Control[[PC]][nchar(df_Control[[PC]])>0],"Type"]<- "p"
  df_res[df_res$hgnc_symbol %in% df_Control[[NC]][nchar(df_Control[[NC]])>0],"Type"]<- "n"
  df_PC<-df_res[df_res$Type=="p",]
  df_NC<-df_res[df_res$Type=="n",]

  # Hits
  GenesDiff<-df_res[!is.na(df_res$padj),]
  GenesDiff<-GenesDiff[GenesDiff$padj<0.1,]
  GenesDiffDep<-GenesDiff[order(GenesDiff$log2FoldChange),]
  GenesDiffDep<-GenesDiffDep[nchar(GenesDiffDep$hgnc_symbol)>0,]
  GenesDiffEnr<-GenesDiff[order(GenesDiff$log2FoldChange, decreasing = T),]
  GenesDiffEnr<-GenesDiffEnr[nchar(GenesDiffEnr$hgnc_symbol)>0,]
  tophits<- c( GenesDiffDep[1:10,2], GenesDiffEnr[1:10,2])
    
  # Write gene table
  
  # Guide Table
  df_res_print<- df_res[,c("hgnc_symbol","ensembl_gene_id","Type","BaseMeanA","BaseMeanB","FoldChange","pvalue","padj")]
  df_res_print<- df_res_print[order(df_res_print$FoldChange),]
  if (RoundNumbers==0){
    df_res_print[,c("BaseMeanA","BaseMeanB")]<-round(df_res_print[,c("BaseMeanA","BaseMeanB")],0)
    df_res_print[,c("FoldChange","pvalue","padj")]<-signif(df_res_print[,c("FoldChange","pvalue","padj")],3)
  }                                                                                       
  write.csv(df_res_print, paste0(dirname,"/DESeq2_RNAseq_",con1, "vs",con2,".csv"), row.names = FALSE)
  
  # Write Plots in PDF 
  pdf(paste0(dirname,"/DESeq2_RNAseq_",con1, "vs",con2,"_Plots.pdf"), width=10, height=10)
  
    # DESEq2 plots
    plotDispEsts(dds, main="Dispersion plot")
    
    rld<-rlog(dds, blind=FALSE)
    sampleDists<- dist(t(assay(rld)))
    sampleDistMatrix <- as.matrix(sampleDists)
    rownames(sampleDistMatrix) <- paste(rld$Group, rld$Rep)
    colnames(sampleDistMatrix) <- paste(rld$Group, rld$Rep)
    colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
    
    pheatmap(sampleDistMatrix,
             clustering_distance_rows=sampleDists,
             clustering_distance_cols=sampleDists,
             col=colors)
    
    print(plotPCA(rld, intgroup=c("Group", "Rep")))

    # MA plots
    xrange<-c(min(df_res$logBaseMeanA, na.rm=TRUE)-0.5, max(df_res$logBaseMeanA, na.rm=TRUE)+0.5)
    yrange<-c(min(df_res$log2FoldChange, na.rm=TRUE)-0.5, max(df_res$log2FoldChange, na.rm=TRUE)+0.5)
    
    den_tot<-density(df_res$log2FoldChange, from=min(yrange[1]), to=max(yrange[2]), na.rm=TRUE)
    den_PC<-density(df_PC$log2FoldChange, from=min(yrange[1]), to=max(yrange[2]), na.rm=TRUE)
    den_NC<-density(df_NC$log2FoldChange, from=min(yrange[1]), to=max(yrange[2]), na.rm=TRUE)
    denMax <- max(c(den_tot$y, den_PC$y, den_NC$y))
    
    
    Gene_of_interest<-c(" ", "Hitlist", if (MA_all_genes==0) {df_res$hgnc_symbol}, if (MA_all_genes==2) {tophits}, if (MA_all_genes==3) {interestingGenes})
    for(Gene in Gene_of_interest){
      df_GOI<-df_res[df_res$hgnc_symbol %in% Gene,]
      
      par(mar=c(4,5,2,1))
      par(fig=c(0.1,0.83,0.1,0.83))
      plot(df_res$logBaseMeanA, df_res$log2FoldChange, type="p", col="gray", cex=.7, pch=16, xlab="Log10 Average Read Counts (Control)", 
           ylab="Log2 Fold Change", cex.lab=1.5, cex.axis=1.3, xlim=xrange, ylim=yrange, xaxp = c(0, 10, 10), 
           yaxp = c(-20, 20, 20))
  
      points(df_PC$logBaseMeanA, df_PC$log2FoldChange, type="p", col="red", cex=1, pch=15)
      points(df_NC$logBaseMeanA, df_NC$log2FoldChange, type="p", col="blue", cex=1, pch=17)
      if (Gene=="Hitlist"&& nrow(GenesDiff)>=1){
        points(GenesDiff$logBaseMeanA, GenesDiff$log2FoldChange, type="p", col=1, cex=0.7, pch=19)
      }
      if (nrow(df_GOI)>=1){
        points(df_GOI$logBaseMeanA, df_GOI$log2FoldChange, type="p", col=1, cex=1.5, pch=19)
      }
      legend("topright",legend=c("Total", "Essential", "Non-Essential", if (nchar(Gene)>1 && !is.na(Gene)) {Gene}), cex=1, pch=c(16,15,17,if (nchar(Gene)>1 && !is.na(Gene)) {19}), col=c("gray","red","blue",if (nchar(Gene)>1 && !is.na(Gene)) {1}))
      abline(median(df_res$log2FoldChange, na.rm=TRUE),0, col=1, lty=3, untf=TRUE)
      
      # Density plot fold change
      par(fig=c(0.71,1,0.1,0.83),new=TRUE)
      plot(den_tot$y, den_tot$x, ylim=range(yrange), xlim=(c(0,denMax)), type='l', axes=FALSE, col="gray", xlab="", 
           ylab="", lwd=2)
      par(new=TRUE)
      lines(den_PC$y, den_PC$x, col="red", lwd=2)
      lines(den_NC$y, den_NC$x, col="blue", lwd=2)
      if (Gene=="Hitlist" && nrow(GenesDiff)>1){
        den_hits_total<-density(GenesDiff$log2FoldChange, from=min(yrange[1]), to=max(yrange[2]), na.rm=TRUE)
        lines(den_hits_total$y, den_hits_total$x, col=1, lwd=2)
      }
      if (nrow(df_GOI)>1){
        par(new=TRUE)
        den_GOI<-density(df_GOI$log2FoldChange, from=min(yrange[1]), to=max(yrange[2]), na.rm=TRUE)
        lines(den_GOI$y, den_GOI$x, col=1, lwd=2)
      }
      # Density plot read count
      denX_tot<-density(df_res$logBaseMeanA, from=min(xrange[1]), to=max(xrange[2]), na.rm=TRUE)
      denX_PC<-density(df_PC$logBaseMeanA, from=min(xrange[1]), to=max(xrange[2]), na.rm=TRUE)
      denX_NC<-density(df_NC$logBaseMeanA, from=min(xrange[1]), to=max(xrange[2]), na.rm=TRUE)
      denXMax <- max(c(denX_tot$y, denX_PC$y, denX_NC$y))
      
      par(fig=c(0.1,0.83,0.71,1),new=TRUE)
      denX_tot<-density(df_res$logBaseMeanA, from=min(xrange[1]), to=max(xrange[2]), na.rm=TRUE)
      plot(denX_tot$x, denX_tot$y, xlim=range(xrange), ylim=c(0,denXMax), type='l', axes=FALSE, col="gray", xlab="", 
           ylab="", lwd=2)
      lines(denX_PC$x, denX_PC$y, col="red", lwd=2)
      lines(denX_NC$x, denX_NC$y, col="blue", lwd=2)
      if (Gene=="Hitlist" && nrow(GenesDiff)>1){
        par(new=TRUE)
        denX_hits_total<-density(GenesDiff$logBaseMeanA, from=min(xrange[1]), to=max(xrange[2]), na.rm=TRUE)
        lines(denX_hits_total$x, denX_hits_total$y, col=1, lwd=2)
      }
      if (nrow(df_GOI)>1){
        par(new=TRUE)
        denX_GOI<-density(df_GOI$logBaseMeanA, from=min(xrange[1]), to=max(xrange[2]), na.rm=TRUE)
        lines(denX_GOI$x, denX_GOI$y, col=1, lwd=2)
      }
    }
  dev.off()
  
  # Kinome(-like) genes taken from Brunello CRISPR library
  kinome<- unique(df_Gene_IDX$GeneSymbol)
  kinome<-kinome[1:764]
  dataPrintKin<-df_res_print[df_res_print$hgnc_symbol %in% kinome,]
  write.csv(dataPrintKin, paste0(dirname,"/DESeq2_RNAseq_kinome",con1, "vs",con2,".csv"), row.names=F, quote=F)
}
