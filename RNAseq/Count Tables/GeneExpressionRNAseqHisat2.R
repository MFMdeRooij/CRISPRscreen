# Follow the lines with hastag signs in R studio
# Author: M.F.M. de Rooij PhD, Amsterdam UMC, Spaargaren Lab, 2019, info: m.f.derooij@amsterdamumc.nl
##################################################################################################################################
#                                                         SETTINGS
# Workdirectory
setwd("~/BioLin/RNAseq/")

# Select for protein-coding genes: 0 Yes, 1 = No
pcg <- 0

# Use effective gene lengths: 0 Yes, 1 = No, if yes, what is the mean fragment size?
effLen <- 1
fragmentSize <- 500

# PCA analysis:  0 Yes, 1 = No 
pca <- 0

# Order samples: 0 = Yes, 1 = No (order MapSamples.txt (together with existing samples) in  AllSamplesSorted.txt)
order <- 1
##################################################################################################################################
# Download, install, and configure the sra-toolkit from NCBI website:
# https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=software

# Install Linux Tools on command line:
# sudo apt-get install seqtk
# sudo apt-get install samtools
# sudo apt-get install hisat2

# Build genome
# Download genome in fasta format from NCBI, UCSC, or Ensembl (hg38.fa.gz), unzip and add in folder ~/HumanGenome
# cd ~/HumanGenome
# hisat2-build hg38.fa hg38

# Download GTF file with gene loci from Ensembl download page (Homo_sapiens.GRCh38.105.gtf.gz), unzip, rename to hg38.105.gtf.gz, and add in folder ~/HumanGenome
# Transform the gtf file into a known splice site text file using the python script delivered with hisat2 (in the bin folder)
# python hisat2_extract_splice_sites.py hg38.105.gtf > ~/HumanGenome/hg38.105_spliceSites.txt

# # Make tables with gene symbols and protein coding genes from the GTF file:
# gtf<- read.table("~/HumanGenome/hg38.105.gtf", fill=T)
# df_gene_id<- gtf[,c(10,15,16,22)]
# df_gene_id<-df_gene_id[df_gene_id$V15=="gene_name",]
# df_gene_id$V15<-NULL
# df_gene_protein<-df_gene_id[df_gene_id$V22=="protein_coding",]
# # Add B and T cell receptor constant genes
# df_gene_protein<-rbind(df_gene_protein,df_gene_id[df_gene_id$V22=="IG_C_gene",])
# df_gene_protein<-rbind(df_gene_protein,df_gene_id[df_gene_id$V22=="TR_C_gene",])
# df_gene_id$V22<-NULL
# df_gene_protein$V22<-NULL
# colnames(df_gene_id)<-c("ensembl_gene_id", "hgnc_symbol")
# colnames(df_gene_protein)<-c("ensembl_gene_id", "hgnc_symbol")
# write.csv(df_gene_id, "~/HumanGenome/hg38.105.GeneAnnotations.csv", row.names = F, quote = F)
# write.csv(df_gene_protein, "~/HumanGenome/hg38.105.ProteinCodingGenes.csv", row.names = F, quote = F)
# # Later you can also save a table with gene lengths

# # R Packages
# install.packages("BiocManager")
# BiocManager::install(c("Rsubread", "DESeq2"), dependencies=T)
# For PCA analysis, you have to install more r-packages, see at the end of this script

# Now you are ready for RNAseq analysis

# The first time, make the bash script excutable (chmod 755 MapRNAseqHumanSRA.sh)
# Search for RNAseq data on https://www.ncbi.nlm.nih.gov/sra/, and fill in the SRR IDs, cell 
# line ID and a P (paired-end) or S (single-end) separated by commas (SRR8615345,NAMALWA,P), and run the 
# script (./MapRNAseqHumanSRA.sh on the command line) 
# You can also run the bash script here by uncommenting and run the next line:
#system("bash MapRNAseqHumanSRA.sh")

# With the sorted bam files, you can check the read mapping in IGV-viewer

# Make a count table (TPM with DESEq2 normalisation (median of ratios method) to make the read counts comparable between samples (nTPX))
# Store always the non-normalizered count table to add more samples
library("Rsubread")
library("DESeq2")
for (pair in c("P","S")){
  if (pair=="P") {
     if (length(list.files("pairedEnd/", ".sam"))>0) {
        fc<- featureCounts(files=Sys.glob("pairedEnd/*.sam"),
                          annot.ext="~/HumanGenome/hg38.105.gtf",
                          isGTFAnnotationFile=TRUE,
                          isPairedEnd=TRUE,
                          countMultiMappingReads=T,
                          strandSpecific=0,
                          allowMultiOverlap=F)
      } else {
          next
        }
  }
  if (pair=="S") {
    if (length(list.files("singleEnd/", ".sam"))>0) {
      fc<- featureCounts(files=Sys.glob("singleEnd/*.sam"),
                        annot.ext="~/HumanGenome/hg38.105.gtf",
                        isGTFAnnotationFile=TRUE,
                        isPairedEnd=FALSE,
                        countMultiMappingReads=T,
                        strandSpecific=0,
                        allowMultiOverlap=F)
    } else {
        break
      }
  } 
  df_annotation<- as.data.frame(fc[2])
  CountTableRaw<- as.data.frame(fc[1])
  colnames(CountTableRaw)<- unlist(lapply(strsplit(colnames(CountTableRaw), ".sam"), "[[", 1))
  colnames(CountTableRaw)<- unlist(lapply(strsplit(colnames(CountTableRaw), "_"), "[[", 2))
  CountTableRaw$ensembl_gene_id<- rownames(CountTableRaw)
  
  if (file.exists("RNAseqCountTableRawAll.csv")) {
    # Add to existing count table
    CountTableOld<- read.csv("RNAseqCountTableRawAll.csv", stringsAsFactors = F)
    CountTableRaw<- merge(CountTableOld, CountTableRaw, by='ensembl_gene_id')
  } else {
      # Add gene symbols
      df_gene_id<- read.csv('~/HumanGenome/hg38.105.GeneAnnotations.csv', stringsAsFactors = F)
      CountTableRaw<- merge(df_gene_id, CountTableRaw, by="ensembl_gene_id", all.y=TRUE)
    }
  if (order==0){
    dfSort<-read.table("AllSamplesSorted.txt", sep=",")
    dfSort<-dfSort[dfSort$V2 %in% colnames(CountTableRaw[-1:-2]),]
    CountTableRaw<- CountTableRaw[,c("ensembl_gene_id","hgnc_symbol",dfSort$V2)]
  }
  # The raw counts can be used for DESeq2 analysis
  write.csv(CountTableRaw, "RNAseqCountTableRawAll.csv", row.names=FALSE)
  # When poly-A captured, remove all other genes than protein-coding genes
  if (pcg==0){
    df_gene_protein<-read.csv("~/HumanGenome/hg38.105.ProteinCodingGenes.csv", header=T, stringsAsFactors = F)
    CountTableRawProteinCoding<-CountTableRaw[CountTableRaw$ensembl_gene_id %in% df_gene_protein$ensembl_gene_id,]
    write.csv(CountTableRawProteinCoding, "RNAseqCountTableRawProteinCoding.csv", row.names=FALSE)
  }
}  

# Transform counts to Transcripts per million
#CountTableRaw<- read.csv("RNAseqCountTableRawAll.csv", stringsAsFactors = F)
CountTableTPM<- CountTableRaw
df_geneLength<-data.frame(ensembl_gene_id=df_annotation$annotation.GeneID, Length=df_annotation$annotation.Length)

# # Save gene lengths
# write.csv(df_geneLength, "~/HumanGenome/hg38.105.GeneLengths.csv", row.names=F)
# df_geneLength<-read.csv("~/HumanGenome/hg38.105.GeneLengths.csv", stringsAsFactors = F)

CountTableTPM<- merge(CountTableTPM,df_geneLength, by="ensembl_gene_id", all.x=T)

# Effective gene length is gene length minus fragment length (because the gene ends contain less reads)
if (effLen==0) {
  CountTableTPM$Length<- CountTableTPM$Length-fragmentSize
  CountTableTPM$Length[CountTableTPM$Length < fragmentSize] <- fragmentSize
}
# When poly-A captured, remove all other genes than protein-coding genes
if (pcg==0){
  df_gene_protein<-read.csv("~/HumanGenome/hg38.105.ProteinCodingGenes.csv", header=T, stringsAsFactors = F)
  CountTableTPM<-CountTableTPM[CountTableTPM$ensembl_gene_id %in% df_gene_protein$ensembl_gene_id,]
}

# # Check uniqueness of gene symbols
# x<-data.frame(g=CountTableTPM$hgnc_symbol, c=1)
# y<-aggregate(x$c, by=list(x$g), FUN=sum)
# y<-y[order(y$x, decreasing = T),]
# z<-CountTableTPM[CountTableTPM$hgnc_symbol %in% y$Group.1[y$x>1],]
# z<-z[order(z$hgnc_symbol),]

# Calculate TPM
CountTableTPM[,-1:-2]<- CountTableTPM[,-1:-2]/(CountTableTPM$Length)
CountTableTPM$Length<-NULL
CountTableTPM[,-1:-2]<- as.data.frame(t(t(CountTableTPM[,-1:-2])/colSums(as.data.frame(CountTableTPM[,-1:-2]))*1000000))
write.csv(CountTableTPM, "RNAseqCountTableTPM.csv", row.names=FALSE)

# Sample normalization (TPM -> between sample normalization -> nTPX)
#CountTableTPM <- read.csv("RNAseqCountTableTPM.csv", stringsAsFactors = F)
sf<- estimateSizeFactorsForMatrix(CountTableTPM[,-1:-2])
CountTableNorTPX<- CountTableTPM
CountTableNorTPX[,-1:-2]<- as.data.frame(round(t(t(CountTableTPM[,-1:-2])/sf),1))
write.csv(CountTableNorTPX, "RNAseqCountTableNorTPX.csv", row.names=FALSE)

if (pca==0) {
  #PCA & UMAP
  
  # Add RNAseqDesign.csv in the workdirectory (this should include a column with Sample (Cell line names)
  # and a culumn with Group (Cell types)
  
  #install.packages(c("corrplot", "ggplot2", "ggfortify", "factoextra", "umap", "basicPlotteR"))
  library("corrplot")
  library("ggplot2")
  library("ggfortify")
  library("factoextra")
  library("umap")
  library("basicPlotteR")
  
  CountTableNor<- read.csv("RNAseqCountTableNorTPX.csv")
  
  # Sample metadata
  df_design<- read.csv("RNAseqDesign.csv")
  df_design$Color <- rainbow(length(unique(as.factor(df_design$Group))))[as.factor(df_design$Group)]
  
  # Select samples from design table
  CountTableNor<- CountTableNor[,c("ensembl_gene_id", "hgnc_symbol", colnames(CountTableNor)[colnames(CountTableNor) %in% df_design$Sample])]
  
  # PCA analysis
  df_pr<-CountTableNor[3:ncol(CountTableNor)]
  rownames(df_pr)<- CountTableNor$ensembl_gene_id
  df_pr<- log2(df_pr+1)
  df_pr$meanExpr<- apply(df_pr, 1, FUN=mean)
  df_prSel <- df_pr[df_pr$meanExpr>1,]
  df_prSel$meanExpr<- NULL
  
  res.pca<- prcomp(t(df_prSel))
  res.ind <- get_pca_ind(res.pca)
  res.eig<- get_eigenvalue(res.pca)
  pcaind<- as.data.frame(res.ind$coord)
  
  # UMAP
  set.seed(100)
  rna.umap<- umap(t(df_prSel))
  
  pdf("PCA_UMAP.pdf", width=15, height=10)
    # Scree plot
    print(fviz_eig(res.pca))
    
    # PCA plot
    par(mar=c(5,5,5,20), xpd=T)
    plot(pcaind$Dim.1, pcaind$Dim.2, pch=16, col=df_design$Color, xlab=paste0("PC1 (",round(res.eig$variance.percent[1],1),"%)"), 
         ylab=paste0("PC2 (",round(res.eig$variance.percent[2],1),"%)"), main = "PCA")
    addTextLabels(pcaind$Dim.1,pcaind$Dim.2,colnames(df_prSel), avoidPoints = TRUE,
                  keepLabelsInside = TRUE, col.label="gray", cex.label=1)
    legend("topright", inset=c(-0.3,0), legend=unique(df_design$Group), col=unique(df_design$Color), pch=16, cex=1, title="Cell Type")
    
    # UMAP plot
    par(mar=c(5,5,5,20), xpd=T)
    plot(rna.umap$layout, pch=16, col=df_design$Color, xlab="UMAP1", ylab="UMAP2", main="UMAP")
    addTextLabels(rna.umap$layout[,1],rna.umap$layout[,2], colnames(df_prSel), avoidPoints = TRUE,
                  keepLabelsInside = TRUE, col.label="gray", cex.label=1)
    legend("topright", inset=c(-0.3,0), legend=unique(df_design$Group), col=unique(df_design$Color), pch=16, cex=1, title="Cell Type")
    
    # Cor plot
    par(mar=c(5,5,5,5))
    cor<- cor(df_prSel, method = "pearson")
    print(corrplot(cor, method='color', order = "hclust"))
    
    # Dendrogram
    # You can adjust k to the number of groups
    res<- hcut(t(df_prSel), k = ncol(df_prSel)-1, stand = TRUE)
    print(fviz_dend(res, rect = TRUE, cex = 0.5))
  dev.off()
} 
