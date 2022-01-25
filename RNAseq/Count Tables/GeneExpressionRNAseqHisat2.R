# Follow the lines with hastag signs in R studio
# Author: M.F.M. de Rooij PhD, Amsterdam UMC, Spaargaren Lab, 2019, info: m.f.derooij@amsterdamumc.nl
##################################################################################################################################
#                                                         SETTINGS
# Workdirectory
setwd("~/BioLin/RNAseq/")

# Select for protein-coding genes: 0 Yes, 1 = No
pcg <- 0

# Use effective gene lengths: 0 Yes, 1 = No, if yes, what is the mean fragment size?
effLen <- 0
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

# biomaRt dependencies:
# sudo apt-get install libcurl4-openssl-dev libssl-dev

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
                          strandSpecific=0)
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
                        strandSpecific=0)
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
  #PCA analysis
  #install.packages(c("corrplot", "ggplot2", "ggfortify", "factoextra"), dependencies=T)
  library("corrplot")
  library("ggplot2")
  library("ggfortify")
  library("factoextra")
  CountTableNor<- read.csv("RNAseqCountTableNorTPX.csv")
  
  # Remove low noisy counts
  df_pr<-CountTableNor[3:ncol(CountTableNor)]
  rownames(df_pr)<- paste(1:nrow(CountTableNor),CountTableNor$hgnc_symbol, sep="_")
  df_pr$meanExpr<- apply(df_pr, 1, FUN=mean)
  df_prSel <- df_pr[df_pr$meanExpr>10,]
  df_prSel$meanExpr<- NULL
  
  # Z scores
  df_prSel<- t(scale(t(df_prSel), center=T, scale = T))
  df_prSel<- df_prSel[which(!is.na(rowMeans(df_prSel))),]
  
  pdf("PCA.pdf", width=10, height=10)
    cor<- cor(df_prSel, method = "pearson")
    print(corrplot(cor, method='color', order = "hclust"))
    res.pca<- prcomp(t(df_prSel))
    print(fviz_eig(res.pca))                 
    print(fviz_pca_ind(res.pca,
                 col.ind = "cos2", # Color by the quality of representation
                 gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                 repel = TRUE     # Avoid text overlapping
    ))
    res.pca<- prcomp((df_prSel))
    print(fviz_pca_var(res.pca,
                 col.var = "contrib", # Color by contributions to the PC
                 gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                 repel = TRUE     # Avoid text overlapping
    ))
    # You can adjust k to the number of groups
    res<- hcut(t(df_prSel), k = ncol(df_prSel)-1, stand = TRUE)
    print(fviz_dend(res, rect = TRUE, cex = 0.5))
  dev.off()
} 
