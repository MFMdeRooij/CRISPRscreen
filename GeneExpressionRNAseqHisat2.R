# Follow the lines with hastag signs in R studio
# Author: M.F.M. de Rooij PhD, Amsterdam UMC, Spaargaren Lab, 2019, info: m.f.derooij@amsterdamumc.nl
##################################################################################################################################
#                                                         SETTINGS
# Workdirectory
setwd("~/BioLin/RNAseq/")

# Order samples: 0 = Yes, 1 = No (order MapSamples.txt (together with existing samples) in  AllSamplesSorted.txt)
order<- 0
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
# Download genome in fasta format from NCBI or UCSC, unzip and add in folder ~/HumanGenome
# cd ~/HumanGenome
# hisat2-build hg38.fa hg38

# Download reference files for RNAseq from NCBI or UCSC:
# Download file with gene loci (hg38.95.gtf)
# Transform the gtf file into a known splice site text file using the python script delivered with hisat2 (to hg38_splicesites.txt)
# Add both files to the ~/HumanGenome folder

# R Packages
#install.packages("BiocManager")
#BiocManager::install(c("Rsubread", "DESeq2", "biomaRt"), dependencies=T)
# if biomaRt is not working, try command line: 'sudo apt-get update' & 'sudo apt-get install r-bioc-biomart'

# You can find the gene annotations with biomaRt, and store it in a file for reuse
#library(biomaRt)
#ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
#id<- rownames(CountTable)
##AllAttributes<- listAttributes(ensembl) # to see all atrributes
#GeneList<- getBM(filters= "ensembl_gene_id", attributes=c("ensembl_gene_id", "hgnc_symbol", "description"), values=id, mart= ensembl)
#write.csv(GeneList, "~/HumanGenome/GeneAnnotations.csv", row.names = F)

# Now you are ready for RNAseq analysis

# Search for RNAseq data on https://www.ncbi.nlm.nih.gov/sra/, and fill in the SRR IDs, cell 
# line ID and a P (paired-end) or S (single-end) separated by commas (SRR8615345,NAMALWA,P), and run the 
# script (./MapRNAseqHumanSRA.sh on the command line) 
# The first time, make the bash script excutable (chmod 755 MapRNAseqHumanSRA.sh)
system("bash MapRNAseqHumanSRA.sh")
# With the sorted bam files, you can check the read mapping in IGV-viewer

# Make a count table (TPM with DESEq2 normalisation (median of ratios method) to make the read counts comparable between samples (nTPX))
# Store always the non-normalizered count table to add more samples
library("Rsubread")
library("DESeq2")
for (pair in c("P","S")){
  if (pair=="P") {
     if (length(list.files("pairedEnd/", ".sam"))>0) {
        fc<- featureCounts(files=Sys.glob("pairedEnd/*.sam"),
                          annot.ext="~/HumanGenome/hg38.95.gtf",
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
                        annot.ext="~/HumanGenome/hg38.95.gtf",
                        isGTFAnnotationFile=TRUE,
                        isPairedEnd=FALSE,
                        strandSpecific=0)
    } else {
        break
      }
  } 
  df_annotation<- as.data.frame(fc[2])
  CountTable3<- as.data.frame(fc[1])
  colnames(CountTable3)<- unlist(lapply(strsplit(colnames(CountTable3), ".sam"), "[[", 1))
  colnames(CountTable3)<- unlist(lapply(strsplit(colnames(CountTable3), "_"), "[[", 2))
  CountTable2<- CountTable3/(df_annotation$annotation.Length)
  CountTable<- t(t(CountTable2)/colSums(CountTable2)*1000000)
  CountTable<- as.data.frame(CountTable)
  CountTable$ensembl_gene_id<- rownames(CountTable)

  if (file.exists("RNAseqCountTable.csv")) {
    # Add to existing count table
    CountTableOld<- read.csv("RNAseqCountTable.csv", stringsAsFactors = F)
    CountTable<- merge(CountTableOld, CountTable, by='ensembl_gene_id')
  } else {
      # Add gene symbols
      GeneList<- read.csv('~/HumanGenome/GeneAnnotations.csv', stringsAsFactors = F)
      CountTable<- merge(GeneList, CountTable, by="ensembl_gene_id", all.x=TRUE)
      CountTable$description<- NULL
    }
  if (order==0){
    dfSort<-read.table("AllSamplesSorted.txt", sep=",")
    CountTable<- CountTable[,c("ensembl_gene_id","hgnc_symbol",dfSort$V2)]
  }
  write.csv(CountTable, "RNAseqCountTable.csv", row.names=FALSE)
}  

# Normalize samples (TPM -> between sample normalization -> nTPX)
CountTable<- read.csv("RNAseqCountTable.csv", stringsAsFactors = F)
sf<- estimateSizeFactorsForMatrix(CountTable[,3:ncol(CountTable)])
CountTableNor<- CountTable
CountTableNor[,3:ncol(CountTable)]<- as.data.frame(round(t(t(CountTable[,3:ncol(CountTable)])/sf),1))
write.csv(CountTableNor, "RNAseqCountTableNorTPX.csv", row.names=FALSE)

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

df_prSel<- t(scale(t(df_prSel), center=T, scale = T))
df_prSel<- df_prSel[which(!is.na(rowMeans(df_prSel))),]

pdf("PCA.pdf", width=10, height=10)
  cor<- cor(df_prSel, method = "pearson")
  corrplot(cor, method='color', order = "hclust")
  res.pca<- prcomp(t(df_prSel))
  fviz_eig(res.pca)                 
  fviz_pca_ind(res.pca,
               col.ind = "cos2", # Color by the quality of representation
               gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
               repel = TRUE     # Avoid text overlapping
  )
  res.pca<- prcomp((df_prSel))
  fviz_pca_var(res.pca,
               col.var = "contrib", # Color by contributions to the PC
               gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
               repel = TRUE     # Avoid text overlapping
  )
  res<- hcut(t(df_prSel), k = 15, stand = TRUE)
  fviz_dend(res, rect = TRUE, cex = 0.5)
dev.off()