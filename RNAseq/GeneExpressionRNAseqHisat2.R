# Follow the lines with hastag signs
# You can run a couple of cell lines together, but run single-end and paired-end data separately
# After the read mapping, you should rename the sam files (include the cell names)  
# (-> Soon I will automate this, so that the entire script can be run overnight)
# Author: M.F.M. de Rooij PhD, Amsterdam UMC, Spaargaren Lab, 2019, info: m.f.derooij@amsterdamumc.nl

# Workdirectory
Workdirectory<- "~/Data/RNAseqCellLines/"

# Paired-end reads 0 = Yes, 1 = No
pairedEnd<-0

q<- function(...) {
  sapply(match.call()[-1], deparse)
}
# Order count table
Sort<- q(ensembl_gene_id,	hgnc_symbol,	NALM6, REH, SEM, X697, MEC1,	GRANTA519,	JEKO1,	MINO,	REC1, Z138, DAUDI, NAMALWA,	RAJI, DOHH2, SUDHL4, SUDHL5, 
                SUDHL6,	OCILY3,	OCILY10,	TMD8,	U2932, EJM, INA6, L363,	LP1, MM1S, NCIH929,	OPM2,	RPMI8226, U266)
##########################################
setwd(Workdirectory)
# NCBI:
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
# Download File with gene loci (hg38.95.gtf)
# Transform the gtf file into a known splice site text file using the python script delivered with hisat2 (to hg38_splicesites.txt)
# add these files to the ~/HumanGenome folder

# R Packages
#install.packages("BiocManager")
#BiocManager::install(c("Rsubread", "biomaRt", "DESeq2"), dependencies=T)
# if biomaRt is not working, try command line: 'sudo apt-get update' & 'sudo apt-get install r-bioc-biomart'
library(Rsubread)
library(biomaRt)
library(DESeq2)

# Search for RNAseq data on https://www.ncbi.nlm.nih.gov/sra/, and fill in the SRR codes in MapSamples.txt
if (pairedEnd==0){
  system("./MapRNAseqHuman_SRA_PE.sh")
}
if (pairedEnd==1){
  system("./MapRNAseqHuman_SRA_U.sh")
}

# With the sorted bam files, you can check the read mapping in IGV-viewer

# Add an underscore and a cell name (which corresponds to a cell name in the Sort vector) to the SAM filenames ("SRR1031032.sam" -> "SRR1031032_MINO.sam")

# Make Count table (TPM with DESEq2 normalisation (median of ratios method) to make it comparable between samples)
if (pairedEnd==0){
  fc<- featureCounts(files=Sys.glob("*.sam"), 
                      annot.ext="~/HumanGenome/hg38.95.gtf", 
                      isGTFAnnotationFile=TRUE,
                      isPairedEnd=TRUE,
                      strandSpecific=0)
}
if (pairedEnd==1){
  fc<- featureCounts(files=Sys.glob("*.sam"), 
                      annot.ext="~/HumanGenome/hg38.95.gtf", 
                      isGTFAnnotationFile=TRUE,
                      isPairedEnd=FALSE,
                      strandSpecific=0)
}
df_annotation<- as.data.frame(fc[2])
CountTable3<- as.data.frame(fc[1])
colnames(CountTable3)<- unlist(lapply(strsplit(colnames(CountTable3), "_"), "[[", 2))
CountTable2<- CountTable3/(df_annotation$annotation.Length)
CountTable<- t(t(CountTable2)/colSums(CountTable2)*1000000)
CountTable<- as.data.frame(CountTable)

# Add gene symbols
CountTable$ensembl_gene_id<- rownames(CountTable)
 # You can flnd gene annotations with biomaRt, and store it in a file for reuse
 #ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
 #id<- rownames(CountTable)
 ##AllAttributes<- listAttributes(ensembl)
 #GeneList<- getBM(filters= "ensembl_gene_id", attributes=c("ensembl_gene_id", "hgnc_symbol", "description"), values=id, mart= ensembl)
 #write.csv(GeneList, "GeneList.csv", row.names = F)
GeneList<- read.csv('GeneList.csv', stringsAsFactors = F)
CountTable<- merge(GeneList, CountTable, by="ensembl_gene_id", all.x=TRUE)
CountTable$description<- NULL
write.csv(CountTable, "RNAseqCountTable2.csv", row.names=FALSE)

# Add new cell lines to an existing table (rename the old RNAseqCountTable.csv to RNAseqCountTable1.csv)
CountTable1<- read.csv("RNAseqCountTable1.csv", stringsAsFactors = F)
CountTable2<- read.csv("RNAseqCountTable2.csv", stringsAsFactors = F)
CountTable<- merge(CountTable1, CountTable2[,-2], by='ensembl_gene_id')
CountTable<- CountTable[,Sort]
write.csv(CountTable, "RNAseqCountTable.csv", row.names=FALSE)
#CountTable<- read.csv("RNAseqCountTable.csv", stringsAsFactors = F)
sf<- estimateSizeFactorsForMatrix(CountTable[,3:ncol(CountTable)])
CountTableNor<- CountTable
CountTableNor[,3:ncol(CountTable)]<- as.data.frame(round(t(t(CountTable[,3:ncol(CountTable)])/sf),1))
CountTableNor<- CountTableNor[,Sort]
write.csv(CountTableNor, "RNAseqCountTableNorTPM.csv", row.names=FALSE)


#PCA analysis
#install.packages(c("corrplot", "ggplot2", "ggfortify", "factoextra"), dependencies=T)
library("corrplot")
library("ggplot2")
library("ggfortify")
library("factoextra")
CountTableNor<- read.csv("RNAseqCountTableNorTPM.csv")
# If there is RNAseq data mapped to another genome reference, you can remove genes which are not available in all samples
#CountTableNor$minExp<- apply(CountTableNor[3:ncol(CountTableNor)], 1, FUN=min)
#CountTableNor<- CountTableNor[CountTableNor$minExp>1,]
#CountTableNor$minExp<- NULL
#sf<- estimateSizeFactorsForMatrix(CountTableNor[,3:ncol(CountTableNor)])
#CountTableNor[,3:ncol(CountTableNor)]<- as.data.frame(round(t(t(CountTableNor[,3:ncol(CountTableNor)])/sf),1))
CountTableNor<- CountTableNor[,Sort]
df_pr<- CountTableNor[,3:ncol(CountTableNor)]
rownames(df_pr)<- paste(1:nrow(CountTableNor),CountTableNor$hgnc_symbol, sep="_")
df_prScale<- t(scale(t(df_pr), center=T, scale = T))
df_pr_NoNA<- which(!is.na(rowMeans(df_prScale)))
df_prSel<- df_prScale[df_pr_NoNA,]
samples<- data.frame(Sample = as.factor(colnames(CountTableNor)[3:ncol(CountTableNor)]))
rownames(samples)<- colnames(CountTableNor[3:ncol(CountTableNor)])

pdf("PCA.pdf", width=10, height=10)
  autoplot(prcomp(t(df_prSel)), data=samples, colour="Sample")
  #pairs(df_pr, cex=0.1, log='xy')
  cor<- cor(df_prSel, method = "pearson")
  corrplot(cor, method='color', order = "hclust")
  res.pca<- prcomp(t(df_prSel), scale.=T, center=T)
  fviz_eig(res.pca)                 
  fviz_pca_ind(res.pca,
               col.ind = "cos2", # Color by the quality of representation
               gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
               repel = TRUE     # Avoid text overlapping
  )
  res.pca<- prcomp((df_prSel), scale.=T, center=T)
  fviz_pca_var(res.pca,
               col.var = "contrib", # Color by contributions to the PC
               gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
               repel = TRUE     # Avoid text overlapping
  )
  my_data<- t(df_prSel)
  #fviz_nbclust(my_data, kmeans, method = "gap_stat")
  res<- hcut(t(df_prSel), k = 15, stand = TRUE)
  fviz_dend(res, rect = TRUE, cex = 0.5)
dev.off()
