# Workdirectory
Workdirectory<- "~/Data/RNAseqCellLines/"

# Remove SRA and intermediate FASTQ files 0 = Yes, 1 = No
removeSRA<- 0

# Cells which are paired-end sequenced)
pairedEnd<- c("SRR1031032_MINO.sra", "SRR8615896_U266.sra", "SRR8616164_X697.sra", "SRR7611497_BJAB.sra","SRR7659713_UPN1.sra")

# Order count table
q <- function(...) {
  sapply(match.call()[-1], deparse)
}
Sort<- q(ensembl_gene_id,	hgnc_symbol,	NALM6, REH, SEM, X697, MEC1,	GRANTA519,	JEKO1,	MINO,	REC1, Z138, DAUDI, NAMALWA,	RAJI, DOHH2, SUDHL4, SUDHL5, 
                SUDHL6,	OCILY3,	OCILY10,	TMD8,	U2932, EJM, INA6, L363,	LP1, MM1S, NCIH929,	OPM2,	RPMI8226, U266)
##########################################
setwd(Workdirectory)
##NCBI:
#Download, install, and configure the sra-toolkit from NCBI website:
#https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=software

##Linux Tools:
#sudo apt-get install seqtk
#sudo apt-get install samtools
#sudo apt-get install bowtie2
#sudo apt install hisat2

#biomaRt dependencies:
#sudo apt-get install libcurl4-openssl-dev libssl-dev

##Build genome index
#setwd("~/HumanGenome")
#system("bowtie2-build hg38.fa hg38") 

##Download reference files for RNAseq:
# File with known splicesites (hg38_splicesites.txt)
# File with gene loci (hg38.95.gtf)

##R Packages
#install.packages("BiocManager")
#BiocManager::install(c("Rsubread", "biomaRt", "DESeq2"), dependencies=T)
#if biomaRt is not working, try command line: 'sudo apt-get update' & 'sudo apt-get install r-bioc-biomart'
library(Rsubread)
library(biomaRt)
library(DESeq2)

##Download SRR files:
##https://www.ncbi.nlm.nih.gov/sra/
##On command line (SRA tools):
#prefetch SRR6047361
##for multiple SRR files:
#for i in {2..4}; do prefetch SRR604736$i; done

##To tranfer files to another location
#for i in *.sra; do sudo cp $i ~/SharedFolder/$1; done


# Convert SRA files to FastQ files
files <- (Sys.glob("*.sra"))
for (cell in files) {
  system(paste0("fastq-dump ",cell))
  #"fastq-dump --split-files FILE": for immediately splitting paired-end fastq files
  if (removeSRA == 0) {
    file.remove(cell)
  }
}

# Count read length
system("perl GeneExpressionRNAseqReadSize.pl")
df_readLength<- read.csv('ReadSize.csv', stringsAsFactors = F)

# Split paired-end FASTQ files
pairedEndfastq<-paste0(pairedEnd, '.fastq')
for (cell in pairedEndfastq) {
    numberCutOff<-df_readLength[df_readLength$cell==cell,"readLength"]/2 
    # -e remove nt for the end, -b remove nt from begin
    system(paste0("seqtk trimfq -e ", numberCutOff, " ", cell, " > 1.", cell))
    system(paste0("seqtk trimfq -b ", numberCutOff, " ", cell, " > 2.", cell))
    if (removeSRA == 0) {
      file.remove(cell)
    }
}    

# Trim FASTQ files (to remove low quality sequences)
files <- (Sys.glob("*.sra.fastq"))
for (cell in files) {    
    system(paste0("seqtk trimfq ", cell, " > ", cell, ".fq"))
    if (removeSRA == 0) {
      file.remove(cell)
    }   
} 

# Run hisat2 (mapping reads to genome)

# first paired-ends
for (PE in pairedEnd){
  dir.create(paste0(PE,".fastq.fq1"))
  file.rename(paste0("1.",PE,".fastq.fq"), paste0(PE,".fastq.fq1/1.",PE,".fastq.fq"))
  file.rename(paste0("2.",PE,".fastq.fq"), paste0(PE,".fastq.fq1/2.",PE, ".fastq.fq"))
  setwd(paste0(PE,".fastq.fq1"))
  system(paste0("hisat2 -x ~/HumanGenome/hg38 --known-splicesite-infile ~/HumanGenome/hg38_splicesites.txt -1 1.", PE, ".fastq.fq -2 2.", PE, ".fastq.fq -S ", PE, ".fastq.fq.sam --summary-file ", PE, ".fastq.fq.Summary.txt")) 
  file.rename(paste0(PE,".fastq.fq.sam"), paste0("../",PE,".fastq.fq.sam"))
  setwd(Workdirectory)
}

# Then the single-end
files <- (Sys.glob("*.fq"))
for (cell in files) {
  dir.create(paste0(cell,1))
  file.rename(cell,paste0(cell,1,"/",cell))
  setwd(paste0(cell,1))
  system(paste0("hisat2 -x ~/HumanGenome/hg38 --known-splicesite-infile ~/HumanGenome/hg38_splicesites.txt -U ", cell, " -S ", cell, ".sam --summary-file ", cell, ".Summary.txt")) 
  file.rename(paste0(cell,".sam"), paste0("../",cell,".sam"))
  setwd(Workdirectory)
}

# Make Count table (TPM with DESEq2 normalisation to make it comparable between samples)
fc <- featureCounts(files=Sys.glob("*.sam"), 
                    annot.ext="~/HumanGenome/hg38.95.gtf", 
                    isGTFAnnotationFile=TRUE,
                    isPairedEnd=FALSE,
                    strandSpecific=0)

df_annotation<-as.data.frame(fc[2])

CountTable3<- as.data.frame(fc[1])
colnames(CountTable3)<- unlist(lapply(strsplit(colnames(CountTable3), ".sra"), "[[", 1))
colnames(CountTable3)<- unlist(lapply(strsplit(colnames(CountTable3), "_"), "[[", 2))

CountTable2<- CountTable3/(df_annotation$annotation.Length)
CountTable<- t(t(CountTable2)/colSums(CountTable2)*1000000)
CountTable<- as.data.frame(CountTable)

# Add gene symbols
CountTable$ensembl_gene_id<- row.names(CountTable)
# ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
# genes <- rownames(CountTable)
# #AllAttributes<-listAttributes(ensembl)
# GeneList <- getBM(filters= "ensembl_gene_id", attributes=c("ensembl_gene_id", "hgnc_symbol", "description"), values=genes,mart= ensembl)
# 
# write.csv(GeneList, "GeneList.csv", row.names = F)
GeneList<-read.csv('GeneList.csv', stringsAsFactors = F)

CountTable<- merge(GeneList, CountTable, by="ensembl_gene_id", all.x=TRUE)
CountTable$description<-NULL


write.csv(CountTable, "RNAseqCountTable2.csv", row.names=FALSE)

# Add new cell lines to an existing table
CountTable1<-read.csv("RNAseqCountTable1.csv", stringsAsFactors = F)
CountTable2<-read.csv("RNAseqCountTable2.csv", stringsAsFactors = F)
CountTable<- merge(CountTable1, CountTable2[,-2], by='ensembl_gene_id')
write.csv(CountTable, "RNAseqCountTable.csv", row.names=FALSE)
#CountTable<-read.csv("RNAseqCountTable.csv", stringsAsFactors = F)
sf <- estimateSizeFactorsForMatrix(CountTable[,3:ncol(CountTable)])
CountTableNor<-CountTable
CountTableNor[,3:ncol(CountTable)]<- as.data.frame(round(t(t(CountTable[,3:ncol(CountTable)])/sf),1))

CountTableNor<-CountTableNor[,Sort]

write.csv(CountTableNor, "RNAseqCountTableNorTPM.csv", row.names=FALSE)

#PCA analysis
#install.packages(c("corrplot", "ggplot2", "ggfortify"), dependencies=T)
library("corrplot")
library("ggplot2")
library("ggfortify")
CountTableNor<-read.csv("RNAseqCountTableNorTPM.csv")

CountTableNor$minExp <- apply(CountTableNor[3:ncol(CountTableNor)], 1, FUN=min)
CountTableNor<- CountTableNor[CountTableNor$minExp>1,]
CountTableNor$minExp <- NULL
sf <- estimateSizeFactorsForMatrix(CountTableNor[,3:ncol(CountTableNor)])
CountTableNor[,3:ncol(CountTableNor)]<- as.data.frame(round(t(t(CountTableNor[,3:ncol(CountTableNor)])/sf),1))
CountTableNor<-CountTableNor[,Sort]

df_pr<- CountTableNor[,3:ncol(CountTableNor)]
rownames(df_pr)<-paste(1:nrow(CountTableNor),CountTableNor$hgnc_symbol, sep="_")

df_prScale<-t(scale(t(df_pr), center=T, scale = T))
df_pr_NoNA<- which(!is.na(rowMeans(df_prScale)))
df_prSel<-df_prScale[df_pr_NoNA,]

samples<-data.frame(Sample = as.factor(colnames(CountTableNor)[3:ncol(CountTableNor)]))
rownames(samples)<-colnames(CountTableNor[3:ncol(CountTableNor)])

pdf("PCA.pdf", width=10, height=10)
  autoplot(prcomp(t(df_prSel)), data=samples, colour="Sample")
  #pairs(df_pr, cex=0.1, log='xy')
  cor<-cor(df_prSel, method = "pearson")
  corrplot(cor, method='color', order = "hclust")
  
  #install.packages("factoextra", dependecies=T)
  library(factoextra)
  res.pca <- prcomp(t(df_prSel), scale.=T, center=T)
  fviz_eig(res.pca)                 
  fviz_pca_ind(res.pca,
               col.ind = "cos2", # Color by the quality of representation
               gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
               repel = TRUE     # Avoid text overlapping
  )
  res.pca <- prcomp((df_prSel), scale.=T, center=T)
  fviz_pca_var(res.pca,
               col.var = "contrib", # Color by contributions to the PC
               gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
               repel = TRUE     # Avoid text overlapping
  )
  my_data <- t(df_prSel)
  #fviz_nbclust(my_data, kmeans, method = "gap_stat")
  res <- hcut(t(df_prSel), k = 15, stand = TRUE)
  # Visualize
  fviz_dend(res, rect = TRUE, cex = 0.5)

dev.off()
