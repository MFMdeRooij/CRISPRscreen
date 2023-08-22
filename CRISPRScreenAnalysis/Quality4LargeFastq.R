## CRISPR Screen quality analysis for large fastq files:
# - We used the screen design of Jastrzebski et al Methods Mol Biol 2016
# - Use this script in Linux
# - Make count tables from FASTQ files with the Perl script FastqToCountTable.pl and a library.csv file
# - Put Quality4LargeFastq.pl in the same folder as the fastq files
# - The output of this script is a table with the number of read per barcode, and how much reads map to the used library, 
#   and shows the sequence of the most abundant read per barcode
# Author: M.F.M. de Rooij PhD, Amsterdam UMC, Spaargaren Lab, 2023, info: m.f.derooij@amsterdamumc.nl
######################################################################################

# Workdirectory (folder in which the count tables are located, use always slash (/) instead of backslash)
setwd("~/BioLin/Screens/")

# Fastq files
FastqFiles <- Sys.glob('*.fastq.gz')
#FastqFiles <- c("test.fastq.gz.csv", "test2.fastq.gz", "test3.fastq.gz")

# The index barcode ID (in the same order as the order of the fastq files)
print(FastqFiles)
IndexID<- c("ATCACG", "CGATGT", "TTAGGC")
# # c("ATCACG", "CGATGT", "TTAGGC", "TGACCA", "ACAGTG", "GCCAAT")

# Reverse-complement of the template part of the reverse primer 
revCompRevPrimerSeq = "GGAATTGGCTCCGGTGCCCGTCAGT"

############################################################################################################
for (file in FastqFiles) {
  system(paste0('zcat ', file, ' | split -l 100000000 - PartOf', file))
}
system('./Quality4LargeFastq.pl PartOf*')

ReadTables <- paste0('ReadTable_PartOf', FastqFiles)
CountTables <- paste0('CountTable_', FastqFiles, '.csv')
df_index <- data.frame(ReadTable=ReadTables, CountTable=CountTables, IndexID=IndexID, stringsAsFactors = F)

for (index in 1:nrow(df_index)) {
  files <- Sys.glob(paste0(df_index[index,1],"*"))

  for (i in files) {
    df_ReadTable<- read.csv(file=i, sep=",", header=TRUE, stringsAsFactors = FALSE)
    df_ReadTable['X']<-NULL
    
    # Total reads per barcode
    df_Summary<- data.frame(Index=df_index[index,3], Barcode=colnames(df_ReadTable[2:(ncol(df_ReadTable))]), stringsAsFactors = F)
    df_Summary$TotalReads<- colSums(df_ReadTable[,2:ncol(df_ReadTable)])
    
    # Number of CRISPR and shRNA sequences
    df_ReadTable$CRISPR <- 0
    df_ReadTable$shRNA <- 0
    df_ReadTable$primerDimer <- 0
    df_ReadTable[grep("CACCG[ACGT]{18,21}GTTTT", df_ReadTable$Read_without_barcode), "CRISPR"]<- 1
    df_ReadTable[grep("ACCGG[ACGT]{19,22}CTCGAG", df_ReadTable$Read_without_barcode), "shRNA"]<- 1
    df_ReadTable[grep(revCompRevPrimerSeq, df_ReadTable$Read_without_barcode), "primerDimer"]<- 1
    df_Summary$CRISPR<- colSums(df_ReadTable[df_ReadTable$CRISPR==1,2:(ncol(df_ReadTable)-3)])
    df_Summary$shRNA<- colSums(df_ReadTable[df_ReadTable$shRNA==1,2:(ncol(df_ReadTable)-3)])  
    df_Summary$primerDimer<- colSums(df_ReadTable[df_ReadTable$primerDimer==1,2:(ncol(df_ReadTable)-3)])  
    
    if (i==files[1]) {
      # Reads with highest counts
      df_Summary$HighestReadSeq100ntFirstPartFq<-"x"
      for (bc in 1:nrow(df_Summary)) {
        df_high<- df_ReadTable[,c(1,bc+1)]
        colnames(df_high)<- c("Seq", "BC")
        df_high$Seq<-substr(df_high$Seq,1,100)
        df_high2 <- aggregate(df_high$BC, by=list(df_high$Seq), FUN=sum)
        colnames(df_high2)<- c("Seq", "BC")
        df_high2 <- df_high2[order(df_high2$BC, decreasing = T),]
        df_Summary[bc,'HighestReadSeq100ntFirstPartFq']<- paste0(df_high2[1,2], " out of " , sum(df_high2$BC), " - ", df_high2[1,1])
      } 
    }
    assign(paste0("df_",i), df_Summary)
  }  

  dfs <- ls(pattern=paste0(df_index[index,1]))
  df_SumSummary <- get(dfs[1])
  for (i in 2:length(dfs)){
    df_temp <- get(dfs[i])
    df_SumSummary[,3:(ncol(df_SumSummary)-1)] <- df_SumSummary[,3:(ncol(df_SumSummary)-1)] + df_temp[,3:ncol(df_temp)]
  }
  
  df_CountTable <- read.csv(df_index[index,2], sep= ",", header=T, stringsAsFactors = FALSE)
  df_CountTable$X<-NULL

  # Usable reads per barcode
  df_SumSummary$UsableReadsCount<- 0
  df_SumSummary[2:nrow(df_SumSummary),"UsableReadsCount"] <- colSums(df_CountTable[,4:(nrow(df_SumSummary)+2)])
 
  assign(paste0(" df_sumSummaryN",index), df_SumSummary)
}


dfsS <- ls(pattern="df_sumSummaryN")
sumSummary<-get(dfsS[1])
if (length(dfsS)>1) {
    for (i in 2:length(dfsS)) {
      sumSummary = rbind(sumSummary, get(dfsS[i]))
    }
}

sumSummary<- sumSummary[order(sumSummary$Index),]
sumSummary[nrow(sumSummary)+1,]<- c(0, 0, colSums(sumSummary[3:(ncol(sumSummary)-2)]), 0, colSums(sumSummary[ncol(sumSummary)]))
sumSummary$RelUsableReads<- round(sumSummary$UsableReadsCount/sumSummary$TotalReads*100,1)
sumSummary[nrow(sumSummary),1:2]<- "Total"
x<-unlist(lapply(strsplit(sumSummary$HighestReadSeq100ntFirstPartFq[-length(sumSummary$HighestReadSeq100ntFirstPartFq)], " - "), "[", 2))
sumSummary[nrow(sumSummary),'HighestReadSeq100ntFirstPartFq'] <- unique(x)[which.max(tabulate(match(x, unique(x))))]
sumSummary<- cbind(sumSummary[-7],sumSummary[7])

write.csv(sumSummary, "Summary.csv", row.names = F)
