## CRISPR Screen quality analysis:
# - We used the screen design of Jastrzebski et al Methods Mol Biol 2016
# - Make count tables from FASTQ files with the Perl script FastqToCountTable.pl and a library.csv file
# - Make read tables from FASTQ files with the Perl script Quality.pl
# - The output of this script is a table with the number of read per barcode, and how much reads map to the used library, 
#   and shows the sequence of the most abundant read per barcode
# Author: M.F.M. de Rooij PhD, Amsterdam UMC, Spaargaren Lab, 2019, info: m.f.derooij@amsterdamumc.nl
######################################################################################
# Workdirectory (folder in which the count tables are located, use always slash (/) instead of backslash)
Workdirectory<- "~/Data/Screens"

# Read table files
ReadTable<- c("ReadTable_test.fastq.gz.csv", "ReadTable_test2.fastq.gz.csv")   

# Count table files (in the same order as the order of the read table files)
CountTable<- c("CountTable_test.fastq.gz.csv", "CountTable_test2.fastq.gz.csv")              

# The index barcode ID (in the same order as the order of the read/count table files)
IndexID<- c("TGACCA", "ATCACG")    
# c("ATCACG", "CGATGT", "TTAGGC", "TGACCA", "ACAGTG", "GCCAAT")

# Reverse-complement of the template part of the reverse primer 
revCompRevPrimerSeq = "GGAATTGGCTCCGGTGCCCGTCAGT"
############################################################################################################
setwd(Workdirectory)
df_index<- data.frame(ReadTable=ReadTable, CountTable=CountTable, IndexID=IndexID, stringsAsFactors = F)

df_SumSummary<- data.frame(Index=numeric(0), Barcode=character(0), TotalReads=numeric(0), UsableReadsCount= numeric(0), 
                            RelUsableReads=numeric(0), shRNA=numeric(0), CRISPR=numeric(0), HighestReadCount=numeric(0), HighestReadSeq=character(0))

for (index in 1:nrow(df_index)) {
  df_ReadTable <- read.csv(df_index[index,1], sep= ",", header=T, stringsAsFactors = FALSE)
  df_ReadTable$X<-NULL
  df_CountTable <- read.csv(df_index[index,2], sep= ",", header=T, stringsAsFactors = FALSE)
  df_CountTable$X<-NULL
  
  # Total reads per barcode
  df_Summary<- data.frame(Index=df_index[index,3], Barcode=colnames(df_ReadTable[2:(ncol(df_ReadTable))]), stringsAsFactors = F)
  df_Summary$TotalReads<- colSums(df_ReadTable[,2:ncol(df_ReadTable)])

  # Usable reads per barcode
  df_Summary$UsableReadsCount<- 0
  df_Summary[2:nrow(df_Summary),4] <- colSums(df_CountTable[,4:(nrow(df_Summary)+2)])
  df_Summary$RelUsableReads<- round(df_Summary$UsableReadsCount/df_Summary$TotalReads*100,1)
  
  # Number of shRNA and CRISPR sequences
  df_ReadTable$CRISPR <- 0
  df_ReadTable$shRNA <- 0
  df_ReadTable$primerDimer <- 0
  df_ReadTable[grep("CACCG[ACGT]{18,21}GTT", df_ReadTable$Read_without_barcode), "CRISPR"]<- 1
  df_ReadTable[grep("ACCGG[ACGT]{19,22}CT", df_ReadTable$Read_without_barcode), "shRNA"]<- 1
  df_ReadTable[grep(revCompRevPrimerSeq, df_ReadTable$Read_without_barcode), "primerDimer"]<- 1
  df_Summary$CRISPR<- colSums(df_ReadTable[df_ReadTable$CRISPR==1,2:(ncol(df_ReadTable)-3)])
  df_Summary$shRNA<- colSums(df_ReadTable[df_ReadTable$shRNA==1,2:(ncol(df_ReadTable)-3)])  
  df_Summary$primerDimer<- colSums(df_ReadTable[df_ReadTable$primerDimer==1,2:(ncol(df_ReadTable)-3)])  
  # Reads with highest counts
  df_Summary$HighestReadCount<-0
  df_Summary$HighestReadSeq<-0
  
  for (bc in 1:nrow(df_Summary)) {
    df_high<- df_ReadTable[,c(1,bc+1)]  
    colnames(df_high)<- c("Seq", "BC")
    df_high <- df_high[order(df_high$BC, decreasing = T),]
    df_Summary[bc,9]<- df_high[1,2]
    df_Summary[bc,10]<- df_high[1,1]
  }
 df_SumSummary<-rbind(df_SumSummary, df_Summary)
}

df_SumSummary<- df_SumSummary[order(df_SumSummary$Index),]
write.csv(df_SumSummary, "Summary.csv", row.names = F)