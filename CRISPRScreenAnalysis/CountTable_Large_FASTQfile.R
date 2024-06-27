## CRISPR Screen FASTQ to countTable for large fastq files:
# - We used the screen design of Jastrzebski et al Methods Mol Biol 2016
# - Nextera PCR primers can be found in 'PCR design Nextera-NovaSeq.xlsx'
# - Use this script in Linux
# - Replace 'zcat' for 'cat' in FastqToCountTable.pl (line 59)
# - Put the FastqToCountTable.pl and library.csv files in the same folder as the fastq.gz file
# This script cuts the fastq file in pieces, make count tables, and counts them up
# Author: M.F.M. de Rooij PhD, Amsterdam UMC, Spaargaren Lab, 2024, info: m.f.derooij@amsterdamumc.nl
######################################################################################

# Workdirectory (folder in which the count tables are located, use always slash (/) instead of backslash)
setwd("~/BioLin/Screens/")

# Fastq files
FastqFile <- "test.fastq.gz"

############################################################################################################
# Cut FASTQ file in pieces
system(paste0('zcat ', file, ' | split -l 100000000 - PartOf', file))

# Make Count tables
system('./FastqToCountTable.pl PartOf*')

# Count them up
CountTables<- Sys.glob('CountTable_PartOf*')
ct<- read.csv(CountTables[1])
for (i in 2:length(CountTables)){
  temp<- read.csv(CountTables[2])
  ct[,-1:-3]<- ct[,-1:-3]+temp[,-1:-3]
}
write.csv(ct, paste0("CountTable_", FastqFile, ".csv"), quote = F, row.names = F)
