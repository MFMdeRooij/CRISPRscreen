## RNAseq gene familiy pie plot:
# - Make count tables with the GeneExpressionRNAseqHisat2.R script
# - Use the RNAseqCountTableNorTPX.csv (or RNAseqCountTableTPM.csv) count table
#   (It is important that the counts are adjusted for gene length!!)
# - Adjust the settings and run the script
# Author: M.F.M. de Rooij PhD, Amsterdam UMC, Spaargaren Lab, 2022, info: m.f.derooij@amsterdamumc.nl
##########################################################################################################
#                                     SETTINGS

# Workdirectory (folder in which the count tables are located, use always slash (/) instead of backslash)
setwd("H:/BioWin/RNAseq/")

# Which count table?
Filename<- read.csv(file="RNAseqCountTableNorTPX.csv", sep=",", header=TRUE, stringsAsFactors = FALSE)

# Title gene family
Title<- "ChemokineReceptors"

# Select genes with regular expression
sel<- CountTableNor[grep("^CXCR[[:digit:]]$", toupper(CountTableNor$hgnc_symbol)),]
sel<- rbind(sel,CountTableNor[grep("^CCR[[:digit:]]$", toupper(CountTableNor$hgnc_symbol)),])

# # Select genes from a file (Make a csv file with gene symbols in the first column with Gene as colname)
# GeneList<- read.csv("GeneList.csv", sep=",", header=TRUE, stringsAsFactors = FALSE)
# Genes<- toupper(GeneList$Gene)
# sel<- CountTableNor[toupper(CountTableNor$hgnc_symbol) %in% Genes,]

##########################################################################################################
# Generate pie plots
rownames(sel)<-sel$hgnc_symbol
sel$hgnc_symbol<-NULL
sel$ensembl_gene_id<-NULL
pdf(paste0("PiePlots_", Title,".pdf"),10,20)
  par(mfrow=c(6,3))
  for (cell in colnames(sel)){
    total<-as.numeric(colSums(sel[cell]))
    pie(sel[[cell]], label=paste0(rownames(sel),  " (", as.numeric(t(sel[cell])),")"), main=cell, sub=paste0("Total counts: ",round(total,0)," nTPX"))
  }
dev.off()
