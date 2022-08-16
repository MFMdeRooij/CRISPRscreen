# With this script you can look at splicing of human CD44, which has constant and variable exon usage, or or CXCL12 which has different isoforms
# Follow the lines with hastag signs in R studio
# Author: M.F.M. de Rooij PhD, Amsterdam UMC, Spaargaren Lab, 2022, info: m.f.derooij@amsterdamumc.nl
##################################################################################################################################
#                                                         SETTINGS

# Workdirectory
wd <- "~/BioLin/RNAseq/"
#wd <- "H:/BioWin/RNAseq/"
GeneStructure<- "hCD44_GeneStructure.csv"
#GeneStructure<- "hCXCL12_GeneStructure.csv"
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

# R Packages
#install.packages("BiocManager")
#BiocManager::install(c("scales","Rsamtools","stringr"))

# For other genes, you can find in Ensembl the exon loci (or in the GTF file), which you can put in hGeneX_GeneStructure.csv, and 
# also adjust the most upstream and downstream loci in hGeneX_Splicing.sh

# Now you are ready to look at the splice variants

# Search for RNAseq data on https://www.ncbi.nlm.nih.gov/sra/, and fill in the SRR IDs, cell 
# line ID and a P (paired-end) or S (single-end) separated by commas (SRR8615282,HCT116,P), and run the 
# script (./MapRNAseqHumanSRA.sh on the command line) 
# The first time, make the bash script excutable (chmod 755 MapRNAseqHumanSRA.sh)
# You can also run the bash script here by uncommenting and run the next line:
#system("bash hCD44_Splicing.sh")

# This script will make a plot with read depth of all reads (blue) and reads with gaps (black), and shows the introns of human CD44. As a control gene the read depth of actin will be shown.
library("scales")
library("Rsamtools")
library("stringr")
# A function to read bam file 
readBAM <- function(bamFile){ 
  bam <- scanBam(bamFile) 
  # A function for collapsing the list of lists into a single list 
  # as per the Rsamtools vignette 
  .unlist <- function (x){ 
    x1 <- x[[1L]] 
    if (is.factor(x1)){ 
      structure(unlist(x), class = "factor", levels = levels(x1)) 
    } else { 
      do.call(c, x) 
    } 
  } 
  bam_field <- names(bam[[1]]) 
  list <- lapply(bam_field, function(y) .unlist(lapply(bam, "[[", y))) 
  bam_df <- do.call("DataFrame", list) 
  names(bam_df) <- bam_field 
  #return a list that can be called as a data frame 
  return(bam_df) 
}
setwd(wd)
dfID<- read.csv("MapSamples.txt", header=F, stringsAsFactors = F)
GeneStructure_GOI<- read.csv(GeneStructure, stringsAsFactors = F)
direction <- "pos"
if (GeneStructure_GOI$Start[1]>GeneStructure_GOI$Start[nrow(GeneStructure_GOI)]){
  direction<- "neg"
}
GeneSymbol<- strsplit(strsplit(GeneStructure, "_")[[1]][1],"h")[[1]][2]

for (i in 1:nrow(dfID)){
  ID<-paste0(dfID[i,1], "_", dfID[i,2])
  setwd(paste0(wd,ID))
  GOI<-data.frame(V1=1,V2=1,V3=1)
  splice_GOI<-data.frame(V1=1,V2=1,V3=1)
  actb<-data.frame(V1=1,V2=1,V3=1)
  splice_actb<-data.frame(V1=1,V2=1,V3=1)
  try(
    GOI<-read.table(paste0(ID, ".", tolower(GeneSymbol), ".bam.txt"), sep="\t", fill=T)
  )  
  try(
    splice_GOI<-read.table(paste0("splice_",ID,".", tolower(GeneSymbol), ".bam.txt"), sep="\t", fill=T)
  ) 
  try(
    actb<-read.table(paste0(ID,".actb.bam.txt"), sep="\t", fill=T)
  )
  try(
    splice_actb<-read.table(paste0("splice_",ID,".actb.bam.txt"), sep="\t", fill=T)
  )
  
  # Count introns
  
  # Load the bam file
  bam<-data.frame(qname=1,flag=1,rname=1,strand=1,pos=1,qwidth=1,mapq=1,cigar=1,mrnm=1,mpos=1,isize=1,seq=1,qual=1)
  
  try(
    bam2 <- as.data.frame(readBAM(paste0("splice_",ID,".", tolower(GeneSymbol), ".bam")))
  )
  if (nrow(bam2)>0){
    bam<-bam2
  }
  # Use only 1-times mapped reads (some reads are inconclusive)
  bam<- bam[bam$mapq==60,]
  bam$numberIntr<- str_count(bam$cigar, pattern = "N")
  bam$up<- unlist(lapply(strsplit(as.character(bam$cigar), "N"), "[", 1)) 
  bam$startintron<-0
  for (i in 1:nrow(bam)){
    bam$startintron[i]<-bam$pos[i]+sum(as.numeric(sub("M", "", str_extract_all(bam$up[i],"(\\d+)(M)",simplify = T))))
  }
  bam$endintron<- bam$startintron+as.numeric(str_extract_all(bam$up, "\\d+$"))

  # Include multiple introns per read
  if (max(bam$numberIntr)>1){
    bamMultInt<- bam[bam$numberIntr>1,]
    for (n in 2:max(bam$numberIntr)){
      bamMultInt<- bamMultInt[bamMultInt$numberIntr>=n,]
      for (i in 1:nrow(bamMultInt)){
        bamMultInt$cigar[i]<-strsplit(bamMultInt$cigar[i], paste0(bamMultInt$up[i], "N"))[[1]][2]
      }
      bamMultInt$up<- unlist(lapply(strsplit(as.character(bamMultInt$cigar), "N"), "[", 1))
      for (i in 1:nrow(bamMultInt)){
        bamMultInt$startintron[i]<-bamMultInt$endintron[i]+sum(as.numeric(sub("M", "", str_extract_all(bamMultInt$up[i],"(\\d+)(M)",simplify = T))))
      }
      bamMultInt$endintron<- bamMultInt$startintron+as.numeric(str_extract_all(bamMultInt$up[i], "\\d+$"))
      bam<-rbind(bam,bamMultInt)
    }
  }
  dfsum<- data.frame(Intron=(paste0(bam$startintron,",",bam$endintron)), count=1)
  dfIntron<-aggregate(dfsum$count, by=list(dfsum$Intron), FUN=sum)
  dfIntron<-dfIntron[order(dfIntron$x, decreasing = T),]
  
  dfIntron$start<-as.numeric(
    unlist(lapply(strsplit(as.character(dfIntron$Group.1), ","), "[", 1))
  )    
  dfIntron$end<-as.numeric(
    unlist(lapply(strsplit(as.character(dfIntron$Group.1), ","), "[", 2))
  ) 
  dfIntron$size<- dfIntron$end-dfIntron$start
  dfIntron<- dfIntron[order(dfIntron$size),]
  dfIntron<- dfIntron[order(dfIntron$x, decreasing = T),]
  if (direction=="pos"){
    dfIntron<- dfIntron[order(dfIntron$start),]
  } else{
    dfIntron<- dfIntron[order(dfIntron$end, decreasing = T),]
  }
  # Use only introns which are abundant > 1% of the most frequent intron
  dfIntron<-dfIntron[dfIntron$x>max(dfIntron$x)/100,]
  
  # Read densities
  df<-GOI
  dfs<-splice_GOI
  xlim<-c(min(GeneStructure_GOI[c("Start","End")]),max(GeneStructure_GOI[c("Start","End")]))
  if (direction=="neg"){
    xlim<-rev(xlim)
  }
  chr<-paste0("Chr",GeneStructure_GOI$Chr[1])
  main=GeneSymbol
  df2<- data.frame(V2=min(xlim):max(xlim), V4=0)
  df<-merge(df,df2,by="V2", all.y=T)
  df$V3[is.na(df$V3)]<-0
  df$V3<-df$V4+df$V3
  dfs2<- data.frame(V2=min(xlim):max(xlim), V4=0)
  dfs<-merge(dfs,dfs2,by="V2", all.y=T)
  dfs$V3[is.na(dfs$V3)]<-0
  dfs$V3<-dfs$V4+dfs$V3
  ylim<-max(df$V3)
  if (ylim<10){
    ylim=10
  }
  pdf(paste0("../", ID,".pdf"),10,5)
    par(mar=c(4,5,2,1))
    
    # Plot read density
    par(fig=c(0.1,0.83,0.5,1))
    plot(df$V2, df$V3, type="l" , cex=0.1, xlim=xlim, ylim=c(0, ylim), xlab="", xaxt="n", ylab="depth", col="white" , main=main)
    polygon(c(min(df$V2)-1,df$V2, max(df$V2)+1), c(0,df$V3,0), col="lightblue", border=F)
    polygon(c(min(dfs$V2)-1,dfs$V2, max(dfs$V2)+1), c(0,dfs$V3,0), col="black", border=F)
    segments(xlim[1], 0, xlim[2],0, col="gray", lwd=0.1)

    # Draw gene structures 
    par(fig=c(0.1,0.83,0.1,0.73),new=TRUE)
    ylimMax<-nrow(dfIntron)*-10*(5/4)
    plot(1, type='n', xlab=chr, ylab="Introns", xlim=xlim, yaxt="n", ylim=c(ylimMax,0), yaxs="i")
    abline(h=(ylimMax/5))
    for (i in 1:nrow(GeneStructure_GOI)){
      if (!grepl("st", GeneStructure_GOI$Exon[i])){
        rect(GeneStructure_GOI$Start[i],	ylimMax-1000, GeneStructure_GOI$End[i], ylimMax/5, col=alpha(GeneStructure_GOI$Color[i],0.1), border=F)
      }
    }
    for (i in 1:nrow(GeneStructure_GOI)){
      if (!grepl("st", GeneStructure_GOI$Exon[i])){
        text(((GeneStructure_GOI$Start[i]+GeneStructure_GOI$End[i])/2), ylimMax/5*0.5, GeneStructure_GOI$Exon[i], cex=0.3, col=GeneStructure_GOI$Color[i],
              adj=if (GeneSymbol=="CD44"){if(GeneStructure_GOI$Exon[i]=="v4"){0.7}else{if(GeneStructure_GOI$Exon[i]=="v5"){0.3}else{0.5}}})
      }
    }
    for (i in 1:nrow(GeneStructure_GOI)){
      segments(GeneStructure_GOI$Start[i],	ylimMax/5*0.7, GeneStructure_GOI$End[i], ylimMax/5*0.7,lwd=GeneStructure_GOI$Size[i], col=GeneStructure_GOI$Color[i], lend=1)
    }
  
    # Draw introns
    for (i in 1:nrow(dfIntron)){
      rect(dfIntron$start[i],	(ylimMax/5)+i*((ylimMax-(ylimMax/5))/(nrow(dfIntron)+1))-ylimMax/100, dfIntron$end[i], (ylimMax/5)+i*((ylimMax-(ylimMax/5))/(nrow(dfIntron)+1))+ylimMax/100, lwd=0.1, col="white", border="black")
      rect(dfIntron$start[i],	(ylimMax/5)+i*((ylimMax-(ylimMax/5))/(nrow(dfIntron)+1))-ylimMax/100, dfIntron$end[i], (ylimMax/5)+i*((ylimMax-(ylimMax/5))/(nrow(dfIntron)+1))+ylimMax/100, lwd=0.1, col=alpha("black",dfIntron$x[i]/max(dfIntron$x)), border="black")
      text(if (direction=="neg"){dfIntron$start[i]} else{dfIntron$end[i]}, (ylimMax/5)+i*((ylimMax-(ylimMax/5))/(nrow(dfIntron)+1)), pos=4,offset=0.1, paste0("(",dfIntron$x[i],")"), srt = 0, xpd = TRUE, cex=0.3)
    }

    # Special features for CD44
    if (GeneSymbol=="CD44"){
      dfFeatures<-data.frame(
        feature=c("ECD","TM","ICD"),
        start=c("1start","17","18stop"),
        end=c("16","17","18stop")
      )
      for (i in 1:nrow(dfFeatures)) {
        text(((GeneStructure_GOI$Start[GeneStructure_GOI$Exon==dfFeatures$start[i]]+GeneStructure_GOI$End[GeneStructure_GOI$Exon==dfFeatures$end[i]])/2), ylimMax/5*0.15, dfFeatures$feature[i], font=2, cex=0.3, adj=0.5)
        segments(GeneStructure_GOI$Start[GeneStructure_GOI$Exon==dfFeatures$start[i]],	ylimMax/5*0.3, GeneStructure_GOI$End[GeneStructure_GOI$Exon==dfFeatures$end[i]], ylimMax/5*0.3, lwd=0.5, col="black", lend=1)
      }
      dfFeatures2<-data.frame(
        feature=c("HS","MET"),
        exon=c("v3","v6")
      )
      for (i in 1:nrow(dfFeatures2)) {
        text(((GeneStructure_GOI$Start[GeneStructure_GOI$Exon==dfFeatures2$exon[i]]+GeneStructure_GOI$End[GeneStructure_GOI$Exon==dfFeatures2$exon[i]])/2), ylimMax/5*0.15, dfFeatures2$feature[i], col="gray", cex=0.3, srt=45)
      }    
    }    
    
    #Actin
    df<-actb
    dfs<-splice_actb
    xlim<-c(5530601,5527147)
    chr<-"chr7"
    main="ACTB"
    df2<- data.frame(V2=min(xlim):max(xlim), V4=0)
    df<-merge(df,df2,by="V2", all.y=T)
    df$V3[is.na(df$V3)]<-0
    df$V3<-df$V4+df$V3
    dfs2<- data.frame(V2=min(xlim):max(xlim), V4=0)
    dfs<-merge(dfs,dfs2,by="V2", all.y=T)
    dfs$V3[is.na(dfs$V3)]<-0
    dfs$V3<-dfs$V4+dfs$V3
    ylim<-max(df$V3)
    if (ylim<10){
      ylim=10
    }
    plot(df$V2, df$V3, type="l" , cex=0.1, xlim=xlim, ylim=c(ylim/-6, ylim), xlab=chr, ylab="depth", col="white" , main=main)
    polygon(c(min(df$V2)-1,df$V2, max(df$V2)+1), c(0,df$V3,0), col="lightblue", border=F)
    polygon(c(min(dfs$V2)-1,dfs$V2, max(dfs$V2)+1), c(0,dfs$V3,0), col="black", border=F)
    abline(h=0)
    segments(5530601,	ylim/-9, 5530542, ylim/-9,lwd=2, col="red", lend=1)
    segments(5529684,	ylim/-9, 5529535, ylim/-9,lwd=2, col="red", lend=1)
    segments(5529657,	ylim/-9, 5529535, ylim/-9,lwd=5, col="red", lend=1)
    segments(5529400,	ylim/-9, 5529161, ylim/-9,lwd=5, col="red", lend=1)
    segments(5528719,	ylim/-9, 5528281, ylim/-9,lwd=5, col="red", lend=1)
    segments(5528185,	ylim/-9, 5528004, ylim/-9,lwd=5, col="red", lend=1)
    segments(5527891,	ylim/-9, 5527751, ylim/-9,lwd=5, col="red", lend=1)
    segments(5527891,	ylim/-9, 5527147, ylim/-9,lwd=2, col="red", lend=1)
  dev.off()
}
