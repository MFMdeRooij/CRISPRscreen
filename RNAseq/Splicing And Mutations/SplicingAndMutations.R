# With this script you can look at splicing and mutations, of human CD44 which has constant and variable exon usage, 
# or CXCL12 which has different isoforms
# Follow the lines with hastag signs in R studio
# Author: M.F.M. de Rooij PhD, Amsterdam UMC, Spaargaren Lab, 2022, info: m.f.derooij@amsterdamumc.nl
##################################################################################################################################
#                                                         SETTINGS

# Workdirectory Linux (This R script works also in Windows, but the mapping should be in Linux)
wd <- "~/BioLin/RNAseq/"

# Gene structure file (the first letter determines the species followed by gene symbol)
GeneStructure<- "hcxcl12.csv"

# Minimal mapping quality: Hisat2: 0 = multiple mapped reads with lost of mismatches/indels, 1 = multiple mapped reads, 60 = unique mapped reads
minMapQ <- 60

# Minimal frequency of intron abundancy (compared to the most abundant intron (%)
minIntronFreq<- 1

# Tresholds mutations: minMutFreq = minimal frequency of mutation (%), minMutCoverage =  minimal coverage at mutation spot (number of reads)
MutminFreq<- 0
MutminCoverage<- 0
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

# For other genes, you can find in Ensembl the exon loci (or in the GTF file), which you can put in hgeneX.csv, and 
# also adjust the most upstream and downstream loci in SplicingAndMutationsHg38.sh
# cxcl12 <- gtf[grep("CXCL12", gtf$V34),]
# write.csv(cxcl12,"cxcl12.csv", row.names=F, quote=F)
# Adjust this table to the example (hgeneX.csv)

# R Packages
#install.packages("BiocManager")
#BiocManager::install(c("scales","Rsamtools","stringr"))

# Now you are ready to look at the splice variants and mutations of your gene of interest

# Search for RNAseq data on https://www.ncbi.nlm.nih.gov/sra/, and fill in the SRR IDs, cell 
# line ID and a P (paired-end) or S (single-end) separated by commas (SRR8615282,HCT116,P), and run the 
# script (./SplicingAndMutationsHg38.sh on the command line) 
# The first time, make the bash script excutable (chmod 755 SplicingAndMutationsHg38.sh)
# You can also run the bash script here by uncommenting and run the next line:

#system("bash SplicingAndMutationsHg38.sh")

# For mouse you do the same for mm39, and use SplicingAndMutationsMm39.sh

#system("bash SplicingAndMutationsMm39.sh")

# To run the this script in Rstudio, put this script file, MapSamples.txt, hgeneX.csv, and hgeneX.fa 
# (hcxcl12.fa/hcd44.fa) and hactb.fa and the folders with BAM files in the workdirectory.
##################################################################################################################################
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

DetectMutations<-function(wd,ID, GeneSymbol, species, xlim, minMapQ, df, MutminFreq, MutminCoverage){
  GOI<- pileup(paste0(ID, ".",species, tolower(GeneSymbol), ".bam"), paste0(ID, ".", species, tolower(GeneSymbol), ".bam.bai"), pileupParam=PileupParam(max_depth=10^9, min_mapq=minMapQ, include_insertions=T))
  ref<-readChar(paste0("../", species, tolower(GeneSymbol),".fa"), file.info(paste0("../", species, tolower(GeneSymbol),".fa"))$size)
  ref<- sub('^>.+\n',"", ref, perl=T)
  ref<- gsub('\n',"", ref)
  poslim<-c(min(xlim),max(xlim))
  df_ref<- data.frame(pos=poslim[1]:poslim[2], RefNucleotide=unlist(strsplit(ref, split = "")), stringsAsFactors = F)
  GOImut<-merge(GOI, df_ref, by="pos")
  GOImut$mutation<-ifelse(toupper(GOImut$nucleotide) == toupper(GOImut$RefNucleotide), 0, 1)
  GOImut<-GOImut[GOImut$mutation==1,]
  coverage<-df[1:2]
  colnames(coverage)<-c("pos", "coverage")
  GOImut<-merge(GOImut, coverage, by="pos", all.x=T)
  GOImut$freq<-GOImut$count/GOImut$coverage*100 
  GOImut<-GOImut[GOImut$freq>=MutminFreq,]
  GOImut<-GOImut[GOImut$coverage>=MutminCoverage,]
  GOImut$color<-ifelse(toupper(GOImut$nucleotide) == "A", "green", 
                    ifelse(toupper(GOImut$nucleotide) == "C", "blue",
                           ifelse(toupper(GOImut$nucleotide) == "G", "black", 
                                  ifelse(toupper(GOImut$nucleotide) == "T", "red", 
                                         ifelse(toupper(GOImut$nucleotide) == "-", "white",
                                                ifelse(toupper(GOImut$nucleotide) == "+", "yellow", "pink"))))))
  return(GOImut)
}

setwd(wd)
dfID<- read.csv("MapSamples.txt", header=F, stringsAsFactors = F)
GeneStructure_GOI<- read.csv(GeneStructure, stringsAsFactors = F)
species<- strsplit(GeneStructure, "")[[1]][1]
direction <- "pos"
if (GeneStructure_GOI$Start[1]>GeneStructure_GOI$Start[nrow(GeneStructure_GOI)]){
  direction<- "neg"
}

for (i in 1:nrow(dfID)){
  ID<-paste0(dfID[i,1], "_", dfID[i,2])
  setwd(paste0(wd,ID))
  GeneSymbol<- sub(".", "", strsplit(GeneStructure, ".csv")[[1]][1])
  GOI<-data.frame(seqnames=0,pos=0,count=0)
  splice_GOI<-data.frame(seqnames=0,pos=0,count=0)
  actb<-data.frame(seqnames=0,pos=0,count=0)
  splice_actb<-data.frame(seqnames=0,pos=0,count=0)
  try(
    GOI<- pileup(paste0(ID, ".", species, tolower(GeneSymbol), ".bam"), paste0(ID, ".", species, tolower(GeneSymbol), ".bam.bai"), pileupParam=PileupParam(max_depth=10^9, min_mapq=minMapQ))
  )  
  try(
    splice_GOI<- pileup(paste0("splice_", ID, ".", species, tolower(GeneSymbol), ".bam"), paste0("splice_", ID, ".", species, tolower(GeneSymbol), ".bam.bai"), pileupParam=PileupParam(max_depth=10^9, min_mapq=minMapQ))
  ) 
  try(
    actb<-pileup(paste0(ID, ".", species, "actb.bam"), paste0(ID, ".", species, "actb.bam.bai"), pileupParam=PileupParam(max_depth=10^9, min_mapq=minMapQ))
  )
  try(
    splice_actb<-pileup(paste0("splice_", ID, ".", species, "actb.bam"), paste0("splice_", ID, ".", species, "actb.bam.bai"), pileupParam=PileupParam(max_depth=10^9, min_mapq=minMapQ))
  )
  
  # Count introns
  
  # Load the bam file
  bam<-data.frame(bam=character())
  
  try(
    bam <- as.data.frame(readBAM(paste0("splice_",ID,".", species, tolower(GeneSymbol), ".bam")))
  )
  dfIntron<-data.frame(Intron=character())
  if (nrow(bam)>0){
    # Use only 1-times mapped reads (some reads are inconclusive)
    bam<- bam[bam$mapq>=minMapQ,]
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
          bamMultInt$endintron[i]<- bamMultInt$startintron[i]+as.numeric(str_extract_all(bamMultInt$up[i], "\\d+$"))
        }
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
    dfIntron$end<-dfIntron$end-1
    # Filter low abundant introns out
    dfIntron<-dfIntron[dfIntron$x>=(max(dfIntron$x)*(minIntronFreq/100)),]
  }
  
  # Read densities
  df<-GOI
  df<-df[,c("pos","count")]
  df<-aggregate(df$count, by=list(df$pos), FUN=sum)
  colnames(df)<-c("pos","count")
  dfs<-splice_GOI
  if (nrow(dfs)>0){
    dfs<-aggregate(dfs$count, by=list(dfs$pos), FUN=sum)
    colnames(dfs)<-c("pos","count")
  }
  xlim<-c(min(GeneStructure_GOI[c("Start","End")]),max(GeneStructure_GOI[c("Start","End")]))
  if (direction=="neg"){
    xlim<-rev(xlim)
  }
  chr<-paste0("Chr",GeneStructure_GOI$Chr[1])
  df2<- data.frame(pos=min(xlim):max(xlim), startCount=0)
  df<-merge(df,df2,by="pos", all.y=T)
  df$count[is.na(df$count)]<-0
  df$count<-df$startCount+df$count
  dfs2<- data.frame(pos=min(xlim):max(xlim), startCount=0)
  if (nrow(dfs)>0){
    dfs<-merge(dfs,dfs2,by="pos", all.y=T)
    dfs$count[is.na(dfs$count)]<-0
    dfs$count<-dfs$startCount+dfs$count
  } else{
      dfs<-dfs2
      colnames(dfs)[2]<-"count"
    }
  ylim<-max(df$count)
  if (ylim<10){
    ylim=10
  }

  # Mutations
  GOImut<-data.frame(mut=character())
  try(
    GOImut<- DetectMutations(wd=wd,ID=ID, GeneSymbol=GeneSymbol, species=species, xlim=xlim, minMapQ=minMapQ, df=df, MutminFreq=MutminFreq, MutminCoverage=MutminCoverage)
  )
  pdf(paste0("../", ID, "_", GeneSymbol,".", minMapQ, ".",minIntronFreq,".pdf"),10,5)
    par(mar=c(4,5,2,1))
    
    # Plot read density
    par(fig=c(0.1,0.83,0.5,1))
    
    if (species=="h"){ 
      main=toupper(GeneSymbol)
    }
    if (species=="m"){ 
      main<-gsub("(^.)(.+)","\\U\\1\\L\\2", GeneSymbol, perl=T)
    }

    plot(df$pos, df$count, type="n" , xlim=xlim, ylim=c(0, ylim), xlab="", xaxt="n", ylab="Depth", main=main)
    polygon(c(min(df$pos)-1,df$pos, max(df$pos)+1), c(0,df$count,0), col="lightblue", border=F)
    polygon(c(min(dfs$pos)-1,dfs$pos, max(dfs$pos)+1), c(0,dfs$count,0), col="azure4", border=F)
    segments(xlim[1], 0, xlim[2],0, col="gray", lwd=0.1)
  
    # Draw mutation
    if(nrow(GOImut)>0){
      for (p in unique(GOImut$pos)){
        temp<-GOImut[GOImut$pos %in% p,]
        temp<- temp[,c("count", "color")]
        temp<- aggregate(temp$count, by=list(temp$color), FUN=sum)
        colnames(temp)<- c("color","count")
        temp<- temp[order(temp$count, decreasing = T),]
        segments(p, 0,p,temp$count[1], col=temp$color[1], lwd=0.25, lend=1)
        if (nrow(temp)>1){
          for (i in 2:nrow(temp)){
            segments(p, sum(temp$count[1:i-1]),p,sum(temp$count[1:i]), col=temp$color[i],lwd=0.25, lend=1)
          }
        }
      }
    }
    
    # Draw gene structures 
    par(fig=c(0.1,0.83,0.1,0.73),new=TRUE)
    if (nrow(dfIntron)>0){
      ylimMax<-nrow(dfIntron)*-10*(5/4)
    } else {ylimMax=-12.5}
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
             adj=if (toupper(GeneSymbol)=="CD44"){if(GeneStructure_GOI$Exon[i]=="v4"){0.7}else{if(GeneStructure_GOI$Exon[i]=="v5"){0.3}else{0.5}}})
      }
    }
    for (i in 1:nrow(GeneStructure_GOI)){
      segments(GeneStructure_GOI$Start[i],	ylimMax/5*0.7, GeneStructure_GOI$End[i], ylimMax/5*0.7,lwd=GeneStructure_GOI$Size[i], col=GeneStructure_GOI$Color[i], lend=1)
    }
    
    # Draw introns
    if(nrow(dfIntron)>0){
      for (i in 1:nrow(dfIntron)){
        rect(dfIntron$start[i],	(ylimMax/5)+i*((ylimMax-(ylimMax/5))/(nrow(dfIntron)+1))-ylimMax/100, dfIntron$end[i], (ylimMax/5)+i*((ylimMax-(ylimMax/5))/(nrow(dfIntron)+1))+ylimMax/100, lwd=0.1, col="white", border="black")
        rect(dfIntron$start[i],	(ylimMax/5)+i*((ylimMax-(ylimMax/5))/(nrow(dfIntron)+1))-ylimMax/100, dfIntron$end[i], (ylimMax/5)+i*((ylimMax-(ylimMax/5))/(nrow(dfIntron)+1))+ylimMax/100, lwd=0.1, col=alpha("black",dfIntron$x[i]/max(dfIntron$x)), border="black")
        text(if (direction=="neg"){dfIntron$start[i]} else{dfIntron$end[i]}, (ylimMax/5)+i*((ylimMax-(ylimMax/5))/(nrow(dfIntron)+1)), pos=4,offset=0.1, paste0("(",dfIntron$x[i],")"), srt = 0, xpd = TRUE, cex=0.3)
      }
    }
    
    # Special features for CD44
    if (toupper(GeneSymbol)=="CD44"){
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
    
    # Actin as control gene
    df<-actb
    df<-df[,c("pos","count")]
    df<-aggregate(df$count, by=list(df$pos), FUN=sum)
    colnames(df)<-c("pos","count")
    dfs<-splice_actb
    dfs<-dfs[,c("pos","count")]
    dfs<-aggregate(dfs$count, by=list(dfs$pos), FUN=sum)
    colnames(dfs)<-c("pos","count")
    if (species=="h"){ 
      xlim<-c(5530601,5527147)
      chr<-"Chr7"
      GeneSymbol<-"ACTB"
    } 
    if (species=="m"){
      xlim<-c(142892509,142888870)
      chr<-"Chr5"
      GeneSymbol<-"Actb"
    }
    df2<- data.frame(pos=min(xlim):max(xlim), startCount=0)
    df<-merge(df,df2,by="pos", all.y=T)
    df$count[is.na(df$count)]<-0
    df$count<-df$startCount+df$count
    dfs2<- data.frame(pos=min(xlim):max(xlim), startCount=0)
    dfs<-merge(dfs,dfs2,by="pos", all.y=T)
    dfs$count[is.na(dfs$count)]<-0
    dfs$count<-dfs$startCount+dfs$count
    ylim<-max(df$count)
    if (ylim<10){
      ylim=10
    }
    
    # Mutations
    GOImut<-data.frame(mut=character())
    try(
      GOImut<- DetectMutations(wd=wd,ID=ID, GeneSymbol=GeneSymbol, species=species, xlim=xlim, minMapQ=minMapQ, df=df, MutminFreq=MutminFreq, MutminCoverage=MutminCoverage)
    )
    
    # Plot read density
    plot(df$pos, df$count, type="n" , xlim=xlim, ylim=c(ylim/-6, ylim), xlab=chr, ylab="Depth", main=GeneSymbol)
    polygon(c(min(df$pos)-1,df$pos, max(df$pos)+1), c(0,df$count,0), col="lightblue", border=F)
    polygon(c(min(dfs$pos)-1,dfs$pos, max(dfs$pos)+1), c(0,dfs$count,0), col="azure4", border=F)
    segments(xlim[1], 0, xlim[2],0, col="gray", lwd=0.1)
    abline(h=ylim/-30)
    
    # Draw mutation
    if(nrow(GOImut)>0){
      for (p in unique(GOImut$pos)){
        temp<-GOImut[GOImut$pos %in% p,]
        temp<- temp[order(temp$count, decreasing = T),]
        segments(p, 0,p,temp$count[1], col=temp$color[1], lwd=0.25, lend=1)
        if (nrow(temp)>1){
          for (i in 2:nrow(temp)){
            segments(p, sum(temp$count[1:i-1]),p,sum(temp$count[1:i]), col=temp$color[i], lwd=0.25, lend=1)
          }
        }
      }
    }
    # Draw actin exons
    if (species=="h"){ 
      # hg38
      segments(5530601,	ylim/-8, 5530542, ylim/-8,lwd=2, col="red", lend=1)
      segments(5529684,	ylim/-8, 5529535, ylim/-8,lwd=2, col="red", lend=1)
      segments(5529657,	ylim/-8, 5529535, ylim/-8,lwd=5, col="red", lend=1)
      segments(5529400,	ylim/-8, 5529161, ylim/-8,lwd=5, col="red", lend=1)
      segments(5528719,	ylim/-8, 5528281, ylim/-8,lwd=5, col="red", lend=1)
      segments(5528185,	ylim/-8, 5528004, ylim/-8,lwd=5, col="red", lend=1)
      segments(5527891,	ylim/-8, 5527751, ylim/-8,lwd=5, col="red", lend=1)
      segments(5527891,	ylim/-8, 5527147, ylim/-8,lwd=2, col="red", lend=1)
    }
    if (species=="m"){ 
      # mm39
      segments(142892509,	ylim/-9, 142892407, ylim/-9,lwd=2, col="red", lend=1)
      segments(142891447,	ylim/-9, 142891319, ylim/-9,lwd=2, col="red", lend=1)
      segments(142891441,	ylim/-9, 142891319, ylim/-9,lwd=5, col="red", lend=1)   
      segments(142891231,	ylim/-9, 142890992, ylim/-9,lwd=5, col="red", lend=1)
      segments(142890537,	ylim/-9, 142890099, ylim/-9,lwd=5, col="red", lend=1)
      segments(142890003,	ylim/-9, 142889822, ylim/-9,lwd=5, col="red", lend=1)
      segments(142889696,	ylim/-9, 142889553, ylim/-9,lwd=5, col="red", lend=1)
      segments(142889696,	ylim/-9, 142888870, ylim/-9,lwd=2, col="red", lend=1)
    }
  dev.off()
}
