# With this script you can look at splicing of human CD44, which has constant and variable exon usage.
# Follow the lines with hastag signs in R studio
# Author: M.F.M. de Rooij PhD, Amsterdam UMC, Spaargaren Lab, 2022, info: m.f.derooij@amsterdamumc.nl
##################################################################################################################################
#                                                         SETTINGS

# Workdirectory
wd <- "~/BioLin/RNAseq/"

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

# For other genes than CD44, you can find in Ensembl the exon loci (or in the GTF file), which you can put in hCD44_GeneStructure.csv, and 
# also adjust the most upstream and downstream loci in hCD44_Splicing.sh

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
GeneStructure_hCD44<- read.csv("hCD44_GeneStructure.csv", stringsAsFactors = F)
for (i in 1:nrow(dfID)){
  ID<-paste0(dfID[i,1], "_", dfID[i,2])
  setwd(paste0(wd,ID))
  cd44<-data.frame(V1=1,V2=1,V3=1)
  splice_cd44<-data.frame(V1=1,V2=1,V3=1)
  actb<-data.frame(V1=1,V2=1,V3=1)
  splice_actb<-data.frame(V1=1,V2=1,V3=1)
  try(
    cd44<-read.table(paste0(ID,".cd44.bam.txt"), sep="\t", fill=T)
  )  
  try(
    splice_cd44<-read.table(paste0("splice_",ID,".cd44.bam.txt"), sep="\t", fill=T)
  ) 
  try(
    actb<-read.table(paste0(ID,".actb.bam.txt"), sep="\t", fill=T)
  )
  try(
    splice_actb<-read.table(paste0("splice_",ID,".actb.bam.txt"), sep="\t", fill=T)
  )
  # Load the bam file
  bam<-data.frame(qname=1,flag=1,rname=1,strand=1,pos=1,qwidth=1,mapq=1,cigar=1,mrnm=1,mpos=1,isize=1,seq=1,qual=1)
  
  try(
    bam2 <- as.data.frame(readBAM(paste0("splice_",ID,".cd44.bam")))
  )
  if (nrow(bam2)>0){
    bam<-bam2
  }
  bam$up<- unlist(lapply(strsplit(as.character(bam$cigar), "N"), "[", 1)) 
  bam$startintron<-0
  for (i in 1:nrow(bam)){
    bam$startintron[i]<-bam$pos[i]+sum(as.numeric(sub("M", "", str_extract_all(bam$up[i],"(\\d+)(M)",simplify = T))))
  }
  bam$endintron<- bam$startintron+as.numeric(str_extract_all(bam$up, "\\d+$"))
  
  dfsum<- data.frame(Intron=(paste0(bam$startintron,",",bam$endintron)), count=1)
  dfIntron<-aggregate(dfsum$count, by=list(dfsum$Intron), FUN=sum)
  dfIntron<-dfIntron[order(dfIntron$x, decreasing = T),]
  
  dfIntron$start<-as.numeric(
    unlist(lapply(strsplit(as.character(dfIntron$Group.1), ","), "[", 1))
  )    
  dfIntron$start<-as.numeric(
    unlist(lapply(strsplit(as.character(dfIntron$Group.1), ","), "[", 1))
  ) 
  dfIntron$end<-as.numeric(
    unlist(lapply(strsplit(as.character(dfIntron$Group.1), ","), "[", 2))
  ) 
  dfIntron$size<- dfIntron$end-dfIntron$start
  dfIntron<- dfIntron[order(dfIntron$size, decreasing = T),]
  dfIntron<- dfIntron[order(dfIntron$start),]
  dfIntron<-dfIntron[dfIntron$x>max(dfIntron$x)/100,]
  ##############################################
  
  pdf(paste0("../", ID,".pdf"),10,5)
    df<-cd44
    dfs<-splice_cd44
    xlim<-c(min(GeneStructure_hCD44[c("Start","End")]),max(GeneStructure_hCD44[c("Start","End")]))
    chr<-"chr11"
    main="CD44"
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
    par(mar=c(4,5,2,1))
    par(fig=c(0.1,0.83,0.5,1))
    plot(df$V2, df$V3, type="l" , cex=0.1, xlim=xlim, ylim=c(0, ylim), xlab="", xaxt="n", ylab="depth", col="white" , main=main)
  
    polygon(c(min(df$V2)-1,df$V2, max(df$V2)+1), c(0,df$V3,0), col="lightblue", border=F)
    polygon(c(min(dfs$V2)-1,dfs$V2, max(dfs$V2)+1), c(0,dfs$V3,0), col="black", border=F)
  
    par(fig=c(0.1,0.83,0.1,0.73),new=TRUE)
    ylimMax<-nrow(dfIntron)*-10*(5/4)
    plot(1, type='n', xlab=chr, ylab="Splicing", xlim=xlim, yaxt="n", ylim=c(ylimMax,0))
    dfFeatures<-data.frame(
      feature=c("ECD","TM","ICD"),
      start=c("1start","17","18stop"),
      end=c("16","17","18stop")
    )
    dfFeatures2<-data.frame(
      feature=c("HS","MET"),
      exon=c("v3","v6")
    )
    for (i in 1:nrow(dfFeatures)) {
      text(((GeneStructure_hCD44$Start[GeneStructure_hCD44$Exon==dfFeatures$start[i]]+GeneStructure_hCD44$End[GeneStructure_hCD44$Exon==dfFeatures$end[i]])/2), ylimMax/5*0.1, dfFeatures$feature[i], font=2, cex=0.3, adj=0.5)
      segments(GeneStructure_hCD44$Start[GeneStructure_hCD44$Exon==dfFeatures$start[i]],	ylimMax/5*0.3, GeneStructure_hCD44$End[GeneStructure_hCD44$Exon==dfFeatures$end[i]], ylimMax/5*0.3, lwd=0.5, col="black", lend=1)
    }
    for (i in 1:nrow(dfFeatures2)) {
      text(((GeneStructure_hCD44$Start[GeneStructure_hCD44$Exon==dfFeatures2$exon[i]]+GeneStructure_hCD44$End[GeneStructure_hCD44$Exon==dfFeatures2$exon[i]])/2), ylimMax/5*0.1, dfFeatures2$feature[i], col="gray", cex=0.3, srt=45)
    }    

    for (i in 1:nrow(GeneStructure_hCD44)){
      if (!grepl("st", GeneStructure_hCD44$Exon[i])){
        rect(GeneStructure_hCD44$Start[i],	ylimMax-1000, GeneStructure_hCD44$End[i], ylimMax/5*0.95, col=alpha(GeneStructure_hCD44$Color[i],0.1), border=F)
      }
    }
    for (i in 1:nrow(GeneStructure_hCD44)){
      if (!grepl("st", GeneStructure_hCD44$Exon[i])){
        text(((GeneStructure_hCD44$Start[i]+GeneStructure_hCD44$End[i])/2), ylimMax/5*0.5, GeneStructure_hCD44$Exon[i], cex=0.3, adj=if(GeneStructure_hCD44$Exon[i]=="v4"){0.7}else{if(GeneStructure_hCD44$Exon[i]=="v5"){0.3}else{0.5}})
      }
    }
    for (i in 1:nrow(GeneStructure_hCD44)){
      segments(GeneStructure_hCD44$Start[i],	ylimMax/5*0.7, GeneStructure_hCD44$End[i], ylimMax/5*0.7,lwd=GeneStructure_hCD44$Size[i], col=GeneStructure_hCD44$Color[i], lend=1)
    }
  
    abline(h=(ylimMax/5*0.95))
    for (i in 1:nrow(dfIntron)){
      segments(dfIntron$start[i],	ylimMax/5-i*10, dfIntron$end[i], ylimMax/5-i*10, lwd=3, col=alpha("black",dfIntron$x[i]/max(dfIntron$x)), lend=1)
    }
    # for (i in 1:nrow(dfIntron)){
    #  segments(dfIntron$start[i],	ylimMax/5-i*10, dfIntron$end[i], ylimMax/5-i*10, lwd=3, col=alpha("red",0.1), lend=1)
    # }
    
    
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
    plot(df$V2, df$V3, type="l" , cex=0.1, xlim=xlim, ylim=c(ylim/-6, ylim), xlab=chr, ylab="depth", col="black" , main=main)
    polygon(c(min(df$V2)-1,df$V2, max(df$V2)+1), c(0,df$V3,0), col="lightblue", border="lightblue")
    polygon(c(min(dfs$V2)-1,dfs$V2, max(dfs$V2)+1), c(0,dfs$V3,0), col="black", border=" black")
    abline(h=0)
    segments(5530601,	ylim/-9, 5530542, ylim/-9,lwd=3, col="red")
    segments(5529684,	ylim/-9, 5529535, ylim/-9,lwd=3, col="red")
    segments(5529400,	ylim/-9, 5529161, ylim/-9,lwd=3, col="red")
    segments(5528719,	ylim/-9, 5528281, ylim/-9,lwd=3, col="red")
    segments(5528185,	ylim/-9, 5528004, ylim/-9,lwd=3, col="red")
    segments(5527891,	ylim/-9, 5527147, ylim/-9,lwd=3, col="red")
    
  dev.off()
}
