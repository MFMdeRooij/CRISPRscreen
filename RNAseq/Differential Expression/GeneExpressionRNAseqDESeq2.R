## RNAseq differential expression:
# - Make count tables with the GeneExpressionRNAseqHisat2.R script
# - Use the RNAseqCountTableRawProteinCoding.csv count table
# - This script also uses the CRISPRScreenAnalysisLibraries.csv file to highlight the essentials/nonessential genes, 
#   you can change these genes by for example certain pathway genes; this script also extract the kinome in the end
# - Fill in the Design table in RNAseqDesign.csv
# - Adjust the settings and run the script
# Author: M.F.M. de Rooij PhD, Amsterdam UMC, Spaargaren Lab, 2022, info: m.f.derooij@amsterdamumc.nl
######################################################################################
## Install DESeq2 and other packages
#install.packages(c("BiocManager", "devtools", "pheatmap", "RColorBrewer", "gtools", "scales"))
#BiocManager::install("DESeq2")
#devtools::install_github("JosephCrispell/basicPlotteR")
library("DESeq2")
library("pheatmap")
library("RColorBrewer")
library("gtools")
library("scales")
library("basicPlotteR")
######################################################################################
#                                     SETTINGS

# Workdirectory (folder in which the count tables are located, use always slash (/) instead of backslash)
Workdirectory<- "C:/BioWin/RNAseq/"

# Which count table?
Filename<- "RNAseqCountTableRawProteinCoding.csv"

# Fill in the table in RNAseqDesign.csv (Group = tumor-subtypes, Rep = replicates (when paired, this should be matched))

# Round numbers in output table: 0 = Yes, 1 = No
RoundNumbers<- 0

# Paired replicates: 0 = Paired, 1 = Unpaired
Paired<- 1

# Minimal fold change of guides to be a hit: 1 = No minimal fold change,  >1: The minimal fold change (linear scale, 2^abs(l2fc))
minimalFoldChange<- 1

# MA plots of all genes: 0 = Yes, 1 = No, 2 = Top 10, 3 = Genes of interest
MA_all_genes<- 2

# Genes of interest
if (MA_all_genes==3){
  q<- function(...) {
    sapply(match.call()[-1], deparse)
  }
  interestingGenes <-q(
    # Which genes?  
    BTK, SYK, PIK3CA, PIK3CD, PIK3R1, LYN
  )
}

# Colors MA plots (All genes, positive and negative controls (for CRISPR screens), hits):
ColAll<- "lightgray"
ColP<- "lightpink1"
#ColP<- "red"
ColN<- "lightskyblue"  
#ColN<- "blue"
ColH<- "black"
######################################################################################
setwd(Workdirectory)

# Make a data folder
dirname2<- paste0(Filename,Sys.time())
dirname1<- gsub("[[:punct:]]", "", dirname2) 
dirname<- gsub("\\s", "", dirname1) 
dir.create(dirname)

# Read count table
df_raw<- read.csv(file=Filename, sep=",", header=TRUE, stringsAsFactors = FALSE)

# Gene annotation
df_Gene_ID<- df_raw[,1:2]  

# Sample metadata
df_design<- read.csv("RNAseqDesign.csv")
df_colData<- data.frame(Group=df_design$Group, Rep=df_design$Rep) 
df_colData$Group<- as.factor(df_colData$Group)
df_colData$Rep<- as.factor(df_colData$Rep)

# Count table for DESEq2
counts<- df_raw[,-1:-2]
rownames(counts)<- df_raw[,1]
counts<-counts[,df_design$Sample]

# DESeq2 pipeline
if (Paired==0) {
  dds<- DESeqDataSetFromMatrix(countData = counts, colData = df_colData, design = ~ Rep + Group)
}
if (Paired==1) {
  dds<- DESeqDataSetFromMatrix(countData = counts, colData = df_colData, design = ~ Group)
}  
dds<- DESeq(dds, betaPrior = TRUE)

# Define comparisons
Groups<- unique(df_design$Group)
combi<- as.data.frame(combinations(length(Groups),2,Groups))
combi<- combi[nrow(combi):1,ncol(combi):1]

# Analyse comparisons
for (i in 1:nrow(combi)){
  con1<- combi[i,1]
  con2<- combi[i,2]
  res<- results(dds, contrast=c("Group",con2,con1), addMLE=TRUE)
  
  # DESeq2 data table
  df_res<- as.data.frame(res)
  
  df_baseMeanPerLvl<- as.data.frame(sapply(levels(dds$Group), function(lvl) rowMeans(counts(dds, normalized=TRUE)[,dds$Group==lvl])))
  df_res$BaseMeanA<- df_baseMeanPerLvl[[con1]]
  df_res$logBaseMeanA<- log(df_baseMeanPerLvl[[con1]]+1)/log(10)
  df_res$BaseMeanB<- df_baseMeanPerLvl[[con2]]
  
  # GeneIDs
  df_res$ensembl_gene_id<- rownames(df_res)
  df_res<- merge(df_Gene_ID, df_res, by='ensembl_gene_id', all.y=T)
  df_res$log2FoldChange[is.na(df_res$log2FoldChange)]<-0
  df_res$FoldChange<- 2^df_res$log2FoldChange
  
  if (minimalFoldChange > 1) {
    df_res$padj[!is.na(df_res$padj) & df_res$FoldChange >= (1/minimalFoldChange) & df_res$FoldChange <= minimalFoldChange]<- 1
  }
  
  # CRISPR screen controls
  df_Gene_IDX<- read.csv("CRISPRScreenAnalysisLibraries.csv", sep=',', header=TRUE, stringsAsFactors = FALSE)
  df_Control<- df_Gene_IDX[1:927,1:3]
  PC<- "Essential"
  NC<- "NonEssential"
  df_res$Type<- "x"
  df_res[df_res$hgnc_symbol %in% df_Control[[PC]][nchar(df_Control[[PC]])>0],"Type"]<- "p"
  df_res[df_res$hgnc_symbol %in% df_Control[[NC]][nchar(df_Control[[NC]])>0],"Type"]<- "n"
  df_PC<- df_res[df_res$Type=="p",]
  df_NC<- df_res[df_res$Type=="n",]
  
  # Hits
  GenesDiff<- df_res[!is.na(df_res$padj),]
  GenesDiff<- GenesDiff[GenesDiff$padj<0.1,]
  GenesDiffDep<- GenesDiff[order(GenesDiff$log2FoldChange),]
  tophitsDep<- GenesDiffDep$hgnc_symbol[nchar(GenesDiffDep$hgnc_symbol)>0]
  GenesDiffEnr<- GenesDiff[order(GenesDiff$log2FoldChange, decreasing = T),]
  tophitsEnr<- GenesDiffEnr$hgnc_symbol[nchar(GenesDiffEnr$hgnc_symbol)>0]
  tophits<- c(tophitsDep[1:10], tophitsEnr[1:10])
  
  # Top hits volcano (5 best FC & padj - Dep & Enr)
  tophitsVolcanoA<- c(tophitsDep[1:5], tophitsEnr[1:5])
  GenesDiffDep<- GenesDiff[GenesDiff$log2FoldChange<0,]
  GenesDiffDep<- GenesDiffDep[order(GenesDiffDep$padj),]
  tophitsDep<- GenesDiffDep$hgnc_symbol[nchar(GenesDiffDep$hgnc_symbol)>0]
  GenesDiffEnr<- GenesDiff[GenesDiff$log2FoldChange>0,]
  GenesDiffEnr<- GenesDiffEnr[order(GenesDiffEnr$padj),] 
  tophitsEnr<- GenesDiffEnr$hgnc_symbol[nchar(GenesDiffEnr$hgnc_symbol)>0]
  tophitsVolcanoB<- c(tophitsDep[1:5], tophitsEnr[1:5])
  tophitsVolcano<- unique(c(tophitsVolcanoA,tophitsVolcanoB))
  
  # Write gene table
  
  # Guide Table
  df_res_print<- df_res[,c("hgnc_symbol","ensembl_gene_id","Type","BaseMeanA","BaseMeanB","FoldChange","pvalue","padj")]
  df_res_print<- df_res_print[order(df_res_print$FoldChange),]
  if (RoundNumbers==0){
    df_res_print[,c("BaseMeanA","BaseMeanB")]<- round(df_res_print[,c("BaseMeanA","BaseMeanB")],0)
    df_res_print[,c("FoldChange","pvalue","padj")]<- signif(df_res_print[,c("FoldChange","pvalue","padj")],4)
  }                                                                                       
  write.csv(df_res_print, paste0(dirname,"/DESeq2_RNAseq_",con1, "vs",con2,".csv"), row.names = FALSE)
  
  # Write Plots in PDF 
  pdf(paste0(dirname,"/DESeq2_RNAseq_",con1, "vs",con2,"_Plots.pdf"), width=7, height=7)
  
  # DESEq2 plots
  plotDispEsts(dds, main="Dispersion plot")
  
  rld<-rlog(dds, blind=FALSE)
  sampleDists<- dist(t(assay(rld)))
  sampleDistMatrix<- as.matrix(sampleDists)
  rownames(sampleDistMatrix)<- paste(rld$Group, rld$Rep)
  colnames(sampleDistMatrix)<- paste(rld$Group, rld$Rep)
  colors<- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
  
  pheatmap(sampleDistMatrix,
           clustering_distance_rows=sampleDists,
           clustering_distance_cols=sampleDists,
           col=colors)
  
  print(plotPCA(rld, intgroup=c("Group", "Rep")))
  
  # Volcano & MA plots
  Genes_of_interest<- c(" ", "Hitlist", if (MA_all_genes==0) {df_hits_A$GeneSymbol}, 
                        if (MA_all_genes==2) {tophits}, if (MA_all_genes==3) {interestingGenes})
  
  # Color and marker information
  df_res$col<- ColAll
  df_res$col[df_res$Type=="p"]<- ColP
  df_res$col[df_res$Type=="n"]<- ColN
  df_res$pch<- 16
  df_res$pch[df_res$padj<0.1 & df_res$log2FoldChange<0]<- 25
  df_res$pch[df_res$padj<0.1 & df_res$log2FoldChange>0]<- 24
  
  # Mix essential and non-essential randomly
  df_resx<- df_res[df_res$Type=="x",]
  df_resc<- df_res[df_res$Type!="x",]
  set.seed(101)
  df_resc<- df_resc[sample(1:nrow(df_resc)),]
  df_res<- rbind(df_resx,df_resc)
  
  volcano<- df_res
  volcano$padj[is.na(volcano$padj)]<- 1
  volcano$l2bm<- log2(volcano$baseMean+1)
  palette <- colorRampPalette(c("cyan","black"))
  volcano$col<-palette(max(volcano$l2bm)*100+1)[volcano$l2bm*100+1]

  # Prevent a log-bug when Rho is 0
  volcano$padj[volcano$padj==0]<- min(volcano$padj[volcano$padj!=0])/10
  
  # Axes limits
  xrange<- c(min(df_res$logBaseMeanA, na.rm=TRUE)-0.5, max(df_res$logBaseMeanA, na.rm=TRUE)+0.5)
  yrange<- c(min(df_res$log2FoldChange, na.rm=TRUE)-0.5, max(df_res$log2FoldChange, na.rm=TRUE)+0.5)
  
  # Density plot fold change
  denY_tot<-density(df_res$log2FoldChange, from=yrange[1], to=yrange[2], na.rm=T)
  denY_tot$y[1]<- 0
  denY_tot$y[length(denY_tot$y)]<- 0
  denY_PC<-density(df_PC$log2FoldChange, from=yrange[1], to=yrange[2], na.rm=T)
  denY_PC$y[1]<- 0
  denY_PC$y[length(denY_PC$y)]<- 0
  denY_NC<-density(df_NC$log2FoldChange, from=yrange[1], to=yrange[2], na.rm=T)
  denY_NC$y[1]<- 0
  denY_NC$y[length(denY_NC$y)]<- 0
  denYMax<- max(c(denY_tot$y, denY_PC$y, denY_NC$y))
  # Density plot read count
  denX_tot<- density(df_res$logBaseMeanA, from=xrange[1], to=xrange[2], na.rm=T)
  denX_tot$y[1]<- 0
  denX_tot$y[length(denX_tot$y)]<- 0
  denX_PC<- density(df_PC$logBaseMeanA, from=xrange[1], to=xrange[2], na.rm=T)
  denX_PC$y[1]<- 0
  denX_PC$y[length(denX_PC$y)]<- 0
  denX_NC<- density(df_NC$logBaseMeanA, from=xrange[1], to=xrange[2], na.rm=T)
  denX_NC$y[1]<- 0
  denX_NC$y[length(denX_NC$y)]<- 0
  denXMax<- max(c(denX_tot$y, denX_PC$y, denX_NC$y))
  
  # Volcano plot
  
  # Margins
  par(mar=c(4,4,4,4))
  par(fig=c(0.1,0.8,0.1,0.8))
  par(bty="l")
  
  # Empty plot
  plot(0, pch = '', 
       main= "Volcano plot",
       xlab= expression("Log"[2]*" fold change"), 
       ylab= expression("-Log"[10]*" FDR"),
       cex.lab=1, cex.axis=1, las=1, xlim=yrange, ylim=c(0,max(-log10(volcano$padj))))
  
  # Vertical and horizontal lines
  abline(v=2.5*(-10:10), lty=3, col="gray")
  abline(h=10*(0:20), lty=3, col="gray")
  
  # Actual plot
  points(volcano$log2FoldChange, -log10(volcano$padj), type="p", col=volcano$col, bg=volcano$col, cex=1, pch=volcano$pch)
  
  # Extra lines
  abline(v=c(-1,1), lty=2)
  abline(h=-log10(0.1), lty=2)
  
  # Show gene symbols of tophits
  df_volcanoGenes<- volcano[volcano$hgnc_symbol %in% tophitsVolcano,]
  addTextLabels(df_volcanoGenes$log2FoldChange,-log10(df_volcanoGenes$padj),df_volcanoGenes$hgnc_symbol, avoidPoints = TRUE,
                keepLabelsInside = TRUE, col.label="black", cex.label=1)
  
  # colorbar
  par(new=T)
  par(fig=c(0.75,1,0.2,0.8))
  x=1
  y=seq(1,max(volcano$l2bm)*100+1,len=100)
  z=matrix(1:100,nrow=1)
  image(x,((y-1)/100),z,col=palette(max(y))[y], xaxt='n', xlab="", 
        ylab=expression("Log"[2]*"(average read counts + 1)"))
  
  # MA plot
  par(mfrow=c(1,1))
  for(Gene in Genes_of_interest){
    df_GOI<- df_res[df_res$hgnc_symbol %in% Gene,]
    # Main MA plot
    par(mar=c(4,4,0,0))
    par(fig=c(0.1,0.7,0.1,0.7))
    plot(df_res$logBaseMeanA, df_res$log2FoldChange, type="p", col=df_res$col, bg=df_res$col, cex=1, pch=df_res$pch, xlab="Log10 Average Read Counts (Control)", 
         ylab="Log2 Fold Change", cex.lab=1, cex.axis=1, xlim=xrange, ylim=yrange)
    if (Gene=="Hitlist"&& nrow(GenesDiff)>=1){
      GenesDiff<- df_res[df_res$ensembl_gene_id %in% GenesDiff$ensembl_gene_id,]
      points(GenesDiff$logBaseMeanA, GenesDiff$log2FoldChange, type="p", col=ColH, bg=ColH, cex=1, pch=GenesDiff$pch)
    }
    if (nrow(df_GOI)>=1){
      points(df_GOI$logBaseMeanA, df_GOI$log2FoldChange, type="p", col=ColH, bg=ColH, cex=1.5, pch=df_GOI$pch)
    }
    legend("bottomleft",legend=c( if (nchar(Gene)>1 && !is.na(Gene)) {Gene}, "All genes", "Essentials","Non-essentials",
                                  'Significantly enriched', 'Significantly depleted'), cex=0.8, pch=c(if (nchar(Gene)>1 && !is.na(Gene)) 
                                  {16},16,16,16, 24,25), col=c(if (nchar(Gene)>1 && !is.na(Gene)) {ColH}, ColAll,ColP,ColN, "black", "black"))
    abline(median(df_res$log2FoldChange, na.rm=TRUE),0, col=1, lty=3, untf=TRUE)
    
    # Density fold change
    par(mar=c(4,0,0,4))
    par(fig=c(0.7,0.9,0.1,0.7),new=TRUE)
    plot(denY_tot$y, denY_tot$x, ylim=yrange, xlim=(c(0,denYMax)), type='l', axes=FALSE, col=ColAll, xlab="", 
         ylab="", lwd=2)
    par(new=TRUE)
    lines(denY_PC$y, denY_PC$x, col=ColP, lwd=2)
    lines(denY_NC$y, denY_NC$x, col=ColN, lwd=2)
    rgb.val<- col2rgb(ColAll)
    polygon(denY_tot$y, denY_tot$x, col=rgb(rgb.val[1]/255,rgb.val[2]/255,rgb.val[3]/255,alpha=0.3), lwd=0.1)
    rgb.val<- col2rgb(ColP)
    polygon(denY_PC$y, denY_PC$x, col=rgb(rgb.val[1]/255,rgb.val[2]/255,rgb.val[3]/255,alpha=0.3), lwd=0.1)
    rgb.val<- col2rgb(ColN)
    polygon(denY_NC$y, denY_NC$x, col=rgb(rgb.val[1]/255,rgb.val[2]/255,rgb.val[3]/255,alpha=0.3), lwd=0.1)
    if (Gene=="Hitlist" && nrow(GenesDiff)>1){
      denY_hits_total<-density(GenesDiff$log2FoldChange, from=yrange[1], to=yrange[2], na.rm=T)
      denY_hits_total$y[1]<- 0
      denY_hits_total$y[length(denY_hits_total$y)]<- 0
      lines(denY_hits_total$y, denY_hits_total$x, col=ColH, lwd=2)
      rgb.val<- col2rgb(ColH)
      polygon(denY_hits_total$y, denY_hits_total$x, col=rgb(rgb.val[1]/255,rgb.val[2]/255,rgb.val[3]/255,alpha=0.3), lwd=0.1)
    }
    if (nrow(df_GOI)>1){
      par(new=TRUE)
      denY_GOI<- density(df_GOI$log2FoldChange, from=yrange[1], to=yrange[2], na.rm=T)
      denY_GOI$y[1]<- 0
      denY_GOI$y[length(denY_GOI$y)]<- 0
      lines(denY_GOI$y, denY_GOI$x, col=ColH, lwd=2)
      rgb.val<- col2rgb(ColH)
      polygon(denY_GOI$y, denY_GOI$x, col=rgb(rgb.val[1]/255,rgb.val[2]/255,rgb.val[3]/255,alpha=0.3), lwd=0.1)
    }
    
    # Density read counts
    par(mar=c(0,4,4,0))
    par(fig=c(0.1,0.7,0.7,0.9),new=TRUE)
    plot(denX_tot$x, denX_tot$y, main="MA plot", cex.main=1.5, xlim=xrange, ylim=c(0,denXMax), type='l', axes=FALSE, col=ColAll, xlab="", 
         ylab="", lwd=2)
    lines(denX_PC, col=ColP, lwd=2)
    lines(denX_NC, col=ColN, lwd=2)
    rgb.val<- col2rgb(ColAll)
    polygon(denX_tot, col=rgb(rgb.val[1]/255,rgb.val[2]/255,rgb.val[3]/255,alpha=0.3), lwd=0.1)
    rgb.val<- col2rgb(ColP)
    polygon(denX_PC, col=rgb(rgb.val[1]/255,rgb.val[2]/255,rgb.val[3]/255,alpha=0.3), lwd=0.1)
    rgb.val<- col2rgb(ColN)
    polygon(denX_NC, col=rgb(rgb.val[1]/255,rgb.val[2]/255,rgb.val[3]/255,alpha=0.3), lwd=0.1)
    if (Gene=="Hitlist" && nrow(GenesDiff)>1){
      denX_hits_total<-density(GenesDiff$logBaseMeanA, from=xrange[1], to=xrange[2], na.rm=T)
      denX_hits_total$y[1]<- 0
      denX_hits_total$y[length(denX_hits_total$y)]<- 0
      lines(denX_hits_total, col=ColH, lwd=2)
      rgb.val<- col2rgb(ColH)
      polygon(denX_hits_total, col=rgb(rgb.val[1]/255,rgb.val[2]/255,rgb.val[3]/255,alpha=0.3), lwd=0.1)
    }
    if (nrow(df_GOI)>1){
      par(new=TRUE)
      denX_GOI<- density(df_GOI$logBaseMeanA, from=xrange[1], to=xrange[2], na.rm=T)
      denX_GOI$y[1]<- 0
      denX_GOI$y[length(denX_GOI$y)]<- 0
      lines(denX_GOI, col=ColH, lwd=2)
      rgb.val<- col2rgb(ColH)
      polygon(denX_GOI, col=rgb(rgb.val[1]/255,rgb.val[2]/255,rgb.val[3]/255,alpha=0.3), lwd=0.1)
    }
  }
  
  dev.off()
  
  # Kinome(-like) genes taken from Brunello CRISPR library
  kinome<- unique(df_Gene_IDX$GeneSymbol)
  kinome<-kinome[1:764]
  dataPrintKin<-df_res_print[df_res_print$hgnc_symbol %in% kinome,]
  write.csv(dataPrintKin, paste0(dirname,"/DESeq2_RNAseq_kinome",con1, "vs",con2,".csv"), row.names=F, quote=F)
}
