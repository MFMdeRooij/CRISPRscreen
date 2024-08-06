# Use the CRISPRScreenAnalysis.R output files, put this script in the same folder, adjust the settings, and run the script in Rstudio.
# This script plots the log2 median fold change vs RRAscore (Dep or Enr are dependent on fold change) of all genes of all comparisons (T1/T0, T2/T0, and T2/T1).
# Author: M.F.M. de Rooij PhD, Amsterdam UMC, Spaargaren Lab, 2024, info: m.f.derooij@amsterdamumc.nl
################################################################################
#install.packages(c("scales", "RColorBrewer","devtools"))
#devtools::install_github("JosephCrispell/basicPlotteR")
library("scales")
library("RColorBrewer")
library("basicPlotteR")
################################################################################
#                              SETTINGS

# Cell line ID and graph titles
CellLine<- "Z138"
Titles<- c("Control (T1/T0)","Venetoclax (T1/T0)","Venetoclax/Control")

# Workdirectory
setwd(paste0("C:/BioWin/Volcano/", CellLine))

# Tophits: 0: Best per volcano, 1: Top5 in Synthetic lethal comparison, 2: Genes of interest
th<- 2
# If 2, Which genes:
Genes_of_interest<- c("CSNK2A1", "CSNK2A2", "CSNK2B")

# Colors (All genes, positive and negative controls):
ColAll<- "lightgray"
ColP<- "red"
ColN<- "blue"
if (th!=0){
  ColP<- "lightpink1"
  ColN<- "lightskyblue"
  # colors of hits are taken by a color pallete
}
################################################################################

files<- list.files(pattern="Genes.csv$")

pdf(paste0(CellLine, "_VolcanoPlots",as.character(th),".pdf"), width=4*length(files), height=4)
par(mfrow=c(1,length(files)))

g<- 1
for (file in files) {
  #file<- files[1]  
  # Read DESeq2 data
  df_genes<- read.csv(file, stringsAsFactors = FALSE)
  df_genes<- df_genes[df_genes$GeneSymbol!="NonTargetingControlGuideForHuman",]
  
  df_genes$l2mfc<- log2(df_genes$MedianFoldChange)
  
  # Give NA a number
  df_genes$l2mfc[is.na(df_genes$l2mfc)]<- 0
  df_genes$rhoDepleted[is.na(df_genes$rhoDepleted)]<- 1
  df_genes$rhoEnriched[is.na(df_genes$rhoEnriched)]<- 1
  df_genes$fdrDepleted[is.na(df_genes$fdrDepleted)]<- 1
  df_genes$fdrEnriched[is.na(df_genes$fdrEnriched)]<- 1

  # RRA score and FDR dependent on depleted/enriched
  df_genes$rho<- df_genes$rhoDepleted
  df_genes$rho[df_genes$l2mfc > 0]<- df_genes$rhoEnriched[df_genes$l2mfc > 0]
  df_genes$fdr<- df_genes$fdrDepleted
  df_genes$fdr[df_genes$l2mfc > 0]<- df_genes$fdrEnriched[df_genes$l2mfc > 0]
  
  # Colors
  df_genes$col<- ColAll
  df_genes$col[df_genes$Type=="p"]<- ColP
  df_genes$col[df_genes$Type=="n"]<- ColN
  df_genes$col2 <- df_genes$col
  df_genes$pch<- 21
  df_genes$pch[df_genes$fdr < 0.1 & df_genes$l2mfc < 0]<- 25
  df_genes$pch[df_genes$fdr < 0.1 & df_genes$l2mfc > 0]<- 24
  df_genes$alpha <- 0.7
  df_genes$cex <- 0.7
  
  # Axes limits
  xrangeVOL<- c(floor(min(df_genes$l2mfc, na.rm=TRUE)), ceiling(max(df_genes$l2mfc, na.rm=TRUE)))
  yrangeVOL<- c(0, max(-log10(df_genes$rho)*1.2, na.rm=TRUE))
  
  if (th==2){
    tophitsVolcano<- Genes_of_interest
  }
  
  if (th==1){
    GenesHits<- read.csv(files[3], stringsAsFactors = FALSE)
    GenesHits<- GenesHits[GenesHits$GeneSymbol!="NonTargetingControlGuideForHuman",]
    GenesHits$l2mfc<- log2(GenesHits$MedianFoldChange)
    # Give NA a number
    GenesHits$l2mfc[is.na(GenesHits$l2mfc)]<- 0
    GenesHits$rhoDepleted[is.na(GenesHits$rhoDepleted)]<- 1
    GenesHits$rhoEnriched[is.na(GenesHits$rhoEnriched)]<- 1
    GenesHits$fdrDepleted[is.na(GenesHits$fdrDepleted)]<- 1
    GenesHits$fdrEnriched[is.na(GenesHits$fdrEnriched)]<- 1
    
    # Select hits (Top 5 by depletion-statistics)
    tophitsVolcano<- GenesHits$GeneSymbol[order(GenesHits$rhoDepleted)][1:5]
  }
  
  if (th==0){
    GenesDep<- df_genes[df_genes$l2mfc < 0,]
    GenesEnr<- df_genes[df_genes$l2mfc > 0,]
  
    GenesDep<- GenesDep[order(GenesDep$rhoDepleted),]
    GenesEnr<- GenesEnr[order(GenesEnr$rhoEnriched),]
    tophitsVolcanoA<- c(GenesDep$GeneSymbol[1:5], GenesEnr$GeneSymbol[1:5])
  
    GenesDep<- GenesDep[order(GenesDep$l2mfc),]
    GenesEnr<- GenesEnr[order(GenesEnr$l2mfc, decreasing = T),]
    tophitsVolcanoB<- c(GenesDep$GeneSymbol[1:5], GenesEnr$GeneSymbol[1:5])
  
    tophitsVolcano<- unique(c(tophitsVolcanoA,tophitsVolcanoB))
  }
  
  # Mix essential and non-essential randomly
  df_genesx<- df_genes[df_genes$Type=="x",]
  df_genesx<- df_genesx[!df_genesx$GeneSymbol %in% tophitsVolcano,]
  df_genesc<- df_genes[df_genes$Type!="x",]
  df_genesc<- df_genesc[!df_genesc$GeneSymbol %in% tophitsVolcano,]
  df_genesh<- df_genes[match(tophitsVolcano,df_genes$GeneSymbol),]
  
  set.seed(101)
  df_genesc<- df_genesc[sample(1:nrow(df_genesc)),]
  df_genes<-rbind(df_genesx,df_genesc)
  df_genes<-rbind(df_genes,df_genesh)
  
  if (th != 0){
    # Color hits
    df_genes$col[df_genes$GeneSymbol %in% tophitsVolcano]<- brewer.pal(length(tophitsVolcano), "Spectral")
    df_genes$alpha[df_genes$GeneSymbol %in% tophitsVolcano]<- 1
    df_genes$cex[df_genes$GeneSymbol %in% tophitsVolcano]<- 1.2
    df_genes$col2[df_genes$GeneSymbol %in% tophitsVolcano]<- "black"
  }
  
  # Write Volcano Plot in PDF 
  par(bty="l")
  # Plot axes
  plot(0, pch = '', 
       main= Titles[g], 
       xlab= "Log2 median fold change",
       ylab= "RRA score",
       cex.lab=1, cex.axis=1, las=1, xlim=xrangeVOL, ylim=yrangeVOL, xaxt = "n", yaxt = "n")
  
  # Vertical and horizontal lines
  xrange<- xrangeVOL[2]-xrangeVOL[1]
  xticks<- 2
  if (xrange < 8){xticks<- 1}
  if (xrange <= 3){xticks<- 0.5}
  
  yrange<- yrangeVOL[2]-yrangeVOL[1]
  yticks<- 10
  if (yrange < 35){yticks<- 5}
  if (yrange <= 14){yticks<- 2}

  axis(1, at = seq(xrangeVOL[1], xrangeVOL[2], by = xticks))
  axis(2, at = seq(yrangeVOL[1], yrangeVOL[2], by = yticks))
  abline(v=xticks*(-100:100), lty=3, col="gray")
  abline(h=yticks*(0:100), lty=3, col="gray")
  
  # Actual points
  points(df_genes$l2mfc, -log10(df_genes$rho), type="p", pch=df_genes$pch, bg=alpha(df_genes$col,df_genes$alpha), col=alpha(df_genes$col2,df_genes$alpha), cex=df_genes$cex, lwd=0.5)

  # 0 line
  abline(v=0, lty=2)

  # Show gene symbols of top hits
  addTextLabels(df_genesh$l2mfc,-log10(df_genesh$rho),df_genesh$GeneSymbol, avoidPoints = TRUE,
                keepLabelsInside = TRUE, col.label="black", cex.label=0.7)
  
  # Extra titles
  text(0,yrangeVOL[2], substitute(paste(italic('Depleted    Enriched'))))
  g<- g+1
}  
dev.off()
