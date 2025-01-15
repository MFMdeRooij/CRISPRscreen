# Use the CRISPRScreenAnalysis.R output files, put this script in the same folder, adjust the settings, and run the script in Rstudio.
# This script plots all the guides or genes of max 11 genes if interest in MA and Volcano plots of all comparisons (T1/T0, T2/T0, and T2/T1).
# Author: M.F.M. de Rooij PhD, Amsterdam UMC, Spaargaren Lab, 2024, info: m.f.derooij@amsterdamumc.nl
################################################################################
## Install package
#install.packages(c("rstudioapi", "scales", "RColorBrewer","devtools"))
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
setwd(paste0("C:/BioWin/MA/", CellLine))

# Tophits: 0: Best per comparison, 1: Top5 in Synthetic lethal comparison, 2: Genes of interest
th<- 1
# If 2, Which genes:
Genes_of_interest<- c("CSNK2A1", "CSNK2A2", "CSNK2B")

# Add guide ID's in MA plot, 0: No, 1: Yes
guideID<- 0

# Colors (All genes, positive and negative controls):
ColAll<- "lightgray"
ColP<- "lightpink1"
ColN<- "lightskyblue"
# colors of hits are taken by a color pallete

################################################################################

filesMA<- list.files(pattern="Guides.csv$")
filesVOL<- list.files(pattern="Genes.csv$")

# Tophits
tophits<- list()
for (f in 1:length(filesVOL)){
  if (th==0){
    df_genes<- read.csv(filesVOL[f], stringsAsFactors = FALSE)
    df_genes<- df_genes[df_genes$GeneSymbol!="NonTargetingControlGuideForHuman",]
    
    df_genes$l2mfc<- log2(df_genes$MedianFoldChange)
    
    # Give NA a number
    df_genes$l2mfc[is.na(df_genes$l2mfc)]<- 0
    df_genes$rhoDepleted[is.na(df_genes$rhoDepleted)]<- 1
    df_genes$rhoEnriched[is.na(df_genes$rhoEnriched)]<- 1
    df_genes$fdrDepleted[is.na(df_genes$fdrDepleted)]<- 1
    df_genes$fdrEnriched[is.na(df_genes$fdrEnriched)]<- 1
    GenesDep<- df_genes[df_genes$l2mfc < 0,]
    GenesEnr<- df_genes[df_genes$l2mfc > 0,]
    
    GenesDep<- GenesDep[order(GenesDep$rhoDepleted),]
    GenesEnr<- GenesEnr[order(GenesEnr$rhoEnriched),]
    tophitsVolcano<- c(GenesEnr$GeneSymbol[1:5], rev(GenesDep$GeneSymbol[1:5]))
    
    # GenesDep<- GenesDep[order(GenesDep$l2mfc),]
    # GenesEnr<- GenesEnr[order(GenesEnr$l2mfc, decreasing = T),]
    # tophitsVolcanoB<- c(GenesDep$GeneSymbol[1:5], GenesEnr$GeneSymbol[1:5])
    # tophitsVolcano<- unique(c(tophitsVolcano,tophitsVolcanoB))
  }
  
  if (th==1){
    GenesHits<- read.csv(filesVOL[3], stringsAsFactors = FALSE)
    GenesHits<- GenesHits[GenesHits$GeneSymbol!="NonTargetingControlGuideForHuman",]
    GenesHits$l2mfc<- log2(GenesHits$MedianFoldChange)
    # Give NA a number
    GenesHits$l2mfc[is.na(GenesHits$l2mfc)]<- 0
    GenesHits$rhoDepleted[is.na(GenesHits$rhoDepleted)]<- 1
    GenesHits$rhoEnriched[is.na(GenesHits$rhoEnriched)]<- 1
    GenesHits$fdrDepleted[is.na(GenesHits$fdrDepleted)]<- 1
    GenesHits$fdrEnriched[is.na(GenesHits$fdrEnriched)]<- 1
    
    # Select hits (Top 5 by depletion-statistics)
    tophitsVolcano<- rev(GenesHits$GeneSymbol[order(GenesHits$rhoDepleted)][1:5])
  }
  
  if (th==2){
    tophitsVolcano<- rev(Genes_of_interest)
  }
  
  tophits[[f]]<- tophitsVolcano
}

pdf(paste0(CellLine, "_MA_VolcanoPlots_",as.character(th),".pdf"), width=4*length(filesVOL), height=8)
par(mfrow=c(2,length(filesVOL)))

# MA plot
f<- 1
for (file in filesMA) {
  df_res<- read.csv(file, stringsAsFactors = F)
  df_res$logBaseMeanA<- log(df_res$BaseMeanA)/log(10)
  df_res$log2FoldChange<- log(df_res$FoldChange)/log(2)
  
  # Give NA a number
  df_res$logBaseMeanA[is.na(df_res$logBaseMeanA)]<- 0
  df_res$log2FoldChange[is.na(df_res$log2FoldChange)]<- 0
  df_res$padj[is.na(df_res$padj)]<- 1
  
  # Colors
  df_res$col<- ColAll
  df_res$col[df_res$Type=="p"]<- ColP
  df_res$col[df_res$Type=="n"]<- ColN
  df_res$col2 <- df_res$col
  df_res$pch<- 21
  df_res$pch[df_res$padj<0.1 & df_res$log2FoldChange<0] <- 25
  df_res$pch[df_res$padj<0.1 & df_res$log2FoldChange>0] <- 24
  df_res$alpha <- 0.7
  df_res$cex <- 1
  
  # Axes limits
  xrangeMA<- c(0, ceiling(max(df_res$logBaseMeanA, na.rm=TRUE)))
  yrangeMA<- c(floor(min(df_res$log2FoldChange, na.rm=TRUE)), ceiling(max(df_res$log2FoldChange, na.rm=TRUE)))
  
  # Mix essential and non-essential randomly
  df_resx<- df_res[df_res$Type=="x",]
  df_resx<- df_resx[!df_resx$GeneSymbol %in% tophits[[f]],]
  df_resc<- df_res[df_res$Type!="x",]
  df_resc<- df_resc[!df_resc$GeneSymbol %in% tophits[[f]],]
  df_resh<- df_res[df_res$GeneSymbol %in% tophits[[f]],]
  df_resh<- df_resh[order(match(df_resh$GeneSymbol,tophits[[f]])),]
  set.seed(101)
  df_resc<- df_resc[sample(1:nrow(df_resc)),]
  df_res<-rbind(df_resx,df_resc)
  df_res<-rbind(df_res,df_resh)
  
  # Color hits
  if (length(tophits[[f]])>2){
    df_colorhits<- data.frame(GeneSymbol=tophits[[f]], Color= rev(brewer.pal(length(tophits[[f]]), "Spectral")))
  } else{
    df_colorhits<- data.frame(GeneSymbol=tophits[[f]], Color= c("yellow", "lightgreen")[1:length(tophits[[f]])])
  }
  df_res$col[df_res$GeneSymbol %in% df_colorhits$GeneSymbol]<- apply(df_res[df_res$GeneSymbol %in% df_colorhits$GeneSymbol,], 1, function(x) df_colorhits$Color[df_colorhits$GeneSymbol==x[1]]) 
  df_res$alpha[df_res$GeneSymbol %in% tophits[[f]]]<- 1
  df_res$cex[df_res$GeneSymbol %in% tophits[[f]]]<- 1.5
  df_res$col2[df_res$GeneSymbol %in% tophits[[f]]]<- "black"
  
  # Draw plot
  par(bty="l")
  # Plot axes
  plot(0, pch = '', 
       main= Titles[f], 
       xlab= "Log10 avergae read counts",
       ylab= "Log2 fold change",
       cex.lab=1, cex.axis=1, las=1, xlim=xrangeMA, ylim=yrangeMA, xaxt = "n", yaxt = "n")
  
  # Vertical and horizontal lines
  xrange<- xrangeMA[2]-xrangeMA[1]
  xticks<- 1
  if (xrange <= 3){xticks<- 0.5}
  
  yrange<- yrangeMA[2]-yrangeMA[1]
  yticks<- 2
  if (yrange < 8){yticks<- 1}
  if (yrange <= 3){yticks<- 0.5}
  
  axis(1, at = seq(xrangeMA[1], xrangeMA[2], by = xticks))
  axis(2, at = seq(yrangeMA[1], yrangeMA[2], by = yticks))
  abline(v=xticks*(-100:100), lty=3, col="gray")
  abline(h=yticks*(-100:100), lty=3, col="gray")
  
  # Actual points
  points(df_res$logBaseMeanA, df_res$log2FoldChange, type="p", pch=df_res$pch, bg=alpha(df_res$col,df_res$alpha), col=alpha(df_res$col2,df_res$alpha), cex=df_res$cex, lwd=0.5)
  
  if (guideID==1){
    addTextLabels(df_res$logBaseMeanA[df_res$GeneSymbol %in% tophitsVolcano], df_res$log2FoldChange[df_res$GeneSymbol %in% tophitsVolcano],
                  unlist(lapply(strsplit(as.character(df_res$Guide[df_res$GeneSymbol %in% tophitsVolcano]),"-"), "[", 2)), avoidPoints = TRUE,
                  keepLabelsInside = TRUE, col.label="black", cex.label=1)
  }
  
  # 0 line
  abline(h=0, lty=2)
  legend("bottomleft", inset=0.03,legend=df_colorhits$GeneSymbol, cex=1, pch=21,pt.bg=df_colorhits$Color, box.lty=0)
  f<- f+1
}

# Volcano plot
f<- 1
for (file in filesVOL) {
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
  yrangeVOL<- c(0, ceiling(max(-log10(df_genes$rho)*1.2, na.rm=TRUE)))
  
  # Mix essential and non-essential randomly
  df_genesx<- df_genes[df_genes$Type=="x",]
  df_genesx<- df_genesx[!df_genesx$GeneSymbol %in% tophits[[f]],]
  df_genesc<- df_genes[df_genes$Type!="x",]
  df_genesc<- df_genesc[!df_genesc$GeneSymbol %in% tophits[[f]],]
  df_genesh<- df_genes[match(tophits[[f]],df_genes$GeneSymbol),]
  set.seed(101)
  df_genesc<- df_genesc[sample(1:nrow(df_genesc)),]
  df_genes<-rbind(df_genesx,df_genesc)
  df_genes<-rbind(df_genes,df_genesh)
  
  # Color hits
  if (length(tophits[[f]])>2){
    df_colorhits<- data.frame(GeneSymbol=tophits[[f]], Color= rev(brewer.pal(length(tophits[[f]]), "Spectral")))
  } else{
    df_colorhits<- data.frame(GeneSymbol=tophits[[f]], Color= c("yellow", "lightgreen")[1:length(tophits[[f]])])
  }
  df_genes$col[df_genes$GeneSymbol %in% df_colorhits$GeneSymbol]<- apply(df_genes[df_genes$GeneSymbol %in% df_colorhits$GeneSymbol,], 1, function(x) df_colorhits$Color[df_colorhits$GeneSymbol==x[1]]) 
  df_genes$alpha[df_genes$GeneSymbol %in% tophits[[f]]]<- 1
  df_genes$cex[df_genes$GeneSymbol %in% tophits[[f]]]<- 1.5
  df_genes$col2[df_genes$GeneSymbol %in% tophits[[f]]]<- "black"
  
  # Draw Volcano plot 
  par(bty="l")
  # Plot axes
  plot(0, pch = '', 
       main= Titles[f], 
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
  f<- f+1
}
dev.off()
