# Use the CRISPRScreenAnalysis.R output files of a synthetic lethality screen, adjust the setting, and run the script in R studio.
# This script normalizes the median log2 fold change to the essential and non-essential genes of a synthetic lethality screen, and plots T1control/T0 against T1treated/T0. 
# This normalization can improve the comparison treated - control if the treated arm did not have equal cell divisions, however the separation between the essentials and 
# non-essentials will not be improved. Synthetic lethal genes will be on the lower half of the vertical 0 axis.
# Author: M.F.M. de Rooij PhD, Amsterdam UMC, Spaargaren Lab, 2019, info: m.f.derooij@amsterdamumc.nl

# Screen data
Workdirectory<- "H:/BioWin/RPCIWM1/"

# Cell line ID
cell<- "RPCIWM1"

# Drug  (x-axis is control)
drug<- "Ruxolitinib"

# Genes to emphasize
GenesOfInterest<- c("JAK1", "JAK2")

# X and Y axis limits (Relative median log2 fold change)
min<- -2
max<- 1

# Show gene symbols (0=yes, 1=no)
GeneSymbol<- 0
##############################################################
setwd(Workdirectory)

Control<-read.csv("DESeq2 T0vsT1 Genes.csv", stringsAsFactors=F)
medianCP<-log2(median(Control$MedianFoldChange[Control$Type=='p']))
medianCN<-log2(median(Control$MedianFoldChange[Control$Type=='n']))
Control$Nmfc<- (log2(Control$MedianFoldChange)-medianCN)/abs(medianCP-medianCN)
ControlG<- Control[,c(1,2,15)]

Treated <-read.csv("DESeq2 T0vsT2 Genes.csv", stringsAsFactors=F)
medianTP<-log2(median(Treated$MedianFoldChange[Treated$Type=='p']))
medianTN<- log2(median(Treated$MedianFoldChange[Treated$Type=='n']))
Treated$Nmfc<- (log2(Treated$MedianFoldChange)-medianTN)/abs(medianTP-medianTN)
TreatedG<- Treated[,c(1,15)]

TC <-read.csv("DESeq2 T1vsT2 Genes.csv", stringsAsFactors=F)
TC$Stat<- apply(TC[,c(9,12)], 1, FUN=min)
TCG<- TC[,c(1,15)]

Combi<- merge(ControlG, TreatedG, by="GeneSymbol")
Combi<- merge(Combi, TCG, by="GeneSymbol")

pos<-Combi[Combi$Type=="p",]
neg<-Combi[Combi$Type=="n",]
hit<-Combi[Combi$GeneSymbol %in% GenesOfInterest,]

Combi<- Combi[!(Combi$GeneSymbol %in% GenesOfInterest),]

denY_PC<-density(pos$Nmfc.y, from=min, to=max, na.rm=TRUE)
denY_NC<-density(neg$Nmfc.y, from=min, to=max, na.rm=TRUE)
denYMax <- max(c(denY_PC$y, denY_NC$y))
denX_PC<-density(pos$Nmfc.x, from=min, to=max, na.rm=TRUE)
denX_NC<-density(neg$Nmfc.x, from=min, to=max, na.rm=TRUE)
denXMax <- max(c(denX_PC$y, denX_NC$y))

pdf(paste0("CRISPR_SL_",cell,"_",drug,".pdf"),10,10)
  par(mar=c(4,5,2,1))
  par(fig=c(0.1,0.87,0.1,0.87))
  plot(Combi$Nmfc.x, Combi$Nmfc.y, xlab="Control (Relative median log2 fold change)", ylab=paste0(drug, " (Relative median log2 fold change)"), xlim=c(min,max), ylim=c(min,max), pch=16, col=3, cex=-0.2*log10(Combi$Stat))
  if (GeneSymbol == 0){
    text(Combi$Nmfc.x, Combi$Nmfc.y, labels=Combi$GeneSymbol, cex=0.5, col=8, adj = -0.3)
    points(Combi$Nmfc.x, Combi$Nmfc.y,pch=16, col=3, cex=-0.2*log10(Combi$Stat))
  }
  points(pos$Nmfc.x, pos$Nmfc.y, pch=16, col=2, cex=0.7)
  points(neg$Nmfc.x, neg$Nmfc.y, pch=16, col=4, cex=0.7)
  points(hit$Nmfc.x, hit$Nmfc.y, pch=16, col=1, cex=-0.2*log10(hit$Stat))
  text(hit$Nmfc.x, hit$Nmfc.y, labels=hit$GeneSymbol, cex=1, col=1, adj = -0.3)
  abline(v=0, col=1)
  abline(h=0, col=1)
  abline(0,1)
  
  par(fig=c(0.75,1,0.1,0.87),new=TRUE)
  plot(denY_PC$y, denY_PC$x, ylim=c(min,max), xlim=(c(0,denYMax)), type='l', axes=FALSE, col=2, xlab="", ylab="", lwd=2)
  lines(denY_NC$y, denY_NC$x, col=4, lwd=2)
  
  par(fig=c(0.1,0.87,0.75,1),new=TRUE)
  plot(denX_PC$x, denX_PC$y, xlim=c(min,max), ylim=c(0,denXMax), type='l', axes=FALSE, col=2, xlab="", ylab="", lwd=2, main=cell)
  lines(denX_NC$x, denX_NC$y, col=4, lwd=2)
  legend(min,denXMax,legend=c("Essential","Non-essential"), pch=16, col=c(2,4))
dev.off()
