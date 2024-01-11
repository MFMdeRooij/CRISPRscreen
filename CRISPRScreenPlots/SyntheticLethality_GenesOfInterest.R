# Use the CRISPRScreenAnalysis.R output files of a synthetic lethality screen, adjust the settings, and run the script in R studio.
# This script can normalize the median log2 fold change to the essential and non-essential genes of a synthetic lethality screen, and plots T1control/T0 against T1treated/T0. 
# This normalization can improve the comparison treated - control if the treated arm did not have equal cell divisions, however the separation between the essentials and 
# non-essentials will not be improved. Synthetic lethal genes will be located around the lower half of the vertical 0 axis.
# Author: M.F.M. de Rooij PhD, Amsterdam UMC, Spaargaren Lab, 2023, info: m.f.derooij@amsterdamumc.nl
##################################################################################################################################
## Install required pacakges once
#install.packages("devtools)
#devtools::install_github("JosephCrispell/basicPlotteR")
library("basicPlotteR")
##################################################################################################################################
#                              SETTINGS

# Put this script in the folder where the count tables are located
folder<- dirname(rstudioapi::getActiveDocumentContext()$path)
## Fill in workdirectory (folder in which the count tables are located, use always slash (/) instead of backslash)
#folder<- "H:/BioWin/CRISPRscreen/Namalwa/"

# Cell line:
cellID<- "Z138"

# Is there a T1drug/T1control comparison, 0: no, 1: yes
t2t1com<- 1

# Size graph
size<- 7

# Genes to emphasize, 0: all significant (from T1drug/T1control comparison), 1: specific genes  
allsignif<- 0

GenesOfInterest<- NULL
if (allsignif==1){
  # If specific genes, Which ones?
  GenesOfInterest<- c("BTK", "SYK", "PIK3R1")
}

# Show all gene symbols, 0: no, 1: yes
GeneSymbol<- 0

# Axes labels:
xlab<- "Control (Relative log2 median fold change)"
ylab<- "Venetoclax (Relative log2 median fold change)"

# # BCR-controlled adhesion screens:   
# xlab<- "PMA (Log2 median fold change)"
# ylab<- expression(alpha*"IgM (Log2 median fold change)")

# Normalize to essentials and non-essentials (only for lethality), 0: no, 1: yes
NormalizeX<- 1
NormalizeY<- 1

# Normalize to (log)mean or median of the (non-)essentials, 0: mean, 1: median
meanOrmedian<- 1

# Additionally normalize y by linear regression (so y=0+1x of the controls), 0: no, 1: yes
NormalizeLR<- 1

#Axes limit, 0 = automatic, 1: custom
Axlim<- 0
# If automatic, Equal X and Y axes, 0 = no, 1: yes
XYequal<- 1

if (Axlim==1){
  # Custom axes limits: 
  xmin<- -3
  xmax<- 1
  xticks<- 1
  
  ymin<- -3
  ymax<- 1
  yticks<- 1
}

# Colors:
call<- 'lightgray'
cpos<- 'red'
cneg<- 'blue'
chit<- 'black'
##################################################################################################################################
setwd(folder)

Control<- read.csv("DESeq2 T0vsT1 Genes.csv", stringsAsFactors=F)
if (NormalizeX == 0){
  Control$Nmfc<- log2(Control$MedianFoldChange)
}
if (NormalizeX == 1){
  if (meanOrmedian==0){
    mCP<- mean(log2(Control$MedianFoldChange[Control$Type=='p']))
    mCN<- mean(log2(Control$MedianFoldChange[Control$Type=='n']))
  }
  if (meanOrmedian==1){
    mCP<- median(log2(Control$MedianFoldChange[Control$Type=='p']))
    mCN<- median(log2(Control$MedianFoldChange[Control$Type=='n']))
  }
  Control$Nmfc<- (log2(Control$MedianFoldChange)-mCN)/abs(mCP-mCN)
}
ControlG<- Control[,c("GeneSymbol","Type","Nmfc")]

Treated <- read.csv("DESeq2 T0vsT2 Genes.csv", stringsAsFactors=F)
if (NormalizeY == 0){
  Treated$Nmfc<- log2(Treated$MedianFoldChange)
}
if (NormalizeY == 1){
  if (meanOrmedian==0){
    mTP<- mean(log2(Treated$MedianFoldChange[Treated$Type=='p']))
    mTN<- mean(log2(Treated$MedianFoldChange[Treated$Type=='n']))
  }
  if (meanOrmedian==1){
    mTP<- median(log2(Treated$MedianFoldChange[Treated$Type=='p']))
    mTN<- median(log2(Treated$MedianFoldChange[Treated$Type=='n']))
  }
  Treated$Nmfc<- (log2(Treated$MedianFoldChange)-mTN)/abs(mTP-mTN)
}
TreatedG<- Treated[,c("GeneSymbol","Nmfc")]

if (t2t1com==1){
  TC<- read.csv("DESeq2 T1vsT2 Genes.csv", stringsAsFactors=F)
  TC$Stat<- apply(TC[,c("rhoDepleted","rhoEnriched")], 1, FUN=min)
  TC$fdr<- apply(TC[,c("fdrDepleted","fdrEnriched")], 1, FUN=min)
  TCG<- TC[,c("GeneSymbol","Stat", "fdr")]
  if (allsignif==0){
    GenesOfInterest<- TCG$GeneSymbol[TCG$fdr<0.1]
  }
}
Combi<- merge(ControlG, TreatedG, by="GeneSymbol")
if (t2t1com==1){
  Combi<- merge(Combi, TCG, by="GeneSymbol")
}

if (NormalizeLR==1){
  con<- Combi[Combi$Type!="x",]
  modelC<- lm(Nmfc.y ~ Nmfc.x, data=con)
  Combi$Nmfc.y<- (Combi$Nmfc.y-modelC$coefficients[1])/modelC$coefficients[2]
}

pos<- Combi[Combi$Type=="p",]
neg<- Combi[Combi$Type=="n",]
hit<- Combi[Combi$GeneSymbol %in% GenesOfInterest,]

if (Axlim==0){
  # Calculate axis limits:
  xmin<- round(min(Combi$Nmfc.x),2)-0.3
  xmax<- round(max(Combi$Nmfc.x),2)+0.3
  ymin<- round(min(Combi$Nmfc.y),2)-0.3
  ymax<- round(max(Combi$Nmfc.y),2)+0.3
  if (XYequal==1){
    xmin<- min(xmin,ymin)
    ymin<- min(xmin,ymin)
    xmax<- max(xmax,ymax)
    ymax<- max(xmax,ymax)
  }
  xticks<- round((xmax-xmin)/5.1,2)
  yticks<- round((ymax-ymin)/5.1,2)
}

pdf(paste0("CRISPR_SL_",cellID,"_R.pdf"),size,size)
  par(mar=c(4,4,0,0))
  par(fig=c(0.1,0.7,0.1,0.7))
  plot(Combi$Nmfc.x, Combi$Nmfc.y, xlab=xlab, ylab=ylab, cex.lab=1, cex.axis=1, 
       xlim=c(xmin,xmax), ylim=c(ymin,ymax), 
       xaxp=c(xmin,xmin+xticks*floor((xmax-xmin)/xticks),floor((xmax-xmin)/xticks)),
       yaxp=c(ymin,ymin+yticks*floor((ymax-ymin)/yticks),floor((ymax-ymin)/yticks)),    
       pch=16, col=call, cex=if(t2t1com==1){0.3+-0.1*log10(Combi$Stat)}else{0.5})

  if (GeneSymbol == 1){
    text(Combi$Nmfc.x, Combi$Nmfc.y, labels=Combi$GeneSymbol, cex=0.8, col="gray", adj = c(-0.2,0.5), srt=22.5)
  }

  points(pos$Nmfc.x, pos$Nmfc.y, pch=16, col=cpos, cex=if(t2t1com==1){0.3+-0.1*log10(pos$Stat)}else{0.5})
  points(neg$Nmfc.x, neg$Nmfc.y, pch=16, col=cneg, cex=if(t2t1com==1){0.3+-0.1*log10(neg$Stat)}else{0.5})
  if (length(GenesOfInterest)>0){
    points(hit$Nmfc.x, hit$Nmfc.y, pch=16, col=chit, cex=if(t2t1com==1){0.3+-0.1*log10(hit$Stat)}else{0.5})
    addTextLabels(hit$Nmfc.x, hit$Nmfc.y, hit$GeneSymbol, avoidPoints = TRUE,
                        keepLabelsInside = TRUE, col.label="black", cex.label=1)
  }

  # Linreg prediction interval from controls
  con<- rbind(pos,neg)
  modelC<- lm(Nmfc.y ~ Nmfc.x, data=con)
  
  xvalues <- seq(xmin-1, xmax+1, length.out = 10)
  predictions <- predict(modelC, newdata = data.frame(Nmfc.x = xvalues), 
                         interval = "prediction", level=0.95)
  
  lower_values <- predictions[, "lwr"]
  upper_values <- predictions[, "upr"]
  modal_values <- predictions[, "fit"]
  lines(xvalues, modal_values, col = "purple", lwd = 1)
  lines(xvalues, lower_values, col = "purple", lty = 2)
  lines(xvalues, upper_values, col = "purple", lty = 2)
  
  # 0 lines
  abline(v=0, col=cneg, lty=3)
  abline(h=0, col=cneg, lty=3)
  if (NormalizeX == 1){
    abline(v=-1, col=cpos, lty=3)
  }
  if (NormalizeY == 1){
    abline(h=-1, col=cpos, lty=3)
  }
  abline(0,1, col="black", lty=2)
  legend(xmin,ymax,legend=c("All genes", "Essentials","Non-essentials"), pch=16, cex=0.8, col=c(call,cpos,cneg))

  # Density y axis
  denY_PC<- density(pos$Nmfc.y, from=ymin, to=ymax, na.rm=TRUE)
  denY_PC$y[1]<- 0
  denY_PC$y[length(denY_PC$y)]<- 0
  denY_NC<- density(neg$Nmfc.y, from=ymin, to=ymax, na.rm=TRUE)
  denY_NC$y[1]<- 0
  denY_NC$y[length(denY_NC$y)]<- 0
  denYMax <- max(c(denY_PC$y, denY_NC$y))
  par(mar=c(4,0,0,4))
  par(fig=c(0.7,0.9,0.1,0.7),new=TRUE)
  plot(denY_PC$y, denY_PC$x, ylim=c(ymin,ymax), xlim=(c(0,denYMax)), type='l', axes=FALSE, col=2, xlab="", ylab="", lwd=2)
  lines(denY_NC$y, denY_NC$x, col=4, lwd=2)
  rgb.val<- col2rgb(cpos)
  polygon(denY_PC$y, denY_PC$x, col=rgb(rgb.val[1]/255,rgb.val[2]/255,rgb.val[3]/255,alpha=0.3), lwd=0.1)
  rgb.val<- col2rgb(cneg)
  polygon(denY_NC$y, denY_NC$x, col=rgb(rgb.val[1]/255,rgb.val[2]/255,rgb.val[3]/255,alpha=0.3), lwd=0.1)
  
  # Density x axis
  denX_PC<- density(pos$Nmfc.x, from=xmin, to=xmax, na.rm=TRUE)
  denX_PC$y[1]<- 0
  denX_PC$y[length(denX_PC$y)]<- 0
  denX_NC<- density(neg$Nmfc.x, from=xmin, to=xmax, na.rm=TRUE)
  denX_NC$y[1]<- 0
  denX_NC$y[length(denX_NC$y)]<- 0
  denXMax <- max(c(denX_PC$y, denX_NC$y))
  
  par(mar=c(0,4,4,0))
  par(fig=c(0.1,0.7,0.7,0.9),new=TRUE)
  plot(denX_PC$x, denX_PC$y, xlim=c(xmin,xmax), ylim=c(0,denXMax), type='l', axes=FALSE, col=2, xlab="", ylab="", lwd=2, main=cellID)
  lines(denX_NC$x, denX_NC$y, col=4, lwd=2)
  rgb.val<- col2rgb(cpos)
  polygon(denX_PC, col=rgb(rgb.val[1]/255,rgb.val[2]/255,rgb.val[3]/255,alpha=0.3), lwd=0.1)
  rgb.val<- col2rgb(cneg)
  polygon(denX_NC, col=rgb(rgb.val[1]/255,rgb.val[2]/255,rgb.val[3]/255,alpha=0.3), lwd=0.1)
dev.off()