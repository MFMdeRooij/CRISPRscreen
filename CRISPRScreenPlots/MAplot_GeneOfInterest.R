# Use the CRISPRScreenAnalysis.R output files, put this script in the same folder, adjust the settings, and run the script in Rstudio.
# This script plots all the guides of a particular gene in MA plots of all comparisons (T1/T0, T2/T0, and T2/T1).
# Author: M.F.M. de Rooij PhD, Amsterdam UMC, Spaargaren Lab, 2023, info: m.f.derooij@amsterdamumc.nl
################################################################################
## Install package
#install.packages("rstudioapi")
################################################################################
#                              SETTINGS

# Cell line:
maintitle<- 'Namalwa'

# Workdirectory:
# Put this script in the folder where the files of guide and gene statistics are located
folder<- dirname(rstudioapi::getActiveDocumentContext()$path)
# # Fill in workdirectory (folder in which the count tables are located, use always slash (/) instead of backslash)
# folder<- paste0("C:/BioWin/CRISPRscreen/", CellLine)

# Which genes:
Genes_of_interest<- c(" ", "BTK", "SYK", "PIK3R1")

# Axes labels for IgM/PMA-induced Adhesion:
x1title<- "Log10 average read counts (input)"
y1title<- "Log2 fold change (PMA/input)"

x2title<- "Log10 average read counts (input)"
y2title<- expression("Log2 fold change ("*alpha*"IgM/input)")

x3title<- "Log10 average read counts (PMA)"
y3title<- expression("Log2 fold change ("*alpha*"IgM/PMA)")

# # Axes labels for (synthetic) lethality:
# x1title<- expression("Log10 average read counts (T"[0]*")")
# y1title<- expression("Log2 fold change (T"["1(DMSO)"]*"/T"["0"]*")")
# 
# x2title<- expression("Log10 average read counts (T"[0]*")")
# y2title<- expression("Log2 fold change (T"["1(DRUG)"]*"/T"["0"]*")")
# 
# x3title<- expression("Log10 average read counts (T"["1DMSO"]*")")
# y3title<- expression("Log2 fold change (T"["1(DRUG)"]*"/T"["1(DMSO)"]*")")

#Axes limit, 0 = automatic, 1: custom
Axlim<- 0
if (Axlim==1){
  # Custom axes limits:
  # Log10 average read count:
  xmin<- 0
  xmax<- 4.2
  xticks<- 1
  
  # Log2 fold change:
  ymin<- -1.25
  ymax<- 1
  yticks<- 0.25
}

# Colors (All guides, positive and negative controls, hits):
ColAll<- "lightgray"
ColP<- "lightpink1"
#ColP<- "red"
ColN<- "lightskyblue"  
#ColN<- "blue"
ColH<- "black"

################################################################################
setwd(folder) 
files<- list.files(pattern="Guides.csv$")
n=1
for (file in files) {
  df_res<- read.csv(file, stringsAsFactors = F)
  df_resx<- df_res[df_res$Type=="x",]
  df_resc<- df_res[df_res$Type!="x",]
  set.seed(101)
  df_resc<- df_resc[sample(1:nrow(df_resc)),]
  df_res<-rbind(df_resx,df_resc)
  
  df_res$logBaseMeanA<- log(df_res$BaseMeanA)/log(10)
  df_res$log2FoldChange<- log(df_res$FoldChange)/log(2)
  
  df_res$col<- ColAll
  df_res$col[df_res$Type=="p"]<- ColP
  df_res$col[df_res$Type=="n"]<- ColN
  df_res$pch<- 16
  df_res$pch[df_res$padj<0.1 & df_res$log2FoldChange<0] <- 25
  df_res$pch[df_res$padj<0.1 & df_res$log2FoldChange>0] <- 24
  
  df_PC<- df_res[df_res$Type=="p",]
  df_NC<- df_res[df_res$Type=="n",]
  
  if (Axlim==0){
    # Calculate axis limits:
    xmin<- 0
    xmax<- round(max(df_res$logBaseMeanA, na.rm=TRUE),2)+0.3
    ymin<- round(min(df_res$log2FoldChange, na.rm=TRUE),2)-0.3
    ymax<- round(max(df_res$log2FoldChange, na.rm=TRUE),2)+0.3
    xticks<- round((xmax-xmin)/5.1,2)
    yticks<- round((ymax-ymin)/5.1,2)
  }
    
  if (n==1){
    xlab<- x1title
    ylab<- y1title
  }
  if (n==2){
    xlab<- x2title
    ylab<- y2title
  }
  if (n==3){
    xlab<- x3title
    ylab<- y3title
  }
  pdf(paste0(file,paste(substr(Genes_of_interest, 1,1),collapse=""),".pdf"), 7,7)
  for(Gene in Genes_of_interest){
    df_GOI<-df_res[df_res$GeneSymbol %in% Gene,]
    
    par(mar=c(4,4,0,0))
    par(fig=c(0.1,0.7,0.1,0.7))
    plot(df_res$logBaseMeanA, df_res$log2FoldChange, type="p", col=df_res$col, bg=df_res$col, cex=1, pch=df_res$pch, xlab=xlab, 
         ylab=ylab, cex.lab=1, cex.axis=1, xlim=c(xmin, xmax), ylim=c(ymin, ymax), xaxp=c(xmin,xmin+xticks*floor((xmax-xmin)/xticks),
              floor((xmax-xmin)/xticks)), yaxp=c(ymin,ymin+yticks*floor((ymax-ymin)/yticks),floor((ymax-ymin)/yticks)))
    if (nrow(df_GOI)>=1){
      points(df_GOI$logBaseMeanA, df_GOI$log2FoldChange, type="p", col=ColH, bg=ColH, cex=1.5, pch=df_GOI$pch)
    }
    legend("bottomleft",legend=c( if (nchar(Gene)>1 && !is.na(Gene)) {Gene}, "All guides", "Essentials","Non-essentials", 
          'Significantly enriched', 'Significantly depleted'), cex=0.8, pch=c(if (nchar(Gene)>1 && !is.na(Gene)) {16},16,16,16, 24,25), 
          col=c(if (nchar(Gene)>1 && !is.na(Gene)) {ColH}, ColAll,ColP,ColN, "black", "black"))
    abline(median(df_res$log2FoldChange, na.rm=TRUE),0, col=1, lty=3, untf=TRUE)
    
    # Density plot fold change
    denY_tot<- density(df_res$log2FoldChange, from=ymin, to=ymax, na.rm=TRUE)
    denY_tot$y[1]<- 0
    denY_tot$y[length(denY_tot$y)]<- 0
    denY_PC<- density(df_PC$log2FoldChange, from=ymin, to=ymax, na.rm=TRUE)
    denY_PC$y[1]<- 0
    denY_PC$y[length(denY_PC$y)]<- 0
    denY_NC<- density(df_NC$log2FoldChange, from=ymin, to=ymax, na.rm=TRUE)
    denY_NC$y[1]<- 0
    denY_NC$y[length(denY_NC$y)]<- 0
    denYMax<- max(c(denY_tot$y, denY_PC$y, denY_NC$y))
    
    par(mar=c(4,0,0,4))
    par(fig=c(0.7,0.9,0.1,0.7),new=TRUE)
    plot(denY_tot$y, denY_tot$x, ylim=c(ymin,ymax), xlim=(c(0,denYMax)), type='l', axes=FALSE, col=ColAll, xlab="", 
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
    if (nrow(df_GOI)>1){
      par(new=TRUE)
      denY_GOI<- density(df_GOI$log2FoldChange, from=ymin, to=ymax, na.rm=TRUE)
      denY_GOI$y[1]<- 0
      denY_GOI$y[length(denY_GOI$y)]<- 0
      lines(denY_GOI$y, denY_GOI$x, col=ColH, lwd=2)
      rgb.val<- col2rgb(ColH)
      polygon(denY_GOI$y, denY_GOI$x, col=rgb(rgb.val[1]/255,rgb.val[2]/255,rgb.val[3]/255,alpha=0.3), lwd=0.1)
    }
    
    # Density plot read count
    denX_tot<- density(df_res$logBaseMeanA, from=xmin, to=xmax, na.rm=TRUE)
    denX_tot$y[1]<- 0
    denX_tot$y[length(denX_tot$y)]<- 0
    denX_PC<- density(df_PC$logBaseMeanA, from=xmin, to=xmax, na.rm=TRUE)
    denX_PC$y[1]<- 0
    denX_PC$y[length(denX_PC$y)]<- 0
    denX_NC<- density(df_NC$logBaseMeanA, from=xmin, to=xmax, na.rm=TRUE)
    denX_NC$y[1]<- 0
    denX_NC$y[length(denX_NC$y)]<- 0
    denXMax<- max(c(denX_tot$y, denX_PC$y, denX_NC$y))
    
    par(mar=c(0,4,4,0))
    par(fig=c(0.1,0.7,0.7,0.9),new=TRUE)
    plot(denX_tot$x, denX_tot$y, main = maintitle, cex.main=1.5, xlim=c(xmin, xmax), ylim=c(0,denXMax), type='l', axes=FALSE, col=ColAll, xlab="", 
         ylab="", lwd=2)
    lines(denX_PC, col=ColP, lwd=2)
    lines(denX_NC, col=ColN, lwd=2)
    rgb.val<- col2rgb(ColAll)
    polygon(denX_tot, col=rgb(rgb.val[1]/255,rgb.val[2]/255,rgb.val[3]/255,alpha=0.3), lwd=0.1)
    rgb.val<- col2rgb(ColP)
    polygon(denX_PC, col=rgb(rgb.val[1]/255,rgb.val[2]/255,rgb.val[3]/255,alpha=0.3), lwd=0.1)
    rgb.val<- col2rgb(ColN)
    polygon(denX_NC, col=rgb(rgb.val[1]/255,rgb.val[2]/255,rgb.val[3]/255,alpha=0.3), lwd=0.1)
    if (nrow(df_GOI)>1){
      par(new=TRUE)
      denX_GOI<- density(df_GOI$logBaseMeanA, from=xmin, to=xmax, na.rm=TRUE)
      denX_GOI$y[1]<- 0
      denX_GOI$y[length(denX_GOI$y)]<- 0
      lines(denX_GOI, col=1, lwd=2)
      rgb.val<- col2rgb(ColH)
      polygon(denX_GOI, col=rgb(rgb.val[1]/255,rgb.val[2]/255,rgb.val[3]/255,alpha=0.3), lwd=0.1)
    }
  }
  dev.off()
  n<- n+1
}