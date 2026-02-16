# If one or both drug are not toxic you cannot use Chou-Talalay, but you can use Bliss. However, due to that most drug-concentration curves are S-curves, 
# you can also find Biss-synergy with combinations of the same drug.
# Combine 2 drugs in different concentrations in a matrix format, normalize the control (no drugs) to 100%, and save it to a CSV file. 
# Omit the drug names and concentration from the CSV file, but fill them in the variable here below. The number of concentrations of drug 1 and 2 should 
# match with the number of columns and row in the CSV file. 
# Author: M.F.M. de Rooij PhD, Amsterdam UMC, Spaargaren Lab, 2023, info: m.f.derooij@amsterdamumc.nl
#################################################################################
# Packages:
#install.packages("ggplot2")
#install.packages("scales")
#install.packages("cowplot")
#install.packages("reshape2")
#install.packages("readxl")
library("ggplot2")
library("scales")
library("cowplot")
library("reshape2")
library("readxl")
#################################################################################
#                                   VARIABLES

# Workdirectory (put in this folder the excel files (.xlsx or csv) from the cell lines)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
#setwd("H:/BioWin/Bliss/")

# Format data file 0: .csv, 1: .xlsx (Data matrix in upperleft corner)
format<- 1

# Find files automatically (Name the files after the cells)
if (format==0){
  cells<- unlist(strsplit(list.files(pattern=".csv"), ".csv"))
} else{
  if (format==1) {
    cells<- unlist(strsplit(list.files(pattern=".xlsx"), ".xlsx"))
  } 
}
# Or manual (and this will be the order):
#cells <- c("Cell1", "Cell2")

# Drug names and concentrations:
# Drug1 = Columns
xname <- expression("Drug1 ("*mu*"M)")
xlab <- c(0,0.5,1,2,4)

# Drug2 = Rows
yname <- expression("Drug2 ("*mu*"M)")
ylab <- c(0,1,2,4,8)

# Viability or cell number
ViabilityOrCellNumber<- "Relative Viability"
#################################################################################
cellNumber<- 1
for (cell in cells){
  # Read files
  if (format==0){
    AB <- read.csv(paste0(cell,".csv"), sep=",", header=FALSE)/100
  } else {
    if (format==1){
      AB <- read_excel(paste0(cell,".xlsx"), col_names=F)/100
    }
  }
  
  # Normalize control to 100%
  AB<- AB/AB[1,1]
  
  # Caluculate expected, delta Bliss, and relative Bliss values
  expectedMat<- as.matrix(AB)[,1] %*% t(as.matrix(AB)[1,])
  dblissMat<- expectedMat-as.matrix(AB)
  rblissMat<- as.matrix(AB)/expectedMat
  
  matInhX<- as.data.frame(t(AB/AB[,1]))
  matInhX$c<- 1:nrow(matInhX)
  colnames(matInhX)<-c(ylab,"c")
  
  matInhY<- as.data.frame(t(t(AB)/t(AB)[,1]))
  matInhY$c<- 1:nrow(matInhY)
  colnames(matInhY)<-c(xlab,"c")
  
  rmatInhX <- as.data.frame(rblissMat)
  rmatInhX$c<- 1:nrow(rmatInhX)
  colnames(rmatInhX)<-c(xlab,"c")
  
  rmatInhY <- as.data.frame(t(rblissMat))
  rmatInhY$c<- 1:nrow(rmatInhY)
  colnames(rmatInhY)<-c(ylab,"c")   
  
  # Collect the data
  Combo <- data.frame(matrix(nrow=length(xlab)*length(ylab),ncol=6))
  colnames(Combo) <- c("c", "r", "pgrowth", "pr.inh", "dbliss", "rbliss")  
  # Fill the columns with data
  k=1
  for (j in 1:length(xlab)){
    for (i in 1:length(ylab)){
      Combo$c[k] <- j
      Combo$r[k] <- length(ylab)+1-i
      Combo$pgrowth[k] <- AB[i,j]
      Combo$pr.inh[k] <- expectedMat[i,j]
      Combo$dbliss[k] <- dblissMat[i,j]
      Combo$rbliss[k] <- rblissMat[i,j]
      k<-k+1
    }
  }
  
  # Make graphs
  g1<-ggplot(Combo, aes(x=c, y=r, label=(round(pgrowth*100, digits=0)))) +
    geom_tile(aes(fill=pgrowth, height=1, width=1)) +
    scale_fill_gradient2(label=percent,
                         name=ViabilityOrCellNumber,
                         midpoint=0.5,
                         limits=c(0,1),
                         low="red", mid="darkred", high="midnightblue", na.value="midnightblue") +
    geom_text(size=9, colour="white") +
    scale_x_discrete(xname, limits=as.character(xlab), position = "top" ) +
    scale_y_discrete(yname, limits=rev(as.character(ylab))) +
    ggtitle(cell) +
    theme(text=element_text(size=25),
          axis.text=element_text(size=25, colour="black"),
          plot.title=element_text(size=30, face="bold", hjust = 0.5),
          plot.margin = unit(c(1, 2, 0, 1), "cm"),
          legend.title = element_text(vjust = 2),
          legend.position = c(1.03, 0.5),
          legend.justification = c(0, 0.5)
    ) + guides(fill = guide_colourbar(label.hjust = 1))
  
  g2<-ggplot(Combo, aes(x=c, y=r, label=(round(dbliss*100, digits=0)))) +
    geom_tile(aes(fill=dbliss, height=1, width=1)) +
    scale_fill_gradient2(label=percent,
                         name= expression(Delta*"Bliss"),
                         midpoint=0,
                         limits=c(-0.5,0.5),
                         low="red4", mid="lightgray", high="chartreuse", na.value="chartreuse") +
    geom_text(size=9, colour="black") +
    scale_x_discrete(xname, limits=as.character(xlab), position = "bottom" ) +
    scale_y_discrete(yname, limits=rev(as.character(ylab))) +
    theme(text=element_text(size=25),
          axis.text=element_text(size=25, colour="black"),
          plot.title=element_text(size=30, face="bold", hjust = 0.5),
          plot.margin = unit(c(0, 2, 2, 1), "cm"),
          legend.title = element_text(vjust = 2),
          legend.position = c(1.03, 0.5),
          legend.justification = c(0, 0.5)
    ) + guides(fill = guide_colourbar(label.hjust = 1))
  
  g3<-ggplot(Combo, aes(x=c, y=r, label=(round(rbliss*100, digits=0)))) +
    geom_tile(aes(fill=rbliss, height=1, width=1)) +
    scale_fill_gradient2(label=percent,
                         name=ViabilityOrCellNumber,
                         midpoint=1,
                         limits=c(0,2),
                         low="red4", mid="lightgray", high="green3", na.value="green3") +
    geom_text(size=9, colour="black") +
    scale_x_discrete(xname, limits=as.character(xlab), position = "bottom" ) +
    scale_y_discrete(yname, limits=rev(as.character(ylab))) +
    ggtitle(cell) +
    theme(text=element_text(size=25),
          axis.text=element_text(size=25, colour="black"),
          plot.title=element_text(size=30, face="bold", hjust = 0.5),
          plot.margin = unit(c(0, 2, 1, 1), "cm"),
          legend.title = element_text(vjust = 2),
          legend.position = c(1.03, 0.5),
          legend.justification = c(0, 0.5)
    ) + guides(fill = guide_colourbar(label.hjust = 1))
  
  
  fc <- colorRampPalette(c("lightblue", "darkblue"))
  coly<-c("red", fc(length(ylab)-1))
  colx<-c("red", fc(length(xlab)-1))
  
  l1<- ggplot(melt(matInhX,id = "c"), aes(x=c, y=value*100, group=variable, colour=variable)) +
    geom_line(lwd=1.5) + ylab(paste0(ViabilityOrCellNumber," (%)")) +
    scale_x_discrete(xname, limits=as.character(xlab)) +
    theme(text=element_text(size=25),
          axis.text=element_text(size=25, colour="black"),
          plot.title=element_text(size=30, face="bold", hjust = 0.5),
          plot.margin = unit(c(0, 1, 2, 1), "cm"),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          legend.justification = c(0, 0.5)
    ) + scale_color_manual(values=coly, name = yname) + expand_limits(y = 0)
  
  l2<- ggplot(melt(matInhY,id = "c"), aes(x=c, y=value*100, group=variable, colour=variable)) +
    geom_line(lwd=1.5) + ylab(paste0(ViabilityOrCellNumber," (%)")) +
    scale_x_discrete(yname, limits=as.character(ylab)) +
    theme(text=element_text(size=25),
          axis.text=element_text(size=25, colour="black"),
          plot.title=element_text(size=30, face="bold", hjust = 0.5),
          plot.margin = unit(c(0, 1, 2, 1), "cm"),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          legend.justification = c(0, 0.5)
    ) + scale_color_manual(values=colx, name = xname) + expand_limits(y = 0)
  
  l3<- ggplot(melt(rmatInhY,id = "c"), aes(x=c, y=value*100, group=variable, colour=variable)) +
    geom_line(lwd=1.5) + ylab(paste0(ViabilityOrCellNumber," (%)")) +
    scale_x_discrete(xname, limits=as.character(xlab)) +
    theme(text=element_text(size=25),
          axis.text=element_text(size=25, colour="black"),
          plot.title=element_text(size=30, face="bold", hjust = 0.5),
          plot.margin = unit(c(0, 1, 2, 1), "cm"),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          legend.justification = c(0, 0.5)
    ) + scale_color_manual(values=coly, name = yname) + expand_limits(y = 0)
  
  l4<- ggplot(melt(rmatInhX,id = "c"), aes(x=c, y=value*100, group=variable, colour=variable)) +
    geom_line(lwd=1.5) + ylab(paste0(ViabilityOrCellNumber," (%)")) +
    scale_x_discrete(yname, limits=as.character(ylab)) +
    theme(text=element_text(size=25),
          axis.text=element_text(size=25, colour="black"),
          plot.title=element_text(size=30, face="bold", hjust = 0.5),
          plot.margin = unit(c(0, 1, 2, 1), "cm"),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          legend.justification = c(0, 0.5)
    ) + scale_color_manual(values=colx, name = xname) + expand_limits(y = 0)
  
  # Store all graphs in a variable
  assign(paste0("plotgraph1", cellNumber, cell), g1)
  assign(paste0("plotgraph2", cellNumber, cell), g2)
  assign(paste0("plotgraph3", cellNumber, cell), g3)
  assign(paste0("plotgraph4", cellNumber, cell), l1)
  assign(paste0("plotgraph5", cellNumber, cell), l2)
  assign(paste0("plotgraph6", cellNumber, cell), l3)
  assign(paste0("plotgraph7", cellNumber, cell), l4)
  
  cellNumber<- cellNumber+1
}

# List all graphs
plotlist1 <- lapply(ls(pattern="plotgraph1"), function(x) {get(x)})
plotlist2 <- lapply(ls(pattern="plotgraph2"), function(x) {get(x)})
plotlist3 <- lapply(ls(pattern="plotgraph3"), function(x) {get(x)})
plotlist4 <- lapply(ls(pattern="plotgraph4"), function(x) {get(x)})
plotlist5 <- lapply(ls(pattern="plotgraph5"), function(x) {get(x)})
plotlist6 <- lapply(ls(pattern="plotgraph6"), function(x) {get(x)})
plotlist7 <- lapply(ls(pattern="plotgraph7"), function(x) {get(x)})
plotlist<- append(plotlist1,  plotlist2)
plotlist<- append(plotlist,  plotlist4)
plotlist<- append(plotlist,  plotlist5)
plotlist<- append(plotlist,  plotlist3)
plotlist<- append(plotlist,  plotlist6)
plotlist<- append(plotlist,  plotlist7)

# Plot all graphs in PDF file
pdf(file=paste0("BlissSynergyR.pdf"), width=10*length(cells), height=48)
  print(plot_grid(plotlist=plotlist, nrow=7, ncol=length(cells), align="v"))
dev.off()

# Remove all variables (To prevent plotting older graphs)
rm(list = ls())