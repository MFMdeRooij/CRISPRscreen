# Check the quality of lethality screens, of which only CRISPR scores are published.
# This script compares the distribution of the essential and non-essential genes.
# Create a CSV file with the gene symbols in the first column and the CRISPR scores in the following columns,and with 
# the cell line names in the first row.Put Controls.csv in the same folder.
# rZ': robust Z factor, NP50: CRISPR score in which the essential and non-essential curves cross, so at this point
# there is 50% chance to be essential. F1: harmonic mean of precision and recall: https://en.wikipedia.org/wiki/F-score
# Probability of essentiality: chance to be like an essential (ratio essential/nonessential curve at CRISPR score)
# Author: M.F.M. de Rooij, PhD, Amsterdam UMC, Spaargaren Lab, 2023, info: m.f.derooij@amsterdamumc.nl
#########################################################################################
#                                   VARIABLES
# Folder
Workdirectory<- "H:/BioWin/CRISPRess/"

# Data file (layout as described above)
Filename<- "CRISPRscores.csv"      

#########################################################################################
# Data
setwd(Workdirectory)
df_data <- read.csv(file=Filename, sep=",", header=TRUE, stringsAsFactors = FALSE)
colnames(df_data)[1]<-"GeneSymbol"
NumCells<-ncol(df_data)-1

# Guide ID
df_Controls <- read.csv("CRISPRscreenControls.csv", sep=',', header=TRUE, stringsAsFactors = FALSE)

df_data$Type<- "x"
df_data[df_data$GeneSymbol %in% df_Controls$Essential[nchar(df_Controls$Essential)>0],"Type"]<- "p"
df_data[df_data$GeneSymbol %in% df_Controls$Non_Essential[nchar(df_Controls$Non_Essential)>0],"Type"]<- "n"
df_PC<-df_data[df_data$Type=="p",]
df_NC<-df_data[df_data$Type=="n",]

df_int<- data.frame(cell=numeric(0), int=numeric(0))

pdf(paste0(Filename,".pdf"), width=10, height=10)
par(mfrow=c(2,2))
for (n in 1:NumCells) { 
  xMin <- min(unlist(df_data[n+1])) 
  xMax <- max(unlist(df_data[n+1]))
  den_cell1To<- density(unlist(df_data[n+1]), from=xMin, to=xMax)
  den_cell1PC<- density(unlist(df_PC[n+1]), from=xMin, to=xMax)
  den_cell1NC<- density(unlist(df_NC[n+1]), from=xMin, to=xMax)
  yMax<- max(c(den_cell1To$y,den_cell1PC$y, den_cell1NC$y))

  # Robust Z prime (median and median absolute deviation)
  Mp<- median(unlist(df_PC[n+1]))
  Mn<- median(unlist(df_NC[n+1]))
  SDp<-mad(unlist(df_PC[n+1]))
  SDn<-mad(unlist(df_NC[n+1]))
  Zprime <- 1-((3*(SDp+SDn))/(abs(Mp-Mn))) 
  
  # Intersection positive and negative controls
  poi <- which(diff(den_cell1PC$y > den_cell1NC$y) != 0) 
  intersection<- den_cell1PC$x[poi][den_cell1PC$x[poi]<0][which.min(abs(den_cell1PC$x[poi][den_cell1PC$x[poi]<0]))]
  df_int_temp<- data.frame(cell=n, int=intersection)
  df_int<- rbind(df_int, df_int_temp)

  # F measure (with cutoff the intersection), harmonic mean precision and recall
  truePos<- length(df_PC[n+1][df_PC[n+1]<intersection])
  posPredict<- length(df_PC[n+1][df_PC[n+1]<intersection])+length(df_NC[n+1][df_NC[n+1]<intersection])*nrow(df_PC)/nrow(df_NC)
  posControls<- nrow(df_PC)
  precision<- truePos/posPredict
  recall<- truePos/posControls
  F1<- 2/(1/precision + 1/recall)

  plot(den_cell1To$x, den_cell1To$y, xlim=(c(xMin,xMax)), ylim=c(0,yMax), type='l', col=3, xlab=~"CRISPR score", 
       ylab="Density", lwd=2, main=colnames(df_data[n+1]), sub = paste0("(rZ': ", format(round(Zprime,2), nsmall=2), 
                                                                        " NP50: ", format(round(intersection,2), nsmall=2), " F1: ", format(round(F1,2),nsmall = 2), " (Pre: ", format(round(precision,2),nsmall = 2), 
                                                                        " Rec: ", format(round(recall,2),nsmall = 2), "))"))
  par(new=TRUE)
  plot(den_cell1PC$x, den_cell1PC$y, xlim=(c(xMin,xMax)), ylim=c(0,yMax), type='l', axes=FALSE, col=2, xlab="", 
       ylab="", lwd=2)
  par(new=TRUE)
  plot(den_cell1NC$x, den_cell1NC$y, xlim=(c(xMin,xMax)), ylim=c(0,yMax), type='l', axes=FALSE, col=4, xlab="", 
       ylab="", lwd=2)
  legend(xMin,yMax,legend=c("Total", "Essential", "Non-Essential"), cex=1, pch=c(20,20,20), col=c(3,2,4))
  
  
  ddp<-approxfun(den_cell1PC$x, den_cell1PC$y)
  ddn<-approxfun(den_cell1NC$x, den_cell1NC$y)
  funChance<- function(x){
    ddp(x)/(ddp(x)+ddn(x)+0.00001) 
  }
  chance <- apply(df_data[n+1], MARGIN=2, FUN = funChance)
  plot(unlist(df_data[n+1]),chance, cex=0.1, xlim=c(xMin,xMax), ylim=c(0,1), xlab="CRISPR score", ylab="Probability of Essentiality", main=colnames(df_data[n+1]))
  
  df_fit<- data.frame(x=df_data[n+1], y=chance)
  df_fit<- df_fit[order(df_fit[,1]),]
  df_fit<- df_fit[complete.cases(df_fit),]
  x<- unlist(df_fit[1])
  y<- unlist(df_fit[2])
  m.s <- nls(y ~ (1/(1 + exp(-c * (x - d)))), start = list(c = -1, d = round(median(x))), trace = TRUE)
  lines(x, fitted(m.s), lwd=2, col = 2)
  
  probability<- round(predict(m.s, list(x=unlist(df_data[n+1]))),3)
  df_data<-cbind(df_data,probability)
  colnames(df_data)[ncol(df_data)]<- paste0(colnames(df_data[n+1]), "_ProbabilityOfEssentiality")
}  
dev.off()

# Export essentiality scores
df_data2<- df_data[,c("GeneSymbol", "Type", colnames(df_data)[grep("Probability", colnames(df_data))])]
write.csv(df_data2, paste0(Filename, "_EssScores.csv"), row.names = F)
