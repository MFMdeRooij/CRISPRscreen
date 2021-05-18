# Use a constant drugs ratio dependent on its IC50's
# For 2 drugs the best is to use dilution factors of 2
# Normalize the viability values/number of cells between 0 (0%) and 1 (100%), adjust te settings, and run the script in R studio
# Author: M.F.M. de Rooij PhD, Amsterdam UMC, Spaargaren Lab, 2018, info: m.f.derooij@amsterdamumc.nl
##################################################################################################################################
#                                                         SETTINGS

# Put NPI normalized data in SynergyTestData.csv template
Workdirectory<- "H:/BioWin/Synergy"
Filename<- "SynergyTestData.csv"

# Names of the drugs
DrugA<- "Drug A"
DrugB<- "Drug B"

# Name output file
output <- "SynergyCubic.pdf"

##################################################################################################################################

setwd(Workdirectory)
df_raw <- read.csv(file=Filename, sep=",", header=TRUE, stringsAsFactors = FALSE)

# Data table
df_data<-data.frame(RelDoseLog2=df_raw$RelDoseLog2, stringsAsFactors = FALSE)
df_data$DrugA<- rowMeans(df_raw[2:4]) 
df_data$DrugB<- rowMeans(df_raw[5:7])
df_data$Combination<- rowMeans(df_raw[8:10])
df_data$SDDrugA <- apply(df_raw[2:4], 1, sd)
df_data$SDDrugB <- apply(df_raw[5:7], 1, sd)
df_data$SDCombination <- apply(df_raw[8:10], 1, sd)

# Synergy curves
dosepoints<- seq(min(df_data$RelDoseLog2), max(df_data$RelDoseLog2), (max(df_data$RelDoseLog2)-min(df_data$RelDoseLog2))/100)

df_dosepoints<- data.frame(x1=dosepoints, x2=dosepoints^2, x3=dosepoints^3)
# Creating Graphs
pdf(output, width=10, height=10)
par(fig=c(0,0.55,0.39,0.94))
plot(df_data$RelDoseLog2, df_data$DrugA, col=2, cex=0.5, ylim=c(0,1), xlab="", ylab="Viability", xaxt='n')
arrows(df_data$RelDoseLog2, df_data$DrugA-df_data$SDDrugA, df_data$RelDoseLog2, df_data$DrugA+df_data$SDDrugA, 
                                                                            length=0.01, angle=90, code=3, col=2)


x <- df_data$RelDoseLog2
y <- df_data$DrugA
xy<- data.frame(y=y, x1=df_data$RelDoseLog2, x2=df_data$RelDoseLog2^2, x3=df_data$RelDoseLog2^3)

m.sA <- lm(y~x1+x2+x3, data=xy)
dosepointDrugAFit <-predict(m.sA, df_dosepoints)
lines(dosepoints, dosepointDrugAFit, lty = 1, lwd = 1, col = 2)
df_sigmoid<- data.frame(coeffiecient=c("a","b","c","d"), stringsAsFactors = FALSE)
df_sigmoid$logDoseDrugA<- coef(m.sA)
par(new=TRUE)
plot(df_data$RelDoseLog2, df_data$DrugB, col=3, cex=0.5, ylim=c(0,1), xlab="", ylab="", axes=FALSE, xaxt='n', 
                                                                                                          yaxt='n')
arrows(df_data$RelDoseLog2, df_data$DrugB-df_data$SDDrugB, df_data$RelDoseLog2, df_data$DrugB+df_data$SDDrugB, 
                                                                              length=0.01, angle=90, code=3, col=3)
y <- df_data$DrugB
xy<- data.frame(y=y, x1=df_data$RelDoseLog2, x2=df_data$RelDoseLog2^2, x3=df_data$RelDoseLog2^3)

m.sB <- lm(y~x1+x2+x3, data=xy)

dosepointDrugBFit <-predict(m.sB, df_dosepoints)
lines(dosepoints, dosepointDrugBFit, lty = 1, lwd = 1, col = 3)
df_sigmoid$logDoseDrugB<- coef(m.sB)
par(new=TRUE)
plot(df_data$RelDoseLog2, df_data$Combination, col=1, cex=0.5, ylim=c(0,1), xlab="", ylab="", axes=FALSE, 
                                                                                              xaxt='n', yaxt='n')
arrows(df_data$RelDoseLog2, df_data$Combination-df_data$SDCombination, df_data$RelDoseLog2, 
                      df_data$Combination+df_data$SDCombination, length=0.01, angle=90, code=3, col=1)
y <- df_data$Combination
xy<- data.frame(y=y, x1=df_data$RelDoseLog2, x2=df_data$RelDoseLog2^2, x3=df_data$RelDoseLog2^3)

m.sC <- lm(y~x1+x2+x3, data=xy)
dosepointCombinationFit <-predict(m.sC, df_dosepoints)
lines(dosepoints, dosepointCombinationFit, lty = 1, lwd = 1, col = 1)
df_sigmoid$logDoseCombination<- coef(m.sC)
legend(1,0.5,legend=c(DrugA, DrugB, "Combination"), cex=1, pch=1, col=c(2,3,1))

df_sigmoid[is.na(df_sigmoid)]<-0.000001
#install.packages("RConics")
library(RConics)
blissDosepointsIC50s <- c(df_sigmoid$logDoseCombination[4],df_sigmoid$logDoseCombination[3], 
       df_sigmoid$logDoseCombination[2], df_sigmoid$logDoseCombination[1]-0.5)


blissDosepointsIC50<-cubic(blissDosepointsIC50s)
blissDosepointsIC50<-blissDosepointsIC50[which.min(abs(blissDosepointsIC50 - median(df_data$RelDoseLog2)))]

blissDosepointsIC50 <- as.numeric(blissDosepointsIC50)

df_blissDosepointsIC50<- data.frame(x1=blissDosepointsIC50, x2=blissDosepointsIC50^2, x3=blissDosepointsIC50^3)


# Chou-Talalay
df_data$CombinationFit<-predict(m.sC, xy[-1])
df_data$DrugAFit<-predict(m.sA, xy[-1])
df_data$DrugBFit<-predict(m.sB, xy[-1])
  

df_data$logDoseDrugAFit<-0
for (i in 1:nrow(df_data)){
xvals<- cubic(c(df_sigmoid$logDoseDrugA[4],df_sigmoid$logDoseDrugA[3], 
  df_sigmoid$logDoseDrugA[2], df_sigmoid$logDoseDrugA[1]-df_data$CombinationFit[i]))
df_data$logDoseDrugAFit[i]<-xvals[which.min(abs(xvals - median(df_data$RelDoseLog2)))]
print(xvals)
}
df_data$logDoseDrugBFit<-0
for (i in 1:nrow(df_data)){
  xvals<- cubic(c(df_sigmoid$logDoseDrugB[4],df_sigmoid$logDoseDrugB[3], 
                  df_sigmoid$logDoseDrugB[2], df_sigmoid$logDoseDrugB[1]-df_data$CombinationFit[i]))
  df_data$logDoseDrugBFit[i]<-xvals[which.min(abs(xvals - median(df_data$RelDoseLog2)))]
  print(xvals)
}
df_data$CI <- (2^(df_data$RelDoseLog2-df_data$logDoseDrugAFit))+(2^(df_data$RelDoseLog2-df_data$logDoseDrugBFit))
par(new=TRUE)
par(fig=c(0.45,1,0.39,0.94))
plot(df_data$CI, df_data$CombinationFit, col=1, cex=0.5, ylim=c(0,1), xlim=c(0,2), xlab="", ylab="",xaxt='n', 
                                                                                                            yaxt='n')
dosepointlogDoseDrugAFit <-rep(0,length(dosepointCombinationFit)) 
for (i in 1:length(dosepointCombinationFit)){
  xvals<- cubic(c(df_sigmoid$logDoseDrugA[4],df_sigmoid$logDoseDrugA[3], 
                  df_sigmoid$logDoseDrugA[2], df_sigmoid$logDoseDrugA[1]-dosepointCombinationFit[i]))
  dosepointlogDoseDrugAFit[i]<-xvals[which.min(abs(xvals - median(df_data$RelDoseLog2)))]
  print(xvals)
}  

dosepointlogDoseDrugBFit <-rep(0,length(dosepointCombinationFit)) 
for (i in 1:length(dosepointCombinationFit)){
  xvals<- cubic(c(df_sigmoid$logDoseDrugB[4],df_sigmoid$logDoseDrugB[3], 
                  df_sigmoid$logDoseDrugB[2], df_sigmoid$logDoseDrugB[1]-dosepointCombinationFit[i]))
  dosepointlogDoseDrugBFit[i]<-xvals[which.min(abs(xvals - median(df_data$RelDoseLog2)))]
  print(xvals)
}   
dosepointCI <- (2^(dosepoints-dosepointlogDoseDrugAFit))+(2^(dosepoints-dosepointlogDoseDrugBFit))
lines(dosepointCI, dosepointCombinationFit)
abline(v=1)

# Calculation of CI50
Combination50CIs <- cubic(c(df_sigmoid$logDoseCombination[4],df_sigmoid$logDoseCombination[3], 
                           df_sigmoid$logDoseCombination[2], df_sigmoid$logDoseCombination[1]-0.5))
Combination50CI<-Combination50CIs[which.min(abs(Combination50CIs - median(df_data$RelDoseLog2)))]

DrugA50CIs <- cubic(c(df_sigmoid$logDoseDrugA[4],df_sigmoid$logDoseDrugA[3], 
                     df_sigmoid$logDoseDrugA[2], df_sigmoid$logDoseDrugA[1]-0.5))
DrugA50CI<-DrugA50CIs[which.min(abs(DrugA50CIs - median(df_data$RelDoseLog2)))]

DrugB50CIs <- cubic(c(df_sigmoid$logDoseDrugB[4],df_sigmoid$logDoseDrugB[3], 
                     df_sigmoid$logDoseDrugB[2], df_sigmoid$logDoseDrugB[1]-0.5))
DrugB50CI<-DrugB50CIs[which.min(abs(DrugB50CIs - median(df_data$RelDoseLog2)))]

CI50<- as.numeric(round((2^(Combination50CI-DrugA50CI))+(2^(Combination50CI-DrugB50CI)),2))
text(1.9,0.5, bquote(~""*CI[50] == .(CI50)*""), cex=0.9, srt=270) 

# Bliss
df_data$BlissExpected<- df_data$DrugAFit*df_data$DrugBFit
df_data$BlissScore <- df_data$CombinationFit/df_data$BlissExpected
par(new=TRUE)
par(fig=c(0,0.55,0,0.55))
plot(df_data$RelDoseLog2, df_data$BlissScore, col=1, cex=0.5, ylim=c(0,2), xlab="Relative dose (log2)", 
     ylab="Relative Bliss score")
blissDosepoints<-dosepointCombinationFit/(dosepointDrugAFit*dosepointDrugBFit)
lines(dosepoints,blissDosepoints)
abline(h=1)

# Calculation of BS50
CombinationBI <- predict(m.sC, df_blissDosepointsIC50)



DrugABI<- predict(m.sA, df_blissDosepointsIC50)
DrugBBI<- predict(m.sB, df_blissDosepointsIC50)
BS50<- round(0.5/(DrugABI*DrugBBI),2)
text(median(df_data$RelDoseLog2),1.9, bquote(~""*BS[50] == .(BS50)*""), cex=0.9) 

# Bliss vs CI
par(new=TRUE)
par(fig=c(0.45,1,0,0.55))
plot(df_data$CI, df_data$BlissScore, type="p", col=1, cex=0.5, ylim=c(0,2), xlim=c(0,2), 
                                                    xlab="Combination index", ylab="", yaxt='n')
lines(dosepointCI, blissDosepoints)
abline(v=1)
abline(h=1)

# Main title
par(new=TRUE)
par(fig=c(0,0.95,0,0.95))
plot(0:10, type = "n", xaxt="n", yaxt="n", bty="n", xlab = "", ylab = "", main="Chou-Talalay vs Bliss")

dev.off()
