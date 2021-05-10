# Use a constant drugs ratio dependent on its IC50's
# For 2 drugs the best is to use dilution factors of 2

# Author: M.F.M. de Rooij PhD, Amsterdam UMC, Spaargaren Lab, 2019, info: m.f.derooij@amsterdamumc.nl



# Put NPI normalized data in SynergyTestData.csv template
Workdirectory<- "H:/BioWin/Synergy"
Filename<- "SynergyTestData.csv"

# Names of the drugs
DrugA<- "Drug A"
DrugB<- "Drug B"

# Name output file
output <- "Synergy4comp.pdf"

####################################################################################

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

# Creating Graphs
pdf(output, width=10, height=10)
par(fig=c(0,0.55,0.39,0.94))
plot(df_data$RelDoseLog2, df_data$DrugA, col=2, cex=0.5, ylim=c(0,1), xlab="", ylab="Viability", xaxt='n')
arrows(df_data$RelDoseLog2, df_data$DrugA-df_data$SDDrugA, df_data$RelDoseLog2, df_data$DrugA+df_data$SDDrugA, 
                                                                            length=0.01, angle=90, code=3, col=2)
x <- df_data$RelDoseLog2
y <- df_data$DrugA
m.sA <- nls(y ~ a + ((b - a)/(1 + exp(-c * (x - d)))), start = list(a = min(y), b = max(y), c = 1, 
                                                                   d = round(median(x))), trace = TRUE)
dosepointDrugAFit <-predict(m.sA, list(x=dosepoints))
lines(dosepoints, dosepointDrugAFit, lty = 1, lwd = 1, col = 2)
df_sigmoid<- data.frame(coeffiecient=c("a","b","c","d"), stringsAsFactors = FALSE)
df_sigmoid$logDoseDrugA<- coef(m.sA)
par(new=TRUE)
plot(df_data$RelDoseLog2, df_data$DrugB, col=3, cex=0.5, ylim=c(0,1), xlab="", ylab="", axes=FALSE, xaxt='n', 
                                                                                                          yaxt='n')
arrows(df_data$RelDoseLog2, df_data$DrugB-df_data$SDDrugB, df_data$RelDoseLog2, df_data$DrugB+df_data$SDDrugB, 
                                                                              length=0.01, angle=90, code=3, col=3)
y <- df_data$DrugB
m.sB <- nls(y ~ a + ((b - a)/(1 + exp(-c * (x - d)))), start = list(a = min(y), b = max(y), c = 1, 
                                                                   d = round(median(x))), trace = TRUE)
dosepointDrugBFit <-predict(m.sB, list(x=dosepoints))

lines(dosepoints, dosepointDrugBFit, lty = 1, lwd = 1, col = 3)
df_sigmoid$logDoseDrugB<- coef(m.sB)
par(new=TRUE)
plot(df_data$RelDoseLog2, df_data$Combination, col=1, cex=0.5, ylim=c(0,1), xlab="", ylab="", axes=FALSE, 
                                                                                              xaxt='n', yaxt='n')
arrows(df_data$RelDoseLog2, df_data$Combination-df_data$SDCombination, df_data$RelDoseLog2, 
                      df_data$Combination+df_data$SDCombination, length=0.01, angle=90, code=3, col=1)
y <- df_data$Combination
m.sC <- nls(y ~ a + ((b - a)/(1 + exp(-c * (x - d)))), start = list(a = min(y), b = max(y), c = 1, 
                                                                   d = round(median(x))), trace = TRUE)
dosepointCombinationFit <-predict(m.sC, list(x=dosepoints))
lines(dosepoints, dosepointCombinationFit, lty = 1, lwd = 1, col = 1)
df_sigmoid$logDoseCombination<- coef(m.sC)
legend(1,0.25,legend=c(DrugA, DrugB, "Combination"), cex=1, pch=1, col=c(2,3,1))

# Chou-Talalay
df_data$CombinationFit<- predict(m.sC, list(x=df_data$RelDoseLog2))
df_data$DrugAFit<- predict(m.sA, list(x=df_data$RelDoseLog2))
df_data$DrugBFit<- predict(m.sB, list(x=df_data$RelDoseLog2))
# y = a + ((b - a)/(1 + exp(-c * (x - d)))) --> x = ((log(((b-a)/(y-a))-1))/-c)+d
df_data$logDoseDrugAFit <- ((log(((df_sigmoid$logDoseDrugA[2]-df_sigmoid$logDoseDrugA[1])/
                  (df_data$CombinationFit-df_sigmoid$logDoseDrugA[1]))-1))/-df_sigmoid$logDoseDrugA[3])+
                                                                                  df_sigmoid$logDoseDrugA[4]
df_data$logDoseDrugBFit <- ((log(((df_sigmoid$logDoseDrugB[2]-df_sigmoid$logDoseDrugB[1])/
                (df_data$CombinationFit-df_sigmoid$logDoseDrugB[1]))-1))/-df_sigmoid$logDoseDrugB[3])+
                                                                                  df_sigmoid$logDoseDrugB[4]
df_data$CI <- (2^(df_data$RelDoseLog2-df_data$logDoseDrugAFit))+(2^(df_data$RelDoseLog2-df_data$logDoseDrugBFit))
par(new=TRUE)
par(fig=c(0.45,1,0.39,0.94))
plot(df_data$CI, df_data$CombinationFit, col=1, cex=0.5, ylim=c(0,1), xlim=c(0,2), xlab="", ylab="",xaxt='n', 
                                                                                                            yaxt='n')
# y = a + ((b - a)/(1 + exp(-c * (x - d)))) --> x = ((log(((b-a)/(y-a))-1))/-c)+d
dosepointlogDoseDrugAFit <- ((log(((df_sigmoid$logDoseDrugA[2]-df_sigmoid$logDoseDrugA[1])/
                                    (dosepointCombinationFit-df_sigmoid$logDoseDrugA[1]))-1))/-df_sigmoid$logDoseDrugA[3])+
                                                                                                          df_sigmoid$logDoseDrugA[4]
dosepointlogDoseDrugBFit <- ((log(((df_sigmoid$logDoseDrugB[2]-df_sigmoid$logDoseDrugB[1])/
                                    (dosepointCombinationFit-df_sigmoid$logDoseDrugB[1]))-1))/-df_sigmoid$logDoseDrugB[3])+
                                                                                                          df_sigmoid$logDoseDrugB[4]
dosepointCI <- (2^(dosepoints-dosepointlogDoseDrugAFit))+(2^(dosepoints-dosepointlogDoseDrugBFit))
lines(dosepointCI, dosepointCombinationFit)
abline(v=1)

# Calculation of CI50
Combination50CI <- ((log(((df_sigmoid$logDoseCombination[2]-df_sigmoid$logDoseCombination[1])/
          (0.5-df_sigmoid$logDoseCombination[1]))-1))/-df_sigmoid$logDoseCombination[3])+df_sigmoid$logDoseCombination[4]
DrugA50CI <- ((log(((df_sigmoid$logDoseDrugA[2]-df_sigmoid$logDoseDrugA[1])/(0.5-df_sigmoid$logDoseDrugA[1]))-1))/
                                                                  -df_sigmoid$logDoseDrugA[3])+df_sigmoid$logDoseDrugA[4]
DrugB50CI <- ((log(((df_sigmoid$logDoseDrugB[2]-df_sigmoid$logDoseDrugB[1])/(0.5-df_sigmoid$logDoseDrugB[1]))-1))/
                                                                  -df_sigmoid$logDoseDrugB[3])+df_sigmoid$logDoseDrugB[4]
CI50<- round((2^(Combination50CI-DrugA50CI))+(2^(Combination50CI-DrugB50CI)),2)
text(1.9,0.5, bquote(~""*CI[50] == .(CI50)*""), cex=0.9, srt=270) 

# Bliss
df_data$BlissExpected<- df_data$DrugAFit*df_data$DrugBFit
df_data$BlissIndex <- df_data$CombinationFit/df_data$BlissExpected
par(new=TRUE)
par(fig=c(0,0.55,0,0.55))
plot(df_data$RelDoseLog2, df_data$BlissIndex, col=1, cex=0.5, ylim=c(0,2), xlab="Relative dose (log2)", 
     ylab="Relative Bliss index")
blissDosepoints<-dosepointCombinationFit/(dosepointDrugAFit*dosepointDrugBFit)
lines(dosepoints,blissDosepoints)
abline(h=1)

# Calculation of BS50
CombinationBI <- ((log(((df_sigmoid$logDoseCombination[2]-df_sigmoid$logDoseCombination[1])/
        (0.5-df_sigmoid$logDoseCombination[1]))-1))/-df_sigmoid$logDoseCombination[3])+df_sigmoid$logDoseCombination[4]
DrugABI<- predict(m.sA, list(x=CombinationBI))
DrugBBI<-predict(m.sB, list(x=CombinationBI))
BS50<- round(0.5/(DrugABI*DrugBBI),2)
text(median(df_data$RelDoseLog2),1.9, bquote(~""*BS[50] == .(BS50)*""), cex=0.9)

# Bliss vs CI
par(new=TRUE)
par(fig=c(0.45,1,0,0.55))
plot(df_data$CI, df_data$BlissIndex, type="p", col=1, cex=0.5, ylim=c(0,2), xlim=c(0,2), 
                                                    xlab="Combination index", ylab="", yaxt='n')
lines(dosepointCI, blissDosepoints)
abline(v=1)
abline(h=1)

# Main title
par(new=TRUE)
par(fig=c(0,0.95,0,0.95))
plot(0:10, type = "n", xaxt="n", yaxt="n", bty="n", xlab = "", ylab = "", main="Chou-Talalay vs Bliss")

dev.off()
