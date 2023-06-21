## CRISPR Screen analysis:
# - We used the screen design of Jastrzebski et al Methods Mol Biol 2016
# - Make count tables from FASTQ files with the Perl script FastqToCountTable.pl and a libraryX.csv file
# - Analyze the count tables with this script and the CRISPRScreenAnalysisLibraries.csv file
# Author: M.F.M. de Rooij PhD, Amsterdam UMC, Spaargaren Lab, 2019, info: m.f.derooij@amsterdamumc.nl
# A small part (aRRA) is adapted from Thomas Kuilman PhD, NKI Amsterdam, Peeper Lab (https://github.com/PeeperLab?tab=repositories)
# The random seeds are different in Windows and Linux R versions, so the p-values on gene level can be different
######################################################################################
## Install DESeq2 and other packages

## Dependencies for devtools in Ubuntu, command line: 
## sudo apt-get install libzmq3-dev libharfbuzz-dev libfribidi-dev libfreetype6-dev libpng-dev libtiff5-dev\ 
## libjpeg-dev build-essential libcurl4-openssl-dev libxml2-dev libssl-dev libfontconfig1-dev

#install.packages(c("BiocManager", "devtools", "rstudioapi", "pheatmap", "RColorBrewer", "scales"))
#BiocManager::install("DESeq2")
#devtools::install_github("JosephCrispell/basicPlotteR")
library("DESeq2")
library("pheatmap")
library("RColorBrewer")
library("scales")
library("basicPlotteR")
######################################################################################
#                                     SETTINGS

# Put this script in the folder where the count tables are located
Workdirectory<- dirname(rstudioapi::getActiveDocumentContext()$path)
## Fill in workdirectory (folder in which the count tables are located, use always slash (/) instead of backslash)
#Workdirectory<- "H:/BioWin/Screens/"

# Count table files: 0 = custom, 1 = all count tables (csv files) in the workdirectory
Data<- 1
# If custom, which count tables?
Filenames<- c("CountTable_test.fastq.gz.csv") # or "CountTable_test.fastq.csv"

# Round numbers in output tables: 0 = Yes, 1 = No
RoundNumbers<- 0

# Data selection (T0/T1/T2 - replicate 1/2/3) (If a replicate is not present, fill in 
# an existing one, and exclude it in the next option)
# If you have 1 replicate per arm, DESeq2 will not run, but aRRA will
# Fill in the column names of the count table
Guide<- "sgRNA"
T0Rep1<- "CGTGAT"
T0Rep2<- "ACATCG"
T0Rep3<- "GCCTAA"
T1Rep1<- "GATCTG"
T1Rep2<- "TCAAGT"
T1Rep3<- "CTGATC"
T2Rep1<- "TGGTCA"
T2Rep2<- "CACTGT"
T2Rep3<- "ATTGGC"

# Exclude replicates: 0 = Include, 1 = Exclude (you cannot exclude T0 or T1 completely -> change your design)
T0_Rep1<- 0
T0_Rep2<- 0
T0_Rep3<- 0
T1_Rep1<- 0
T1_Rep2<- 0
T1_Rep3<- 0
T2_Rep1<- 0
T2_Rep2<- 0
T2_Rep3<- 0

# If there are no replicates:
# Number of pseudocounts to calculate fold change (to reduce noise at low counts, minimum = 1)
pseudocounts <- 10

# Paired replicates: 0 = Paired, 1 = Unpaired
Paired<- 0

# Library: 0: In countTable (1st column = sgRNA, 2nd = Gene), 1 = Brunello-Kinome, 2 = Brunello-WholeGenome, 
# 3 = New1, 4 = New2 (add to CRISPRScreenAnalysisLibraries.csv file)
Library<- 0

# Positive controls: 0 = Essential (Hart et al 2017 G3)
ControlsP<- 0
# Negative controls: 0 = Non-Essential (Hart et al 2014 Mol Sys Biol), 1 = Non-Targeting
ControlsN<- 0
# Even it is not a lethality screen, it is still worth to look how cell viability 
# affect your screen (separation of Essential-Non-Essential)

# Statistics of Essential/Non-Essential genes (only for lethality screens): 0 = With, 1 Without
ConStat<- 0

# Use shrinkage of fold changes noninformative guides: 0 = Yes, 1 = No
Shrinkage<- 0

# Exclude noninformative guides for DESEq2-FDR & aRRA: 0 = Yes, 1 = No
ExclNonInf<- 0

# Type of Screen: 0 = Drop out, 1 = Resistance (Just for ranking of the genes)
Type_of_Screen<- 0

# Minimal fold change of guides to be a hit: 1 = No minimal fold change,  >1: The minimal fold change (linear scale, 2^abs(l2fc))
minimalFoldChange<- 1

# Number of permutation for aRRA (0 = 250 times the number of genes, 1 = Custom)
numPer<- 0
# If number of permutation is custom, how much? (The best p value from alphaRRA is 1/numOfPer)
numOfPer<- 1000

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

# Colors Correlation & MA plots (All guides, positive and negative controls, hits):
ColAll<- "lightgray"
ColP<- "lightpink1"
#ColP<- "red"
ColN<- "lightskyblue"  
#ColN<- "blue"
ColH<- "black"

# Correlation plots of all samples: 0 = Yes, 1 = No (Only recommended when the replicates are messed up) 
pairs<- 1

# MAGeCK: 0 = Yes, 1 = No (this does not work in Windows)
mageck<- 1

# MAGeCK's aRRA with DESeq2 data: 0 = Yes, 1 = No (this does not work in Windows)
mageckRRA<- 1

# MAGeCK should be installed in linux with all its dependencies, see website: https://sourceforge.net/p/mageck/wiki/Home/
######################################################################################
setwd(Workdirectory)
if (Data==0){
  Filenames<- Filenames
} 
if (Data==1){
  Filenames<- Sys.glob("*.csv")
  Filenames<- Filenames["CRISPRScreenAnalysisLibraries.csv" != Filenames]
}

for (Filename in Filenames) {
  # Make a data folder
  dirname2<- paste0(Filename,Sys.time())
  dirname1<-gsub("[[:punct:]]", "", dirname2) 
  dirname<- gsub("\\s", "", dirname1) 
  dir.create(dirname)
  
  # Read count table
  df_raw<- read.csv(file=Filename, sep=",", header=TRUE, stringsAsFactors = FALSE)
  
  # Select replicates
  if(T2_Rep1+T2_Rep2+T2_Rep3 < 3) {
    df_sel<-data.frame(Guide=df_raw[[Guide]], T0R1=df_raw[[T0Rep1]], T0R2=df_raw[[T0Rep2]],T0R3=df_raw[[T0Rep3]], 
                       T1R1=df_raw[[T1Rep1]], T1R2=df_raw[[T1Rep2]], T1R3=df_raw[[T1Rep3]],
                       T2R1=df_raw[[T2Rep1]], T2R2=df_raw[[T2Rep2]], T2R3=df_raw[[T2Rep3]], stringsAsFactors = F)
  } else {
    df_sel<-data.frame(Guide=df_raw[[Guide]], T0R1=df_raw[[T0Rep1]], T0R2=df_raw[[T0Rep2]],T0R3=df_raw[[T0Rep3]], 
                       T1R1=df_raw[[T1Rep1]], T1R2=df_raw[[T1Rep2]], T1R3=df_raw[[T1Rep3]], stringsAsFactors = F)
  }
  
  # Load library info
  df_Gene_IDX<- read.csv("CRISPRScreenAnalysisLibraries.csv", sep=',', header=TRUE, stringsAsFactors = FALSE)
  if (Library == 0) {
    df_Gene_ID<- df_raw[,1:2]
  } 
  if (Library == 1) {
    df_Gene_ID<- df_Gene_IDX[,c(5,6)]
  }
  if (Library == 2) {
    df_Gene_ID<- df_Gene_IDX[,c(8,9)]
  }
  if (Library == 3) {
    df_Gene_ID<- df_Gene_IDX[,c(11,12)]
  }
  if (Library == 4) {
    df_Gene_ID<- df_Gene_IDX[,c(14,15)]
  }
  colnames(df_Gene_ID)<- c("Guide", "GeneSymbol")
  
  # Controls
  df_Control<- df_Gene_IDX[1:927,1:3]
  if (ControlsP==0){
    PC<- "Essential"  
  }
  if (ControlsN==0){
    NC<- "NonEssential"
  }
  if (ControlsN==1){
    NC<- "NonTargeting"
  }
  
  if (T0_Rep2+T0_Rep3+T1_Rep2+T1_Rep3 < 4){
    replicates<- 0
  }else{
    replicates<- 1
  }
  
  # Count table for DESEq2
  counts<- df_sel[,c(
    if (T0_Rep1==0) {2},
    if (T0_Rep2==0) {3},
    if (T0_Rep3==0) {4},
    if (T1_Rep1==0) {5},
    if (T1_Rep2==0) {6},
    if (T1_Rep3==0) {7},
    if (T2_Rep1==0) {8},
    if (T2_Rep2==0) {9},
    if (T2_Rep3==0) {10})]
  rownames(counts)<- df_sel[,1]
  
  # Sample metadata
  time<- c(if (T0_Rep1==0) {"T0"},
           if (T0_Rep2==0) {"T0"},
           if (T0_Rep3==0) {"T0"},
           if (T1_Rep1==0) {"T1"},
           if (T1_Rep2==0) {"T1"},
           if (T1_Rep3==0) {"T1"},
           if (T2_Rep1==0) {"T2"},
           if (T2_Rep2==0) {"T2"},
           if (T2_Rep3==0) {"T2"})
  
  rep<- c(if (T0_Rep1==0) {"R1"},
          if (T0_Rep2==0) {"R2"},
          if (T0_Rep3==0) {"R3"},
          if (T1_Rep1==0) {"R1"},
          if (T1_Rep2==0) {"R2"},
          if (T1_Rep3==0) {"R3"},
          if (T2_Rep1==0) {"R1"},
          if (T2_Rep2==0) {"R2"},
          if (T2_Rep3==0) {"R3"})
  
  df_colData<- data.frame(Time=time, Rep=rep) 
  df_colData$Time<- as.factor(df_colData$Time)
  df_colData$Rep<- as.factor(df_colData$Rep)
  
  if (replicates==0){
    # DESeq2 pipeline
    if (Paired==0) {
      dds<- DESeqDataSetFromMatrix(countData = counts, colData = df_colData, design = ~ Rep + Time)
    }
    if (Paired==1) {
      dds<- DESeqDataSetFromMatrix(countData = counts, colData = df_colData, design = ~ Time)
    }  
    dds$Time<- relevel(dds$Time, "T0")
    dds<- DESeq(dds, fitType = 'local', betaPrior = TRUE)
    # fitType, default is parametric, for CRISPR screens local is better (parametric is in most cases not possible; 
    # betaPrior=TRUE is for shrinkage of noninformative fold changes)
  }
  
  if (replicates==1){
    sf<- estimateSizeFactorsForMatrix(counts)
    CountTableNor<- as.data.frame(t(t(counts)/sf))
  }
  
  # Loop over time point/condition comparisons
  nr<-1
  if ("T2" %in% time) {nr=3}
  
  for (r in 1:nr) {
    if (replicates==0) {
      if (ExclNonInf == 0) {
        if (r==1){
          res<- results(dds, contrast=c("Time","T1","T0"), addMLE=TRUE)
        }
        if (r==2){
          res<- results(dds, contrast=c("Time","T2","T0"), addMLE=TRUE)
        }
        if (r==3){
          res<- results(dds, contrast=c("Time","T2","T1"), addMLE=TRUE)
        }
      } else {
        if (ExclNonInf == 1) {
          # To don't have NA in padj: add independentFiltering=FALSE:
          if (r==1){
            res<- results(dds, contrast=c("Time","T1","T0"), addMLE=TRUE, independentFiltering=FALSE)
          }
          if (r==2){
            res<- results(dds, contrast=c("Time","T2","T0"), addMLE=TRUE, independentFiltering=FALSE)
          }
          if (r==3){
            res<- results(dds, contrast=c("Time","T2","T1"), addMLE=TRUE, independentFiltering=FALSE)
          }
        }
      } 
      
      # DESeq2 data table
      df_res<- as.data.frame(res)
      
      if (Shrinkage==1) {
        df_res$log2FoldChange<- df_res$lfcMLE
      }
      df_baseMeanPerLvl<- as.data.frame(sapply(levels(dds$Time), 
                                               function(lvl) rowMeans(counts(dds, normalized=TRUE)[,dds$Time==lvl])))
      if (r==1){
        df_res$BaseMeanA<- df_baseMeanPerLvl$T0
        df_res$logBaseMeanA<- log(df_baseMeanPerLvl$T0+1)/log(10)
        df_res$BaseMeanB<- df_baseMeanPerLvl$T1
      } 
      if (r==2){
        df_res$BaseMeanA<- df_baseMeanPerLvl$T0
        df_res$logBaseMeanA<- log(df_baseMeanPerLvl$T0+1)/log(10)
        df_res$BaseMeanB<- df_baseMeanPerLvl$T2
      } 
      if (r==3){
        df_res$BaseMeanA<- df_baseMeanPerLvl$T1
        df_res$logBaseMeanA<- log(df_baseMeanPerLvl$T1+1)/log(10)
        df_res$BaseMeanB<- df_baseMeanPerLvl$T2
      } 
    }
    
    if (replicates==1){
      df_res<- data.frame(Guide=rownames(CountTableNor), stringsAsFactors = F)
      if (r==1){
        df_res$BaseMeanA<- CountTableNor$T0R1
        df_res$logBaseMeanA<- log(df_res$BaseMeanA+1)/log(10)
        df_res$BaseMeanB<- CountTableNor$T1R1
      } 
      if (r==2){
        df_res$BaseMeanA<- CountTableNor$T0R1
        df_res$logBaseMeanA<- log(df_res$BaseMeanA+1)/log(10)
        df_res$BaseMeanB<- CountTableNor$T2R1
      } 
      if (r==3){
        df_res$BaseMeanA<- CountTableNor$T1R1
        df_res$logBaseMeanA<- log(df_res$BaseMeanA+1)/log(10)
        df_res$BaseMeanB<- CountTableNor$T2R1
      } 
      df_res$FoldChange<- (df_res$BaseMeanB+pseudocounts)/(df_res$BaseMeanA+pseudocounts)
      df_res$log2FoldChange<- log(df_res$FoldChange)/log(2)
      df_res$log2FoldChange[is.na(df_res$log2FoldChange)]<- 0
      df_res$pvalue<- 1
      df_res$padj<- 1
      rownames(df_res)<- df_res$Guide
    }
    
    # Guide IDs
    df_res$Guide<-rownames(df_res)
    df_res<- merge(df_Gene_ID, df_res, by='Guide', all.y=T)
    df_res$log2FoldChange[is.na(df_res$log2FoldChange)]<-0
    df_res$FoldChange<- 2^df_res$log2FoldChange
    
    # Exclude noninformative guides from aRRA
    df_res2<-df_res
    df_res2<-df_res2[!is.na(df_res2$padj),]
    
    if (minimalFoldChange > 1) {
      df_res$padj[!is.na(df_res$padj) & df_res$FoldChange >= (1/minimalFoldChange) & df_res$FoldChange <= minimalFoldChange]<- 1
    }
    
    # Control genes
    df_res$Type<- "x"
    df_res[toupper(df_res$GeneSymbol) %in% df_Control[[PC]][nchar(df_Control[[PC]])>0],"Type"]<- "p"
    df_res[toupper(df_res$GeneSymbol) %in% df_Control[[NC]][nchar(df_Control[[NC]])>0],"Type"]<- "n"
    df_PC<-df_res[df_res$Type=="p",]
    df_NC<-df_res[df_res$Type=="n",]
    
    # Determine hits
    df_hits_total<- df_res[which(df_res$padj < 0.1),]
    df_res$count<- 1
    df_hits_A<-aggregate(df_res$count, by=list(df_res$GeneSymbol), FUN=sum)
    colnames(df_hits_A)<- c("GeneSymbol", "TotalGuides")
    if (nrow(df_hits_total)>0) {
      df_hits_down<-df_hits_total[df_hits_total$log2FoldChange<median(df_res$log2FoldChange,na.rm=TRUE),]
      df_hits_up<-df_hits_total[df_hits_total$log2FoldChange>median(df_res$log2FoldChange,na.rm=TRUE),]
      df_hits_down$GeneSymbol[is.na(df_hits_down$GeneSymbol)]<- 0
      df_hits_up$GeneSymbol[is.na(df_hits_up$GeneSymbol)]<- 0
      if (nrow(df_hits_down)>0) {
        df_hits_down$count<- 1
        df_hits_B<-aggregate(df_hits_down$count, by=list(df_hits_down$GeneSymbol), FUN=sum)
        colnames(df_hits_B)<- c("GeneSymbol", "HitsDown")
        df_hits_A<-merge(df_hits_A, df_hits_B, by="GeneSymbol",all.x=TRUE)
        df_hits_A$HitsDown[is.na(df_hits_A$HitsDown)]<-0
      }
      if (nrow(df_hits_down)==0) {
        df_hits_A$HitsDown<- 0
      }
      if (nrow(df_hits_up)>0) {
        df_hits_up$count<- 1
        df_hits_C<-aggregate(df_hits_up$count, by=list(df_hits_up$GeneSymbol), FUN=sum)
        colnames(df_hits_C)<- c("GeneSymbol", "HitsUp")
        df_hits_A<-merge(df_hits_A, df_hits_C, by="GeneSymbol",all.x=TRUE)
        df_hits_A$HitsUp[is.na(df_hits_A$HitsUp)]<-0
      }
      if (nrow(df_hits_up)==0) {
        df_hits_A$HitsUp<- 0
      }
    }
    if (nrow(df_hits_total)==0) {
      df_hits_A$HitsDown<- 0
      df_hits_A$HitsUp<- 0
    }
    df_res_type<- data.frame(GeneSymbol=df_res$GeneSymbol, Type=df_res$Type)
    df_res_type<-unique(df_res_type)
    df_hits_A<- merge(df_hits_A, df_res_type, by="GeneSymbol")
    
    # Fold change info
    df_hits_D<- aggregate(df_res$FoldChange, by=list(df_res$GeneSymbol), FUN=min, na.rm=T)
    colnames(df_hits_D)<-c("GeneSymbol", "MinFoldChange")
    df_hits_E<- aggregate(df_res$FoldChange, by=list(df_res$GeneSymbol), FUN=max, na.rm=T)
    colnames(df_hits_E)<-c("GeneSymbol", "MaxFoldChange")
    df_hits_F<- aggregate(df_res$FoldChange, by=list(df_res$GeneSymbol), FUN=median, na.rm=T)
    colnames(df_hits_F)<-c("GeneSymbol", "MedianFoldChange")
    df_hits_A<-merge(df_hits_A, df_hits_D, by="GeneSymbol")
    df_hits_A<-merge(df_hits_A, df_hits_E, by="GeneSymbol")
    df_hits_A<-merge(df_hits_A, df_hits_F, by="GeneSymbol")
    
    if (replicates==0){
      # AlphaRRA (from guide to gene statistics)
      if (r==1){
        resDepleted<- results(dds, contrast=c("Time","T1","T0"), altHypothesis = "less", addMLE=TRUE)
        resEnriched<- results(dds, contrast=c("Time","T1","T0"), altHypothesis = "greater", addMLE=TRUE)
        con<- "T0vsT1"
      }
      if (r==2){
        resDepleted<- results(dds, contrast=c("Time","T2","T0"), altHypothesis = "less", addMLE=TRUE)
        resEnriched<- results(dds, contrast=c("Time","T2","T0"), altHypothesis = "greater", addMLE=TRUE)
        con<- "T0vsT2"
      }
      if (r==3){
        resDepleted<- results(dds, contrast=c("Time","T2","T1"), altHypothesis = "less", addMLE=TRUE)
        resEnriched<- results(dds, contrast=c("Time","T2","T1"), altHypothesis = "greater", addMLE=TRUE)
        con<- "T1vsT2"
      }
      
      resDepleted<- resDepleted[order(rownames(resDepleted)), ]
      resEnriched<- resEnriched[order(rownames(resEnriched)), ]
      if(Shrinkage==0) {
        df_RRA<- data.frame(log2fc = resDepleted$log2FoldChange,
                            pvalueDepleted = resDepleted$pvalue,
                            pvalueEnriched = resEnriched$pvalue,
                            row.names = rownames(resDepleted))
      }
      if(Shrinkage==1) {
        df_RRA<- data.frame(log2fc = resDepleted$lfcMLE,
                            pvalueDepleted = resDepleted$pvalue,
                            pvalueEnriched = resEnriched$pvalue,
                            row.names = rownames(resDepleted))
      }
      df_RRA$Guide<-rownames(df_RRA)
      df_RRA<- merge(df_Gene_ID, df_RRA, by='Guide', all.y=T)
      df_RRA<- df_RRA[df_RRA$Guide %in% df_res2$Guide,]
      
      # MAGeCK's aRRA with DESeq2 data (in Linux)
      if (mageckRRA==0){
        df_RRAmageck<- df_RRA
        perDep<- nrow(df_RRAmageck[df_RRAmageck$pvalueDepleted < 0.25,])/nrow(df_RRAmageck) 
        perEnr<- nrow(df_RRAmageck[df_RRAmageck$pvalueEnriched < 0.25,])/nrow(df_RRAmageck) 
        
        df_RRAmageckDep<- df_RRAmageck[order(df_RRAmageck$log2fc),]
        df_RRAmageckDep$listID<- "x"
        df_RRAmageckDep$pval<-df_RRAmageckDep$pvalueDepleted 
        df_RRAmageckDep$pvalueDepleted<- NULL
        df_RRAmageckDep$pvalueEnriched<-NULL
        df_RRAmageckDep$log2fc<-NULL
        write.table(df_RRAmageckDep, "/tmp/Dep.txt", row.names=F, sep="\t")
        system(paste0("~/anaconda3/bin/RRA -i /tmp/Dep.txt -o ", paste0(dirname,"/DESeq2_aRRAmageck_",con,"_Dep.txt"), " -p ",perDep))
        
        df_RRAmageckEnr<- df_RRAmageck[order(df_RRAmageck$log2fc, decreasing = T),]
        df_RRAmageckEnr$listID<- "x"
        df_RRAmageckEnr$pval<-df_RRAmageckEnr$pvalueEnriched
        df_RRAmageckEnr$pvalueDepleted<- NULL
        df_RRAmageckEnr$pvalueEnriched<-NULL
        df_RRAmageckEnr$log2fc<-NULL           
        write.table(df_RRAmageckEnr, "/tmp/Enr.txt", row.names=F, sep="\t")
        system(paste0("~/anaconda3/bin/RRA -i /tmp/Enr.txt -o ", paste0(dirname,"/DESeq2_aRRAmageck_",con,"_Enr.txt"), " -p ",perEnr))
      }
    }
    if (replicates==1){
      df_RRA <- data.frame(GeneSymbol=df_res$GeneSymbol,
                           log2fc = df_res$log2FoldChange,
                           pvalueDepleted = 1,
                           pvalueEnriched = 1,
                           row.names = df_res$Guide)
      
      if (r==1){
        con<- "T0vsT1"
      }
      if (r==2){
        con<- "T0vsT2"
      }
      if (r==3){
        con<- "T1vsT2"
      }
    }
    
    # aRRA in R
    
    # Rank by fold change
    df_RRA$scoreDepleted<- rank(df_RRA$log2fc) / nrow(df_RRA)
    df_RRA$scoreEnriched<- rank(-df_RRA$log2fc) / nrow(df_RRA)
    
    if (replicates==0){
      # Apply alpha criterion based on pvalues
      df_RRA$scoreDepleted[df_RRA$pvalueDepleted > 0.25]<- 1
      df_RRA$scoreEnriched[df_RRA$pvalueEnriched > 0.25]<- 1
    }
    if (replicates==1){
      # Apply alpha criterion based on rank
      df_RRA$scoreDepleted[df_RRA$scoreDepleted > 0.25]<- 1
      df_RRA$scoreEnriched[df_RRA$scoreEnriched > 0.25]<- 1
    }
    
    if (minimalFoldChange > 1) {
      df_RRA$scoreDepleted[df_RRA$log2fc >= (log(1/minimalFoldChange)/log(2))]<- 1
      df_RRA$scoreEnriched[df_RRA$log2fc <= (log(minimalFoldChange)/log(2))]<- 1
    }
    
    # Perform RRA
    alphaBeta<- function(p.in) {
      p.in<- sort(p.in)
      n<- length(p.in)
      return(min(pbeta(p.in, 1:n, n - (1:n) + 1)))
    }
    
    # Calculate rho per gene
    df_RRA$rhoDepleted<- unsplit(sapply(split(df_RRA$scoreDepleted,
                                              df_RRA$GeneSymbol), alphaBeta), df_RRA$GeneSymbol)
    df_RRA$rhoEnriched<- unsplit(sapply(split(df_RRA$scoreEnriched,
                                              df_RRA$GeneSymbol), alphaBeta), df_RRA$GeneSymbol)
    
    # Make a null distribution and calculate pvalues per gene
    n.guides<- sort(unique(table(df_RRA$GeneSymbol)))
    makeRhoNull<- function(n, p, nperm) {
      sapply(1:nperm, function(x) {
        p.in<- sort.int(sample(p, n, replace = FALSE))
        alphaBeta(p.in)
      })
    }
    if (numPer==0){
      permutations<- 250*nrow(df_hits_A)
    }
    if (numPer==1){ 
      permutations<- numOfPer
    }
    
    # Depletion
    set.seed(12345)
    rho.nulls<- lapply(n.guides, makeRhoNull, df_RRA$scoreDepleted, permutations)
    names(rho.nulls)<- as.character(n.guides)
    # Prevent pvalues of 0
    for (i in 1:length(rho.nulls)){
      rho.nulls[[i]][[1]]<-0
    }
    df_RRA$pvalueDepleted<- unsplit(sapply(split(df_RRA$rhoDepleted, df_RRA$GeneSymbol), function(x) {
      mean(rho.nulls[[as.character(length(x))]] <= x[1])
    }), df_RRA$GeneSymbol)
    # Enrichment
    set.seed(12345)
    rho.nulls<- lapply(n.guides, makeRhoNull, df_RRA$scoreEnriched, permutations)
    names(rho.nulls)<- as.character(n.guides)
    # Prevent pvalues of 0
    for (i in 1:length(rho.nulls)){
      rho.nulls[[i]][[1]]<-0
    }
    df_RRA$pvalueEnriched<- unsplit(sapply(split(df_RRA$rhoEnriched,df_RRA$GeneSymbol), function(x) {
      mean(rho.nulls[[as.character(length(x))]] <= x[1])
    }), df_RRA$GeneSymbol)
    
    # Select aRRA data and calculate fdr
    df_resultsRRA<- data.frame(GeneSymbol=as.character(df_RRA$GeneSymbol), rhoDepleted=df_RRA$rhoDepleted, pvalueDepleted=df_RRA$pvalueDepleted, 
                               rhoEnriched=df_RRA$rhoEnriched, pvalueEnriched=df_RRA$pvalueEnriched, stringsAsFactors = F)
    df_resultsRRA<-unique(df_resultsRRA)
    df_resultsRRA$fdrDepleted<-p.adjust(df_resultsRRA$pvalueDepleted, method='fdr')
    df_resultsRRA$fdrEnriched<-p.adjust(df_resultsRRA$pvalueEnriched, method='fdr')
    
    # Merge to DESEq2 data 
    df_geneRRA<-merge(df_resultsRRA, df_hits_A, by="GeneSymbol", all.y=T)
    
    # Order genes
    if (Type_of_Screen==0) {
      df_geneRRA<-df_geneRRA[order(df_geneRRA$rhoDepleted),]
    }
    if (Type_of_Screen==1) {
      df_geneRRA<-df_geneRRA[order(df_geneRRA$rhoEnriched),]
    }  
    
    # Write guide and gene tables
    
    # Guide Table
    df_res_print<- df_res[,c("GeneSymbol","Guide","Type","BaseMeanA","BaseMeanB","FoldChange","pvalue","padj")]
    df_res_print<- df_res_print[order(df_res_print$Type),]
    if (RoundNumbers==0){
      df_res_print[,c("BaseMeanA","BaseMeanB")]<-round(df_res_print[,c("BaseMeanA","BaseMeanB")],0)
      df_res_print[,c("FoldChange","pvalue","padj")]<-signif(df_res_print[,c("FoldChange","pvalue","padj")],4)
    }
    write.csv(df_res_print, paste0(dirname,"/DESeq2 ",con," Guides.csv"), row.names = FALSE, quote = F)
    # Gene Table
    df_geneRRA_print<-df_geneRRA[,c("GeneSymbol", "Type", "TotalGuides", "HitsDown","HitsUp","MinFoldChange",
                                    "MedianFoldChange","MaxFoldChange", "rhoDepleted",	"pvalueDepleted",	"fdrDepleted",	"rhoEnriched",	"pvalueEnriched", "fdrEnriched")] 
    if (RoundNumbers==0){
      df_geneRRA_print[,c("MinFoldChange","MedianFoldChange","MaxFoldChange", "rhoDepleted",	"pvalueDepleted",	"fdrDepleted",	"rhoEnriched",	"pvalueEnriched", 
                          "fdrEnriched")]<-signif(df_geneRRA_print[,c("MinFoldChange","MedianFoldChange","MaxFoldChange", "rhoDepleted", "pvalueDepleted",	"fdrDepleted",	
                                                                      "rhoEnriched",	"pvalueEnriched", "fdrEnriched")],4)
    }
    write.csv(df_geneRRA_print, paste0(dirname,"/DESeq2 ",con," Genes.csv"), row.names = FALSE, quote = F)
    
    
    # Prepare plot data
    
    if (replicates==0){
      # Correlation plots
      df_normCounts<- as.data.frame(counts(dds, normalized=TRUE))
      xyrangeMA<-c(0, max(df_res$BaseMeanA, na.rm=TRUE))
      repx<-c(if (T0_Rep1==0 && T0_Rep2==0) {"T0R1"},
              if (T0_Rep1==0 && T0_Rep3==0) {"T0R1"},
              if (T0_Rep2==0 && T0_Rep3==0) {"T0R2"}, 
              if (T1_Rep1==0 && T1_Rep2==0) {"T1R1"},
              if (T1_Rep1==0 && T1_Rep3==0) {"T1R1"},
              if (T1_Rep2==0 && T1_Rep3==0) {"T1R2"},
              if (T2_Rep1==0 && T2_Rep2==0) {"T2R1"},
              if (T2_Rep1==0 && T2_Rep3==0) {"T2R1"},
              if (T2_Rep2==0 && T2_Rep3==0) {"T2R2"})
      
      repy<-c(if (T0_Rep1==0 && T0_Rep2==0) {"T0R2"}, 
              if (T0_Rep1==0 && T0_Rep3==0) {"T0R3"},
              if (T0_Rep2==0 && T0_Rep3==0) {"T0R3"},
              if (T1_Rep1==0 && T1_Rep2==0) {"T1R2"}, 
              if (T1_Rep1==0 && T1_Rep3==0) {"T1R3"},
              if (T1_Rep2==0 && T1_Rep3==0) {"T1R3"},
              if (T2_Rep1==0 && T2_Rep2==0) {"T2R2"}, 
              if (T2_Rep1==0 && T2_Rep3==0) {"T2R3"},
              if (T2_Rep2==0 && T2_Rep3==0) {"T2R3"}) 
      df_rep<-data.frame(repx, repy, stringsAsFactors = FALSE)
      
      df_normCounts$Guide<-rownames(df_normCounts)
      if (r==1){
        # Write count table with the normalized counts
        write.csv(df_normCounts[,c(ncol(df_normCounts),1:ncol(df_normCounts)-1)], paste0(dirname,"/NormalizedCounts.csv"), row.names = F, quote = F)
      }
      df_normCounts<-merge(df_Gene_ID, df_normCounts, by='Guide', all.y=T)
      df_normCounts$Guide<-NULL
      df_normCountsP<- df_normCounts[df_normCounts$GeneSymbol %in% df_Control[[PC]][nchar(df_Control[[PC]])>0],]
      df_normCountsN<- df_normCounts[df_normCounts$GeneSymbol %in% df_Control[[NC]][nchar(df_Control[[NC]])>0],]
      df_normCounts$GeneSymbol<-NULL
      df_normCountsP$GeneSymbol<-NULL
      df_normCountsN$GeneSymbol<-NULL
    }
    if (replicates==1 & r==1){
      # Write count table with the normalized counts
      df_normCounts<- CountTableNor
      df_normCounts$Guide<-rownames(df_normCounts)
      write.csv(df_normCounts[,c(ncol(df_normCounts),1:ncol(df_normCounts)-1)], paste0(dirname,"/NormalizedCounts.csv"), row.names = F, quote = F)
    }
    # Tophits
    RRAdep<-df_geneRRA[order(df_geneRRA$rhoDepleted),]
    RRAenr<-df_geneRRA[order(df_geneRRA$rhoEnriched),]
    tophits<- c(RRAdep[1:10,1], RRAenr[1:10,1])
    
    # Volcano plot
    
    # Median log2 fold change
    df_geneRRA$ml2fc<- log2(df_geneRRA$MedianFoldChange)
    df_geneRRA$ml2fc[is.na(df_geneRRA$ml2fc)]<-0
    
    # Take rho/fdr depleted/enriched dependent on fold change
    df_geneRRA$rhoDepleted[is.na(df_geneRRA$rhoDepleted)]<-1
    df_geneRRA$rhoEnriched[is.na(df_geneRRA$rhoEnriched)]<-1
    df_geneRRA$rho<- df_geneRRA$rhoDepleted
    df_geneRRA$rho[df_geneRRA$ml2fc > 0]<- df_geneRRA$rhoEnriched[df_geneRRA$ml2fc > 0]
    
    df_geneRRA$fdrDepleted[is.na(df_geneRRA$fdrDepleted)]<-1
    df_geneRRA$fdrEnriched[is.na(df_geneRRA$fdrEnriched)]<-1
    df_geneRRA$fdr<- df_geneRRA$fdrDepleted
    df_geneRRA$fdr[df_geneRRA$ml2fc > 0]<- df_geneRRA$fdrEnriched[df_geneRRA$ml2fc > 0]
    
    # Colors
    df_geneRRA$col<- ColAll
    df_geneRRA$col[df_geneRRA$Type=="p"]<- ColP
    df_geneRRA$col[df_geneRRA$Type=="n"]<- ColN
    df_geneRRA$col[df_geneRRA$GeneSymbol %in% tophits]<- ColH
    df_geneRRA$pch<- 16
    df_geneRRA$pch[df_geneRRA$fdr < 0.1 & df_geneRRA$ml2fc < 0]<- 25
    df_geneRRA$pch[df_geneRRA$fdr < 0.1 & df_geneRRA$ml2fc > 0]<- 24
    
    # Mix essential and non-essential randomly
    df_geneRRAx<- df_geneRRA[df_geneRRA$Type=="x",]
    df_geneRRAc<- df_geneRRA[df_geneRRA$Type!="x",]
    set.seed(101)
    df_geneRRAc<- df_geneRRAc[sample(1:nrow(df_geneRRAc)),]
    df_geneRRA<-rbind(df_geneRRAx,df_geneRRAc)
    
    if (ConStat==0) {
      # F measure (with cutoff the intersection), harmonic mean of precision and recall
      truePosGene<- nrow(df_geneRRA[!is.na(df_geneRRA$fdrDepleted) & df_geneRRA$fdrDepleted<0.1 & df_geneRRA$Type=='p',])
      posPredictGene<- nrow(df_geneRRA[!is.na(df_geneRRA$fdrDepleted) & df_geneRRA$fdrDepleted<0.1 & df_geneRRA$Type=='p',])+
        nrow(df_geneRRA[!is.na(df_geneRRA$fdrDepleted) & df_geneRRA$fdrDepleted<0.1 & df_geneRRA$Type=='n',])*sum(df_geneRRA$Type=='p')/sum(df_geneRRA$Type=='n')
      posControlsGene<- sum(df_geneRRA$Type=='p')
      precisionGene<- truePosGene/posPredictGene
      recallGene<- truePosGene/posControlsGene
      F1Gene<- 2/(1/precisionGene + 1/recallGene)
    }
    
    # Axes limits
    xrangeVOL<- c(min(df_geneRRA$ml2fc, na.rm=TRUE), max(df_geneRRA$ml2fc, na.rm=TRUE))
    yrangeVOL<- c(0, max(-log10(df_geneRRA$rho), na.rm=TRUE))
    
    # MA plots
    Genes_of_interest<-c(" ", "Hitlist", if (MA_all_genes==0) {df_geneRRA$GeneSymbol}, 
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
    df_res<-rbind(df_resx,df_resc)
    
    # Axes limits
    xrangeMA<-c(min(df_res$logBaseMeanA, na.rm=TRUE)-0.5, max(df_res$logBaseMeanA, na.rm=TRUE)+0.5)
    yrangeMA<-c(min(df_res$log2FoldChange, na.rm=TRUE)-0.5, max(df_res$log2FoldChange, na.rm=TRUE)+0.5)
    
    # Density plot fold change
    denY_tot<-density(df_res$log2FoldChange, from=yrangeMA[1], to=yrangeMA[2], na.rm=T)
    denY_tot$y[1]<- 0
    denY_tot$y[length(denY_tot$y)]<- 0
    denY_PC<-density(df_PC$log2FoldChange, from=yrangeMA[1], to=yrangeMA[2], na.rm=T)
    denY_PC$y[1]<- 0
    denY_PC$y[length(denY_PC$y)]<- 0
    denY_NC<-density(df_NC$log2FoldChange, from=yrangeMA[1], to=yrangeMA[2], na.rm=T)
    denY_NC$y[1]<- 0
    denY_NC$y[length(denY_NC$y)]<- 0
    denYMax<- max(c(denY_tot$y, denY_PC$y, denY_NC$y))
    # Density plot read count
    denX_tot<- density(df_res$logBaseMeanA, from=xrangeMA[1], to=xrangeMA[2], na.rm=T)
    denX_tot$y[1]<- 0
    denX_tot$y[length(denX_tot$y)]<- 0
    denX_PC<- density(df_PC$logBaseMeanA, from=xrangeMA[1], to=xrangeMA[2], na.rm=T)
    denX_PC$y[1]<- 0
    denX_PC$y[length(denX_PC$y)]<- 0
    denX_NC<- density(df_NC$logBaseMeanA, from=xrangeMA[1], to=xrangeMA[2], na.rm=T)
    denX_NC$y[1]<- 0
    denX_NC$y[length(denX_NC$y)]<- 0
    denXMax<- max(c(denX_tot$y, denX_PC$y, denX_NC$y))
    
    if (ConStat==0) {
      # Robust Z prime (median and median absolute deviation)
      Mp<- median(df_PC$log2FoldChange,na.rm=TRUE)
      Mn<- median(df_NC$log2FoldChange, na.rm=TRUE)
      SDp<-mad(df_PC$log2FoldChange, na.rm=TRUE)
      SDn<-mad(df_NC$log2FoldChange, na.rm=TRUE)
      Zprime<- 1-((3*(SDp+SDn))/(abs(Mp-Mn))) 
      
      # Intersection postive and negative controls
      poi<- which(diff(denY_PC$y > denY_NC$y) != 0) 
      intersection<- denY_PC$x[poi][denY_PC$x[poi]<0][which.min(abs(denY_PC$x[poi][denY_PC$x[poi]<0]))]
      
      # F measure (with cutoff the intersection), harmonic mean of precision and recall
      if (replicates==0){
        truePos<- nrow(df_PC[!is.na(df_PC$padj) & df_PC$padj<0.1 & df_PC$log2FoldChange<0,])
        posPredict<- nrow(df_PC[!is.na(df_PC$padj) & df_PC$padj<0.1 & df_PC$log2FoldChange<0,])+nrow(df_NC[!is.na(df_NC$padj) & df_NC$padj<0.1 & df_NC$log2FoldChange<0,])*nrow(df_PC)/nrow(df_NC)
        posControls<- nrow(df_PC)
        precision<- truePos/posPredict
        recall<- truePos/posControls
        F1<- 2/(1/precision + 1/recall)
      }
      if (replicates==1){
        truePos<- nrow(df_PC[df_PC$log2FoldChange < intersection,])
        posPredict<- nrow(df_PC[df_PC$log2FoldChange < intersection,])+nrow(df_NC[df_NC$log2FoldChange < intersection,])*nrow(df_PC)/nrow(df_NC)
        posControls<- nrow(df_PC)
        precision<- truePos/posPredict
        recall<- truePos/posControls
        F1<- 2/(1/precision + 1/recall)
      }
    }
    
    # Draw Plots in PDF 
    pdf(paste0(dirname,"/DESeq2 ",con," Plots.pdf"), width=7, height=7)
    if (replicates==0){
      # DESEq2 plots
      plotDispEsts(dds, main="Dispersion plot")
      
      rld<-rlog(dds, blind=FALSE)
      sampleDists<- dist(t(assay(rld)))
      sampleDistMatrix<- as.matrix(sampleDists)
      rownames(sampleDistMatrix)<- paste(rld$Time, rld$Rep)
      colnames(sampleDistMatrix)<- paste(rld$Time, rld$Rep)
      colors<- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
      
      pheatmap(sampleDistMatrix,
               clustering_distance_rows=sampleDists,
               clustering_distance_cols=sampleDists,
               col=colors)
      
      print(plotPCA(rld, intgroup=c("Time", "Rep")))
    }  
    if (pairs==0) {
      pairs(df_sel[-1], cex=0.1)
      pairs(df_normCounts, cex=0.1)
    }
    
    # Home made plots
    
    # Correlation plots
    if (replicates==0){
      par(mfrow=c(3,3))
      for (i in c(1:length(repx))) {
        e1<- df_rep[i,1]
        e2<- df_rep[i,2]
        reg1<- lm(df_normCounts[,e2]~df_normCounts[,e1])
        coefs<- coef(reg1)
        intercept<- round(coefs[1], 3)
        slope<- round(coefs[2],3)
        r2<- round(summary(reg1)$r.squared, 3)
        eq1<- bquote(italic(y) == .(slope)*italic(x) + .(intercept) ~"("*r^2 == .(r2)*")")
        
        plot(df_normCounts[,e1], df_normCounts[,e2], type="p", col=ColAll, cex=.1, main="Correlation plot",
             xlab=df_rep[i,1], ylab=df_rep[i,2], xlim=xyrangeMA, ylim=xyrangeMA)
        points(df_normCountsP[,e1], df_normCountsP[,e2], type="p", col=ColP, cex=.3)
        points(df_normCountsN[,e1], df_normCountsN[,e2], type="p", col=ColN, cex=.3)
        abline(reg1, col=6)
        abline(0,1, col=1)
        text(x = xyrangeMA[2]/2, y = xyrangeMA[2]*0.95, labels= eq1, cex=0.8)
      }
    }
    # Volcano plot
    
    # Margins
    par(mfrow=c(1,1))
    par(mar=c(4,4,4,4))
    par(fig=c(0.1,0.9,0.1,0.9))
    par(bty="l")
    # Plot axes
    plot(0, pch = '', 
         main= "Volcano plot",
         xlab= "Log2 median fold change",
         ylab= "RRA score",
         cex.lab=1, cex.axis=1, las=1, xlim=xrangeVOL, ylim=yrangeVOL)
    if (ConStat==0) {
      mtext(side=3, line=0.5, at=mean(xrangeVOL), cex=0.7, 
            paste0("(F1: ", format(round(F1Gene,2),nsmall = 2), " (Pre: ", format(round(precisionGene,2),nsmall = 2), 
                   " Rec: ", format(round(recallGene,2),nsmall = 2), "))"))
    }
    
    # Vertical and horizontal lines
    abline(v=0.2*(-100:100), lty=3, col="gray")
    abline(h=2*(0:100), lty=3, col="gray")
    
    # Actual plot
    points(df_geneRRA$ml2fc, -log10(df_geneRRA$rho), type="p", pch=16, col=alpha(df_geneRRA$col,0.7), cex=0.7)
    
    # lines
    abline(v=0, lty=2)
    
    # Show gene symbols of tophits
    df_volcanoGenes<- df_geneRRA[df_geneRRA$GeneSymbol %in% tophits,]
    addTextLabels(df_volcanoGenes$ml2fc,-log10(df_volcanoGenes$rho),df_volcanoGenes$GeneSymbol, avoidPoints = TRUE,
                  keepLabelsInside = TRUE, col.label="black", cex.label=0.7)
    
    text(0,yrangeVOL[2], substitute(paste(italic('Depleted    Enriched'))))
    
    # MA plots
    par(mfrow=c(1,1))
    par(bty="o")
    for(Gene in Genes_of_interest){
      df_GOI<-df_res[df_res$GeneSymbol %in% Gene,]
      
      # Main MA plot
      par(mar=c(4,4,0,0))
      par(fig=c(0.1,0.7,0.1,0.7))
      plot(df_res$logBaseMeanA, df_res$log2FoldChange, type="p", col=df_res$col, bg=df_res$col, cex=1, pch=df_res$pch, xlab="Log10 Average Read Counts (Control)", 
           ylab="Log2 Fold Change", cex.lab=1, cex.axis=1, xlim=xrangeMA, ylim=yrangeMA)
      if (Gene=="Hitlist"&& nrow(df_hits_total)>=1){
        df_hits_total<-df_res[df_res$Guide %in% df_hits_total$Guide,]
        points(df_hits_total$logBaseMeanA, df_hits_total$log2FoldChange, type="p", col=ColH, bg=ColH, cex=1, pch=df_hits_total$pch)
      }
      if (nrow(df_GOI)>=1){
        points(df_GOI$logBaseMeanA, df_GOI$log2FoldChange, type="p", col=ColH, bg=ColH, cex=1.5, pch=df_GOI$pch)
      }
      legend("bottomleft",legend=c( if (nchar(Gene)>1 && !is.na(Gene)) {Gene}, "All guides", "Essentials",if (ControlsN==0 | ControlsN==1){"Non-essentials"}, 
                                    if (ControlsN==2){"Non-targeting"}, 'Significantly enriched', 'Significantly depleted'), cex=0.8, pch=c(if (nchar(Gene)>1 && !is.na(Gene)) 
                                    {16},16,16,16, 24,25), col=c(if (nchar(Gene)>1 && !is.na(Gene)) {ColH}, ColAll,ColP,ColN, "black", "black"))
      abline(median(df_res$log2FoldChange, na.rm=TRUE),0, col=1, lty=3, untf=TRUE)
      if (ConStat==0) {
        text(median(xrangeMA), yrangeMA[2],cex=0.7, labels=paste0("(rZ': ", format(round(Zprime,2), nsmall=2), " NP50: ", format(round(intersection,2), nsmall=2), 
                                                                  " F1: ", format(round(F1,2),nsmall = 2), " (Pre: ", format(round(precision,2),nsmall = 2), " Rec: ", format(round(recall,2),nsmall = 2), "))"))
      }
      
      # Density fold change
      par(mar=c(4,0,0,4))
      par(fig=c(0.7,0.9,0.1,0.7),new=TRUE)
      plot(denY_tot$y, denY_tot$x, ylim=yrangeMA, xlim=(c(0,denYMax)), type='l', axes=FALSE, col=ColAll, xlab="", 
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
      if (Gene=="Hitlist" && nrow(df_hits_total)>1){
        denY_hits_total<- density(df_hits_total$log2FoldChange, from=yrangeMA[1], to=yrangeMA[2], na.rm=T)
        denY_hits_total$y[1]<- 0
        denY_hits_total$y[length(denY_hits_total$y)]<- 0
        lines(denY_hits_total$y, denY_hits_total$x, col=ColH, lwd=2)
        rgb.val<- col2rgb(ColH)
        polygon(denY_hits_total$y, denY_hits_total$x, col=rgb(rgb.val[1]/255,rgb.val[2]/255,rgb.val[3]/255,alpha=0.3), lwd=0.1)
      }
      if (nrow(df_GOI)>1){
        par(new=TRUE)
        denY_GOI<- density(df_GOI$log2FoldChange, from=yrangeMA[1], to=yrangeMA[2], na.rm=T)
        denY_GOI$y[1]<- 0
        denY_GOI$y[length(denY_GOI$y)]<- 0
        lines(denY_GOI$y, denY_GOI$x, col=ColH, lwd=2)
        rgb.val<- col2rgb(ColH)
        polygon(denY_GOI$y, denY_GOI$x, col=rgb(rgb.val[1]/255,rgb.val[2]/255,rgb.val[3]/255,alpha=0.3), lwd=0.1)
      }
      
      # Density read counts
      par(mar=c(0,4,4,0))
      par(fig=c(0.1,0.7,0.7,0.9),new=TRUE)
      plot(denX_tot$x, denX_tot$y, main="MA plot", cex.main=1.5, xlim=xrangeMA, ylim=c(0,denXMax), type='l', axes=FALSE, col=ColAll, xlab="", 
           ylab="", lwd=2)
      lines(denX_PC, col=ColP, lwd=2)
      lines(denX_NC, col=ColN, lwd=2)
      rgb.val<- col2rgb(ColAll)
      polygon(denX_tot, col=rgb(rgb.val[1]/255,rgb.val[2]/255,rgb.val[3]/255,alpha=0.3), lwd=0.1)
      rgb.val<- col2rgb(ColP)
      polygon(denX_PC, col=rgb(rgb.val[1]/255,rgb.val[2]/255,rgb.val[3]/255,alpha=0.3), lwd=0.1)
      rgb.val<- col2rgb(ColN)
      polygon(denX_NC, col=rgb(rgb.val[1]/255,rgb.val[2]/255,rgb.val[3]/255,alpha=0.3), lwd=0.1)
      if (Gene=="Hitlist" && nrow(df_hits_total)>1){
        denX_hits_total<- density(df_hits_total$logBaseMeanA, from=xrangeMA[1], to=xrangeMA[2], na.rm=T)
        denX_hits_total$y[1]<- 0
        denX_hits_total$y[length(denX_hits_total$y)]<- 0
        lines(denX_hits_total, col=ColH, lwd=2)
        rgb.val<- col2rgb(ColH)
        polygon(denX_hits_total, col=rgb(rgb.val[1]/255,rgb.val[2]/255,rgb.val[3]/255,alpha=0.3), lwd=0.1)
      }
      if (nrow(df_GOI)>1){
        par(new=TRUE)
        denX_GOI<- density(df_GOI$logBaseMeanA, from=xrangeMA[1], to=xrangeMA[2], na.rm=T)
        denX_GOI$y[1]<- 0
        denX_GOI$y[length(denX_GOI$y)]<- 0
        lines(denX_GOI, col=ColH, lwd=2)
        rgb.val<- col2rgb(ColH)
        polygon(denX_GOI, col=rgb(rgb.val[1]/255,rgb.val[2]/255,rgb.val[3]/255,alpha=0.3), lwd=0.1)
      }
    }
    dev.off()
    
    if (mageck==0){
      if (r==1) {
        df_mageck<- read.csv(file=Filename, sep=",", header=TRUE, stringsAsFactors = FALSE)
        df_mageck$Sequence<-NULL
        df_mageck$X<-NULL
        dir.create(paste0(dirname,"/MAGeCK"))
        setwd(paste0(Workdirectory,"/", dirname, "/MAGeCK"))
        write.csv(df_mageck, "MAGeCKCountTable.csv", row.names = FALSE, quote = F)
        system(paste0("mageck test -k MAGeCKCountTable.csv -t ", paste0(if (T1_Rep1==0) {T1Rep1}, if (T1_Rep1==0 & (T1_Rep2==0|T1_Rep3==0)) {","},
                                                                        if (T1_Rep2==0) {T1Rep2}, if (T1_Rep2==0 & T1_Rep3==0) {","},
                                                                        if (T1_Rep3==0) {T1Rep3}),
                      " -c ", paste0(if (T0_Rep1==0) {T0Rep1}, if (T0_Rep1==0 & (T0_Rep2==0|T0_Rep3==0)) {","},
                                     if (T0_Rep2==0) {T0Rep2}, if (T0_Rep2==0 & T0_Rep3==0) {","},
                                     if (T0_Rep3==0) {T0Rep3}),
                      " -n MAGeCK_T0vsT1"))
        source("MAGeCK_T0vsT1.R")
        setwd(Workdirectory)
      }
      if (r==2) {
        setwd(paste0(Workdirectory,"/", dirname, "/MAGeCK"))
        system(paste0("mageck test -k MAGeCKCountTable.csv -t ", paste0(if (T2_Rep1==0) {T2Rep1}, if (T2_Rep1==0 & (T2_Rep2==0|T2_Rep3==0)) {","},
                                                                        if (T2_Rep2==0) {T2Rep2}, if (T2_Rep2==0 & T2_Rep3==0) {","},
                                                                        if (T2_Rep3==0) {T2Rep3}),
                      " -c ", paste0(if (T0_Rep1==0) {T0Rep1}, if (T0_Rep1==0 & (T0_Rep2==0|T0_Rep3==0)) {","},
                                     if (T0_Rep2==0) {T0Rep2}, if (T0_Rep2==0 & T0_Rep3==0) {","},
                                     if (T0_Rep3==0) {T0Rep3}),
                      " -n MAGeCK_T0vsT2"))
        source("MAGeCK_T0vsT2.R")
        setwd(Workdirectory)
      }
      if (r==3) {
        setwd(paste0(Workdirectory,"/", dirname, "/MAGeCK"))
        system(paste0("mageck test -k MAGeCKCountTable.csv -t ", paste0(if (T2_Rep1==0) {T2Rep1}, if (T2_Rep1==0 & (T2_Rep2==0|T2_Rep3==0)) {","},
                                                                        if (T2_Rep2==0) {T2Rep2}, if (T2_Rep2==0 & T2_Rep3==0) {","},
                                                                        if (T2_Rep3==0) {T2Rep3}),
                      " -c ", paste0(if (T1_Rep1==0) {T1Rep1}, if (T1_Rep1==0 & (T1_Rep2==0|T1_Rep3==0)) {","},
                                     if (T1_Rep2==0) {T1Rep2}, if (T1_Rep2==0 & T1_Rep3==0) {","},
                                     if (T1_Rep3==0) {T1Rep3}),
                      " -n MAGeCK_T1vsT2"))
        source("MAGeCK_T1vsT2.R")
        setwd(Workdirectory)
      }
    }
  }
}