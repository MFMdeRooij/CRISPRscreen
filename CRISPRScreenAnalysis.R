## CRISPR Screen analysis:
# - We used the screen design of Jastrzebski et al Methods Mol Biol 2016
# - Make count tables from FASTQ files with the Perl script FastqToCountTable.pl and a libraryX.csv file
# - Analyze the count tables with this script and the CRISPRScreenAnalysisLibraries.csv file
# Author: M.F.M. de Rooij PhD, Amsterdam UMC, Spaargaren Lab, 2019, info: m.f.derooij@amsterdamumc.nl
# A small part (aRRA) is adapted from Thomas Kuilman PhD, NKI Amsterdam, Peeper Lab (https://github.com/PeeperLab?tab=repositories)
# The random seeds are different in Windows and Linux R versions, so the p-values on gene level can be different
######################################################################################
## Install DESeq2
#install.packages("BiocManager")
#BiocManager::install(c("DESeq2", "pheatmap", "RColorBrewer"), dependencies=TRUE)
library("DESeq2")
library("pheatmap")
library("RColorBrewer")
######################################################################################
#                                     SETTINGS

# Workdirectory (folder in which the count tables are located, use always slash (/) instead of backslash)
Workdirectory <- "H:/BioWin/Screens/"

# Count table files: 0 = custom, 1 = all count tables (csv files) in the workdirectory
Data <- 0
# If custom, which count tables?
Filenames<- c("CountTable_test.fastq.gz.csv")

# Data selection (T0/T1/T2 - replicate 1/2/3) (If a replicate is not present, fill in 
# an existing one, and exclude it in the next option)
# Fill in the column names of the count table
Guide <- "sgRNA"
T0Rep1 <- "BC1"
T0Rep2 <- "BC2"
T0Rep3 <- "BC3"
T1Rep1 <- "BC7"
T1Rep2 <- "BC8"
T1Rep3 <- "BC9"
T2Rep1 <- "BC4"
T2Rep2 <- "BC5"
T2Rep3 <- "BC6"

# Exclude replicates: 0 = Include, 1 = Exclude (you cannot exclude T0 or T1 completely -> change your design)
T0_Rep1 <- 0
T0_Rep2 <- 0
T0_Rep3 <- 0
T1_Rep1 <- 0
T1_Rep2 <- 0
T1_Rep3 <- 0
T2_Rep1 <- 0
T2_Rep2 <- 0
T2_Rep3 <- 0

# Paired replicates: 0 = Paired, 1 = Unpaired
Paired <- 0

# Library: 0: In countTable (1st column = sgRNA, 2nd = Gene), 1 = Brunello-Kinome, 2 = Brunello-FullGenome, 3 = Avana-FullGenome, 4 = New (add to CRISPRScreenAnalysisLibraries.csv file)
Library <- 1

# Positive controls: 0 = Essential (Hart et al 2017 G3)
ControlsP <- 0
# Negative controls: 0 = Non-Essential (Hart et al 2014 Mol Sys Biol), 1 = Non-Targeting
ControlsN <- 0
# Even it is not a lethality screen, it is still worth to look how cell viability affect your screen (separation of Essential-Non-Essential)

# Statistics of Essential/Non-Essential genes (only for lethality screens): 0 = With, 1 Without
ConStat<- 1

# Use shrinkage of fold changes noninformative guides: 0 = Yes, 1 = No
Shrinkage<- 0

# Exclude noninformative guides for DESEq2-FDR & aRRA: 0 = Yes, 1 = No
ExclNonInf = 0

# Type of Screen: 0 = Drop out, 1 = Resistance (Just for ranking of the genes)
Type_of_Screen <- 0

# Minimal fold change of guides to be a hit: 1 = No minimal fold change,  >1: The minimal fold change
minimalFoldChange <- 1

# Number of permutation for aRRA (0 = 250 times the number of genes, 1 = Custom)
numPer <- 0
# If number of permutation is custom, how much? (The best p value form aRRA is 1/numOfPer)
numOfPer <- 1000

# MA plots of all genes: 0 = Yes, 1 = No, 2 = Top 10, 3 = Genes of interest
MA_all_genes<- 2
# Genes of interest
if (MA_all_genes==3){
  q <- function(...) {
    sapply(match.call()[-1], deparse)
  }
  interestingGenes <-q(
    # Which genes?  
    BTK, SYK, PIK3CA, PIK3CD, PIK3R1, LYN
  )
}

# Correlation plots of all samples: 0 = Yes, 1 = No (Only recommended when the replicates are messed up) 
pairs<- 1

# MAGeCK: 0 = Yes, 1 = No (this does not work in Windows)
mageck<- 1
# MAGeCK should be installed in linux with all its dependencies, see website: https://sourceforge.net/p/mageck/wiki/Home/
######################################################################################
setwd(Workdirectory)
if (Data==0){
  Filenames <- Filenames
} 
if (Data==1){
  Filenames <- Sys.glob("*.csv")
  Filenames<- Filenames["CRISPRScreenAnalysisLibraries.csv" != Filenames]
}

for (Filename in Filenames) {
  # Make a data folder
  dirname2<- paste0(Filename,Sys.time())
  dirname1<-gsub("[[:punct:]]", "", dirname2) 
  dirname<- gsub("\\s", "", dirname1) 
  dir.create(dirname)
  
  # Read count table
  df_raw <- read.csv(file=Filename, sep=",", header=TRUE, stringsAsFactors = FALSE)
  
  # Select replicates
  if(T2_Rep1+T2_Rep2+T2_Rep3 <3) {
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
  
  # Count table for DESEq2
  counts <- df_sel[,c(
    if (T0_Rep1==0) {2},
    if (T0_Rep2==0) {3},
    if (T0_Rep3==0) {4},
    if (T1_Rep1==0) {5},
    if (T1_Rep2==0) {6},
    if (T1_Rep3==0) {7},
    if (T2_Rep1==0) {8},
    if (T2_Rep2==0) {9},
    if (T2_Rep3==0) {10})]
  rownames(counts) <- df_sel[,1]
  
  # Sample metadata
  time <- c(if (T0_Rep1==0) {"T0"},
            if (T0_Rep2==0) {"T0"},
            if (T0_Rep3==0) {"T0"},
            if (T1_Rep1==0) {"T1"},
            if (T1_Rep2==0) {"T1"},
            if (T1_Rep3==0) {"T1"},
            if (T2_Rep1==0) {"T2"},
            if (T2_Rep2==0) {"T2"},
            if (T2_Rep3==0) {"T2"})
  
  rep <- c(if (T0_Rep1==0) {"R1"},
           if (T0_Rep2==0) {"R2"},
           if (T0_Rep3==0) {"R3"},
           if (T1_Rep1==0) {"R1"},
           if (T1_Rep2==0) {"R2"},
           if (T1_Rep3==0) {"R3"},
           if (T2_Rep1==0) {"R1"},
           if (T2_Rep2==0) {"R2"},
           if (T2_Rep3==0) {"R3"})
  
  df_colData <- data.frame(Time=time, Rep=rep) 
  df_colData$Time <- as.factor(df_colData$Time)
  df_colData$Rep <- as.factor(df_colData$Rep)
  
  # DESeq2 pipeline
  if (Paired==0) {
    dds <- DESeqDataSetFromMatrix(countData = counts, colData = df_colData, design = ~ Rep + Time)
  }
  if (Paired==1) {
    dds <- DESeqDataSetFromMatrix(countData = counts, colData = df_colData, design = ~ Time)
  }  
  dds$time <- relevel(dds$Time, "T0")
  dds <- DESeq(dds, fitType = 'local', betaPrior = TRUE)
  #fitType, default is parametric, for CRISPR screens local is better (parametric is in most cases not possible; betaPrior=TRUEis for shrinkage of noninformative fold changes)
  
  # Loop over time point/condition comparisons
  nr<-1
  if ("T2" %in% time) {nr=3}
  for (r in 1:nr) {
    if (ExclNonInf == 0) {
      if (r==1){
        res <- results(dds, contrast=c("Time","T1","T0"), addMLE=TRUE)
      }
      if (r==2){
        res <- results(dds, contrast=c("Time","T2","T0"), addMLE=TRUE)
      }
      if (r==3){
        res <- results(dds, contrast=c("Time","T2","T1"), addMLE=TRUE)
      }
    } else {
      if (ExclNonInf == 1) {
        # To don't have NA in padj: add independentFiltering=FALSE:
        if (r==1){
          res <- results(dds, contrast=c("Time","T1","T0"), addMLE=TRUE, independentFiltering=FALSE)
        }
        if (r==2){
          res <- results(dds, contrast=c("Time","T2","T0"), addMLE=TRUE, independentFiltering=FALSE)
        }
        if (r==3){
          res <- results(dds, contrast=c("Time","T2","T1"), addMLE=TRUE, independentFiltering=FALSE)
        }
      }
    } 
    
    # DESeq2 data table
    df_res <- as.data.frame(res)
    
    if (Shrinkage==1) {
      df_res$log2FoldChange <- df_res$lfcMLE
    }
    df_baseMeanPerLvl <- as.data.frame(sapply(levels(dds$Time), function(lvl) rowMeans(counts(dds, normalized=TRUE)
                                                                                       [,dds$time==lvl])))
    if (r==1){
      df_res$BaseMeanA <-round(df_baseMeanPerLvl$T0,0)
      df_res$logBaseMeanA <-log(df_baseMeanPerLvl$T0+1)/log(10)
      df_res$BaseMeanB <-round(df_baseMeanPerLvl$T1,0)
    } 
    if (r==2){
      df_res$BaseMeanA <-round(df_baseMeanPerLvl$T0,0)
      df_res$logBaseMeanA <-log(df_baseMeanPerLvl$T0+1)/log(10)
      df_res$BaseMeanB <-round(df_baseMeanPerLvl$T2,0)
    } 
    if (r==3){
      df_res$BaseMeanA <-round(df_baseMeanPerLvl$T1,0)
      df_res$logBaseMeanA <-log(df_baseMeanPerLvl$T1+1)/log(10)
      df_res$BaseMeanB <-round(df_baseMeanPerLvl$T2,0)
    } 
    
    # Guide IDs
    df_res$Guide<-rownames(df_res)
    df_res<- merge(df_Gene_ID, df_res, by='Guide', all.y=T)
    df_res$log2FoldChange[is.na(df_res$log2FoldChange)]<-0
    df_res$FoldChange<- round(2^df_res$log2FoldChange,3)
    
    # Exclude noninformative guides from aRRA
    df_res2<-df_res
    df_res2<-df_res2[!is.na(df_res2$padj),]
    
    # Control genes
    df_res$Type<- "x"
    df_res[df_res$GeneSymbol %in% df_Control[[PC]][nchar(df_Control[[PC]])>0],"Type"]<- "p"
    df_res[df_res$GeneSymbol %in% df_Control[[NC]][nchar(df_Control[[NC]])>0],"Type"]<- "n"
    df_PC<-df_res[df_res$Type=="p",]
    df_NC<-df_res[df_res$Type=="n",]
    
    # Determine hits
    df_hits_total <- df_res[which(df_res$padj < 0.1),]
    if (minimalFoldChange > 1) {
      df_hits_total <- df_hits_total[which(df_hits_total$FoldChange <= 1/minimalFoldChange | df_hits_total$FoldChange >= minimalFoldChange),]
    }
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
    
    # AlphaRRA (from guide to gene statistics)
    if (r==1){
      resDepleted <- results(dds, contrast=c("Time","T1","T0"), altHypothesis = "less", addMLE=TRUE)
      resEnriched <- results(dds, contrast=c("Time","T1","T0"), altHypothesis = "greater", addMLE=TRUE)
    }
    if (r==2){
      resDepleted <- results(dds, contrast=c("Time","T2","T0"), altHypothesis = "less", addMLE=TRUE)
      resEnriched <- results(dds, contrast=c("Time","T2","T0"), altHypothesis = "greater", addMLE=TRUE)
    }
    if (r==3){
      resDepleted <- results(dds, contrast=c("Time","T2","T1"), altHypothesis = "less", addMLE=TRUE)
      resEnriched <- results(dds, contrast=c("Time","T2","T1"), altHypothesis = "greater", addMLE=TRUE)
    }
    resDepleted <- resDepleted[order(rownames(resDepleted)), ]
    resEnriched <- resEnriched[order(rownames(resEnriched)), ]
    if(Shrinkage==0) {
      df_RRA <- data.frame(log2fc = resDepleted$log2FoldChange,
                           pvalueDepleted = resDepleted$pvalue,
                           pvalueEnriched = resEnriched$pvalue,
                           row.names = rownames(resDepleted))
    }
    if(Shrinkage==1) {
      df_RRA <- data.frame(log2fc = resDepleted$lfcMLE,
                           pvalueDepleted = resDepleted$pvalue,
                           pvalueEnriched = resEnriched$pvalue,
                           row.names = rownames(resDepleted))
    }
    df_RRA$Guide<-rownames(df_RRA)
    df_RRA<- merge(df_Gene_ID, df_RRA, by='Guide', all.y=T)
    df_RRA<- df_RRA[df_RRA$Guide %in% df_res2$Guide,]
    
    # Rank by fold change
    df_RRA$scoreDepleted <- rank(df_RRA$log2fc) / nrow(df_RRA)
    df_RRA$scoreEnriched <- rank(-df_RRA$log2fc) / nrow(df_RRA)
    
    # Apply alpha criterion based on pvalues
    df_RRA$scoreDepleted[df_RRA$pvalueDepleted > 0.25] <- 1
    df_RRA$scoreEnriched[df_RRA$pvalueEnriched > 0.25] <- 1
    
    # Perform RRA
    alphaBeta <- function(p.in) {
      p.in <- sort(p.in)
      n <- length(p.in)
      return(min(pbeta(p.in, 1:n, n - (1:n) + 1)))
    }
    
    # Calculate rho per gene
    df_RRA$rhoDepleted <- unsplit(sapply(split(df_RRA$scoreDepleted,
                                               df_RRA$GeneSymbol), alphaBeta), df_RRA$GeneSymbol)
    df_RRA$rhoEnriched <- unsplit(sapply(split(df_RRA$scoreEnriched,
                                               df_RRA$GeneSymbol), alphaBeta), df_RRA$GeneSymbol)
    
    # Make a null distribution and calculate pvalues per gene
    n.guides <- sort(unique(table(df_RRA$GeneSymbol)))
    makeRhoNull <- function(n, p, nperm) {
      sapply(1:nperm, function(x) {
        p.in <- sort.int(sample(p, n, replace = FALSE))
        alphaBeta(p.in)
      })
    }
    if (numPer==0){
      permutations <- 250*nrow(df_hits_A)
    }
    if (numPer==1){ 
      permutations <- numOfPer
    }
    
    # Depletion
    set.seed(12345)
    rho.nulls <- lapply(n.guides, makeRhoNull, df_RRA$scoreDepleted, permutations)
    names(rho.nulls) <- as.character(n.guides)
    # Prevent pvalues of 0
    for (i in 1:length(rho.nulls)){
      rho.nulls[[i]][[1]]<-0
    }
    df_RRA$pvalueDepleted <- unsplit(sapply(split(df_RRA$rhoDepleted, df_RRA$GeneSymbol), function(x) {
      mean(rho.nulls[[as.character(length(x))]] <= x[1])
    }), df_RRA$GeneSymbol)
    # Enrichment
    set.seed(12345)
    rho.nulls <- lapply(n.guides, makeRhoNull, df_RRA$scoreEnriched, permutations)
    names(rho.nulls) <- as.character(n.guides)
    # Prevent pvalues of 0
    for (i in 1:length(rho.nulls)){
      rho.nulls[[i]][[1]]<-0
    }
    df_RRA$pvalueEnriched <- unsplit(sapply(split(df_RRA$rhoEnriched,df_RRA$GeneSymbol), function(x) {
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
      df_genes<-df_geneRRA[order(df_geneRRA$rhoDepleted),]
    }
    if (Type_of_Screen==1) {
      df_genes<-df_geneRRA[order(df_geneRRA$rhoEnriched),]
    }  
    
    # Write guide and gene tables
    if (r==1){
      titel<-"T0vsT1"
    }
    if (r==2){
      titel<-"T0vsT2"
    }
    if (r==3){
      titel<-"T1vsT2"
    }
    # Guide Table
    df_res_print<- df_res[,c("GeneSymbol","Guide","Type","BaseMeanA","BaseMeanB","FoldChange","pvalue","padj")]
    df_res_print<- df_res_print[order(df_res_print$Type),]
    write.csv(df_res_print, paste0(dirname,"/DESeq2 ",titel," Guides.csv"), row.names = FALSE)
    # Gene Table
    df_genes_print<-df_genes[,c("GeneSymbol", "Type", "TotalGuides", "HitsDown","HitsUp","MinFoldChange","MedianFoldChange","MaxFoldChange", "rhoDepleted",	"pvalueDepleted",	"fdrDepleted",	"rhoEnriched",	"pvalueEnriched", "fdrEnriched")] 
    write.csv(df_genes_print, paste0(dirname,"/DESeq2 ",titel," Genes.csv"), row.names = FALSE)
    
    # Correlation plots
    df_normCounts <- as.data.frame(counts(dds, normalized=TRUE))
    xyrange<-c(0, max(df_res$BaseMeanA, na.rm=TRUE))
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
    df_normCounts<-merge(df_Gene_ID, df_normCounts, by='Guide', all.y=T)
    df_normCounts$Guide<-NULL
    df_normCountsP<- df_normCounts[df_normCounts$GeneSymbol %in% df_Control[[PC]][nchar(df_Control[[PC]])>0],]
    df_normCountsN<- df_normCounts[df_normCounts$GeneSymbol %in% df_Control[[NC]][nchar(df_Control[[NC]])>0],]
    df_normCounts$GeneSymbol<-NULL
    df_normCountsP$GeneSymbol<-NULL
    df_normCountsN$GeneSymbol<-NULL
    
    # Write Plots in PDF 
    pdf(paste0(dirname,"/DESeq2 ",titel," Plots.pdf"), width=10, height=10)
    
    # DESEq2 plots
    plotDispEsts(dds, main="Dispersion plot")
    
    rld<-rlog(dds, blind=FALSE)
    sampleDists<- dist(t(assay(rld)))
    sampleDistMatrix <- as.matrix(sampleDists)
    rownames(sampleDistMatrix) <- paste(rld$Time, rld$Rep)
    colnames(sampleDistMatrix) <- paste(rld$Time, rld$Rep)
    colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
    
    pheatmap(sampleDistMatrix,
             clustering_distance_rows=sampleDists,
             clustering_distance_cols=sampleDists,
             col=colors)
    
    print(plotPCA(rld, intgroup=c("Time", "Rep")))
    
    if (pairs==0) {
      pairs(df_sel[-1], cex=0.1)
      pairs(df_normCounts, cex=0.1)
    }
    
    # Home made plots
    par(mfrow=c(3,3))
    for (i in c(1:length(repx))) {
      e1 <- df_rep[i,1]
      e2 <- df_rep[i,2]
      reg1 <- lm(df_normCounts[,e2]~df_normCounts[,e1])
      coefs <- coef(reg1)
      intercept <- round(coefs[1], 3)
      slope <- round(coefs[2],3)
      r2 <- round(summary(reg1)$r.squared, 3)
      eq1 <- bquote(italic(y) == .(slope)*italic(x) + .(intercept) ~"("*r^2 == .(r2)*")")
      
      plot(df_normCounts[,e1], df_normCounts[,e2], type="p", col=3, cex=.1, main="Correlation plot", 
           xlab=df_rep[i,1], ylab=df_rep[i,2], xlim=xyrange, ylim=xyrange)
      points(df_normCountsP[,e1], df_normCountsP[,e2], type="p", col=2, cex=.1)
      points(df_normCountsN[,e1], df_normCountsN[,e2], type="p", col=4, cex=.1)
      abline(reg1, col=6)
      abline(0,1, col=1)
      text(x = xyrange[2]/2, y = xyrange[2]*0.95, labels= eq1)
    }
    
    if (MA_all_genes==2) {
      RRAdep<-df_geneRRA[order(df_geneRRA$rhoDepleted),]
      RRAenr<-df_geneRRA[order(df_geneRRA$rhoEnriched),]
      tophits<- c(RRAdep[1:10,1], RRAenr[1:10,1])
    }
    
    # MA plots
    xrange<-c(min(df_res$logBaseMeanA, na.rm=TRUE)-0.5, max(df_res$logBaseMeanA, na.rm=TRUE)+0.5)
    yrange<-c(min(df_res$log2FoldChange, na.rm=TRUE)-0.5, max(df_res$log2FoldChange, na.rm=TRUE)+0.5)
    
    den_tot<-density(df_res$log2FoldChange, from=min(yrange[1]), to=max(yrange[2]), na.rm=TRUE)
    den_PC<-density(df_PC$log2FoldChange, from=min(yrange[1]), to=max(yrange[2]), na.rm=TRUE)
    den_NC<-density(df_NC$log2FoldChange, from=min(yrange[1]), to=max(yrange[2]), na.rm=TRUE)
    denMax <- max(c(den_tot$y, den_PC$y, den_NC$y))
    
    if (ConStat==0) {
      # Robust Z prime (median and median absolute deviation)
      Mp<- median(df_PC$log2FoldChange,na.rm=TRUE)
      Mn<- median(df_NC$log2FoldChange, na.rm=TRUE)
      SDp<-mad(df_PC$log2FoldChange, na.rm=TRUE)
      SDn<-mad(df_NC$log2FoldChange, na.rm=TRUE)
      Zprime <- 1-((3*(SDp+SDn))/(abs(Mp-Mn))) 
      
      # Intersection postive and negative controls
      poi <- which(diff(den_PC$y > den_NC$y) != 0) 
      intersection<- den_PC$x[poi][den_PC$x[poi]<0][which.min(abs(den_PC$x[poi][den_PC$x[poi]<0]))]
      
      # F measure (with cutoff the intersection), harmonic mean of precision and recall
      truePos<- nrow(df_PC[!is.na(df_PC$padj) & df_PC$padj<0.1 & df_PC$log2FoldChange < (log(1/minimalFoldChange)/log(2)),])
      posPredict<- nrow(df_PC[!is.na(df_PC$padj) & df_PC$padj<0.1 & df_PC$log2FoldChange < (log(1/minimalFoldChange)/log(2)),])+nrow(df_NC[!is.na(df_NC$padj) & df_NC$padj<0.1 & df_NC$log2FoldChange < (log(1/minimalFoldChange)/log(2)),])*nrow(df_PC)/nrow(df_NC)
      posControls<- nrow(df_PC)
      precision<- truePos/posPredict
      recall<- truePos/posControls
      F1<- 2/(1/precision + 1/recall)
    }
    
    Gene_of_interest<-c(" ", "Hitlist", if (MA_all_genes==0) {df_hits_A$GeneSymbol}, if (MA_all_genes==2) {tophits}, if (MA_all_genes==3) {interestingGenes})
    for(Gene in Gene_of_interest){
      df_GOI<-df_res[df_res$GeneSymbol %in% Gene,]
      
      par(mar=c(4,5,2,1))
      par(fig=c(0.1,0.83,0.1,0.83))
      plot(df_res$logBaseMeanA, df_res$log2FoldChange, type="p", col=3, cex=.7, pch=16, xlab=~""*""^10*"Log BaseMeanA", 
           ylab=~""*""^2*"Log Fold Change", cex.lab=1.5, cex.axis=1.3, xlim=xrange, ylim=yrange, xaxp = c(0, 10, 10), 
           yaxp = c(-10, 10, 20))
      if (ConStat==0) {
        text(median(xrange), yrange[2],labels=paste0("(rZ': ", format(round(Zprime,2), nsmall=2), 
                                                     " NP50: ", format(round(intersection,2), nsmall=2), " F1: ", format(round(F1,2),nsmall = 2), 
                                                     " (Pre: ", format(round(precision,2),nsmall = 2), " Rec: ", format(round(recall,2),nsmall = 2), "))"))
      }
      points(df_PC$logBaseMeanA, df_PC$log2FoldChange, type="p", col=2, cex=1, pch=15)
      points(df_NC$logBaseMeanA, df_NC$log2FoldChange, type="p", col=4, cex=1, pch=17)
      if (Gene=="Hitlist"&& nrow(df_hits_total)>1){
        points(df_hits_total$logBaseMeanA, df_hits_total$log2FoldChange, type="p", col=1, cex=0.7, pch=19)
      }
      if (nrow(df_GOI)>1){
        points(df_GOI$logBaseMeanA, df_GOI$log2FoldChange, type="p", col=1, cex=1.5, pch=19)
      }
      legend(xrange[1],yrange[2],legend=c("Total", "Essential", if (ControlsN==0 | ControlsN==1){"Non-Essential"}, if (ControlsN==2)
      {"Non-Targeting"}, if (nchar(Gene)>1 && !is.na(Gene)) {Gene}), cex=1, pch=c(16,15,17,if (nchar(Gene)>1 && !is.na(Gene)) {19}), col=c(3,2,4,if (nchar(Gene)>1 && !is.na(Gene)) {1}))
      abline(median(df_res$log2FoldChange, na.rm=TRUE),0, col=1, lty=3, untf=TRUE)
      
      # Density plot fold change
      par(fig=c(0.75,1,0.1,0.83),new=TRUE)
      plot(den_tot$y, den_tot$x, ylim=range(yrange), xlim=(c(0,denMax)), type='l', axes=FALSE, col=3, xlab="", 
           ylab="", lwd=2)
      par(new=TRUE)
      lines(den_PC$y, den_PC$x, col=2, lwd=2)
      lines(den_NC$y, den_NC$x, col=4, lwd=2)
      if (Gene=="Hitlist" && nrow(df_hits_total)>1){
        den_hits_total<-density(df_hits_total$log2FoldChange, from=min(yrange[1]), to=max(yrange[2]), na.rm=TRUE)
        lines(den_hits_total$y, den_hits_total$x, col=1, lwd=2)
      }
      if (nrow(df_GOI)>1){
        par(new=TRUE)
        den_GOI<-density(df_GOI$log2FoldChange, from=min(yrange[1]), to=max(yrange[2]), na.rm=TRUE)
        lines(den_GOI$y, den_GOI$x, col=1, lwd=2)
      }
      # Density plot read count
      denX_tot<-density(df_res$logBaseMeanA, from=min(xrange[1]), to=max(xrange[2]), na.rm=TRUE)
      denX_PC<-density(df_PC$logBaseMeanA, from=min(xrange[1]), to=max(xrange[2]), na.rm=TRUE)
      denX_NC<-density(df_NC$logBaseMeanA, from=min(xrange[1]), to=max(xrange[2]), na.rm=TRUE)
      denXMax <- max(c(denX_tot$y, denX_PC$y, denX_NC$y))
      
      par(fig=c(0.1,0.83,0.75,1),new=TRUE)
      denX_tot<-density(df_res$logBaseMeanA, from=min(xrange[1]), to=max(xrange[2]), na.rm=TRUE)
      plot(denX_tot$x, denX_tot$y, xlim=range(xrange), ylim=c(0,denXMax), type='l', axes=FALSE, col=3, xlab="", 
           ylab="", lwd=2)
      lines(denX_PC$x, denX_PC$y, col=2, lwd=2)
      lines(denX_NC$x, denX_NC$y, col=4, lwd=2)
      if (Gene=="Hitlist" && nrow(df_hits_total)>1){
        par(new=TRUE)
        denX_hits_total<-density(df_hits_total$logBaseMeanA, from=min(xrange[1]), to=max(xrange[2]), na.rm=TRUE)
        lines(denX_hits_total$x, denX_hits_total$y, col=1, lwd=2)
      }
      if (nrow(df_GOI)>1){
        par(new=TRUE)
        denX_GOI<-density(df_GOI$logBaseMeanA, from=min(xrange[1]), to=max(xrange[2]), na.rm=TRUE)
        lines(denX_GOI$x, denX_GOI$y, col=1, lwd=2)
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
