## Cell-cell communication from bulk RNAseq data using communication scores:
# A heatmap of communications scores is created (expression level ligand on cell A x expression level receptor in cell B).
# - Make count tables with the GeneExpressionRNAseqHisat2.R script
# - This script also uses the LRtransmembrane.csv for ligand-receptor interactions (in which partnerB contains a transmembrane domain), 
#   and complexGeneSymbols.csv for "protein"-complexes.  Both files derived from CellPhoneDB v4 database. 
# - Adjust the settings and run the script
# Author: M.F.M. de Rooij PhD, Amsterdam UMC, Spaargaren Lab, 2022, info: m.f.derooij@amsterdamumc.nl
#######################################################################################################################################
#install.packages("circlize")
#install.packages("BiocManager")
#BiocManager::install(c("ComplexHeatmap","DESeq2"))
library(circlize)
library(ComplexHeatmap)
library(DESeq2)
#######################################################################################################################################
#                                     SETTINGS

wd<-"H:/BioWin/CCC/"
setwd(wd)

# Autocrine signaling: 0 = No (paracrine), 1 = Yes
autocrine<- 0

# Producer cells (ligands)
a<-read.csv("RNAseqCountTableTPMstr.csv", stringsAsFactors = F)
# Which cell-columns
colAa<-c(3:ncol(a))
# General name producer cells
nameA<- "Stroma"

# Receiver cells (receptors) (For autocrine signaling, use the same count table as producer cells)
b<-read.csv("RNAseqCountTableTPMbcl.csv", stringsAsFactors = F)
# Which cell-columns
colBb<-c(3:ncol(b))
# General name receiver cells
nameB<- "Lymphoma"

## Merge different count tables
# b1<- read.csv("RNAseqCountTableTPMbcl.csv", stringsAsFactors = F)
# b2<- read.csv("RNAseqCountTableTPMwm.csv", stringsAsFactors = F)
# b1<-b1[,c(1:2,9,27)]
# b2<-b2[,c(1,9,12)]
# b<- merge(b1,b2, by="ensembl_gene_id")

# LRpairs
LigandReceptorFile<- "LRtransmembrane.csv"

#minimal expression (nTPX)
minExp<- 10
# Minimal present (# of cell types)
minCell<- 1
#######################################################################################################################################

# Load RNAseq data
colnames(a)[-1:-2]<-gsub("_","",colnames(a)[-1:-2])
colnames(b)[-1:-2]<-gsub("_","",colnames(b)[-1:-2])
a<-a[,c(1,2,colAa)]
colnames(a)[-1:-2]<-paste0(colnames(a)[-1:-2],"producerCell")
b<-b[,c(1,2,colBb)]
colnames(b)[-1:-2]<-paste0(colnames(b)[-1:-2],"receiverCell")
# Merge and normalize RNAseq counts
c<-merge(a,b,by='ensembl_gene_id')
c$hgnc_symbol.y<-NULL
colnames(c)[2]<-"hgnc_symbol"
sf<- estimateSizeFactorsForMatrix(c[,-1:-2])
d<- c
d[,-1:-2]<- t(t(c[,-1:-2])/sf)

# Load ligand-receptor pairs
LRpairs<- read.csv(LigandReceptorFile)

# Calculate communication scores
data<-list()
colA <- 2+1:length(colAa)
colB <- 2+length(colAa)+1:length(colBb)
n<-1
for (i in colA){
  for (j in colB){
    type1<-colnames(d)[i]
    type2<-colnames(d)[j]

    if (autocrine==0 | (autocrine==1 & as.character(strsplit(type1,"producerCell"))==as.character(strsplit(type2, "receiverCell")))){
      Type1<-d[,c(2,i)]
      colnames(Type1)[2]<-"Type1"
      Type1ligands<- Type1[Type1$hgnc_symbol %in% LRpairs$partner_a,]
      Type1ligands<- Type1ligands[Type1ligands$Type1>=minExp,]
      colnames(Type1ligands)[1]<-"partner_a"
      
      Type2<-d[,c(2,j)]
      colnames(Type2)[2]<-"Type2"
      Type2receptors<-Type2[Type2$hgnc_symbol %in% LRpairs$partner_b,]
      Type2receptors<- Type2receptors[Type2receptors$Type2>=minExp,] 
      colnames(Type2receptors)[1]<-"partner_b"
      
      LRCom<- LRpairs
      LRCom<- merge(LRCom,Type1ligands, by="partner_a")
      LRCom<- merge(LRCom,Type2receptors, by="partner_b")
      LRCom$com<-LRCom$Type1 * LRCom$Type2
      LRCom<-LRCom[order(LRCom$com,decreasing = T),]
      LRCom<-LRCom[,c("partner_a", "Type1", "partner_b", "Type2", "com")]
      colnames(LRCom)<- c("partner_a", type1, "partner_b", type2, "ComScore")
      LRCom$ComScore<-round(LRCom$ComScore,0)
      data[[n]]<-LRCom
      n<-n+1
    }
  }
}

# Merge all communication scores
df_heatmap<- LRpairs
df_heatmap$pair<- paste(df_heatmap$partner_a, df_heatmap$partner_b,sep=" - ")
for (i in 1:length(data)){
  temp<- data[[i]]
  temp$pair<- paste(temp$partner_a, temp$partner_b,sep=" - ")
  df_heatmap<-merge(df_heatmap,temp[5:6], by="pair", all=T)
  colnames(df_heatmap)[ncol(df_heatmap)]<-paste(colnames(temp)[c(2,4)], collapse="_")
}
df_heatmap[is.na(df_heatmap)]<-0
rownames(df_heatmap)<- df_heatmap$pair
df_heatmap2<-df_heatmap
df_heatmap<-df_heatmap[-1:-3]

# Include Complexes
complexGeneSymbols<-read.csv("complexGeneSymbols.csv", stringsAsFactors = F)
complexes<-NULL
for(i in 2:ncol(complexGeneSymbols)) {
  complexes <- c(complexes,complexGeneSymbols[,i])
}
complexes <- unique(complexes) 
complexes<-complexes[complexes!=""]
df_heatmapComplex<-unique(rbind(df_heatmap2[df_heatmap2$partner_a %in% complexes,],df_heatmap2[df_heatmap2$partner_b %in% complexes,]))
df_heatmapInclComplex<-df_heatmap2

# Complex in partner a
for (i in 1:nrow(complexGeneSymbols)){
  temp<- df_heatmapComplex[df_heatmapComplex$partner_a %in% complexGeneSymbols[i,2:ncol(complexGeneSymbols)],]
  if (all(complexGeneSymbols[i,2:ncol(complexGeneSymbols)][complexGeneSymbols[i,2:ncol(complexGeneSymbols)]!=""] %in% temp$partner_a)){
    for (pb in unique(temp$partner_b)){
      temp2<-temp[temp$partner_b==pb,]
      if (all(complexGeneSymbols[i,2:ncol(complexGeneSymbols)][complexGeneSymbols[i,2:ncol(complexGeneSymbols)]!=""] %in% temp2$partner_a)){
        temp3<-temp2[temp2$partner_a %in% complexGeneSymbols[i,2:ncol(complexGeneSymbols)][complexGeneSymbols[i,2:ncol(complexGeneSymbols)]!=""],]
        ligandComplex<-paste0(complexGeneSymbols$complex_name[i]," - ",pb)
        df_heatmapInclComplex<-rbind(df_heatmapInclComplex,data.frame(t(c(pair=ligandComplex,
          partner_a=complexGeneSymbols$complex_name[i],partner_b=pb,apply(temp3[-1:-3],2,min))), row.names=ligandComplex))
      }
    }
  }
}
df_heatmapInclComplex<-df_heatmapInclComplex[!df_heatmapInclComplex$partner_a %in% complexes,]

# Complex in partner b
for (i in 1:nrow(complexGeneSymbols)){
  temp<- df_heatmapInclComplex[df_heatmapInclComplex$partner_b %in% complexGeneSymbols[i,2:ncol(complexGeneSymbols)],]
  if (all(complexGeneSymbols[i,2:ncol(complexGeneSymbols)][complexGeneSymbols[i,2:ncol(complexGeneSymbols)]!=""] %in% temp$partner_b)){
    for (pa in unique(temp$partner_a)){
      #pa<-unique(temp$partner_a)[6]
      temp2<-temp[temp$partner_a==pa,]
      if (all(complexGeneSymbols[i,2:ncol(complexGeneSymbols)][complexGeneSymbols[i,2:ncol(complexGeneSymbols)]!=""] %in% temp2$partner_b)){
        temp3<-temp2[temp2$partner_b %in% complexGeneSymbols[i,2:ncol(complexGeneSymbols)][complexGeneSymbols[i,2:ncol(complexGeneSymbols)]!=""],]
        receptorComplex<-paste0(pa," - ",complexGeneSymbols$complex_name[i])
        df_heatmapInclComplex<-rbind(df_heatmapInclComplex,data.frame(t(c(pair=receptorComplex,
          partner_a=pa,partner_b=complexGeneSymbols$complex_name[i],apply(temp3[-1:-3],2,min))), row.names=receptorComplex))
        
      }
    }
  }
}
df_heatmapInclComplex<-df_heatmapInclComplex[!df_heatmapInclComplex$partner_b %in% complexes,]
df_heatmapInclComplex<-data.frame(sapply(df_heatmapInclComplex[-1:-3], function(x) as.numeric(as.character(x))), row.names=rownames(df_heatmapInclComplex))

# Order tables
order<- NULL
if (autocrine==0){
  for (c in colnames(a)[-1:-2]){
    order<- c(order,paste(c, colnames(b)[-1:-2], sep="_"))
  }
} else {
  i<-3
  for (c in colnames(a)[-1:-2]){
    order<- c(order,paste(c, colnames(b)[i], sep="_"))
    i<-i+1
  } 
}
df_heatmap<- df_heatmap[match(order,colnames(df_heatmap))]
df_heatmapInclComplex<- df_heatmapInclComplex[match(order,colnames(df_heatmapInclComplex))]
colnamesDf<- colnames(df_heatmap)
if (ncol(df_heatmap)>1){
  df_heatmap<- df_heatmap[rowSums(df_heatmap)!=0,]
  df_heatmapInclComplex<- df_heatmapInclComplex[rowSums(df_heatmapInclComplex)!=0,]
} else {
  colnamesDf<- colnames(df_heatmap)
  df_heatmap<- data.frame(df_heatmap[df_heatmap[[1]]!=0,], row.names= rownames(df_heatmap)[df_heatmap[[1]]!=0])
  df_heatmapInclComplex<- data.frame(df_heatmapInclComplex[df_heatmapInclComplex[[1]]!=0,], row.names= rownames(df_heatmapInclComplex)[df_heatmapInclComplex[[1]]!=0])
  colnames(df_heatmap)<-colnamesDf
  colnames(df_heatmapInclComplex)<-colnamesDf
}
#write.csv(df_heatmap, "df_heatmap.csv")
#write.csv(df_heatmapInclComplex, "df_heatmapComplex.csv")

# Generate heatmaps
if (ncol(df_heatmap)>1){
  df_heatmap<- df_heatmap[rowSums(df_heatmap==0)<=(ncol(df_heatmap)-minCell),]
  df_heatmapInclComplex<- df_heatmapInclComplex[rowSums(df_heatmapInclComplex==0)<=(ncol(df_heatmapInclComplex)-minCell),]
}
partner_a_cell <- as.numeric(factor(unlist(lapply(strsplit(as.character(colnames(df_heatmap)),"_"), "[", 1)),levels=colnames(a)[-1:-2]))
partner_b_cell <- as.numeric(factor(unlist(lapply(strsplit(as.character(colnames(df_heatmap)),"_"), "[", 2)),levels=colnames(b)[-1:-2]))
if (length(colAa)>1){
  col_funProd = colorRamp2(c(1,max(partner_a_cell)), c("darkolivegreen1","cornflowerblue"))
  } else{
    col_funProd = colorRamp2(c(1,2), c("darkolivegreen1","cornflowerblue"))
}
if (length(colBb)>1){
  col_funRec = colorRamp2(c(1,max(partner_b_cell)), c("khaki","magenta1"))
  } else{
    col_funRec = colorRamp2(c(1,2), c("khaki","magenta1"))
}
col<- list(`producer cell`=col_funProd, `receiver cell`=col_funRec)
ha = HeatmapAnnotation(`producer cell`=partner_a_cell, `receiver cell`=partner_b_cell, col=col, gp = gpar(col = "black"), border=T, show_legend=F)

col_fun = colorRamp2(c(0,2,6), c("white","blue","red"))
matlist<-list(df_heatmapInclComplex, df_heatmap)
pdf(paste0("CCC",strsplit(strsplit(LigandReceptorFile,".csv")[[1]][1], "LR")[[1]][2],minExp,".",minCell,nameA,"-",nameB,".pdf"), width=10+ncol(df_heatmap)/10, height=10+nrow(df_heatmap)/10)
    for (mat in matlist){
      mat<-as.data.frame(mat)
      colnames(mat)<- gsub("_", " - ",colnames(mat))
      colnames(mat)<- gsub("producerCell", "",colnames(mat))
      colnames(mat)<- gsub("receiverCell", "",colnames(mat))
      ht<-Heatmap(as.matrix((log(mat+1)/log(10))), cluster_columns = F, show_row_dend = F, col=col_fun,
                    row_names_side = "left", column_names_side = "top", 
              heatmap_legend_param = list(title = "Log10 communication score", 
                                      title_position = "leftcenter-rot",
                                      title_gp = gpar(fontsize = 15),
                                      legend_height = unit(25, "cm"),
                                      grid_width = unit(1, "cm"),
                                      at = 2:6,
                                      labels_gp = gpar(fontsize = 15)
                                    ),
          top_annotation = if (ncol(df_heatmap)>1){ha} else{ha[1]},
          rect_gp = gpar(col = "white", lwd = 1)
          )
      draw(ht, padding = unit(c(10, 50, 10, 10), "mm"))
  }
dev.off()

############################################################################################################
# # Generate LRlists from CellPhoneDB database
# #
# # Linux: pip install -U cellphonedb
# #        cellphonedb database download
# 
# setwd("H:/BioWin/CCC/LRlist/")
# 
# complex<-read.csv("import/complex_input.csv", stringsAsFactors = F)
# gene<-read.csv("import/gene_input.csv", stringsAsFactors = F)
# interaction<-read.csv("import/interaction_input.csv", stringsAsFactors = F)
# protein<-read.csv("import/protein_input.csv", stringsAsFactors = F)
# 
# LR<- interaction[,1:2]
# for (c in complex$complex_name){
#   temp1<-LR[LR$partner_a %in% c,]
#   if (nrow(temp1)>0){
#     temp2<-complex[complex$complex_name %in% c,]
#     temp2<- temp2[,grep("uniprot",colnames(temp2))]
#     temp2<-as.character(temp2)
#     temp2<- temp2[temp2!="NA"]
#     temp2<- temp2[temp2!=""]
#     for (row in 1:nrow(temp1)){
#       for (i in temp2){
#         LR<-rbind(LR,c(i,temp1$partner_b[row]))
#       }
#     }
#     LR<-LR[!LR$partner_a %in% c,]
#   }
# }
# for (c in complex$complex_name){
#   temp1<-LR[LR$partner_b %in% c,]
#   if (nrow(temp1)>0){
#     temp2<-complex[complex$complex_name %in% c,]
#     temp2<- temp2[,grep("uniprot",colnames(temp2))]
#     temp2<-as.character(temp2)
#     temp2<- temp2[temp2!="NA"]
#     temp2<- temp2[temp2!=""]
#     for (row in 1:nrow(temp1)){
#       for (i in temp2){
#         LR<-rbind(LR,c(temp1$partner_a[row],i))
#       }
#     }
#     LR<-LR[!LR$partner_b %in% c,]
#   }
# }
# 
# gs<-unique(gene[2:3])
# gs<-gs[gs$hgnc_symbol!="",]
# 
# for (i in 1:nrow(LR)){
#   gsi<- gs[gs$uniprot==LR[i,1],"hgnc_symbol"]
#   if (length(gsi)==1){
#     LR[i,1]<-gsi
#   } else{
#     LR[i,1]<-paste(gsi,collapse=":")
#   }
#   gsi<- gs[gs$uniprot==LR[i,2],"hgnc_symbol"]
#   if (length(gsi)==1){
#     LR[i,2]<-gsi
#   } else{
#     LR[i,2]<-paste(gsi,collapse=":")
#   }
# }
# 
# doubles<- grep(":", LR$partner_a)
# for (i in doubles){
#   temp1<- LR[i,]
#   LR<-rbind(LR,c(strsplit(temp1[1,1],":")[[1]][1],temp1[1,2]))
#   LR<-rbind(LR,c(strsplit(temp1[1,1],":")[[1]][2],temp1[1,2]))
# }
# LR<-LR[-doubles,]
# doubles<- grep(":", LR$partner_b)
# for (i in doubles){
#   temp1<- LR[i,]
#   LR<-rbind(LR,c(temp1[1,1],strsplit(temp1[1,2],":")[[1]][1]))
#   LR<-rbind(LR,c(temp1[1,1],strsplit(temp1[1,2],":")[[1]][2]))
# }  
# LR<-LR[-doubles,]
# LR<-unique(LR)  
# 
# 
# LR2<-LR[,c(2,1)]
# colnames(LR2)<-colnames(LR)
# LR<-rbind(LR,LR2)
# LR<-unique(LR)
# write.csv(LR,"LRall.csv", row.names = F, quote = F)
# #########################################################
# # select membrane proteins in partner_b
# 
# transmembrane<-protein
# for (i in 1:nrow(transmembrane)){
#   gsi<- gs[gs$uniprot==transmembrane[i,1],"hgnc_symbol"]
#   if (length(gsi)==1){
#     transmembrane[i,1]<-gsi
#   } else{
#     transmembrane[i,1]<-paste(gsi,collapse=":")
#   }
# }
# colnames(transmembrane)[1]<-"hgnc_symbol"
# doubles<- grep(":", transmembrane$hgnc_symbol)
# for (i in doubles){
#   temp1<- transmembrane[i,]
#   tempG<- temp1[1,1]
#   temp1[1,1]<- strsplit(tempG,":")[[1]][1]
#   transmembrane<-rbind(transmembrane,temp1)
#   temp1[1,1]<- strsplit(tempG,":")[[1]][2]
#   transmembrane<-rbind(transmembrane,temp1)
# }
# transmembrane<-transmembrane[-doubles,]
# transmembrane<-unique(transmembrane) 
# transmembranesOnly<-transmembrane[transmembrane$transmembrane=="True",]
# LR$transmembraneA<-0
# LR$transmembraneB<-0
# LR[LR$partner_a %in% transmembranesOnly$hgnc_symbol,"transmembraneA"]<-1
# LR[LR$partner_b %in% transmembranesOnly$hgnc_symbol,"transmembraneB"]<-1
# 
# LRtransmembraneB<- LR[LR$transmembraneB==1,]
# #LRtransmembraneB<- LRtransmembraneB[LRtransmembraneB$transmembraneA==1,]
# LRtransmembraneB<- LRtransmembraneB[1:2]
# write.csv(LRtransmembraneB,"LRtransmembrane.csv", row.names = F, quote = F)
# #################################################################
# # Complexes
# complexGeneSymbols<-complex
# gs<-unique(gene[2:3])
# gs<-gs[gs$hgnc_symbol!="",]
# 
# for (j in 2:5){
#   for (i in 1:nrow(complexGeneSymbols)){
#     gsi<- gs[gs$uniprot==complexGeneSymbols[i,j],"hgnc_symbol"]
#     if (length(gsi)==1){
#       complexGeneSymbols[i,j]<-gsi
#     } else{
#       complexGeneSymbols[i,j]<-paste(gsi,collapse=":")
#     }
#   }
# }
# complexGeneSymbols<-complexGeneSymbols[1:5]
# colnames(complexGeneSymbols)<- gsub("uniprot","GeneSymbol",colnames(complexGeneSymbols))
# 
# complexGeneSymbols$complex_name<- paste0(complexGeneSymbols$complex_name, " complex")
# lipids<-grep("by", complexGeneSymbols$complex_name)
# complexGeneSymbols$complex_name[lipids]<- gsub(" complex","",complexGeneSymbols$complex_name[lipids])
# complexGeneSymbols$complex_name[lipids]<- gsub("_by"," by ",complexGeneSymbols$complex_name[lipids])
# complexGeneSymbols$complex_name<- gsub("_receptor"," receptor ",complexGeneSymbols$complex_name)
# complexGeneSymbols$complex_name<- gsub("_","-",complexGeneSymbols$complex_name)
# complexGeneSymbols$complex_name<- gsub(":","-",complexGeneSymbols$complex_name)
# complexGeneSymbols$complex_name<- gsub("  "," ",complexGeneSymbols$complex_name)
# complexGeneSymbols$complex_name<- gsub("complex complex","complex",complexGeneSymbols$complex_name)
# complexGeneSymbols$complex_name<- gsub("-complex"," complex",complexGeneSymbols$complex_name)
# complexGeneSymbols$complex_name<- gsub(" -","-",complexGeneSymbols$complex_name)
# write.csv(complexGeneSymbols,"complexGeneSymbols.csv", row.names = F,quote = F)
