# Use the GeneExpressionRNAseqDESeq2.R to generate a table with rlog normalized counts
# Use the GSEA_RNAseq.R script to run GSEA analysis
# Put this script in the DESeq2 folder, adjust the settings, and run the script in Rstudio
# This script makes heatmaps of the rlog values of the Hallmarks or KEGG pathway genes
# Author: M.F.M. de Rooij PhD, Amsterdam UMC, Spaargaren Lab, 2026, info: m.f.derooij@amsterdamumc.nl
###########################################################################################################################
#install.packages("BiocManager")
#BiocManager::install(c("ComplexHeatmap", "clusterProfiler","msigdbr","org.Hs.eg.db","HGNChelper"))
#install.packages("readxl", "circlize")
library("readxl")
library("ComplexHeatmap")
library("circlize")
library("clusterProfiler")
library("msigdbr")
library("org.Hs.eg.db")
library("HGNChelper")
###########################################################################################################################
#                                     SETTINGS

# Add this file to the DESeq2 folder
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# KEGG of Hallmarks: 0 = Hallmarks, 1 = KEGG
PathwaySet<- 0

if (PathwaySet==0){
  PathSet<- "GSEA_MSIGDB_H"
}
if (PathwaySet==1){
  PathSet<- "GSEA_KEGG"
}
# Cutoff of GSEA's Padjust (For all pathways, use 1)
PadjCutoff <- 0.1

# Load RNAseq data, sample design, and GSEA data
dataAll<- read.csv("../NormalizedCounts_rlog.csv")
RNAseqDesign<- read.csv("../RNAseqDesign.csv")
if (PadjCutoff < 1){
  GSEA<- read_excel("GSEA_STRING/GSEA_ClusterProfiler_DESeq2_RNAsehjsdlksdhkq_L363_shSIK3_cl3vsL363_shSIK3_cl3_dox.csv.xlsx", sheet=PathSet)
}

###########################################################################################################################
dir.create("GSEA_STRING/Heatmaps")

data<- dataAll
data[is.na(data)]<- 0

# Adjust gene symbols to recent hugo symbols
Hugo <- checkGeneSymbols(data$hgnc_symbol, unmapped.as.na=T)
data$Hugo <- Hugo$Suggested.Symbol

# Remove duplicate Hugo symbols, take the one with highest average expression if the original is not there 
duplicates<- data.frame(h=data$Hugo, n=1)
duplicates<- aggregate(duplicates$n, by=list(duplicates$h), FUN=sum)
duplicates<- duplicates[duplicates$x>1,]
dataDupl<- data[data$Hugo %in% duplicates$Group.1,]

data<- data[!data$Hugo %in% duplicates$Group.1,]
for (g in unique(dataDupl$Hugo)){
  temp<- dataDupl[dataDupl$Hugo==g,]
  temp$counts<- rowSums(temp[c(-1:-2,-(ncol(temp)-1):-ncol(temp))])
  idmax<- temp$ensembl_gene_id[temp$count==max(temp$count)][1]
  
  if (sum(temp$hgnc_symbol==g)==1){
    data<- rbind(data, dataDupl[dataDupl$hgnc_symbol==g,])
  } else{
    data<- rbind(data, dataDupl[dataDupl$ensembl_gene_id==idmax,])
  }
}

# ColIDs
colLabels<- paste0(RNAseqDesign$Group, "_",RNAseqDesign$Rep)  

Group <- as.numeric(factor(RNAseqDesign$Group,levels=unique(RNAseqDesign$Group)))
Rep <- as.numeric(factor(RNAseqDesign$Rep,levels=unique(RNAseqDesign$Rep)))
col_funGroup = colorRamp2(c(1,max(Group)), c("oldlace","dodgerblue3"))
col_funRep = colorRamp2(c(1,max(Rep)), c("wheat","indianred"))
col<- list(`Group`=col_funGroup, `Replicate`=col_funRep)
ha = HeatmapAnnotation(`Group`=Group, `Replicate`=Rep, col=col, gp = gpar(col = "black"), border=T, show_legend=F)

# Heatmap colors
col_fun = colorRamp2(c(-2,0,2), c("royalblue","snow","firebrick2"))

# Pathways
if (PadjCutoff < 1){
  GSEApathways<- GSEA$Description[GSEA$p.adjust < PadjCutoff]
}
if (PathwaySet==0){
  h_gene_sets <- msigdbr(species = "Homo sapiens", collection = "H")
  if (PadjCutoff==1){
    GSEApathways<- unique(h_gene_sets$gs_name)
  }
}
if (PathwaySet==1){
  kegg_data <- download_KEGG("hsa")
  pathway_genes <- kegg_data$KEGGPATHID2EXTID
  pathway_names <- kegg_data$KEGGPATHID2NAME
  if (PadjCutoff==1){
    GSEApathways<- unique(pathway_names$to)
  }
}

for (Pathway in GSEApathways){
  if (PathwaySet==0){
    pathwayGeneSymbols<- h_gene_sets$gene_symbol[h_gene_sets$gs_name==Pathway]
  }
  if (PathwaySet==1){
    pathway<- pathway_names$from[pathway_names$to==Pathway]
    pathwayGenes<- pathway_genes$to[pathway_genes$from==pathway]
    pathwayGeneSymbols <- bitr(pathwayGenes, fromType = "ENTREZID", toType = "SYMBOL", OrgDb = org.Hs.eg.db)
    pathwayGeneSymbols<- pathwayGeneSymbols$SYMBOL
  }
  df_heatmap<- data[data$Hugo %in% pathwayGeneSymbols,]
  rownames(df_heatmap)<- df_heatmap$Hugo
  df_heatmap$ensembl_gene_id<- NULL
  df_heatmap$hgnc_symbol<- NULL
  df_heatmap$Hugo<- NULL
  df_heatmap<-df_heatmap[,RNAseqDesign$Sample]
  
  mat_scaled = as.data.frame(t(scale(t(df_heatmap))))
  mat_scaled<- mat_scaled[!is.na(mat_scaled[[1]]),]
 
  pdf(paste0("GSEA_STRING/Heatmaps/Heatmap_", PathSet, "_", sub("/", "-",Pathway),".pdf"), width=2+ncol(mat_scaled)/3, height=5+nrow(mat_scaled)/8)
        ht<-Heatmap(as.matrix(mat_scaled), column_labels = colLabels, cluster_columns = F, show_row_dend = T, row_dend_side = "right", col = col_fun,
                      row_names_side = "left", column_names_side = "top", clustering_distance_rows = "euclidean", 
                heatmap_legend_param = list(title = "rlog normalized counts (z-scores)", 
                                        title_position = "leftcenter-rot",
                                        title_gp = gpar(fontsize = 12),
                                        legend_height = unit(5, "cm"),
                                        grid_width = unit(5, "mm"),
                                        at = -2:2,
                                        labels_gp = gpar(fontsize = 12)
                                      ),
            top_annotation = ha,
            #left_annotation = va,
            rect_gp = gpar(col = "black", lwd = 1)
            )
        draw(ht, padding = unit(c(10, 50, 10, 10), "mm"))
  dev.off()
}  
