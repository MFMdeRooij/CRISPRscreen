## OncoPrint:
# To summarize your CRISPR screen hits (or originally genetic aberrations) for multiple cell lines 
# Add your data to OncoPrintData.csv, and adjust all settings, and run the script in R studio
# Author: M.F.M. de Rooij PhD, Amsterdam UMC, Spaargaren Lab, 2019, info: m.f.derooij@amsterdamumc.nl
##################################################################################################################################
#                                                          SETTINGS

# Folder and data file
setwd("H://BioWin/OncoPrint")
matOriginal<-read.csv("OncoPrintData.csv", sep=",", header=T, row.names=1, stringsAsFactors = F)

# Tumortypes: (if less or more than 3 groups, adjust the SampleOrder lines)
tumortype<-c(rep("MCL",10), rep("DLBCL",12), rep("MM",14))


# Title:
column_title = "Vulnerabilities Cell Lines"

#Colors: (if less or more than 12 hits, adjust the alter_fun list)
# https://htmlcolorcodes.com/
col = c('POS1' = "#ff3333", 
        'POS2' = "#ff9933",
        'POS3' = "#ffec08",
        'POS4' = "#71ff08",
        'POS5' = "#08fff0",
        'POS6' = "#053efd",
        'POS7' = "#7205fd",
        'POS8' = "#f505fd",
        'POS9' = "#117509",
        'POS10' = "#970404",
        'POS11' = "#47adb9",
        'POS12' = "#9fca64",
        "MISSING" = "#f3f3f3")

backgroundColor <- "#c8c8c8"
# Size of the white areas between the blocks
width<-1
heigth<-1
##################################################################################################################################
#install.packages("BiocManager")
#BiocManager::install("ComplexHeatmap")
library(ComplexHeatmap)

mat<-matOriginal


percMut<-rowSums(matOriginal=="POS")/(rowSums(matOriginal=="POS")+rowSums(matOriginal=="NEG"))*100
round(percMut,1)

mat[mat=='NEG']<-""
mat[mat=='MISSING']<-""
for (i in 1:nrow(mat)){
  mat[i,][mat[i,]=="POS"] <- paste0("POS", i)
}
mat <- as.matrix(mat)

alter_fun = list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w-unit(width, "mm"), h-unit(heigth, "mm"),
              gp = gpar(fill = backgroundColor, col = NA))
  },
  POS1 = function(x, y, w, h) {
    grid.rect(x, y, w-unit(width, "mm"), h-unit(heigth, "mm"),
              gp = gpar(fill = col['POS1'], col = NA))
  },
  POS2 = function(x, y, w, h) {
    grid.rect(x, y, w-unit(width, "mm"), h-unit(heigth, "mm"),
              gp = gpar(fill = col['POS2'], col = NA))
  },
  POS3 = function(x, y, w, h) {
    grid.rect(x, y, w-unit(width, "mm"), h-unit(heigth, "mm"),
              gp = gpar(fill = col['POS3'], col = NA))
  },
  POS4 = function(x, y, w, h) {
    grid.rect(x, y, w-unit(width, "mm"), h-unit(heigth, "mm"),
              gp = gpar(fill = col['POS4'], col = NA))
  },
  POS5 = function(x, y, w, h) {
    grid.rect(x, y, w-unit(width, "mm"), h-unit(heigth, "mm"),
              gp = gpar(fill = col['POS5'], col = NA))
  },
  POS6 = function(x, y, w, h) {
    grid.rect(x, y, w-unit(width, "mm"), h-unit(heigth, "mm"),
              gp = gpar(fill = col['POS6'], col = NA))
  },
  POS7 = function(x, y, w, h) {
    grid.rect(x, y, w-unit(width, "mm"), h-unit(heigth, "mm"),
              gp = gpar(fill = col['POS7'], col = NA))
  },
  POS8 = function(x, y, w, h) {
    grid.rect(x, y, w-unit(width, "mm"), h-unit(heigth, "mm"),
              gp = gpar(fill = col['POS8'], col = NA))
  },
  POS9 = function(x, y, w, h) {
    grid.rect(x, y, w-unit(width, "mm"), h-unit(heigth, "mm"),
              gp = gpar(fill = col['POS9'], col = NA))
  },
  POS10 = function(x, y, w, h) {
    grid.rect(x, y, w-unit(width, "mm"), h-unit(heigth, "mm"),
              gp = gpar(fill = col['POS10'], col = NA))
  },
  POS11 = function(x, y, w, h) {
    grid.rect(x, y, w-unit(width, "mm"), h-unit(heigth, "mm"),
              gp = gpar(fill = col['POS11'], col = NA))
  },
  POS12 = function(x, y, w, h) {
    grid.rect(x, y, w-unit(width, "mm"), h-unit(heigth, "mm"),
              gp = gpar(fill = col['POS12'], col = NA))
  },
  MISSING = function(x, y, w, h) {
    grid.rect(x, y, w-unit(width, "mm"), h-unit(heigth, "mm"),
              gp = gpar(fill = col['MISSING'], col = NA))
  }
)


heatmap_legend_param = list(title = "Features", at = names(col), 
                            labels = c(rownames(mat), "Data not available"))
onc<-oncoPrint(mat,
               alter_fun = alter_fun, col = col,
               column_title = column_title, heatmap_legend_param = heatmap_legend_param, 
               top_annotation = HeatmapAnnotation(cbar = anno_oncoprint_barplot(),
                                                  'Tumor type' = tumortype), alter_fun_is_vectorized = FALSE)

rowOrder<- onc@row_names_param$labels[onc@row_order]
#rowOrder<- rownames(mat)
SampleOrder<- onc@column_names_param$labels[onc@column_order]
SampleOrder1<-SampleOrder[SampleOrder%in%colnames(mat)[which(tumortype==unique(tumortype)[1])]]
SampleOrder2<-SampleOrder[SampleOrder%in%colnames(mat)[which(tumortype==unique(tumortype)[2])]]
SampleOrder3<-SampleOrder[SampleOrder%in%colnames(mat)[which(tumortype==unique(tumortype)[3])]]

mat<-matOriginal
mat[mat=='NEG']<-""
for (i in 1:nrow(mat)){
  mat[i,][mat[i,]=="POS"] <- paste0("POS", i)
}
mat <- as.matrix(mat)

pdf("OncoPrint.pdf", 10,5)
set.seed(19)
print(oncoPrint(mat,
                alter_fun = alter_fun, col = col, row_order = rowOrder, column_order = c(SampleOrder1, SampleOrder2, SampleOrder3),
                column_title = column_title, heatmap_legend_param = heatmap_legend_param, 
                top_annotation = HeatmapAnnotation(cbar = anno_oncoprint_barplot(),
                                                   'Type' = tumortype), show_pct = F, left_annotation = rowAnnotation(Frequency = percMut), alter_fun_is_vectorized = FALSE))
dev.off()

