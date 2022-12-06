library(Seurat)
library(dplyr)
library(Matrix)
library(ggplot2)
library(patchwork)

##rm(list = ls())

#data_path <- "/Users/sbombin/Desktop/analysis_mm39.comb/Seurat_Analysis/BEC/"
#fig_path <- "/Users/sbombin/Desktop/analysis_mm39.comb/Seurat_Analysis/BEC/"

# Convenience functions
SaveFigure <- function(plots, name, type = "png", width, height, res){
  if(type == "png") {
    png(paste0(fig_path, name, ".", type),
        width = width, height = height, units = "in", res = 300)
  } else {
    pdf(paste0(fig_path, name, ".", type),
        width = width, height = height)
  }
  print(plots)
  dev.off()
}

SaveObject <- function(object, name){
  saveRDS(object, paste0(data_path, name, ".RDS"))
}

ReadObject <- function(name){
  readRDS(paste0(data_path, name, ".RDS"))
}

#### Find automatically Optimal Number PCAs to explain variation
MIN_CUM_SD <- 90 ## Percent of variance to be captured by PCs.
## Get number of PCs required to capture xx% of variability                 
fgetNPCs <- function(obj,MIN_CUM_SD=90){
  sds <- Stdev(obj,reduction="pca")
  cumsds <- 0
  for(i in 1:length(sds)){
    cumsds <- cumsds + sds[i]
    if(cumsds >= MIN_CUM_SD){
      return(i)
    }
  }
}
