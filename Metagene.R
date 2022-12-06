library(Seurat)

setwd("/Users/sbombin/Desktop/analysis_mm39.comb/Seurat_Analysis/BDL/")
data_path <- "/Users/sbombin/Desktop/analysis_mm39.comb/Seurat_Analysis/BDL/"
fig_path <- "/Users/sbombin/Desktop/analysis_mm39.comb/Seurat_Analysis/BDL/"

#object<- ReadObject("seuratobj_BEC_cKO.scale_all")
object<- ReadObject("seuratobj.Sham-BDL.scale_all")

#a <- as.data.frame(my.data)
#mg0.mean.expr <- a[mg0,] %>% as.data.frame()
b <- my.data[rownames(my.data) %in% mg4, ] %>% as.data.frame()
#c <- as.data.frame(rownames(my.data))

## metagene - Gene Expression
###
mg0 <- c("Apoe", "Mgst1", "Txnip", "Ttr", "Cuta")
mg1 <-c('Cd14', 'Krt19', 'Efhd2', 'Lgals3', 'Neurl3')
mg2 <- c('Fam107b', 'Nop58', 'Eif1a','Cdkn1a','Selenop') 
mg3 <-c('Jund', 'Spag9', 'Fosb', 'Klf6', 'Pdgfa')
mg4 <- c('Slc2a2', 'Sall1', 'Cftr', 'Arid5a', 'Usp16')


my.data=GetAssayData(object, slot = "data") ## normalized counts data
Metagene0 <- colMeans(my.data[mg0,],na.rm=T) 
object@meta.data$Metagene0 <- Metagene0
FeaturePlot(object,features="Metagene0")

Metagene1 <- colMeans(my.data[mg1,],na.rm=T)
object@meta.data$Metagene1 <- Metagene1
FeaturePlot(object,features="Metagene1")

Metagene2 <- colMeans(my.data[mg2,],na.rm=T)
object@meta.data$Metagene2 <- Metagene2
FeaturePlot(object,features="Metagene2")

Metagene3 <- colMeans(my.data[mg3,],na.rm=T)
object@meta.data$Metagene3 <- Metagene3
FeaturePlot(object,features="Metagene3")

Metagene4 <- colMeans(my.data[mg4,],na.rm=T)
object@meta.data$Metagene4 <- Metagene4
FeaturePlot(object,features="Metagene4")

p1 <- DotPlot(object, features = c("Metagene0","Metagene1", "Metagene2", "Metagene3","Metagene4", 'Sox9', 'Ccn1', 'Alb'), 
              dot.scale = 8, cols = c("lightgrey", "red")) & coord_flip() 
p1

SaveFigure(p1, "Dotplot_Metagenes", width = 15, height = 8)

p1 <- FeaturePlot(object, features = c("Metagene0","Metagene1", "Metagene2", "Metagene3","Metagene4"), split.by = "group", cols = c("grey", "red"))
SaveFigure(p1, "FeaturePlot_Metagenes_1", width = 15, height = 18)

p2 <- FeaturePlot(object, features = c("Metagene0","Metagene1", "Metagene2", "Metagene3","Metagene4","Sox9"), split.by = "group", cols = c("grey", "red"), combine = FALSE)
p2 <- wrap_plots(plots = p2, ncol=2) 
SaveFigure(p2, "FeaturePlot_Metagenes_2", width = 15, height = 18)

p1 <- VlnPlot(object, features = c("Metagene0","Metagene1", "Metagene2", "Metagene3","Metagene4","Sox9"))
SaveFigure(p1, "VlnPlot_Metagenes_1", width = 15, height = 10)
p2 <- VlnPlot(object, features = c("Metagene0","Metagene1", "Metagene2", "Metagene3","Metagene4","Sox9"),split.by = "group", combine = FALSE)
p2 <- wrap_plots(plots = p2, ncol=2) 
SaveFigure(p2, "VlnPlot_Metagenes_2", width = 20, height = 15)

p <- DimPlot(object, label = FALSE, split.by = "group", reduction = "umap", repel = FALSE)
SaveFigure(p, "UMAP_GroupSplit_2", width = 12, height = 8)

VlnPlot(object, features = "Metagene4", split.by = "group")

## Only Control Cells
#control <- subset(x = object, subset = group == "Sox9_CTRL") # cKO
control <- subset(x = object, subset = group == "Sham") # BDL
mast2 <- FindAllMarkers(control, only.pos = FALSE, min.pct = 0.1, logfc.threshold = 0.25, test.use = "MAST")
write.csv(mast2,"DE-MAST_AllMarkers_Control.csv", row.names = TRUE)

