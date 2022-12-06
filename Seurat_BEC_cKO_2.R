library(MAST)

setwd("/Users/sbombin/Desktop/analysis_mm39.comb/Seurat_Analysis/BEC/")
data_path <- "/Users/sbombin/Desktop/analysis_mm39.comb/Seurat_Analysis/BEC/"
fig_path <- "/Users/sbombin/Desktop/analysis_mm39.comb/Seurat_Analysis/BEC/"

object <- ReadObject("seuratobj_BEC_cKO.scale_all")

p1 <- DimPlot(object, reduction = "umap", label = TRUE)
p1
object.markers <- FindAllMarkers(object, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.8, test.use = "MAST")
object.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC) -> top20

m0 <- dplyr::filter(top20, cluster == '0')
m0 <- m0$gene
p0 <- DotPlot(object, features = m0, dot.scale = 8, cols = c("lightgrey", "red")) & coord_flip() 
p0
SaveFigure(p0, "Dotplot_Cluster0", width = 10, height = 8)

m1 <- dplyr::filter(top20, cluster == '1')
m1 <- m1$gene
p1 <- DotPlot(object, features = m1, dot.scale = 8, cols = c("lightgrey", "red")) & coord_flip() 
p1
SaveFigure(p1, "Dotplot_Cluster1", width = 10, height = 8)

m2 <- dplyr::filter(top20, cluster == '2')
m2 <- m2$gene
p2 <- DotPlot(object, features = m2, dot.scale = 8, cols = c("lightgrey", "red")) & coord_flip() 
p2
SaveFigure(p2, "Dotplot_Cluster2", width = 10, height = 8)

m3 <- dplyr::filter(top20, cluster == '3')
m3 <- m3$gene
p3 <- DotPlot(object, features = m3, dot.scale = 8, cols = c("lightgrey", "red")) & coord_flip() 
p3
SaveFigure(p3, "Dotplot_Cluster3", width = 10, height = 8)

m4 <- dplyr::filter(top20, cluster == '4')
m4 <- m4$gene
p4 <- DotPlot(object, features = m4, dot.scale = 8, cols = c("lightgrey", "red")) & coord_flip() 
p4
SaveFigure(p4, "Dotplot_Cluster4", width = 10, height = 8)

m5 <- dplyr::filter(top20, cluster == '5')
m5 <- m5$gene
p5 <- DotPlot(object, features = m5, dot.scale = 8, cols = c("lightgrey", "red")) & coord_flip() 
p5
SaveFigure(p5, "Dotplot_Cluster5", width = 10, height = 8)

m6 <- dplyr::filter(top20, cluster == '6')
m6 <- m6$gene
p6 <- DotPlot(object, features = m6, dot.scale = 8, cols = c("lightgrey", "red")) & coord_flip() 
p6
SaveFigure(p6, "Dotplot_Cluster6", width = 10, height = 8)

m7 <- dplyr::filter(top20, cluster == '7')
m7 <- m7$gene
p7 <- DotPlot(object, features = m7, dot.scale = 8, cols = c("lightgrey", "red")) & coord_flip() 
p7
SaveFigure(p7, "Dotplot_Cluster7", width = 10, height = 8)

m8 <- dplyr::filter(top20, cluster == '8')
m8 <- m8$gene
p8 <- DotPlot(object, features = m8, dot.scale = 8, cols = c("lightgrey", "red")) & coord_flip() 
p8
SaveFigure(p8, "Dotplot_Cluster8", width = 10, height = 8)

m9 <- dplyr::filter(top20, cluster == '9')
m9 <- m9$gene
p9 <- DotPlot(object, features = m9, dot.scale = 8, cols = c("lightgrey", "red")) & coord_flip() 
p9
SaveFigure(p9, "Dotplot_Cluster9", width = 10, height = 8)

m10 <- dplyr::filter(top20, cluster == '10')
m10 <- m10$gene
p10 <- DotPlot(object, features = m10, dot.scale = 8, cols = c("lightgrey", "red")) & coord_flip() 
p10
SaveFigure(p10, "Dotplot_Cluster10", width = 10, height = 8)

m11 <- dplyr::filter(top20, cluster == '11')
m11 <- m11$gene
p11 <- DotPlot(object, features = m11, dot.scale = 8, cols = c("lightgrey", "red")) & coord_flip() 
p11
SaveFigure(p11, "Dotplot_Cluster11", width = 10, height = 8)

write.csv(top20,paste0('Top20_Markers.csv'),row.names = FALSE)





















