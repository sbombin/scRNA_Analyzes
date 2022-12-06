ibrary(Seurat)

setwd("/Users/sbombin/Desktop/analysis_mm39.comb/Seurat_Analysis/BDL/")
data_path <- "/Users/sbombin/Desktop/analysis_mm39.comb/Seurat_Analysis/BDL/"
fig_path <- "/Users/sbombin/Desktop/analysis_mm39.comb/Seurat_Analysis/BDL/"

object <- ReadObject("seuratobj.Sham-BDL.scale_all")
### Rename Clusters
object <- RenameIdents(object, `0` = "Sham1", `1` = "BDL1", `2` = "BDL2",`3` = "Sham2", `4` = "Sham3", `5` = "BDL3", `6` = "BDL4", `7` = "Sham4", `8` = "Sham5",
                       `9` = "BDL5", `10` = "BDL6", `11` = "Sham6", `12` = "BDL7", `13` = "Sham7", `14` = "Sham8", `15` = "BDL8")

p1 <- DimPlot(object, reduction = "umap", label = TRUE)
p1
object.markers <- FindAllMarkers(object, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5, test.use = "MAST")
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

m12 <- dplyr::filter(top20, cluster == '12')
m12 <- m12$gene
p <- DotPlot(object, features = m12, dot.scale = 8, cols = c("lightgrey", "red")) & coord_flip() 
p
SaveFigure(p, "Dotplot_Cluster12", width = 10, height = 8)

m13 <- dplyr::filter(top20, cluster == '13')
m13 <- m13$gene
p <- DotPlot(object, features = m13, dot.scale = 8, cols = c("lightgrey", "red")) & coord_flip() 
p
SaveFigure(p, "Dotplot_Cluster13", width = 10, height = 8)

m14 <- dplyr::filter(top20, cluster == '14')
m14 <- m14$gene
p <- DotPlot(object, features = m14, dot.scale = 8, cols = c("lightgrey", "red")) & coord_flip() 
p
SaveFigure(p, "Dotplot_Cluster14", width = 10, height = 8)

m15 <- dplyr::filter(top20, cluster == '15')
m15 <- m15$gene
p <- DotPlot(object, features = m15, dot.scale = 8, cols = c("lightgrey", "red")) & coord_flip() 
p
SaveFigure(p, "Dotplot_Cluster15", width = 10, height = 8)

write.csv(top20,paste0('Top20_Markers.csv'),row.names = FALSE)


p <- DimPlot(object, reduction = "umap", label = TRUE)
p
SaveFigure(p, "UMAP_Clusters_renamed1", width = 10, height = 8)

p <- DimPlot(object, reduction = "umap", label = FALSE)
p
SaveFigure(p, "UMAP_Clusters_renamed2", width = 10, height = 8)

p <- VlnPlot(object, features = c("Mki67", "Pcna", "Ccna2", "Ccl2", "Cxcl1", "Cxcl2", "Sparc", "Tgfb2",
                                               "Cd14", "Fn1", "Cdk1", "Itgb6", "Ccn1", "Ccn2"), split.by = NULL, group.by = "sample", pt.size = 0)

#wrap_plots(plots = p, ncol = 1)
#SaveFigure(wrap_plots(plots = p, ncol = 1), "VlnPlots4_Intgrtn_Sham-BDL", width = 10, height = 15)
SaveFigure(p, "VlnPlot_renamed2", width = 15, height = 15)

write.csv(markers,paste0('DE-MAST_All_renamed.csv'), row.names = TRUE)

markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC) -> top10
p <- DoHeatmap(object, features = top10$gene, size = 3.5) + NoLegend()
SaveFigure(p, "Heatmap_top10_all_renamed", width = 20, height = 15)

#top <- dplyr::filter(top10, cluster == 'Sham5' | cluster == 'BDL7' | cluster == 'Sham7')
#p1 <- DoHeatmap(object, features = top$gene, size = 3.5) + NoLegend()
#p1
#SaveFigure(p1, "Heatmap_Clusters8-12-13_renamed", width = 10, height = 10)



p <- VlnPlot(clean3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)

p1 <- FeaturePlot(object, features = "nCount_RNA")
p2 <- FeaturePlot(object, features = "nFeature_RNA")

SaveFigure(p1, "UMAP_Transcript_Count", width = 6, height = 6)
SaveFigure(p2, "UMAP_Genes_Count", width = 6, height = 6)



p3 <- VlnPlot(object, features = "nCount_RNA", split.by = "group")
p4 <- VlnPlot(object, features = "nFeature_RNA", split.by = "group")

SaveFigure(p3, "VlnPlot_Transcript_Count", width = 6, height = 6)
SaveFigure(p4, "VlnPlot_Genes_Count", width = 6, height = 6)



p <- FeatureScatter(object, feature1 = "nCount_RNA", feature2 = "read_count")

SaveFigure(p, "Scatter_Reads_to_UMI", width = 6, height = 6)

p1 <- object@meta.data %>%
  ggplot(aes(x=(nFeature_RNA/nCount_RNA), color = group, fill=group)) + geom_density(alpha = 0.2) +
  theme_classic() + geom_vline(xintercept = 0.8)
p1

p2 <- object@meta.data %>%
  ggplot(aes(x=log10(nFeature_RNA/nCount_RNA), color = group, fill=group)) + geom_density(alpha = 0.2) +
  theme_classic() + geom_vline(xintercept = 0.8)
p2

p3 <- object@meta.data %>%
  ggplot(aes(x=log10(nFeature_RNA)/log10(nCount_RNA), color = group, fill=group)) + geom_density(alpha = 0.2) +
  theme_classic() + geom_vline(xintercept = 0.8)
p3

p4 <- object@meta.data %>%
  ggplot(aes(x=log(nFeature_RNA/nCount_RNA), color = group, fill=group)) + geom_density(alpha = 0.2) +
  theme_classic() + geom_vline(xintercept = 0.8)
p4

p1+p2+p3+p4

#object$log10GenesPerUMI <- log10(object$nFeature_RNA) / log10(object$nCount_RNA)
#object$nlogGenesPerUMI <- -log(object$nFeature_RNA/object$nCount_RNA)
object$GenesPerUMI <- object$nFeature_RNA/object$nCount_RNA
object$GenesPerRead <- object$nFeature_RNA/object$read_count
object$UMIperRead <- object$nCount_RNA/object$read_count


p1 <- VlnPlot(object, features = "GenesPerUMI", group.by = "group")
p2 <- VlnPlot(object, features = "GenesPerRead", group.by = "group")
p3 <- VlnPlot(object, features = "UMIperRead", group.by = "group")
#p1+p2+p3+p4
SaveFigure(p1, "VlnPlot_GenesPerUMI_Group", width = 8, height = 8)
SaveFigure(p2, "VlnPlot_GenesPerRead_Group", width = 8, height = 8)
SaveFigure(p3, "VlnPlot_UMIperRead_Group", width = 8, height = 8)

p1 <- FeaturePlot(object, features = "GenesPerUMI", split.by = "group")
p2 <- FeaturePlot(object, features = "GenesPerRead", split.by = "group")
p3 <- FeaturePlot(object, features = "UMIperRead", split.by = "group")
SaveFigure(p1, "FeaturePlot_GenesPerUMI_Group", width = 10, height = 6)
SaveFigure(p2, "FeaturePlot_GenesPerRead_Group", width = 10, height = 6)
SaveFigure(p3, "FeaturePlot_UMIperRead_Group", width = 10, height = 6)

p1 <- FeaturePlot(object, features = "GenesPerUMI",label = TRUE,repel = TRUE,label.size = 3)
p2 <- FeaturePlot(object, features = "GenesPerRead",label = TRUE,repel = TRUE,label.size = 3)
p3 <- FeaturePlot(object, features = "UMIperRead",label = TRUE,repel = TRUE,label.size = 3)
SaveFigure(p1, "FeaturePlot_GenesPerUMI_label", width = 8, height = 6)
SaveFigure(p2, "FeaturePlot_GenesPerRead_label", width = 8, height = 6)
SaveFigure(p3, "FeaturePlot_UMIperRead_label", width = 8, height = 6)

p1 <- DotPlot(object, features = "GenesPerUMI")
p2 <- DotPlot(object, features = "GenesPerRead", split.by = "group")
p3 <- DotPlot(object, features = "UMIperRead", split.by = "group")
#p1+p2+p3+p4
SaveFigure(p1, "DotPlot_GenesPerUMI_Group", width = 8, height = 6)
SaveFigure(p2, "DotPlot_GenesPerRead_Group", width = 8, height = 6)
SaveFigure(p3, "DotPlot_UMIperRead_Group", width = 8, height = 6)

a <-as.data.frame(p1$data)
color.by <- ifelse(test = split.colors, yes = 'colors', no = 'avg.exp.scaled')
plot <- ggplot(data = a, mapping = aes_string(x = 'features.plot', y = 'id')) +
  geom_point(mapping = aes_string(size = 'pct.exp', color = color.by)) +
  scale.func(range = c(0, dot.scale), limits = c(scale.min, scale.max)) +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
  guides(size = guide_legend(title = 'Percent Expressed')) +
  labs(
    x = 'Features',
    y = ifelse(test = is.null(x = split.by), yes = 'Identity', no = 'Split Identity')
  )