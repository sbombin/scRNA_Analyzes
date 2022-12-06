
sample1_path <- "/Users/sbombin/Desktop/analysis_mm39.comb/Ms1405_BEC_Sox9con/DGE_filtered"
sample1 <- ReadParseBio(sample1_path)
# Check to see if empty gene names are present
table(rownames(sample1) == "")
rownames(sample1)[rownames(sample1) == ""] <- "unknown"
cell_meta1405 <- read.csv(paste0(sample1_path, "/cell_metadata.csv"), row.names = 1)
obj1405 <- CreateSeuratObject(sample1, min.genes = 100, min.cells = 2, meta.data = cell_meta1405)
obj1405@meta.data[, "group"] <- "Sox9_CTRL" ## Add new metadata variable called group
obj1405@meta.data[, "sample"] <- "Ms1405_CTRL" 
obj1405@meta.data$orig.ident <- factor(rep("obj1405", nrow(obj1405@meta.data)))
Idents(obj1405) <- obj1405@meta.data$orig.ident
SaveObject(obj1405, "obj1405.seurat")
obj1405 <- ReadObject("obj1405.seurat")

sample2_path <- "/Users/sbombin/Desktop/analysis_mm39.comb/Ms1402_BEC_Sox9cKO/DGE_filtered"
sample2 <- ReadParseBio(sample2_path)
# Check to see if empty gene names are present
table(rownames(sample2) == "")
rownames(sample2)[rownames(sample2) == ""] <- "unknown"
cell_meta1402 <- read.csv(paste0(sample2_path, "/cell_metadata.csv"), row.names = 1)
obj1402 <- CreateSeuratObject(sample2, min.genes = 100, min.cells = 2, meta.data = cell_meta1402)
obj1402@meta.data[, "group"] <- "Sox9_cKO"
obj1402@meta.data[, "sample"] <- "Ms1402_cKO" ## Add new metadata variable called sample
obj1402@meta.data$orig.ident <- factor(rep("obj1402", nrow(obj1402@meta.data)))
Idents(obj1402) <- obj1402@meta.data$orig.ident
SaveObject(obj1402, "obj1402.seurat")
obj1402 <- ReadObject("obj1402.seurat")

pbmc <- merge(obj1405, y = obj1402, add.cell.ids = c("Ms1405", "Ms1402"), project = "BEC")

pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 5000)
## Scaling the data
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
#pbmc <- ScaleData(pbmc)
## LDR
pbmc <- RunPCA(pbmc)
##Find optimal PCA dimensions
npcs <- fgetNPCs(pbmc,MIN_CUM_SD)
pbmc <- JackStraw(pbmc, num.replicate = 100, dims = npcs+5)
pbmc <- ScoreJackStraw(pbmc, dims = 1:(npcs+5))
print(JackStrawPlot(pbmc, dims = 1:(npcs+5)))
#Clustering
pbmc <- FindNeighbors(pbmc, dims = 1:npcs)
pbmc <- FindClusters(pbmc, resolution = 0.6)
pbmc <- RunTSNE(pbmc, dims = 1:npcs)
pbmc <- RunUMAP(pbmc, dims = 1:npcs)

p1 <- DimPlot(pbmc, reduction = "umap", label = TRUE)
p2 <- DimPlot(pbmc, label = TRUE, group.by = "sample", reduction = "umap", repel = TRUE)
p3 <- DimPlot(pbmc, label = FALSE, group.by = "group", reduction = "umap")

SaveFigure(p1, "Umap_Clusters", width = 10, height = 8)
SaveFigure(p2, "Umap_Sample", width = 10, height = 8)
SaveFigure(p3, "Umap_Group", width = 10, height = 8)

##MArkers and DE
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.8, test.use = "MAST")
pbmc.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

write.csv(pbmc.markers,"/Users/sbombin/Desktop/analysis_mm39.comb/Seurat_Analysis/BEC//All_Clusters_DE.csv", row.names = TRUE)

#cluster0.markers <- FindMarkers(pbmc, ident.1 = 0, logfc.threshold = 0.8, test.use = "roc", only.pos = TRUE)
pbmc.markers2 <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.8, test.use = "roc")
write.csv(pbmc.markers2,"/Users/sbombin/Desktop/analysis_mm39.comb/Seurat_Analysis/BEC//Markers_Power_AllClusters.csv", row.names = TRUE)

pbmc.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
p4 <- DoHeatmap(pbmc, features = top10$gene) + NoLegend()
SaveFigure(p4, "Heatmap_top10_all", width = 15, height = 15)

##Save Seurat
SaveObject(pbmc, "seuratobj_BEC.scale_all")
pbmc <- ReadObject("seuratobj_BEC")

##Convert Seurat to AnnData
library(SeuratDisk)
setwd('/Users/sbombin/Desktop/analysis_mm39.comb/Seurat_Analysis/BEC/')
SaveH5Seurat(pbmc, filename = "BEC-cKO_all.h5Seurat")
Convert("BEC-cKO_all.h5Seurat", dest = "h5ad")

p5<- FeaturePlot(object = pbmc, features = "Krt19")
p6 <- FeaturePlot(object = pbmc, features = "Anxa2")
SaveFigure((p5 + p6), "Krt19_Anxa2_expression", width = 8, height = 8)

##Subset Filtering
p0 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
SaveFigure(p0, "Scatter", width = 10, height = 8)

pbmc2 <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 3000)
pbmc2 <- NormalizeData(pbmc2, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc2 <- FindVariableFeatures(pbmc2, selection.method = "vst", nfeatures = 2000)
all.genes2 <- rownames(pbmc2)
pbmc2 <- ScaleData(pbmc2, features = all.genes2)
pbmc2 <- RunPCA(pbmc2)
##Find optimal PCA dimensions
npcs2 <- fgetNPCs(pbmc2,MIN_CUM_SD)
#pbmc2 <- JackStraw(pbmc2, num.replicate = 100, dims = npcs2+5)
#pbmc2 <- ScoreJackStraw(pbmc2, dims = 1:(npcs2+5))
#print(JackStrawPlot(pbmc2, dims = 1:(npcs2+5)))
#Clustering
pbmc2 <- FindNeighbors(pbmc2, dims = 1:npcs2)
pbmc2 <- FindClusters(pbmc2, resolution = 0.6)
pbmc2 <- RunTSNE(pbmc2, dims = 1:npcs2)
pbmc2 <- RunUMAP(pbmc2, dims = 1:npcs2)

p2.1 <- DimPlot(pbmc2, reduction = "umap", label = TRUE)
p2.2 <- DimPlot(pbmc2, label = TRUE, group.by = "sample", reduction = "umap", repel = TRUE)
p2.3 <- DimPlot(pbmc2, label = FALSE, group.by = "group", reduction = "umap")

SaveFigure((p1+p2+p2.1+p2.2), "Extra_Filtered_Clusters", width = 8, height = 8)
p1+p2+p2.1+p2.2