
library(Seurat)
library(ggplot2)
library(sctransform)

obj1325 <- ReadObject("obj1325.seurat")
obj1373 <- ReadObject("obj1373.seurat")
obj1339 <- ReadObject("obj1339.seurat_obj")
obj1361 <- ReadObject("obj1361.seurat_obj")
pbmc <- merge(obj1325, y = c(obj1373, obj1339, obj1361), add.cell.ids = c("Ms1325", "Ms1373", "Ms1339", "Ms1361"), project = "IEC")

# store mitochondrial percentage in object meta data
pbmc <- PercentageFeatureSet(pbmc, pattern = "^MT-", col.name = "percent.mt")

# run sctransform
pbmc <- SCTransform(pbmc, vars.to.regress = "percent.mt", verbose = FALSE)

# These are now standard steps in the Seurat workflow for visualization and clustering
pbmc <- RunPCA(pbmc, verbose = FALSE)
pbmc <- RunUMAP(pbmc, dims = 1:30, verbose = FALSE)

pbmc <- FindNeighbors(pbmc, dims = 1:30, verbose = FALSE)
pbmc <- FindClusters(pbmc, verbose = FALSE)
p1 <- DimPlot(pbmc, label = TRUE) 
p2 <- DimPlot(pbmc, label = TRUE, group.by = "sample", reduction = "umap", repel = TRUE)
p3 <- DimPlot(pbmc, label = FALSE, group.by = "group", reduction = "umap")
SaveFigure(p1, "Umap_Clusters_SCT", width = 10, height = 8)
SaveFigure(p2, "Umap_Sample_SCT", width = 10, height = 8)
SaveFigure(p3, "Umap_Group_SCT", width = 10, height = 8)