## Integration

# split the dataset into a list of two seurat objects 
seurat.list <- SplitObject(pbmc, split.by = "group")

# normalize and identify variable features for each dataset independently
seurat.list <- lapply(X = seurat.list, FUN = function(x) {
##  x <- NormalizeData(x, normalization.method = "CLR", margin = 2)
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 5000)
})

# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = seurat.list)

seurat.anchors <- FindIntegrationAnchors(object.list = seurat.list, anchor.features = features)
# this command creates an 'integrated' data assay
seurat.combined <- IntegrateData(anchorset = seurat.anchors)

# specify that we will perform downstream analysis on the corrected data note that the
# original unmodified data still resides in the 'RNA' assay
DefaultAssay(seurat.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
seurat.combined <- ScaleData(seurat.combined, verbose = FALSE)
##seurat.combined <- RunPCA(seurat.combined, npcs = 30, verbose = FALSE)
seurat.combined <- RunPCA(seurat.combined)
#Find PCA dimensions
npcs <- fgetNPCs(seurat.combined,MIN_CUM_SD)

seurat.combined <- RunUMAP(seurat.combined, reduction = "pca", dims = 1:npcs)
seurat.combined <- FindNeighbors(seurat.combined, reduction = "pca", dims = 1:npcs)
seurat.combined <- FindClusters(seurat.combined, resolution = 1.2)
SaveObject(seurat.combined, "seuratobj.Sham-BDL.intgrtn")
seurat.combined <- ReadObject("seuratobj_intgrtn_CLR2_EIC")

p1 <- DimPlot(seurat.combined, reduction = "umap", group.by = "group")
p2 <- DimPlot(seurat.combined, reduction = "umap", label = TRUE, repel = TRUE)
SaveFigure((p1 + p2), "Dimplot.Intgrtn_EIC_CLR2_group", width = 12, height = 8)

p3 <- DimPlot(seurat.combined, reduction = "umap", group.by = "sample")
p4 <- DimPlot(seurat.combined, reduction = "umap", label = TRUE, repel = TRUE)
SaveFigure((p3 + p4), "Dimplot.Intgrtn_EIC_CLR2_sample", width = 12, height = 8)

p5 <- FeaturePlot(object = seurat.combined, features = 'read_count')
p6 <- FeaturePlot(object = seurat.combined, features = 'tscp_count')
SaveFigure((p5 + p6), "FeaturePlot.Intgrtn_EIC_CLR2", width = 12, height = 8)

SaveFigure((p3+p4+p5+p6), "Plot.Intgrtn_EIC_CLR2_sample", width = 12, height = 12)