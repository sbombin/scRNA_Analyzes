pbmc2 <- merge(obj1325, y = obj1373, add.cell.ids = c("Ms1325", "Ms1373"), project = "EIC")

pbmc2 <- NormalizeData(pbmc2, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc2 <- FindVariableFeatures(pbmc2, selection.method = "vst", nfeatures = 5000)
## Scaling the data
pbmc2 <- ScaleData(pbmc2)
## LDR
pbmc2 <- RunPCA(pbmc2)
##Find optimal PCA dimensions
npcs2 <- fgetNPCs(pbmc2,MIN_CUM_SD)
#Clustering
pbmc2 <- FindNeighbors(pbmc2, dims = 1:npcs2)
pbmc2 <- FindClusters(pbmc2, resolution = 0.7)
