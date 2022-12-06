setwd("/Users/sbombin/Desktop/EICC/Projects/Veronika_AMA23057/")
data_path <- "/Users/sbombin/Desktop/EICC/Projects/Veronika_AMA23057/"
fig_path <- "/Users/sbombin/Desktop/EICC/Projects/Veronika_AMA23057/"

s1 <- Read10X(data.dir = "matrix_raw_norm")
s2 <- Read10X(data.dir = "matrix_filtr_norm")
s3 <- Read10X(data.dir = "matrix_raw_hypox")
s4 <- Read10X(data.dir = "matrix_filtr_hypox")

#norm1 <- CreateSeuratObject(counts = s1, project = "s1",min.cells = 3, min.features = 200)
norm2 <- CreateSeuratObject(counts = s2, project = "norm",min.cells = 3, min.features = 200)
#hypox1 <- CreateSeuratObject(counts = s3, project = "s3",min.cells = 3, min.features = 200)
hypox2 <- CreateSeuratObject(counts = s4, project = "hypox",min.cells = 3, min.features = 200)

#norm1[["percent.mt"]] <- PercentageFeatureSet(norm1, pattern = "^MT-")
#p1 <- VlnPlot(norm1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

norm2[["percent.mt"]] <- PercentageFeatureSet(norm2, pattern = "^MT-")
p2 <- VlnPlot(norm2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
p2
SaveFigure(p2, "QC_Norm_Plot1", width = 10, height = 8)

#p1 <- FeatureScatter(norm1, feature1 = "nCount_RNA", feature2 = "percent.mt")
#p2 <- FeatureScatter(norm1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
#p1 + p2

p3 <- FeatureScatter(norm2, feature1 = "nCount_RNA", feature2 = "percent.mt")
SaveFigure(p3, "QC_Norm_Plot2", width = 10, height = 8)
p4 <- FeatureScatter(norm2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
SaveFigure(p4, "QC_Norm_Plot3", width = 10, height = 8)
p3 + p4

p1 + p2 + p3 + p4

summary(norm1$nCount_RNA)
summary(norm2$nCount_RNA)

summary(norm1$nFeature_RNA)
summary(norm2$nFeature_RNA)

#norm1 <- subset(norm1, subset = nFeature_RNA > 200 & nFeature_RNA < 10000 & percent.mt < 5 & nCount_RNA >700 & nCount_RNA < 150000)
norm2 <- subset(norm2, subset = nFeature_RNA > 200 & nFeature_RNA < 10000 & percent.mt < 5 & nCount_RNA >700 & nCount_RNA < 150000)

hypox2[["percent.mt"]] <- PercentageFeatureSet(hypox2, pattern = "^MT-")
p2 <- VlnPlot(hypox2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
p2
SaveFigure(p2, "QC_Hypox_Plot1", width = 10, height = 8)

p3 <- FeatureScatter(hypox2, feature1 = "nCount_RNA", feature2 = "percent.mt")
SaveFigure(p3, "QC_Hypox_Plot2", width = 10, height = 8)
p4 <- FeatureScatter(hypox2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
SaveFigure(p4, "QC_Hypox_Plot3", width = 10, height = 8)
p3 + p4

summary(hypox2$nCount_RNA)
summary(hypox2$nFeature_RNA)

hypox2 <- subset(hypox2, subset = nFeature_RNA > 200 & nFeature_RNA < 10000 & percent.mt < 5 & nCount_RNA >700 & nCount_RNA < 150000)

norm2[["group"]] <- "norm"
hypox2[["group"]] <- "hypox"

object.list <- c(norm2,hypox2)
names(object.list) <- c("norm","hypox")

# normalize and identify variable features for each dataset independently
object.list <- lapply(X = object.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = object.list)

### Perform Integration
object.anchors <- FindIntegrationAnchors(object.list = object.list, anchor.features = features)
object <- IntegrateData(anchorset = object.anchors)

DefaultAssay(object) <- "integrated"

object <- ScaleData(object, verbose = FALSE)
object <- RunPCA(object, verbose = FALSE)

object <- JackStraw(object, num.replicate = 100, dims = 50)
object <- ScoreJackStraw(object, dims = 1:50)
p <- JackStrawPlot(object, dims = 1:50)
p
SaveFigure(p, "JackStrawPlot_PC50", width = 10, height = 8)

p <- ElbowPlot(object, ndims = 50)
p
SaveFigure(p, "ElbowPlot_PC50", width = 10, height = 8)

### Optimal PCAs ~50
object <- RunUMAP(object, dims = 1:50)

object <- FindNeighbors(object, dims = 1:50)
## Clustree to find resolution
library(clustree)
clustree <- FindClusters(object, resolution = c(0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.2,1.4))
p <- clustree(clustree)
p
SaveFigure(p, "Clustree", width = 8, height = 10)
## Cluster 
object <- FindClusters(object, resolution = 0.4)

object <- RunTSNE(object, dims = 1:50)

p <- DimPlot(object, reduction = "umap", label = TRUE)
p
SaveFigure(p, "UMAP_Clusters_04", width = 10, height = 8)

p <- DimPlot(object, reduction = "tsne", label = TRUE)
p
SaveFigure(p, "TSNE_Clusters_04", width = 10, height = 8)

p <- DimPlot(object, reduction = "tsne", label = FALSE,group.by = "group")
p
SaveFigure(p, "tSNE_Group", width = 10, height = 8)

### 
DefaultAssay(object) <- "RNA"
all.genes <- as.data.frame(row.names(object))

mk1 <- c("HTATIP2","MYO10","IL13RA2")
mk2 <- c("GPI","STRADA","SLC2A1","VEGFA","CA9","CTCFL")
mk0 <- c("HTATIP2","MYO10","IL13RA2","GPI","STRADA","SLC2A1","VEGFA","CA9","CTCFL")

p1 <- FeaturePlot(object,features = mk1, reduction = "tsne")
p1
SaveFigure(p1, "tSNE_Markers1", width = 10, height = 8)

p2 <- FeaturePlot(object,features = mk2, reduction = "tsne")
p2
SaveFigure(p2, "tSNE_Markers2", width = 10, height = 10)

p1 <- DotPlot(object, features = mk0, dot.scale = 8) + RotatedAxis()
p1
SaveFigure(p1, "DotPlot_Markers", width = 10, height = 8)

p2 <- VlnPlot(object, features = mk0) 
p2
SaveFigure(p2, "VlnPlot_Markers", width = 15, height = 10)

#### DE
summary(object$percent.mt)
FeaturePlot(object,features = "percent.mt", reduction = "tsne")

library(purrr)
library(tibble)

mast.0 <- FindMarkers(object = object,ident.1 = "hypox", group.by = "group", min.pct = 0.1, logfc.threshold = 0.25,
                           test.use = "MAST", only.pos = FALSE, subset.ident = "0")
wilcox.0 <- FindMarkers(object = object,ident.1 = "hypox", group.by = "group", min.pct = 0.1, logfc.threshold = 0.25,
                      test.use = "wilcox", only.pos = FALSE, subset.ident = "0")
wilcox.1 <- FindMarkers(object = object,ident.1 = "hypox", group.by = "group", min.pct = 0.1, logfc.threshold = 0.25,
                        test.use = "wilcox", only.pos = FALSE, subset.ident = "1")


get_de <- function(cluster){
  FindMarkers(object = object, ident.1 = "hypox", min.pct = 0.1, logfc.threshold = 0.25,
              group.by = "group", test.use = "MAST", only.pos = FALSE, subset.ident = cluster) %>%
    rownames_to_column(var = "gene") %>% cbind(cluster_id = cluster, .)
}

# Iterate function across desired clusters
mast.de <- map_dfr(c(0,1,2,3,4,5,6,7,8), get_de)
mast.all <- FindAllMarkers(object = object, min.pct = 0.1, logfc.threshold = 0.25, test.use = "MAST", only.pos = FALSE)

write.csv(mast.de,"DE_MAST_Hypox_vs_Norm.csv", row.names = FALSE)


### Save Seurat Object as RDS
SaveObject(object, "SeuratObject_Processed")

object <- ReadObject("SeuratObject_Processed")