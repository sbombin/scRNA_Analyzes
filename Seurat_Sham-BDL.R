
sample1_path <- "/Users/sbombin/Desktop/analysis_mm39.comb/Ms802_BEC_sham/DGE_filtered"
sample1 <- ReadParseBio(sample1_path)
# Check to see if empty gene names are present
table(rownames(sample1) == "")
rownames(sample1)[rownames(sample1) == ""] <- "unknown"
cell_meta <- read.csv(paste0(sample1_path, "/cell_metadata.csv"), row.names = 1)
obj802 <- CreateSeuratObject(sample1, min.genes = 100, min.cells = 2, meta.data = cell_meta)
obj802@meta.data[, "group"] <- "Sham" ## Add new metadata variable called group
obj802@meta.data[, "sample"] <- "Ms802_Sham" 
obj802@meta.data$orig.ident <- factor(rep("obj802", nrow(obj802@meta.data)))
Idents(obj802) <- obj802@meta.data$orig.ident
SaveObject(obj802, "obj802.seurat")
obj802 <- ReadObject("obj802.seurat")

sample2_path <- "/Users/sbombin/Desktop/analysis_mm39.comb/Ms840_BEC_sham/DGE_filtered"
sample2 <- ReadParseBio(sample2_path)
# Check to see if empty gene names are present
table(rownames(sample2) == "")
rownames(sample2)[rownames(sample2) == ""] <- "unknown"
cell_meta <- read.csv(paste0(sample2_path, "/cell_metadata.csv"), row.names = 1)
obj840 <- CreateSeuratObject(sample2, min.genes = 100, min.cells = 2, meta.data = cell_meta)
obj840@meta.data[, "group"] <- "BDL"
obj840@meta.data[, "sample"] <- "Ms840_BDL" ## Add new metadata variable called sample
obj840@meta.data$orig.ident <- factor(rep("obj840", nrow(obj840@meta.data)))
Idents(obj840) <- obj840@meta.data$orig.ident
SaveObject(obj840, "obj840.seurat")
obj840 <- ReadObject("obj840.seurat")

sample3_path <- "/Users/sbombin/Desktop/analysis_mm39.comb/Ms805_BEC_BDL/DGE_filtered"
sample3 <- ReadParseBio(sample3_path)
# Check to see if empty gene names are present
table(rownames(sample3) == "")
rownames(sample3)[rownames(sample3) == ""] <- "unknown"
cell_meta <- read.csv(paste0(sample3_path, "/cell_metadata.csv"), row.names = 1)
obj805 <- CreateSeuratObject(sample3, min.genes = 100, min.cells = 2, meta.data = cell_meta)
obj805@meta.data[, "group"] <- "BDL"
obj805@meta.data[, "sample"] <- "Ms805_BDL" ## Add new metadata variable called sample
obj805@meta.data$orig.ident <- factor(rep("obj805", nrow(obj805@meta.data)))
Idents(obj805) <- obj805@meta.data$orig.ident
SaveObject(obj805, "obj805.seurat")
obj805 <- ReadObject("obj805.seurat")

sample4_path <- "/Users/sbombin/Desktop/analysis_mm39.comb/Ms877_BEC_BDL/DGE_filtered"
sample4 <- ReadParseBio(sample4_path)
# Check to see if empty gene names are present
table(rownames(sample4) == "")
rownames(sample4)[rownames(sample4) == ""] <- "unknown"
cell_meta <- read.csv(paste0(sample4_path, "/cell_metadata.csv"), row.names = 1)
obj877 <- CreateSeuratObject(sample4, min.genes = 100, min.cells = 2, meta.data = cell_meta)
obj877@meta.data[, "group"] <- "Sham"
obj877@meta.data[, "sample"] <- "Ms877_Sham"
obj877@meta.data$orig.ident <- factor(rep("obj877", nrow(obj877@meta.data)))
Idents(obj877) <- obj877@meta.data$orig.ident
SaveObject(obj877, "obj877.seurat")
obj877 <- ReadObject("obj877.seurat")

## Merge Seurat objects
pbmc <- merge(obj802, y = c(obj840, obj805, obj877), add.cell.ids = c("Ms802", "Ms840", "Ms805", "Ms877"), project = "SHAM-BDL")
#SaveObject(pbmc, "seuratobj.Sham-BDL")
##Option to downsample to the minimal number of cells
##down <- subset(x = immune.combined, downsample = 2086)

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

SaveFigure(p1, "UMAP_Clusters", width = 10, height = 8)
SaveFigure(p2, "UMAP_Sample", width = 10, height = 8)
SaveFigure(p3, "UMAP_Group", width = 10, height = 8)

## MArkers and DE
markers <- FindAllMarkers(pbmc, only.pos = FALSE, min.pct = 0.1, logfc.threshold = 0.25, test.use = "MAST")
pbmc.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

write.csv(pbmc.markers,"/Users/sbombin/Desktop/analysis_mm39.comb/Seurat_Analysis/BDL//All_Clusters_DE.csv", row.names = TRUE)

## DESeq2
pbmc[["RNA"]]@counts <- as.matrix(pbmc[["RNA"]]@counts)+1
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.5, test.use = "DESeq2")
pbmc.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

write.csv(pbmc.markers,"/Users/sbombin/Desktop/analysis_mm39.comb/Seurat_Analysis/BDL//DESeq2_All_Clusters.csv", row.names = TRUE)

#cluster0.markers <- FindMarkers(pbmc, ident.1 = 0, logfc.threshold = 0.8, test.use = "roc", only.pos = TRUE)
pbmc.markers2 <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.8, test.use = "roc")
write.csv(pbmc.markers2,"/Users/sbombin/Desktop/analysis_mm39.comb/Seurat_Analysis/BDL//Markers_Power_AllClusters.csv", row.names = TRUE)

pbmc.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
p4 <- DoHeatmap(pbmc, features = top10$gene) + NoLegend()
SaveFigure(p4, "Heatmap_top10_all", width = 15, height = 15)

##Save Seurat
SaveObject(pbmc, "seuratobj.Sham-BDL")
pbmc <- ReadObject("seuratobj.Sham-BDL")

##Convert Seurat to AnnData
library(SeuratDisk)
setwd('/Users/sbombin/Desktop/analysis_mm39.comb/Seurat_Analysis/BDL/')
SaveH5Seurat(pbmc, filename = "BDL_all.h5Seurat")
Convert("BDL_all.h5Seurat", dest = "h5ad")

##Subset Filtering
#p0 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
#SaveFigure(p0, "Scatter", width = 10, height = 8)

#pbmc2 <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 4500)
#pbmc2 <- NormalizeData(pbmc2, normalization.method = "LogNormalize", scale.factor = 10000)
#pbmc2 <- FindVariableFeatures(pbmc2, selection.method = "vst", nfeatures = 2000)
#all.genes2 <- rownames(pbmc2)
#pbmc2 <- ScaleData(pbmc2, features = all.genes2)
#pbmc2 <- RunPCA(pbmc2)
##Find optimal PCA dimensions
#n3pcs2 <- fgetNPCs(pbmc2,MIN_CUM_SD)
#pbmc2 <- JackStraw(pbmc2, num.replicate = 100, dims = npcs2+5)
#pbmc2 <- ScoreJackStraw(pbmc2, dims = 1:(npcs2+5))
#print(JackStrawPlot(pbmc2, dims = 1:(npcs2+5)))
#Clustering
#pbmc2 <- FindNeighbors(pbmc2, dims = 1:npcs2)
#pbmc2 <- FindClusters(pbmc2, resolution = 0.6)
#pbmc2 <- RunTSNE(pbmc2, dims = 1:npcs2)
#pbmc2 <- RunUMAP(pbmc2, dims = 1:npcs2)

#p2.1 <- DimPlot(pbmc2, reduction = "umap", label = TRUE)
#p2.2 <- DimPlot(pbmc2, label = TRUE, group.by = "sample", reduction = "umap", repel = TRUE)
#p2.3 <- DimPlot(pbmc2, label = FALSE, group.by = "group", reduction = "umap")

#SaveFigure((p1+p2+p2.1+p2.2), "Extra_Filtered_Clusters", width = 8, height = 8)
#p1+p2+p2.1+p2.2

## DE Between groups
#de <-  FindConservedMarkers(pbmc, ident.1 = 1, min.pct = 0.25, logfc.threshold = 0.25, test.use = "MAST", 
#                            grouping.var = "group", verbose = TRUE)
#de2 <-  FindConservedMarkers(pbmc, ident.1 =1, min.pct = 0.25, logfc.threshold = 0.8, 
#                             test.use = "MAST", grouping.var = "group", verbose = FALSE, only.pos = TRUE)

# Create function to get conserved markers for any given cluster
library(purrr)
library(tibble)
library(multtest)
library(metap)

get_conserved <- function(cluster){
  FindConservedMarkers(pbmc,
                       ident.1 = cluster, min.pct = 0.1, logfc.threshold = 0.8,
                       grouping.var = "group", test.use = "DESeq2", 
                       only.pos = FALSE) %>%
    rownames_to_column(var = "gene") %>%
    cbind(cluster_id = cluster, .)
}

# Iterate function across desired clusters
conserved_markers <- map_dfr(c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15), get_conserved)
write.csv(conserved_markers,"/Users/sbombin/Desktop/analysis_mm39.comb/Seurat_Analysis/BDL//DESeq2_conserved_markers.csv", row.names = FALSE)

