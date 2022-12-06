
sample1_path <- "/Users/sbombin/Desktop/analysis_mm39.comb/Ms1325_Tet1con/DGE_filtered"
sample1 <- ReadParseBio(sample1_path)
# Check to see if empty gene names are present
table(rownames(sample1) == "")
rownames(sample1)[rownames(sample1) == ""] <- "unknown"
cell_meta1 <- read.csv(paste0(sample1_path, "/cell_metadata.csv"), row.names = 1)
obj1325 <- CreateSeuratObject(sample1, min.genes = 100, min.cells = 2, meta.data = cell_meta1)
obj1325@meta.data[, "group"] <- "Tet1con" ## Add new metadata variable called group
obj1325@meta.data[, "sample"] <- "Ms1325_Tet1con"
obj1325@meta.data$orig.ident <- factor(rep("obj1325", nrow(obj1325@meta.data)))
Idents(obj1325) <- obj1325@meta.data$orig.ident
SaveObject(obj1325, "obj1325.seurat")
obj1325 <- ReadObject("obj1325.seurat")


sample2_path <- "/Users/sbombin/Desktop/analysis_mm39.comb/Ms1373_IEC_Tet1con/DGE_filtered"
sample2 <- ReadParseBio(sample2_path)
# Check to see if empty gene names are present
table(rownames(sample2) == "")
rownames(sample2)[rownames(sample2) == ""] <- "unknown"
cell_meta2 <- read.csv(paste0(sample2_path, "/cell_metadata.csv"), row.names = 1)
obj1373 <- CreateSeuratObject(sample2, min.genes = 100, min.cells = 2, meta.data = cell_meta2)
obj1373@meta.data[, "group"] <- "Tet1con"
obj1373@meta.data[, "sample"] <- "Ms1373_Tet1con"
obj1373@meta.data$orig.ident <- factor(rep("obj1373", nrow(obj1373@meta.data)))
Idents(obj1373) <- obj1373@meta.data$orig.ident
SaveObject(obj1373, "obj1373.seurat")
obj1373 <- ReadObject("obj1373.seurat")

sample3_path <- "/Users/sbombin/Desktop/analysis_mm39.comb/Ms1339_IEC_Tet1iKO/DGE_filtered"
sample3 <- ReadParseBio(sample3_path)
# Check to see if empty gene names are present
table(rownames(sample3) == "")
rownames(sample3)[rownames(sample3) == ""] <- "unknown"
cell_meta3 <- read.csv(paste0(sample3_path, "/cell_metadata.csv"), row.names = 1)
obj1339 <- CreateSeuratObject(sample3, min.genes = 100, min.cells = 2, meta.data = cell_meta3)
obj1339@meta.data[, "group"] <- "Tet1iKO"
obj1339@meta.data[, "sample"] <- "Ms1339_Tet1iKO"
obj1339@meta.data$orig.ident <- factor(rep("obj1339", nrow(obj1339@meta.data)))
Idents(obj1339) <- obj1339@meta.data$orig.ident
SaveObject(obj1339, "obj1339.seurat_obj")
obj1339 <- ReadObject("obj1339.seurat_obj")

sample4_path <- "/Users/sbombin/Desktop/analysis_mm39.comb/Ms1361_IEC_Tet1iKO/DGE_filtered"
sample4 <- ReadParseBio(sample4_path)
# Check to see if empty gene names are present
table(rownames(sample4) == "")
rownames(sample4)[rownames(sample4) == ""] <- "unknown"
cell_meta4 <- read.csv(paste0(sample4_path, "/cell_metadata.csv"), row.names = 1)
obj1361 <- CreateSeuratObject(sample4, min.genes = 100, min.cells = 2, meta.data = cell_meta4)
obj1361@meta.data[, "group"] <- "Tet1iKO"
obj1361@meta.data[, "sample"] <- "Ms1361_Tet1iKO"
obj1361@meta.data$orig.ident <- factor(rep("obj1361", nrow(obj1361@meta.data)))
Idents(obj1361) <- obj1361@meta.data$orig.ident
SaveObject(obj1361, "obj1361.seurat_obj")
obj1361 <- ReadObject("obj1361.seurat_obj")

## Merge Seurat objects
pbmc <- merge(obj1325, y = c(obj1373, obj1339, obj1361), add.cell.ids = c("Ms1325", "Ms1373", "Ms1339", "Ms1361"), project = "IEC")

pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 5000)
## Scaling the data
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
## LDR
pbmc <- RunPCA(pbmc)
##Find optimal PCA dimensions
npcs <- fgetNPCs(pbmc,MIN_CUM_SD)

p01 <- ElbowPlot(pbmc)
pbmc <- JackStraw(pbmc, num.replicate = 100, dims = npcs+5)
pbmc <- ScoreJackStraw(pbmc, dims = 1:(npcs+5))
p02 <- JackStrawPlot(pbmc, dims = 1:(npcs+5))

#pbmc <- JackStraw(pbmc, num.replicate = 100, dims = 35)
#pbmc <- ScoreJackStraw(pbmc, dims = 1:35)
#JackStrawPlot(pbmc, dims = 1:35)
SaveFigure(p01, "ElbowPlot", width = 10, height = 8)
SaveFigure(p02, "JackStrawPlot", width = 10, height = 8)

##Clustering
pbmc <- FindNeighbors(pbmc, dims = 1:npcs)
pbmc <- FindClusters(pbmc, resolution = 0.7)

pbmc <- RunUMAP(pbmc, dims = 1:npcs)
pbmc <- RunTSNE(pbmc, dims = 1:npcs)

##Save
SaveObject(pbmc, "seuratobj_IEC")
pbmc <- ReadObject("seuratobj.IEC.scale_all")

##Markers
stem.list <- list(c('Lgr5', 'Ascl2', 'Slc12a2', 'Axin2', 'Olfm4', 'Gkn3'))
stem <- c('Lgr5', 'Ascl2', 'Slc12a2', 'Axin2', 'Olfm4', 'Gkn3')

cell_cycle <- c('Mki67', 'Cdk4', 'Mcm5', 'Mcm6', 'Pcna')
cell_cycle.list <- list(c('Mki67', 'Cdk4', 'Mcm5', 'Mcm6', 'Pcna'))

enterocyte <- c('Alpi', 'Apoa1', 'Apoa4', 'Fabp1')
enterocyte.list <- list(c('Alpi', 'Apoa1', 'Apoa4', 'Fabp1'))

goblet <- c('Muc2', 'Clca3', 'Tff3', 'Agr2')
goblet.list <- list(c('Muc2', 'Clca3', 'Tff3', 'Agr2'))

paneth <- c('Lyz1', 'Defa17', 'Defa22', 'Defa24', 'Ang4')
paneth.list <- list(c('Lyz1', 'Defa17', 'Defa22', 'Defa24', 'Ang4'))

enteroendocrine <- c('Chga', 'Chgb', 'Tac1', 'Tph1', 'Neurog3')
enteroendocrine.list <- list(c('Chga', 'Chgb', 'Tac1', 'Tph1', 'Neurog3'))

tuft  <- c('Dclk1', 'Trpm5', 'Gfi1b', 'Il25')
tuft.list  <- list(c('Dclk1', 'Trpm5', 'Gfi1b', 'Il25'))

sec_pro <-c('Atoh1', 'Dll1', 'Sox4', 'Spdef', 'Gfi1', 'Mycn')
sec_pro.list <-list(c('Atoh1', 'Dll1', 'Sox4', 'Spdef', 'Gfi1', 'Mycn'))

object <- AddModuleScore(object = pbmc, features = stem.list, name = "Stem")
object <- AddModuleScore(object = object, features = cell_cycle.list, name = "Cell_cycle")
object <- AddModuleScore(object = object, features = enterocyte.list, name = "Enterocyte")
object <- AddModuleScore(object = object, features = goblet.list, name = "Goblet")
object <- AddModuleScore(object = object, features = paneth.list, name = "Paneth")
object <- AddModuleScore(object = object, features = enteroendocrine.list, name = "Enteroendocrine")
object <- AddModuleScore(object = object, features = tuft.list, name = "Tuft")
object <- AddModuleScore(object = object, features = sec_pro.list, name = "Sec_pro")

p1 <- DimPlot(pbmc, reduction = "umap", label = TRUE)
p2 <- DimPlot(pbmc, label = TRUE, group.by = "sample", reduction = "umap", repel = TRUE)
p3 <- DimPlot(pbmc, label = FALSE, group.by = "group", reduction = "umap")
p4 <- FeaturePlot(object = object, features = "Stem1")
p5 <- FeaturePlot(object = object, features = "Cell_cycle1")
p6 <- FeaturePlot(object = object, features = "Enterocyte1")
p7 <- FeaturePlot(object = object, features = "Goblet1")
p8 <- FeaturePlot(object = object, features = "Paneth1")
p9 <- FeaturePlot(object = object, features = "Enteroendocrine1")
p10 <- FeaturePlot(object = object, features = "Tuft1")
p11 <- FeaturePlot(object = object, features = "Sec_pro1")

##Separate Plots
p14 <- FeaturePlot(object = pbmc, features = stem)
p15 <- FeaturePlot(object = pbmc, features = cell_cycle)
p16 <- FeaturePlot(object = pbmc, features = enterocyte)
p17 <- FeaturePlot(object = pbmc, features = goblet)
p18 <- FeaturePlot(object = pbmc, features = paneth)
p19 <- FeaturePlot(object = pbmc, features = enteroendocrine)
p20 <- FeaturePlot(object = pbmc, features = tuft)
p21 <- FeaturePlot(object = pbmc, features = sec_pro)

##Save Plots
SaveFigure(p1, "Umap_Clusters", width = 10, height = 8)
SaveFigure(p2, "Umap_Sample", width = 10, height = 8)
SaveFigure(p3, "Umap_Group", width = 10, height = 8)
SaveFigure(p4, "Umap_Stem", width = 10, height = 8)
SaveFigure(p5, "Umap_Cell_cycle", width = 10, height = 8)
SaveFigure(p6, "Umap_Enterocyte", width = 10, height = 8)
SaveFigure(p7, "Umap_Goblet", width = 10, height = 8)
SaveFigure(p8, "Umap_Paneth", width = 10, height = 8)
SaveFigure(p9, "Umap_Enteroendocrine", width = 8, height = 8)
SaveFigure(p10, "Umap_Tuft", width = 10, height = 8)
SaveFigure(p11, "Umap_Sec_pro", width = 10, height = 8)

SaveFigure(p14, "Umap_Stem_Markers", width = 10, height = 8)
SaveFigure(p15, "Umap_Cell_cycle_Markers", width = 10, height = 8)
SaveFigure(p16, "Umap_Enterocyte_Markers", width = 10, height = 8)
SaveFigure(p17, "Umap_Goblet_Markers", width = 10, height = 8)
SaveFigure(p18, "Umap_Paneth_Markers", width = 10, height = 8)
SaveFigure(p19, "Umap_Enteroendocrine_Markers", width = 10, height = 8)
SaveFigure(p20, "Umap_Tuft_Markers", width = 10, height = 8)
SaveFigure(p21, "Umap_Sec_pro_Markers", width = 10, height = 8)

SaveFigure((p1 + p2), "Dimplot.Intgrtn", width = 10, height = 10)


##Find Markers
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.8, test.use = "MAST")
pbmc.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

pbmc.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
p4 <- DoHeatmap(pbmc, features = top10$gene) + NoLegend()
SaveFigure(p4, "Heatmap2_top10_all", width = 20, height = 20)

write.csv(pbmc.markers,"/Users/sbombin/Desktop/analysis_mm39.comb/Seurat_Analysis/IEC//All_Clusters_DE.csv", row.names = TRUE)

#cluster0.markers <- FindMarkers(pbmc, ident.1 = 0, logfc.threshold = 0.8, test.use = "roc", only.pos = TRUE)
pbmc.markers2 <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.8, test.use = "roc")
write.csv(pbmc.markers2,"/Users/sbombin/Desktop/analysis_mm39.comb/Seurat_Analysis/IEC//Markers_Power_AllClusters.csv", row.names = TRUE)

##Sox9 Plots
p30 <- FeaturePlot(object = pbmc, features = "Sox9")
p31 <- VlnPlot(pbmc, features = "Sox9", pt.size = 0, combine = FALSE)
SaveFigure((p30 + p31), "Plot_Sox9", width = 10, height = 10)

new.cluster.ids <- c("Enterocyte Progenitor (late) 1", "Mature Enterocyte (proximal) 1", "Mature Enterocyte (proximal) 2", "TA", "Mature Enterocyte (distal)",
                     "Immature Enterocyte (proximal)", "Mature Enterocyte (proximal) 3", "Enterocyte Progenitor (early) 1", "ISC", "Enterocyte Progenitor (early) 2",
                     "Enterocyte Progenitor (late) 2", "Goblet", "Tuft", "Enteroendocrine", "Mature Enterocyte (proximal) 4", "Paneth")

names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)

p100 <- DimPlot(pbmc, reduction = "umap", label = FALSE, pt.size = 0.5, repel = TRUE) 
SaveFigure(p100, "UMAP_Renamed_Clusters_unlabel", width = 12, height = 8)

## Save AnnDat format
library(SeuratDisk)
setwd('/Users/sbombin/Desktop/analysis_mm39.comb/Seurat_Analysis/IEC/')
SaveH5Seurat(pbmc, filename = "IEC_all_cellassign.h5Seurat")
Convert("IEC_all_cellassign.h5Seurat", dest = "h5ad")

# Create function to get conserved markers for any given cluster
library(purrr)
library(tibble)
library(multtest)
library(metap)

get_conserved <- function(cluster){
  FindConservedMarkers(pbmc,
                       ident.1 = cluster, min.pct = 0.1, logfc.threshold = 0.8,
                       grouping.var = "group", test.use = "MAST", 
                       only.pos = FALSE) %>%
    rownames_to_column(var = "gene") %>%
    cbind(cluster_id = cluster, .)
}

# Iterate function across desired clusters
conserved_markers <- map_dfr(c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15), get_conserved)
write.csv(conserved_markers,"/Users/sbombin/Desktop/analysis_mm39.comb/Seurat_Analysis/IEC//DE-MAST_conserved_markers.csv", row.names = FALSE)

## DE between Control and KO
de_markers1 <- FindMarkers(object = pbmc,ident.1 = "Tet1iKO", group.by = "group", min.pct = 0.1, logfc.threshold = 0.8,
                       test.use = "MAST", only.pos = FALSE,subset.ident ="0")

get_de <- function(cluster){
  FindMarkers(object = pbmc, ident.1 = "Tet1iKO",
                       subset.ident = cluster, min.pct = 0.1, logfc.threshold = 0.8,
                       group.by = "group", test.use = "MAST", 
                       only.pos = FALSE) %>%
    rownames_to_column(var = "gene") %>%
    cbind(cluster_id = cluster, .)
}

# Iterate function across desired clusters
de_markers <- map_dfr(c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15), get_de)

write.csv(de_markers,"/Users/sbombin/Desktop/analysis_mm39.comb/Seurat_Analysis/IEC//DE-MAST_Tet1iKO_vs_Control.csv", row.names = FALSE)

