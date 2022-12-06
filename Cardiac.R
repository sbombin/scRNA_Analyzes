library(readr)
library(dplyr)
library(Seurat)
library(Matrix)
library(ggplot2)
library(MAST)
library(DESeq2)
library(patchwork)

setwd('/Users/sbombin/Desktop/EICC/Projects/Kwong_Jennifer/')
data_path <- "/Users/sbombin/Desktop/EICC/Projects/Kwong_Jennifer/"
fig_path <- "/Users/sbombin/Desktop/EICC/Projects/Kwong_Jennifer/" #Cardiac_Cells_All/

#full <- ReadObject("gene_count_cleaned")
#meta <- read_csv("cell_annotate.csv")
#cardiac <- read_csv("cardiac_cells.csv")
#rownames(meta) <- meta$sample
#a <- meta[,"Main_cell_type", drop=FALSE]
#rownames(a) <- rownames(meta)

### Target Gene ENSMUSG00000003528.14, protein_coding, Slc25a1

object <- ReadObject("seurobj_cardiac_subtrj_raw")

object[["percent.mt"]] <- PercentageFeatureSet(object, pattern = "^MT-")
p <- VlnPlot(object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)
SaveFigure(p, "QC_Plot1", width = 10, height = 8)

p <- FeatureScatter(object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
SaveFigure(p, "QC_Plot2", width = 10, height = 8)

#p <- FeatureScatter(object, feature1 = "nCount_RNA", feature2 = "percent.mt")

## Might be needed to subset  outliers with high counts and/or mitochondrial content
object<- subset(object, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 5)
#object <- subset(object, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & nCount_RNA < 20000 & percent.mt < 5)

## Normalize 
object <- NormalizeData(object, normalization.method = "LogNormalize", scale.factor = 10000)
## Variable Features
object <- FindVariableFeatures(object, selection.method = "vst", nfeatures = 2000)
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(object), 10)
# plot variable features with and without labels
p <- VariableFeaturePlot(object)
p <- LabelPoints(plot = p, points = top10, repel = TRUE)
SaveFigure(p, "Variable_Features", width = 10, height = 8)
## Scale
all.genes <- rownames(object)
object <- ScaleData(object, features = all.genes)
## PCA
object <- RunPCA(object, features = VariableFeatures(object = object))
## Find optimal PCA dimensions
npcs <- fgetNPCs(object,MIN_CUM_SD)
## Find Optimal number of PCs
p <- ElbowPlot(object)
SaveFigure(p, "ElbowPlot", width = 10, height = 8)
object <- JackStraw(object, num.replicate = 100, dims = 50)
object <- ScoreJackStraw(object, dims = 1:50)
p <- JackStrawPlot(object, dims = 1:50)
p
SaveFigure(p, "JackStrawPlot_PC50", width = 10, height = 8)

object <- FindNeighbors(object, dims = 1:39)
## Clustree to find resolution
library(clustree)
clustree <- FindClusters(object, resolution = c(0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.2,1.4))
p <- clustree(clustree)
SaveFigure(p, "Clustree", width = 8, height = 10)
## Cluster 
object <- FindClusters(object, resolution = 0.6)
object <- RunUMAP(object, dims = 1:39)

p <- DimPlot(object, reduction = "umap", label = TRUE)
SaveFigure(p, "UMAP_Clusters", width = 10, height = 8)

p <- DimPlot(object, label = FALSE, group.by = "Sub_trajectory_name", reduction = "umap", repel = FALSE)
p
SaveFigure(p, "UMAP_SubTragectory", width = 10, height = 8)


p <- DimPlot(object, label = FALSE, group.by = "development_stage", reduction = "umap", repel = FALSE)
SaveFigure(p, "UMAP_Develop_Stage", width = 10, height = 8)

p1 <- FeaturePlot(object, features = "Slc25a1", split.by = "development_stage", max.cutoff = 5, cols = c("grey", "red"))
SaveFigure(p1, "Slc25a1_DevelopStage_Cardiac", width = 20, height = 8)

VlnPlot(object, features = "Slc25a1", split.by = "development_stage")
VlnPlot(object, features = "Slc25a1")

p2 <- FeaturePlot(object, features = "Slc25a1", split.by = "development_stage", max.cutoff = 5, cols = c("grey", "red"), combine = FALSE)
p2 <- wrap_plots(plots = p2, ncol = 2)
SaveFigure(p2, "Slc25a1_DevelopStage_Cardiac2", width = 10, height = 10)


p3 <- FeaturePlot(object, features = "Slc25a1", cols = c("grey", "red"))
SaveFigure(p3, "UMAP_Slc25a1_Cardiac", width = 10, height = 8)

SaveObject(object, "seurobj_cardiac_subtrj_raw_processed")

mast <- FindAllMarkers(object, only.pos = TRUE, test.use = "MAST")
mast2 <-  FindAllMarkers(clean, only.pos = FALSE, test.use = "MAST") #clean08

roc <- FindAllMarkers(object, only.pos = TRUE, test.use = "roc") #object06
roc2 <- FindAllMarkers(clean, only.pos = TRUE, test.use = "roc") #clean08
roc3 <- FindAllMarkers(clean, only.pos = TRUE, test.use = "roc") #clean06
roc4 <- FindAllMarkers(clean, only.pos = TRUE, test.use = "roc") #clean09



write.csv(mast2,"Cleaned/DE-MAST_All_Markers_Cardiac.csv", row.names = TRUE)
write.csv(roc2,"Cleaned/ROC_Markers_Cardiac.csv", row.names = TRUE)

### Part2 Labeling Clusters
###########################
my.data=GetAssayData(object, slot = "data")
df <- as.data.frame(my.data)
genes <- unique(row.names(df)) %>% as.data.frame()

object <- ReadObject("seurobj_cardiac_subtrj_raw_processed")

markers <- c('C1qa', 'Cav1','Cdh5','Dcn','Ehd3','Fabp4','Gsn','Myh11', 'Npr3', 'Pecam1', 'Tie1', 'Tnni3') # 'Hbb−b1'
p <- DotPlot(object, features = markers) & coord_flip()
p
VlnPlot(object, features = markers)

fb <- c('Pdgfra', 'Col1a1', 'Tcf21') # Fibroblast
end <- c('Ly6c1', 'Egfl7', 'Cdh5', 'Mgll', 'Emcn', 'Kdr', 'Pecam1') # Endothelial
cm <- c('Tnnt2', 'Tpma', 'Tnnc1', 'Tnni3', 'Actc1', 'Ryr2')  # Cardiomyocytes
ep <- c('Wt1','Sema3d','Tbx18','Scx','Tcf21', 'Aldh1a2','Upk3b','Upk1b','Tmem255a', 'Gpm6a', 'Dmkn')  # Epicardial

VlnPlot(object, features = cm)

# Cluster 1 = Endothelial
# Cluster 3-4 and maybe 8-9 = Cardiomyocytes
# Clusters 0,2,6 = Epicardial (likely)

#######################################################
#### To use different resolutions for downstream analysis
Idents(object) = object$RNA_snn_res.1

DotPlot(object, features = "Myh11") & coord_flip() 


VlnPlot(object, features = "Myh11")


obj13 <- subset(object, development_stage == 13.5)
obj10 <-subset(object, development_stage == 10.5)

p <- DotPlot(obj13, features = markers) & coord_flip()
p

p <- DotPlot(obj10, features = fb) & coord_flip()
p

VlnPlot(obj13, features = end)


#cardiomyocytes
#endocardial
#epicardial

cardiac <- cell_annotate %>% dplyr::filter_all(any_vars(grepl("ardi", .)))

unique(cardiac$Sub_trajectory_name)

heart <- dplyr::filter(cell_annotate, Main_cell_type == "Cardiac muscle lineages")
table(heart$Sub_trajectory_name)