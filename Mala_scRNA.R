library(readr)
library(Seurat)
library(dplyr)
library(RColorBrewer)

#BiocManager::install("tanaylab/metacell")

setwd("/Users/sbombin/Desktop/EICC/Projects/Mala_scRNA/")
data_path <- "/Users/sbombin/Desktop/EICC/Projects/Mala_scRNA/"
fig_path <- "/Users/sbombin/Desktop/EICC/Projects/Mala_scRNA/"

#GSE117156_series_matrix <- read_delim("GSE117156_series_matrix.txt", delim = "\t", escape_double = FALSE, trim_ws = TRUE)

metadata <- read_delim("GSE117156_metadata.txt", delim = "\t", escape_double = FALSE,  trim_ws = TRUE)
metadata <-  metadata[,-8]
#metadata <-dplyr::filter(metadata, Number_of_cells >0)
row.names(metadata) <- metadata$well
#wells <- metadata$well
empty.wells <-dplyr::filter(metadata, Number_of_cells <1)

### Extracting string between first and second underscore
metadata$Experiment_ID <- sub("^[^_]*_([^_]*).*", "\\1", metadata$Experiment_ID)

#mtx <- ReadMtx(mtx = "~/Downloads/GSE117156_series_matrix.txt")

cell <- read_csv("~/Desktop/EICC/Projects/Mala_scRNA/CellInfo.csv") 
cell <-  cell[,-1]
rownames(cell) <- cell$well

#unique(cell$sample_name)

filepath <- read_csv("~/Desktop/EICC/Projects/Mala_scRNA/FilePath.csv")
filepath  <- dplyr::filter(filepath, grepl('MM', sample_name))
filepath$filename <- gsub('./DataFiles', "GSE117156_RAW", filepath$filename)

#table(filepath$sample_name)

ab3058 <- read.delim("~/Desktop/EICC/Projects/Mala_scRNA/GSE117156_RAW/GSM3272410_AB3058.txt.gz")
ab3449 <- read.delim("GSE117156_RAW/GSM3272431_AB3449.txt.gz")

genes <- rownames(ab3449)

#mtx.ab3058 <- ReadMtx(mtx = "AB3058", cells = "cell", features = "genes", cell.column =2)

obj.ab3058 <- CreateSeuratObject(counts = ab3058, project = "AB3058")

obj <- CreateSeuratObject(counts = data[[1]], min.cells = 3, min.features = 200, project = "AB3449", meta.data = cell)

obj.ab3449 <- CreateSeuratObject(counts = data[["AB3449"]], meta.data = cell)

## Cells (Wells) should be in rownames
obj.ab3058 <- AddMetaData(object = obj.ab3058, metadata = cell)

head(obj.ab3058@meta.data)

#data.files <- list.files()
data.files <- filepath$filename

for(i in 1:length(data.files)) {                              # Head of for-loop
  assign(paste0("data", i),                                   # Read and store data frames
         read.delim(paste0("GSE117156_RAW/",data.files[i])))
}


data = lapply(data.files, read.delim)
sample.ids <- filepath$sample_ID
names(data) <- sample.ids

object = lapply(data, CreateSeuratObject)


### All Filtering easier to do before creating Seurat object
## For metagene just make vector of wells with sum < 100 and filter out same as empty wells

## Remove empty wells
# Save if in the vector
a = AB3449[,(colnames(AB3449) %in% wells)]
# Remove if in the vector
a = AB3449[,!(colnames(AB3449) %in% wells)]

## Metagene
img <- read.table("immunoglobulin_genes.txt", quote="\"", comment.char="", stringsAsFactors=FALSE)
img <- img$V1

metagene <- colSums(ab3058[img,],na.rm=T) 
metagene <- metagene[ metagene > 100 ]

a = ab3058[,(colnames(ab3058) %in% names(metagene))]

b = a[,!(colnames(a) %in% rownames(empty.wells))]

### Remove cells/wells if immunoglobulin genes total expression < 100 
img.count <- lapply(data, function(x) colSums(x[img,],na.rm=T))
empty.img <- lapply(img.count, function(x) x[x < 100 ])
#full.img <- lapply(img.count, function(x) x[x > 100 ])
empty.img <- unlist(lapply(names(empty.img), function(x) names(empty.img[[x]])))

data.cleaned <- lapply(data, function(x) x[,!(colnames(x) %in% empty.img)])
data.cleaned <- lapply(data, function(x) x[,!(colnames(x) %in% empty.wells)])

### library size < 500 or > 10000, min cells 100 but in combined matrix

### Create Seurat Objects
project.name <- names(data.cleaned)
object2 = lapply(data, function(x) CreateSeuratObject(counts = x,min.cells = 3,min.features = 200))


### Combine all count tables from list to one dataframe
rn <- rownames(ab3449)
dat <- data.cleaned[[1]]
for(i in 2:length(data.cleaned)) {
  dat <- merge(dat, data.cleaned[[i]],  by= "row.names", all.x= F, all.y= F) [,-1]
  rownames(dat) <- rn
}

### Create Seuratobject from combined table
object <- CreateSeuratObject(counts = dat, min.cells = 100, min.features = 200)
object <- AddMetaData(object = object, metadata = metadata)

SaveObject(object, "SeuratObject_MM")

object[["percent.mt"]] <- PercentageFeatureSet(object, pattern = "^MT-")

p <- VlnPlot(object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1)
SaveFigure(p, "QC_Plot1", width = 10, height = 8)

p <- FeatureScatter(object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
SaveFigure(p, "QC_Plot2", width = 10, height = 8)

p <- FeatureScatter(object, feature1 = "nCount_RNA", feature2 = "percent.mt")
SaveFigure(p, "QC_Plot3", width = 10, height = 8)


object <- subset(object, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & nCount_RNA > 500 & nCount_RNA < 10000 & percent.mt < 5)

FeatureScatter(object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

### Integration 1 (batch correction) ### Out of memmory + too few cells in batches 
object.list <- SplitObject(object, split.by = "Amp_batch_ID")

#object.list[["AB3195"]] <- NULL
### Need to remove cell: WMC0076617 

# normalize and identify variable features for each dataset independently
object.list <- lapply(object.list, FUN = function(x) {
  x <- NormalizeData(x)
})

object.list <- lapply(object.list, FUN = function(x) {
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = object.list)
integ.anchors <- FindIntegrationAnchors(object.list = object.list, anchor.features = features,dims = 1:29,k.score=29,k.filter=200)
integ.combined <- IntegrateData(anchorset = integ.anchors,dims = 1:29,k.weight = 29)

### Integration 2 (batch correction)
object.list <- SplitObject(object, split.by = "Experiment_ID")

# normalize and identify variable features for each dataset independently
object.list <- lapply(object.list, FUN = function(x) {
  x <- NormalizeData(x)
})

object.list <- lapply(object.list, FUN = function(x) {
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = object.list)

integr.anchors <- FindIntegrationAnchors(object.list = object.list, anchor.features = features)
integr.combined <- IntegrateData(anchorset = integr.anchors)

SaveObject(integr.combined, "SeuratObject_Integrated")


integr <- ReadObject("SeuratObject_Integrated")
DefaultAssay(integr) <- "integrated"
integr <- ScaleData(integr, verbose = FALSE)
integr <- RunPCA(integr, verbose = FALSE)

integr <- JackStraw(integr, num.replicate = 100, dims = 50)
integr <- ScoreJackStraw(integr, dims = 1:50)
p <- JackStrawPlot(integr, dims = 1:50)
p
SaveFigure(p, "JackStrawPlot_PC50", width = 10, height = 8)

### Optimal PCAs = 43
integr <- FindNeighbors(integr, dims = 1:43)
## Clustree to find resolution
library(clustree)
clustree <- FindClusters(integr, resolution = c(0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.2,1.4))
p <- clustree(clustree)
SaveFigure(p, "Clustree", width = 8, height = 10)
## Cluster 
integr <- FindClusters(integr, resolution = 0.2)
integr <- RunUMAP(integr, dims = 1:43)
integr <- RunTSNE(integr, dims = 1:43)


Idents(integr) = integr$integrated_snn_res.0.2

p <- DimPlot(integr, reduction = "umap", label = TRUE)
p
SaveFigure(p, "UMAP_Clusters_0.2", width = 10, height = 8)

p <- DimPlot(integr, reduction = "tsne", label = TRUE)
p
SaveFigure(p, "TSNE_Clusters_0.2", width = 10, height = 8)


DefaultAssay(integr) <- "RNA"

etc <- read.table("etc_genes.txt", quote="\"", comment.char="", stringsAsFactors=FALSE)
etc <- etc$V1
#etc <- etc[etc %in% rownames(my.data)]
my.data=GetAssayData(integr, slot = "data", assay = "RNA") ## normalized counts data

etc <- etc[etc %in% rownames(my.data)]
y <- etc[!etc %in% rownames(my.data)]

Metagene0 <- colMeans(my.data[etc,],na.rm=T) 
Metagene1 <- colSums(my.data[etc,],na.rm=T) 
a <- as.data.frame(my.data)
f <- as.data.frame(rownames(my.data))
c <- as.data.frame(etc)
d <- my.data[rownames(my.data) %in% etc, ] %>% as.data.frame()

integr@meta.data$ETC <- Metagene0
integr@meta.data$ETC_Sum <- Metagene1

Idents(integr) = integr$integrated_snn_res.0.6

#FeaturePlot(integr,features="ETC",label=TRUE,cols = c("grey", "red"),reduction = "tsne")
p <- FeaturePlot(integr,features="ETC_Sum",label=TRUE,cols = c("grey", "red"),reduction = "tsne")

Idents(integr) = integr$integrated_snn_res.0.6

p1 <- DotPlot(integr, features = "ETC_Sum",dot.scale = 8, cols = c("lightgrey", "red")) & coord_flip() 
p1
SaveFigure(p1, "DotPlot_ETC_0.6", width = 10, height = 8)
a <- as.data.frame(p1$data)
write.csv(a,"Expression_ETC_0.6.csv", row.names = TRUE)


Idents(integr) = integr$integrated_snn_res.0.2
p2 <- DotPlot(integr, features = "ETC_Sum",dot.scale = 8, cols = c("lightgrey", "red")) & coord_flip() 
p2
b <- as.data.frame(p2$data)

### DE
# MAST
Idents(integr) = integr$integrated_snn_res.0.6
mast <-  FindMarkers(integr, only.pos = FALSE, test.use = "MAST", ident.1 = c(1,11),ident.2 = NULL,logfc.threshold = 0.1)
write.csv(mast,"DE_MAST_0.6.csv", row.names = TRUE)

negbinom <-  FindMarkers(integr, only.pos = FALSE, test.use = "negbinom", ident.1 = c(1,11),ident.2 = NULL,logfc.threshold = 0.25)
write.csv(mast,"DE_MAST_0.6.csv", row.names = TRUE)

c <- as.data.frame(rownames(mast))

heme <- read.table("heme.txt", quote="\"", comment.char="", stringsAsFactors=FALSE)
#heme <- heme$V1
iron <- read.table("iron.txt", quote="\"", comment.char="", stringsAsFactors=FALSE)
#iron <- iron$V1

de.heme <- mast[rownames(mast) %in% heme$V1, ]
de.iron <- mast[rownames(mast) %in% iron$V1, ]

de.heme2 <- negbinom[rownames(negbinom) %in% heme$V1, ]
de.iron2 <- negbinom[rownames(negbinom) %in% iron$V1, ]

library(EnhancedVolcano)
p<-EnhancedVolcano(mast,
                   lab = rownames(mast), x = 'avg_log2FC', y = 'p_val_adj', title = 'ETC_High vs ETC_Low', subtitle = "Clusters 1 and 11 vs All Others",
                   legendLabels=c('Not sig.','Log (base 2) FC','q-value', 'q-value & Log (base 2) FC'),
                   drawConnectors = TRUE,
                   pCutoff = 0.01,
                   FCcutoff = 0.5,
                   pointSize = 3.0,
                   axisLabSize = 12,
                   legendLabSize = 10,
                   legendIconSize = 3.0,
                   labSize = 3.0)
p
ggsave(paste0('VolcanoPlot_ETC_High-Low.png'),plot=p,dpi=300, scale = 1.5)

z <- as.data.frame(rownames(mast))
z2 <- as.data.frame(rownames(wilcox))

y <- table(integr$integrated_snn_res.0.6) %>% as.data.frame()

SaveObject(integr, "SeuratObject_Integr_Processed")

integr <- ReadObject("SeuratObject_Integr_Processed")

p <- VlnPlot(integr, features = "ETC_Score", group.by = "Experiment_ID")
p

SaveFigure(p, "VlnPlot_ETC_Score_Patients", width = 10, height = 8)


integr <- AddModuleScore(integr, features = list(etc), name="ETC_enriched")

FeaturePlot(integr,
            features = "ETC_enriched1", label = TRUE, repel = TRUE,reduction = "tsne") +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

p1 <- FeaturePlot(integr,
            features = "ETC_Sum", label = TRUE, repel = TRUE,reduction = "tsne") +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")))

p2 <- FeaturePlot(integr,
            features = "ETC_Score", label = TRUE, repel = FALSE,reduction = "tsne") +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")))
p2

SaveFigure(p1, "FeaturePlot_ETC_Expression", width = 10, height = 8)

SaveFigure(p2, "FeaturePlot_ETC_Score", width = 10, height = 8)

## Add Z-score
integr@meta.data$ETC_Score <- (integr$ETC_Sum  - mean(integr$ETC_Sum )) / sd(integr$ETC_Sum )


p <- VlnPlot(integr, features = "FTL", group.by = "Experiment_ID")
p
SaveFigure(p, "VlnPlot_FTL_Patients", width = 10, height = 8)

p <- VlnPlot(integr, features = "FTL")
p
SaveFigure(p, "VlnPlot_FTL_Clusters", width = 10, height = 8)

p <- VlnPlot(integr, features = "FTH1", group.by = "Experiment_ID")
p
SaveFigure(p, "VlnPlot_FTH1_Patients", width = 10, height = 8)

p <- VlnPlot(integr, features = "FTH1")
p
SaveFigure(p, "VlnPlot_FTH1_Clusters", width = 10, height = 8)



p1 <- FeaturePlot(integr,
                  features = "FTL", label = TRUE, repel = TRUE,reduction = "tsne") +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")))

p2 <- FeaturePlot(integr,
                  features = "FTH1", label = TRUE, repel = FALSE,reduction = "tsne") +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")))
p2

SaveFigure(p1, "FeaturePlot_FTL", width = 10, height = 8)

SaveFigure(p2, "FeaturePlot_FTH1", width = 10, height = 8)



