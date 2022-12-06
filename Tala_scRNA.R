### IL13ra2

setwd("/Users/sbombin/Desktop/EICC/Projects/Tala_AMA23076/")
data_path <- "/Users/sbombin/Desktop/EICC/Projects/Tala_AMA23076/"
fig_path <- "/Users/sbombin/Desktop/EICC/Projects/Tala_AMA23076/Plots/"

## Import Data
s01 <- Read10X(data.dir = "s01")
s02 <- Read10X(data.dir = "s02")
s03 <- Read10X(data.dir = "s03")
s04 <- Read10X(data.dir = "s04")
s05 <- Read10X(data.dir = "s05")
s06 <- Read10X(data.dir = "s06")
s07 <- Read10X(data.dir = "s07")
s08 <- Read10X(data.dir = "s08")
s09 <- Read10X(data.dir = "s09")
s10 <- Read10X(data.dir = "s10")
s11 <- Read10X(data.dir = "s11")
s12 <- Read10X(data.dir = "s12")

## Create Seurat Objects
obj01 <- CreateSeuratObject(counts = s01, project = "H23",min.cells = 3, min.features = 200)
obj02 <- CreateSeuratObject(counts = s02, project = "H460",min.cells = 3, min.features = 200)
obj03 <- CreateSeuratObject(counts = s03, project = "CK10341",min.cells = 3, min.features = 200)
obj04 <- CreateSeuratObject(counts = s04, project = "CK10151",min.cells = 3, min.features = 200)
obj05 <- CreateSeuratObject(counts = s05, project = "CK5494",min.cells = 3, min.features = 200)
obj06 <- CreateSeuratObject(counts = s06, project = "CK9121",min.cells = 3, min.features = 200)
obj07 <- CreateSeuratObject(counts = s07, project = "H1299",min.cells = 3, min.features = 200)
obj08 <- CreateSeuratObject(counts = s08, project = "H1915",min.cells = 3, min.features = 200)
obj09 <- CreateSeuratObject(counts = s09, project = "H1975",min.cells = 3, min.features = 200)
obj10 <- CreateSeuratObject(counts = s10, project = "A549",min.cells = 3, min.features = 200)
obj11 <- CreateSeuratObject(counts = s11, project = "SKMES1",min.cells = 3, min.features = 200)
obj12 <- CreateSeuratObject(counts = s12, project = "EUH3174",min.cells = 3, min.features = 200)

obj01[["Passage"]] <- "p5"
obj02[["Passage"]] <- "p2"
obj03[["Passage"]] <- "p2"
obj04[["Passage"]] <- "p4"
obj05[["Passage"]] <- "p4"
obj06[["Passage"]] <- "p4"
obj07[["Passage"]] <- "p4"
obj08[["Passage"]] <- "p4"
obj09[["Passage"]] <- "p3"
obj10[["Passage"]] <- "p4"
obj11[["Passage"]] <- "p3"
obj12[["Passage"]] <- "p6"

objects.raw <- c(obj01,obj02,obj03,obj04,obj05,obj06,obj07,obj08,obj09,obj10,obj11,obj12)
names(objects.raw) <- c("H23","H460","CK10341","CK10151","CK5494","CK9121","H1299","H1915","H1975","A549","SKMES1","EUH3174")

### QC
summary(norm1$nCount_RNA)
summary(norm1$nFeature_RNA)
summary(norm1$percent.mt)

## Add mtDNA ratio

for (x in seq_along(objects.raw)) {
  objects.raw[[x]][["percent.mt"]] = PercentageFeatureSet(objects.raw[[x]], pattern = "^MT-")
}

nCount <- lapply(objects.raw, function(x) summary(x$nCount_RNA))
nFeature <- lapply(objects.raw, function(x) summary(x$nFeature_RNA))
nMT <- lapply(objects.raw, function(x) summary(x$percent.mt))

score.count <- lapply(objects.raw, function(x) mean(x$nCount_RNA) + 3*sd(x$nCount_RNA))
#score.count2 <- lapply(objects.raw, function(i) 10^(mean(log10(i$nCount_RNA)) + 2*sd(log10(i$nCount_RNA))))
score.feature <- lapply(objects.raw, function(x) mean(x$nFeature_RNA) + 3*sd(x$nFeature_RNA))

#mad.count <- lapply(objects.raw, function(x) median(x$nCount_RNA) + 3*mad(x$nCount_RNA))
mad.count <- lapply(objects.raw, function(i) 10^(median(log10(i$nCount_RNA)) + 2*mad(log10(i$nCount_RNA))))
#a <- obj01@meta.data$nCount_RNA %>% log10()  
#b <-  10^(median(a) + 3*mad(a))
mad.feature <- lapply(objects.raw, function(i) 10^(median(log10(i$nFeature_RNA)) + 2*mad(log10(i$nFeature_RNA))))

#count.feature.ls <- obj02@meta.data[, c("nCount_RNA", "nFeature_RNA")]
#count.feature.ls %<>% map(log10) %>% map(~c(10^(mean(.x) + 3*sd(.x)), 10^(mean(.x) - 3*sd(.x))))

### Combine Seurat objects, requires too much memory
#sc.list = objects.raw
#allen = sc.list[[1]]
#sc.list[[1]] = NULL
#for (i in seq(length(sc.list))) {
#  allen = merge(allen, sc.list[[1]])
#  sc.list[[1]] = NULL
#}
###############


p3 <- FeatureScatter(objects.raw[[i]], feature1 = "nCount_RNA", feature2 = "percent.mt")
#SaveFigure(p3, "QC_Norm_Plot2", width = 10, height = 8)
ggsave(paste0('QC_2b.png'),plot = p3, dpi=300, scale =1.5)
p4 <- FeatureScatter(norm2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
SaveFigure(p4, "QC_Norm_Plot3", width = 10, height = 8)

### Plot
lapply(seq_along(objects.raw), function(i) {
  ggplot(objects.raw[[i]], aes(x=dpt_pseudotime, y=expression,color = colors)) + 
    geom_point()+ geom_smooth(method = 'loess',formula = y ~ x)+
    ggtitle(names(objects.raw)[i])
  ggsave(paste0('Scatter_',names(objects.raw)[i],'.png'),dpi=300, scale =1.5)})

lapply(seq_along(objects.raw), function(i) {
  p1 <- VlnPlot(objects.raw[[i]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  SaveFigure(p1, paste0('QC_1_',names(objects.raw)[i]), width = 10, height = 8)})

lapply(seq_along(objects.raw), function(i) {
  p1 <- FeatureScatter(objects.raw[[i]], feature1 = "nCount_RNA", feature2 = "percent.mt")
  SaveFigure(p1, paste0('QC_2_',names(objects.raw)[i]), width = 8, height = 6)
  p2 <- FeatureScatter(objects.raw[[i]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  SaveFigure(p2, paste0('QC_3_',names(objects.raw)[i]), width = 8, height = 6)})

## Add groups
group <- c("g1","g3","g3","g2","g3","g2","g1","g1","g2","g3","g1","g1")

for (x in seq_along(objects.raw)) {
  objects.raw[[x]][["group"]] <- group[x]
}

## Filter outliers

for (x in seq_along(objects.raw)) {
  objects.raw[[x]] <- subset(objects.raw[[x]], subset = nFeature_RNA > 200 & nFeature_RNA < score.feature[[x]] & percent.mt < 15 &
                               nCount_RNA >500 & nCount_RNA < score.count[[x]])
}

####
object.list = objects.raw
rm(objects.raw)
# normalize and identify variable features for each dataset independently
object.list <- lapply(X = object.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# select features that are repeatedly variable across datasets for integration
#features <- SelectIntegrationFeatures(object.list = object.list)
#library(future)
#plan()
#plan("multiprocess", workers = 4)
#plan()
#plan("sequential")
Sys.setenv('R_MAX_VSIZE'=64000000000)
### Perform Integration
#object.anchors <- FindIntegrationAnchors(object.list = object.list, anchor.features = features)
SaveObject(object.list, "SeuratObject")

object <- ReadObject("SeuratIntegr")
DefaultAssay(object) <- "integrated"

### Optimal PCs number
p <- JackStrawPlot(object, dims = 1:50)
p
SaveFigure(p, "JackStrawPlot_PC50", width = 10, height = 8)

p <- ElbowPlot(object, ndims = 50)
p
SaveFigure(p, "ElbowPlot_PC50", width = 10, height = 8)

npcs <- fgetNPCs(object,90)

table(object$orig.ident)

object <- FindNeighbors(object, reduction = "pca", dims = 1:50)


## Clustree
library(clustree)
clustree <- FindClusters(object, resolution = c(0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.2,1.4))
p <- clustree(clustree)
p
SaveFigure(p, "Clustree", width = 8, height = 10)

####
DefaultAssay(clustree) <- "RNA"
p1 <- VlnPlot(clustree, features = "IL13RA2", group.by = "integrated_snn_res.0.3")
SaveFigure(p1, "VlnPlot_IL13RA2_res03", width = 10, height = 8)
p2 <- DotPlot(clustree, features = "IL13RA2", group.by = "integrated_snn_res.0.3")
SaveFigure(p2, "DotPlot_IL13RA2_res03", width = 10, height = 8)
p3 <- VlnPlot(clustree, features = "IL13RA2", group.by = "integrated_snn_res.0.8")
SaveFigure(p3, "VlnPlot_IL13RA2_res08", width = 10, height = 8)
p4 <- DotPlot(clustree, features = "IL13RA2", group.by = "integrated_snn_res.0.8")
SaveFigure(p4, "DotPlot_IL13RA2_res08", width = 10, height = 8)
p2+p4
#########

## Cluster 
object <- FindClusters(object, resolution = 0.8)
summary(object$percent.mt)

object <- RunUMAP(object,reduction = "pca", dims = 1:50)
object <- RunTSNE(object, reduction = "pca", dims = 1:50)

p <- DimPlot(object, reduction = "umap", label = TRUE)
p
SaveFigure(p, "UMAP_Clusters_res08", width = 10, height = 8)

p <- DimPlot(object, reduction = "tsne", label = TRUE)
p
SaveFigure(p, "tSNE_Clusters_res08", width = 10, height = 8)

p <- DimPlot(object, reduction = "tsne", label = FALSE, group.by = "orig.ident")
p
SaveFigure(p, "tSNE_Sample", width = 10, height = 8)

p <- DimPlot(object, reduction = "tsne", label = FALSE, split.by = "orig.ident",group.by = "integrated_snn_res.0.8",ncol = 3)
p
SaveFigure(p, "tSNE_SplitSample_res08", width = 15, height = 15)

p <- DimPlot(object, reduction = "tsne", label = FALSE, split.by = "orig.ident",group.by = "integrated_snn_res.0.3", ncol = 3)
p
SaveFigure(p, "tSNE_SplitSample_res03", width = 15, height = 15)


### 
DefaultAssay(object) <- "RNA"
all.genes <- as.data.frame(row.names(object))

p1 <- FeaturePlot(object,features = "IL13RA2", reduction = "tsne")
p1
SaveFigure(p1, "tSNE_IL13RA2", width = 10, height = 8)

#### DE
library(MAST)
library(purrr)
library(tibble)
object2  <- subset(object, subset = group == "g3", invert = TRUE)

get_de <- function(cluster){
  FindMarkers(object = object2, ident.1 = "g1", min.pct = 0.5, logfc.threshold = 0,
              group.by = "group", test.use = "MAST", only.pos = FALSE, subset.ident = cluster) %>%
    rownames_to_column(var = "gene") %>% cbind(cluster_id = cluster, .)
}

get_wilcox <- function(cluster){
  FindMarkers(object = object2, ident.1 = "g1", min.pct = 0.5, logfc.threshold = 0,
              group.by = "group", test.use = "wilcox", only.pos = FALSE, subset.ident = cluster) %>%
    rownames_to_column(var = "gene") %>% cbind(cluster_id = cluster, .)
}

# Iterate function across desired clusters
Idents(object2) = object2$integrated_snn_res.0.8
mast.de.08 <- map_dfr(c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18), get_de)
write.csv(mast.de,"DE_MAST_G1_vs_G2_res08_full.csv", row.names = FALSE)
wilcox.de.08 <- map_dfr(c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18), get_wilcox)
write.csv(wilcox.de.08,"DE_Wilcox_G1_vs_G2_res08_full.csv", row.names = FALSE)
##
Idents(object2) = object2$integrated_snn_res.0.3
mast.de.03 <- map_dfr(c(0,1,2,3,4,5,6,7,8,9,10), get_de)
write.csv(mast.de.03,"DE_MAST_G1_vs_G2_res03_full.csv", row.names = FALSE)
wilcox.de.03 <- map_dfr(c(0,1,2,3,4,5,6,7,8,9,10), get_wilcox)
write.csv(wilcox.de.03,"DE_Wilcox_G1_vs_G2_res03_full.csv", row.names = FALSE)

### Plot DE
library(EnhancedVolcano)
mast.de.03b <- mast.de.03 %>% dplyr::filter(pct.1 > 0.5 & pct.2 > 0.5)

clusters <- unique(mast.de.03b$cluster_id) 

mast.list.03b <- lapply(clusters, FUN = function(x) {
  x <- dplyr::filter(mast.de.03b, cluster_id == x)})
names(mast.list.03b) <- clusters

for(i in 1:11) {
  EnhancedVolcano(mast.list.03b[[i]], lab = mast.list.03b[[i]][["gene"]], x = 'avg_log2FC', y = 'p_val_adj',
                  title = 'Group1 vs Group2', subtitle = paste0('Cluster_',names(mast.list.03b[i])),
                  legendLabels=c('Not sig.','Log (base 2) FC','q-value', 'q-value & Log (base 2) FC'),
                  drawConnectors = TRUE, pCutoff = 0.01, FCcutoff = 0.5,
                  pointSize = 3.0, axisLabSize = 12, legendLabSize = 10, legendIconSize = 3.0, labSize = 3.0)
  ggsave(paste0('Volcano_Cluster_',names(mast.list.03b[i]),'_res03.png'),dpi=300, scale =1.5)}

## Wilcox
wilcox.de.03b <- wilcox.de.03 %>% dplyr::filter(pct.1 > 0.5 & pct.2 > 0.5)
clusters <- unique(wilcox.de.03b$cluster_id) 
  
  wilcox.list.03b <- lapply(clusters, FUN = function(x) {
    x <- dplyr::filter(wilcox.de.03b, cluster_id == x)})
names(wilcox.list.03b) <- clusters

for(i in 1:11) {
  EnhancedVolcano(wilcox.list.03b[[i]], lab = wilcox.list.03b[[i]][["gene"]], x = 'avg_log2FC', y = 'p_val_adj',
                  title = 'Group1 vs Group2', subtitle = paste0('Cluster_',names(wilcox.list.03b[i])),
                  legendLabels=c('Not sig.','Log (base 2) FC','q-value', 'q-value & Log (base 2) FC'),
                  drawConnectors = TRUE, pCutoff = 0.01, FCcutoff = 0.5,
                  pointSize = 3.0, axisLabSize = 12, legendLabSize = 10, legendIconSize = 3.0, labSize = 3.0)
  ggsave(paste0('Volcano_Cluster_',names(wilcox.list.03b[i]),'_res03.png'),dpi=300, scale =1.5)}

##
wilcox.de.08b <- wilcox.de.08 %>% dplyr::filter(pct.1 > 0.5 & pct.2 > 0.5)
clusters <- unique(wilcox.de.08b$cluster_id) 

wilcox.list.08b <- lapply(clusters, FUN = function(x) {
  x <- dplyr::filter(wilcox.de.08b, cluster_id == x)})
names(wilcox.list.08b) <- clusters

for(i in 1:19) {
  EnhancedVolcano(wilcox.list.08b[[i]], lab = wilcox.list.08b[[i]][["gene"]], x = 'avg_log2FC', y = 'p_val_adj',
                  title = 'Group1 vs Group2', subtitle = paste0('Cluster_',names(wilcox.list.08b[i])),
                  legendLabels=c('Not sig.','Log (base 2) FC','q-value', 'q-value & Log (base 2) FC'),
                  drawConnectors = TRUE, pCutoff = 0.01, FCcutoff = 0.5,
                  pointSize = 3.0, axisLabSize = 12, legendLabSize = 10, legendIconSize = 3.0, labSize = 3.0)
  ggsave(paste0('Volcano_Cluster_',names(wilcox.list.08b[i]),'_res08.png'),dpi=300, scale =1.5)}

#ylim = c(0,max(-log10(mast.list.03b[[1]][["p_val_adj"]])+1, na.rm = TRUE) )


SaveObject(object, "SeuratProcessed")