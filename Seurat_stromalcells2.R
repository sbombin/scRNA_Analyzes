library(Seurat)
library(dplyr)
library(MAST)
library(ggplot2)
library(Matrix)
library(DESeq2)
library(patchwork)
library(tibble)
library(multtest)
library(metap)
library(cowplot)

setwd("/Users/sbombin/Desktop/EICC/Projects/Lawrence_Boise/")
data_path <- "/Users/sbombin/Desktop/EICC/Projects/Lawrence_Boise/"
fig_path <- "/Users/sbombin/Desktop/EICC/Projects/Lawrence_Boise/"

#SaveObject(myeloma.integrated, "myeloma.integrated")
myeloma.integrated <- ReadObject("myeloma.integrated")
DefaultAssay(myeloma.integrated) <- "integrated"

p1 <- DimPlot(myeloma.integrated, reduction = "umap", label = TRUE)
SaveFigure(p1, "Umap_Clusters", width = 10, height = 8)

p2 <- DimPlot(myeloma.integrated, label = TRUE, group.by = "source", reduction = "umap", repel = TRUE)

p4 <- FeaturePlot(object = myeloma.integrated, features = "CXCL12")
p5 <- FeaturePlot(object = myeloma.integrated, features = "CD36")
p6 <- FeaturePlot(object = myeloma.integrated, features = "LEPR")
# Main genes
p7 <- FeaturePlot(object = myeloma.integrated, features = "IL6R")
p8 <- FeaturePlot(object = myeloma.integrated, features = "IL6ST")
p9 <- FeaturePlot(object = myeloma.integrated, features = "AK1")
p10 <- FeaturePlot(object = myeloma.integrated, features = "STAT3")
p11 <- FeaturePlot(object = myeloma.integrated, features = "SOCS3")

p7+p8+p9+p10+p11

DefaultAssay(myeloma.integrated) <- "RNA"

## Add module to separate non-hematopoietic based on "CXCL12" and "CD36"
non_hematopoietic  <- c('CXCL12','CD36')
non_hematopoietic.list  <- list(c('CXCL12','CD36'))
object <- AddModuleScore(object = myeloma.integrated, features = non_hematopoietic.list, name = "Non_Hematopoietic")
p12 <- FeaturePlot(object = object, features = "Non_Hematopoietic")

nhc <- subset(myeloma.integrated, subset = CXCL12 > 0 | CD36 > 0  )
dim(nhc)


FeaturePlot(object = nhc, features = "CD34")
p15 <- FeaturePlot(object = nhc, features = "CD45")

#cells.located <- CellSelector(plot = p1)

DefaultAssay(nhc) <- "integrated"

nhc <- ScaleData(object = nhc, verbose = F)
nhc <- RunPCA(object = nhc, verbose = F)
ElbowPlot(object = nhc, ndims = 40)
nhc <- RunUMAP(object = nhc, reduction = "pca", dims = 1:10, return.model = T)
nhc <- FindNeighbors(object = nhc, dims = 1:10)
nhc <- FindClusters(object = nhc, resolution = 0.3)

p1 <- DimPlot(nhc, reduction = "umap", label = TRUE)
SaveFigure(p1, "Umap_Clusters_NHC", width = 10, height = 8)

p12 <- FeaturePlot(object = nhc, features = "CXCL12")
p13 <- FeaturePlot(object = nhc, features = "CD36")
p12+p13

p14 <- DimPlot(nhc, label = TRUE, group.by = "source", reduction = "umap", repel = TRUE)
p14 +p11

table(myeloma.integrated@active.ident)
table(nhc@active.ident)

#hc <- c('Sca-1', 'CD27', 'CD34', 'CD38', 'CD43', 'CD48', 'CD117', 'CD150')
#p18 <- FeaturePlot(object = myeloma.integrated, features = 'CD133')
#p19 <- FeaturePlot(object = nhc, features = "CD133")
#p17 <- FeaturePlot(object = nhc2, features = 'CD133')
#p17+p18+p19
#nhc2 <- subset(myeloma.integrated, subset = CD34 < 0.5 & CD27 <0.5 & CD38 <0.5 & CD48<0.5)
#nhc3 <- subset(myeloma.integrated, subset =  CD38 <0.5)
#dim(nhc2)
#dim(myeloma.integrated)
#dim(nhc3)
#genes <- c('LEPR','CXCL12','CDH5', 'THBD' ,'CD36', 'SELP', 'RUNX2', 'SP7')

nhc4 <- subset(myeloma.integrated, subset =  LEPR> 1|CXCL12>1|CDH5>1| THBD>1 | CD36>1 | SELP>1 | RUNX2>1 | SP7>1)
dim(nhc4)
#p17 <- FeaturePlot(object = nhc4, features = 'CD38')
DefaultAssay(nhc4) <- "integrated"
nhc4 <- ScaleData(object = nhc4, verbose = F)
nhc4 <- RunPCA(object = nhc4, verbose = F)
ElbowPlot(object = nhc4, ndims = 40)
nhc4 <- RunUMAP(object = nhc4, reduction = "pca", dims = 1:10, return.model = T)
nhc4 <- FindNeighbors(object = nhc4, dims = 1:10)
nhc4 <- FindClusters(object = nhc4, resolution = 0.3)

p11 <- DimPlot(nhc4, reduction = "umap", label = TRUE)
SaveFigure(p11, "Umap_Clusters_NHC4", width = 10, height = 8)


SaveObject(nhc, "nhc.integrated")
SaveObject(nhc4, "nhc4.integrated")

nhc <- ReadObject("nhc.integrated")
nhc4 <- ReadObject("nhc4.integrated")


DefaultAssay(nhc) <- "RNA"
DefaultAssay(nhc4) <- "RNA"

get_conserved <- function(cluster){
  FindConservedMarkers(nhc4,
                       ident.1 = cluster, min.pct = 0.1, logfc.threshold = 0.8,
                       grouping.var = "source", test.use = "MAST", 
                       only.pos = FALSE) %>%
    rownames_to_column(var = "gene") %>%
    cbind(cluster_id = cluster, .)
}

# Iterate function across desired clusters
conserved_markers_nhc <- map_dfr(c(0,1,2,3,4,5,6,7), get_conserved)
write.csv(conserved_markers_nhc,"DE-MAST_conserved_markers_nhc.csv", row.names = FALSE)

conserved_markers_nhc4 <- map_dfr(c(0,1,2,3,4,5,6,7), get_conserved)
write.csv(conserved_markers_nhc4,"DE-MAST_conserved_markers_nhc4.csv", row.names = FALSE)

ec4 <- subset(nhc4, idents = 4)
Idents(ec4) <- "source"
de.ec4 <- FindMarkers(ec4, ident.1 = "myeloma", ident.2 = "control", verbose = FALSE,only.pos = FALSE,test.use = "MAST")

get_de <- function(cluster){
  FindMarkers(object = nhc4, ident.1 = "myeloma",
              subset.ident = cluster,
              group.by = "source", test.use = "MAST", 
              only.pos = FALSE) %>%
    rownames_to_column(var = "gene") %>%
    cbind(cluster_id = cluster, .)
}

de.nhc4 <- map_dfr(c(0,1,2,3,4,5,6,7), get_de)
de.nhc <- map_dfr(c(1,2,3,4,5,6,7), get_de)

write.csv(de.nhc4,"DE-MAST_NHC.csv", row.names = FALSE)
#write.csv(de.nhc,"DE-MAST_nhc.csv", row.names = FALSE)

markers.nhc <- FindAllMarkers(nhc, only.pos = FALSE, test.use = "MAST")
markers.nhc %>% group_by(cluster) %>% slice_max(n = 2, order_by = avg_log2FC)

markers.nhc4 <- FindAllMarkers(nhc4, only.pos = FALSE, test.use = "MAST")
markers.nhc4 %>% group_by(cluster) %>% slice_max(n = 2, order_by = avg_log2FC)

write.csv(markers.nhc4,"DE-MAST_All_Markers.csv", row.names = FALSE)

p2.nhc4 <-FeaturePlot(nhc4, features = c("IL6R", "IL6ST", "JAK1", "STAT3", "SOCS3"), split.by = "source", max.cutoff = 3,
            cols = c("grey", "red"))

SaveFigure(p2.nhc4, "UMAP_Genes", width = 8, height = 15)

nhc4 <- RenameIdents(nhc4, `0` = "MSC1", `1` = "MSC2", `2` = "MSC4",`3` = "MSC3", `4` = "EC", `5` = "MSC5", `6` = "SEC", `7` = "OLC")
p3.nhc4 <- DimPlot(nhc4, reduction = "umap", label = TRUE)
SaveFigure(p3.nhc4, "UMAP_cell_rename", width = 10, height = 8)

markers <- c("IL6R", "IL6ST", "JAK1", "STAT3", "SOCS3")

p4.nhc4 <- DotPlot(nhc4, features = markers, cols = c("blue", "red"), dot.scale = 8, split.by = "source") 
SaveFigure(p4.nhc4, "DotPlot_markers", width = 8, height = 8)

nhc4$celltype <- Idents(nhc4)

plots <- VlnPlot(nhc4, features = markers, split.by = "source", group.by = "celltype", pt.size = 0, combine = FALSE)
plots2 <- wrap_plots(plots = plots, ncol = 2)
SaveFigure(plots2, "VlnPlot_Markers", width = 10, height = 10)


de.nhc4 <- de.nhc4 %>% mutate(celltype = cluster_id)

de.nhc4$celltype <- gsub("^0$", "MSC1", de.nhc4$celltype)
de.nhc4$celltype <- gsub("^1$", "MSC2", de.nhc4$celltype)
de.nhc4$celltype <- gsub("^2$", "MSC4", de.nhc4$celltype)
de.nhc4$celltype <- gsub("^3$", "MSC3", de.nhc4$celltype)
de.nhc4$celltype <- gsub("^4$", "Endothelial cells", de.nhc4$celltype)
de.nhc4$celltype <- gsub("^5$", "MSC5", de.nhc4$celltype)
de.nhc4$celltype <- gsub("^6$", "SEC", de.nhc4$celltype)
de.nhc4$celltype <- gsub("^7$", "Osteolineage cells", de.nhc4$celltype)

### Pairwise Comparison

msc1.vs.msc2 <- FindMarkers(nhc4, ident.1 = "MSC1", ident.2 = "MSC2", test.use = "MAST", only.pos = FALSE)
msc1.vs.msc3 <- FindMarkers(nhc4, ident.1 = "MSC1", ident.2 = "MSC3", test.use = "MAST", only.pos = FALSE)
msc1.vs.msc4 <- FindMarkers(nhc4, ident.1 = "MSC1", ident.2 = "MSC4", test.use = "MAST", only.pos = FALSE)
msc1.vs.msc5 <- FindMarkers(nhc4, ident.1 = "MSC1", ident.2 = "MSC5", test.use = "MAST", only.pos = FALSE)

msc2.vs.msc3 <- FindMarkers(nhc4, ident.1 = "MSC2", ident.2 = "MSC3", test.use = "MAST", only.pos = FALSE)
msc2.vs.msc4 <- FindMarkers(nhc4, ident.1 = "MSC2", ident.2 = "MSC4", test.use = "MAST", only.pos = FALSE)
msc2.vs.msc5 <- FindMarkers(nhc4, ident.1 = "MSC2", ident.2 = "MSC5", test.use = "MAST", only.pos = FALSE)

msc3.vs.msc4 <- FindMarkers(nhc4, ident.1 = "MSC3", ident.2 = "MSC4", test.use = "MAST", only.pos = FALSE)
msc3.vs.msc5 <- FindMarkers(nhc4, ident.1 = "MSC3", ident.2 = "MSC5", test.use = "MAST", only.pos = FALSE)

msc4.vs.msc5 <- FindMarkers(nhc4, ident.1 = "MSC4", ident.2 = "MSC5", test.use = "MAST", only.pos = FALSE)


write.csv(msc1.vs.msc2,"DE_msc1.vs.msc2.csv", row.names = TRUE)
write.csv(msc1.vs.msc3,"DE_msc1.vs.msc3.csv", row.names = TRUE)
write.csv(msc1.vs.msc4,"DE_msc1.vs.msc4.csv", row.names = TRUE)
write.csv(msc1.vs.msc5,"DE_msc1.vs.msc5.csv", row.names = TRUE)

write.csv(msc2.vs.msc3,"DE_msc2.vs.msc3.csv", row.names = TRUE)
write.csv(msc2.vs.msc4,"DE_msc2.vs.msc4.csv", row.names = TRUE)
write.csv(msc2.vs.msc5,"DE_msc2.vs.msc5.csv", row.names = TRUE)

write.csv(msc3.vs.msc4,"DE_msc3.vs.msc4.csv", row.names = TRUE)
write.csv(msc3.vs.msc5,"DE_msc3.vs.msc5.csv", row.names = TRUE)
write.csv(msc4.vs.msc5,"DE_msc4.vs.msc5.csv", row.names = TRUE)


msc <- subset(nhc4, idents =  c("MSC1","MSC2","MSC3","MSC4","MSC5"))

markers.msc <- FindAllMarkers(msc, only.pos = FALSE, test.use = "MAST")
write.csv(markers.msc,"DE-MAST_All_Markers_MSC.csv", row.names = FALSE)

table(markers.nhc4$cluster)
table(markers.msc$cluster)

## Count based on q-value and LogFC threshold
length( which( msc1.vs.msc2$p_val_adj < 0.05 & msc1.vs.msc2$avg_log2FC > 1 ) )
length( which( msc1.vs.msc2$p_val_adj < 0.05 & msc1.vs.msc2$avg_log2FC < -1 ) )
length( which( msc1.vs.msc3$p_val_adj < 0.05 & msc1.vs.msc3$avg_log2FC > 1 ) )
length( which( msc1.vs.msc3$p_val_adj < 0.05 & msc1.vs.msc3$avg_log2FC < -1 ) )
length( which( msc1.vs.msc4$p_val_adj < 0.05 & msc1.vs.msc4$avg_log2FC > 1 ) )
length( which( msc1.vs.msc4$p_val_adj < 0.05 & msc1.vs.msc4$avg_log2FC < -1 ) )
length( which( msc1.vs.msc5$p_val_adj < 0.05 & msc1.vs.msc5$avg_log2FC > 1 ) )
length( which( msc1.vs.msc5$p_val_adj < 0.05 & msc1.vs.msc5$avg_log2FC < -1 ) )

length( which( msc2.vs.msc3$p_val_adj < 0.05 & msc2.vs.msc3$avg_log2FC > 1 ) )
length( which( msc2.vs.msc3$p_val_adj < 0.05 & msc2.vs.msc3$avg_log2FC < -1 ) )
length( which( msc2.vs.msc4$p_val_adj < 0.05 & msc2.vs.msc4$avg_log2FC > 1 ) )
length( which( msc2.vs.msc4$p_val_adj < 0.05 & msc2.vs.msc4$avg_log2FC < -1 ) )
length( which( msc2.vs.msc5$p_val_adj < 0.05 & msc2.vs.msc5$avg_log2FC > 1 ) )
length( which( msc2.vs.msc5$p_val_adj < 0.05 & msc2.vs.msc5$avg_log2FC < -1 ) )

length( which( msc3.vs.msc4$p_val_adj < 0.05 & msc3.vs.msc4$avg_log2FC > 1 ) )
length( which( msc3.vs.msc4$p_val_adj < 0.05 & msc3.vs.msc4$avg_log2FC < -1 ) )
length( which( msc3.vs.msc5$p_val_adj < 0.05 & msc3.vs.msc5$avg_log2FC > 1 ) )
length( which( msc3.vs.msc5$p_val_adj < 0.05 & msc3.vs.msc5$avg_log2FC < -1 ) )

length( which( msc4.vs.msc5$p_val_adj < 0.05 & msc4.vs.msc5$avg_log2FC > 1 ) )
length( which( msc4.vs.msc5$p_val_adj < 0.05 & msc4.vs.msc5$avg_log2FC < -1 ) )


########
# Part2#
########
nhc4 <- ReadObject("nhc4.integrated")
DefaultAssay(nhc4) <- "RNA"

## Up Genes Plots
p1 <-FeaturePlot(nhc4, features = c("CXCL8", "IL1B", "IL1A", "CSF2"), split.by = "source", cols = c("grey", "red"))
p1
SaveFigure(p1, "UMAP_UP1", width = 8, height = 15)
p2 <-FeaturePlot(nhc4, features = c("CSF3", "IL24", "CXCL3"), split.by = "source", cols = c("grey", "red"))
p2
SaveFigure(p2, "UMAP_UP2", width = 8, height = 12)

## Down Genes Plots
p1 <-FeaturePlot(nhc4, features = c("CXCL5", "CXCL12", "DKK3", "VEGFA"), split.by = "source", cols = c("grey", "red"))
p1
SaveFigure(p1, "UMAP_Down1", width = 8, height = 15)
p2 <-FeaturePlot(nhc4, features = c("CD70", "ITGA6", "TGFB1"), split.by = "source", cols = c("grey", "red"))
p2
SaveFigure(p2, "UMAP_Down2", width = 8, height = 12)

p2 <-FeaturePlot(nhc4, features = c("IL6", "IGF1"), split.by = "source", cols = c("grey", "red"))
p2
SaveFigure(p2, "UMAP_IL6-IGF1", width = 8, height = 10)


DefaultAssay(nhc4) <- "integrated"

## Subset control samples
control <- subset(x=nhc4, subset = source == "control")

p11 <- DimPlot(control, reduction = "umap", label = TRUE)
SaveFigure(p11, "UMAP_Clusters_Control", width = 10, height = 8)

control <- RenameIdents(control, `0` = "MSC1", `1` = "MSC2", `2` = "MSC4",`3` = "MSC3", `4` = "EC", `5` = "MSC5", `6` = "SEC", `7` = "OLC")
p12 <- DimPlot(control, reduction = "umap", label = TRUE)
SaveFigure(p12, "UMAP_Celltype_Control", width = 10, height = 8)
