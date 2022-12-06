library(sctransform)
library(glmGamPoi)
library(dplyr)

data_path <- "/Users/sbombin/Desktop/analysis_mm39.comb/Seurat_Analysis/IEC/"
fig_path <- "/Users/sbombin/Desktop/analysis_mm39.comb/Seurat_Analysis/IEC/"
setwd("/Users/sbombin/Desktop/analysis_mm39.comb/Seurat_Analysis/IEC/")

iec <- ReadObject("seuratobj.IEC.scale_all")
DefaultAssay(iec) <- "RNA"

isc <- subset(iec, idents = 8)
ta <- subset(iec, idents = 3)
tuft <- subset(iec, idents = 12)
goblet <- subset(iec, idents = 11)
eec <- subset(iec, idents = 13)
paneth <- subset(iec, idents = 15)

#nhc4 <- ScaleData(object = nhc4, verbose = F)
#nhc4 <- RunPCA(object = nhc4, verbose = F)
#ElbowPlot(object = nhc4, ndims = 40)
#nhc4 <- RunUMAP(object = nhc4, reduction = "pca", dims = 1:10, return.model = T)
#nhc4 <- FindNeighbors(object = nhc4, dims = 1:10)
#nhc4 <- FindClusters(object = nhc4, resolution = 0.3)
#p11 <- DimPlot(nhc4, reduction = "umap", label = TRUE)
#SaveFigure(p11, "Umap_Clusters_NHC4", width = 10, height = 8)

## ISC
isc.sct <- SCTransform(isc, verbose = FALSE, vst.flavor = "v2")
isc.sct<- RunPCA(isc.sct, verbose = FALSE)
isc.sct <- RunUMAP(isc.sct, dims = 1:30, verbose = FALSE)
isc.sct <- FindNeighbors(isc.sct, dims = 1:30, verbose = FALSE)
isc.sct <- FindClusters(isc.sct, verbose = FALSE, resolution = 0.6)
p1.isc <- DimPlot(isc.sct, reduction = "umap", label = TRUE)
p2.isc <- DimPlot(isc.sct, reduction = "umap", label = FALSE,group.by = "group")
p1.isc + p2.isc
SaveFigure((p1.isc + p2.isc), "UMAP_ISC_sct", width = 10, height = 10)
#ElbowPlot(object = isc.sct, ndims = 40)
## TA
ta.sct <- SCTransform(ta, verbose = FALSE, vst.flavor = "v2")
ta.sct<- RunPCA(ta.sct, verbose = FALSE)
ta.sct <- RunUMAP(ta.sct, dims = 1:30, verbose = FALSE)
ta.sct <- FindNeighbors(ta.sct, dims = 1:30, verbose = FALSE)
ta.sct <- FindClusters(ta.sct, verbose = FALSE, resolution = 0.4)
p1.ta <- DimPlot(ta.sct, reduction = "umap", label = TRUE)
#ElbowPlot(object = ta.sct, ndims = 40)
p2.ta <- DimPlot(ta.sct, reduction = "umap", label = FALSE,group.by = "group")
p1.ta + p2.ta
SaveFigure((p1.ta + p2.ta), "UMAP_TA_sct", width = 10, height = 10)
## Tuft
tuft.sct <- SCTransform(tuft, verbose = FALSE, vst.flavor = "v2")
tuft.sct<- RunPCA(tuft.sct, verbose = FALSE)
tuft.sct <- RunUMAP(tuft.sct, dims = 1:30, verbose = FALSE)
tuft.sct <- FindNeighbors(tuft.sct, dims = 1:30, verbose = FALSE)
tuft.sct <- FindClusters(tuft.sct, verbose = FALSE, resolution = 1.0)
p1.tuft <- DimPlot(tuft.sct, reduction = "umap", label = TRUE)
p2.tuft <- DimPlot(tuft.sct, reduction = "umap", label = FALSE,group.by = "group")
p1.tuft + p2.tuft
SaveFigure((p1.tuft + p2.tuft), "UMAP_Tuft_sct", width = 10, height = 10)
## Goblet
goblet.sct <- SCTransform(goblet, verbose = FALSE, vst.flavor = "v2")
goblet.sct<- RunPCA(goblet.sct, verbose = FALSE)
goblet.sct <- RunUMAP(goblet.sct, dims = 1:30, verbose = FALSE)
goblet.sct <- FindNeighbors(goblet.sct, dims = 1:30, verbose = FALSE)
goblet.sct <- FindClusters(goblet.sct, verbose = FALSE, resolution = 0.6)
p1.goblet <- DimPlot(goblet.sct, reduction = "umap", label = TRUE)
p2.goblet <- DimPlot(goblet.sct, reduction = "umap", label = FALSE,group.by = "group")
p1.goblet + p2.goblet
SaveFigure((p1.goblet + p2.goblet), "UMAP_Goblet_sct", width = 10, height = 10)
#ElbowPlot(object = goblet.sct, ndims = 40)
## EEC
eec.sct <- SCTransform(eec, verbose = FALSE, vst.flavor = "v2")
eec.sct<- RunPCA(eec.sct, verbose = FALSE)
eec.sct <- RunUMAP(eec.sct, dims = 1:30, verbose = FALSE)
eec.sct <- FindNeighbors(eec.sct, dims = 1:30, verbose = FALSE)
eec.sct <- FindClusters(eec.sct, verbose = FALSE, resolution = 1.4)
p1.eec <- DimPlot(eec.sct, reduction = "umap", label = TRUE)
p2.eec <- DimPlot(eec.sct, reduction = "umap", label = FALSE,group.by = "group")
p1.eec + p2.eec
SaveFigure((p1.eec + p2.eec), "UMAP_EEC_sct", width = 10, height = 10)
#ElbowPlot(object = eec.sct, ndims = 40)
## Paneth
paneth.sct <- SCTransform(paneth, verbose = FALSE, vst.flavor = "v2")
paneth.sct<- RunPCA(paneth.sct, verbose = FALSE)
paneth.sct <- RunUMAP(paneth.sct, dims = 1:15, verbose = FALSE)
paneth.sct <- FindNeighbors(paneth.sct, dims = 1:15, verbose = FALSE)
paneth.sct <- FindClusters(paneth.sct, verbose = FALSE, resolution = 0.6)
p1.paneth <- DimPlot(paneth.sct, reduction = "umap", label = TRUE)
p2.paneth <- DimPlot(paneth.sct, reduction = "umap", label = FALSE,group.by = "group")
p1.paneth + p2.paneth
SaveFigure((p1.paneth + p2.paneth), "UMAP_Paneth_sct", width = 10, height = 10)
#ElbowPlot(object = paneth.sct, ndims = 40)

table(iec@active.ident)

SaveObject(isc.sct, "seuratobj_ISC")
SaveObject(ta.sct, "seuratobj_TA")
SaveObject(tuft.sct, "seuratobj_Tuft")
SaveObject(goblet.sct, "seuratobj_Goblet")
SaveObject(eec.sct, "seuratobj_EEC")
SaveObject(paneth.sct, "seuratobj_paneth")

isc.sct <- PrepSCTFindMarkers(isc.sct)
ta.sct <- PrepSCTFindMarkers(ta.sct)
tuft.sct <- PrepSCTFindMarkers(tuft.sct)
goblet.sct <- PrepSCTFindMarkers(goblet.sct)
eec.sct <- PrepSCTFindMarkers(eec.sct)
paneth.sct <- PrepSCTFindMarkers(paneth.sct)

markers.isc <- FindAllMarkers(isc.sct, only.pos = FALSE, test.use = "MAST", assay = "SCT")
markers.ta <- FindAllMarkers(ta.sct, only.pos = FALSE, test.use = "MAST", assay = "SCT")
markers.tuft <- FindAllMarkers(tuft.sct, only.pos = FALSE, test.use = "MAST", assay = "SCT")
markers.goblet <- FindAllMarkers(goblet.sct, only.pos = FALSE, test.use = "MAST", assay = "SCT")
markers.eec <- FindAllMarkers(eec.sct, only.pos = FALSE, test.use = "MAST", assay = "SCT")
markers.paneth <- FindAllMarkers(paneth.sct, only.pos = FALSE, test.use = "MAST", assay = "SCT")

write.csv(markers.isc,"Markers-ISC_MAST_SCT.csv", row.names = FALSE)
write.csv(markers.ta,"Markers-TA_MAST_SCT.csv", row.names = FALSE)
write.csv(markers.tuft,"Markers-Tuft_MAST_SCT.csv", row.names = FALSE)
write.csv(markers.goblet,"Markers-Goblet_MAST_SCT.csv", row.names = FALSE)
write.csv(markers.eec,"Markers-EEC_MAST_SCT.csv", row.names = FALSE)
write.csv(markers.paneth,"Markers-Paneth_MAST_SCT.csv", row.names = FALSE)

get_de <- function(cluster){
  FindMarkers(object = iec, ident.1 = "Tet1iKO",
              subset.ident = cluster,
              group.by = "group", test.use = "MAST", 
              only.pos = FALSE) %>%
    rownames_to_column(var = "gene") %>%
    cbind(cluster_id = cluster, .)
}

# Iterate function across desired clusters
markers.iec <- map_dfr(c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15), get_de)

markers.iec <- markers.iec %>% mutate(celltype = cluster_id)
markers.iec$celltype <- gsub("^0$", "Enterocyte Progenitor (late) 1", markers.iec$celltype)
markers.iec$celltype <- gsub("^1$", "Mature Enterocyte (proximal) 1", markers.iec$celltype)
markers.iec$celltype <- gsub("^2$", "Mature Enterocyte (proximal) 2", markers.iec$celltype)
markers.iec$celltype <- gsub("^3$", "TA", markers.iec$celltype)
markers.iec$celltype <- gsub("^4$", "Mature Enterocyte (distal)", markers.iec$celltype)
markers.iec$celltype <- gsub("^5$", "Immature Enterocyte (proximal)", markers.iec$celltype)
markers.iec$celltype <- gsub("^6$", "Mature Enterocyte (proximal) 3", markers.iec$celltype)
markers.iec$celltype <- gsub("^7$", "Enterocyte Progenitor (early) 1", markers.iec$celltype)
markers.iec$celltype <- gsub("^8$", "ISC", markers.iec$celltype)
markers.iec$celltype <- gsub("^9$", "Enterocyte Progenitor (early) 2", markers.iec$celltype)
markers.iec$celltype <- gsub("^10$", "Enterocyte Progenitor (late) 2", markers.iec$celltype)
markers.iec$celltype <- gsub("^11$", "Goblet", markers.iec$celltype)
markers.iec$celltype <- gsub("^12$", "Tuft", markers.iec$celltype)
markers.iec$celltype <- gsub("^13$", "Enteroendocrine", markers.iec$celltype)
markers.iec$celltype <- gsub("^14$", "Mature Enterocyte (proximal) 4", markers.iec$celltype)
markers.iec$celltype <- gsub("^15$", "Paneth", markers.iec$celltype)
write.csv(markers.iec,"DE-MAST_Tet1iKO_vs_Control_Full.csv", row.names = FALSE)

de.0 <- markers.iec[markers.iec$cluster_id == '0',]
de.1 <- markers.iec[markers.iec$cluster_id == '1',]
de.2 <- markers.iec[markers.iec$cluster_id == '2',]
de.3 <- markers.iec[markers.iec$cluster_id == '3',]
de.4 <- markers.iec[markers.iec$cluster_id == '4',]
de.5 <- markers.iec[markers.iec$cluster_id == '5',]
de.6 <- markers.iec[markers.iec$cluster_id == '6',]
de.7 <- markers.iec[markers.iec$cluster_id == '7',]
de.8 <- markers.iec[markers.iec$cluster_id == '8',]
de.9 <- markers.iec[markers.iec$cluster_id == '9',]
de.10 <- markers.iec[markers.iec$cluster_id == '10',]
de.11 <- markers.iec[markers.iec$cluster_id == '11',]
de.12 <- markers.iec[markers.iec$cluster_id == '12',]
de.13 <- markers.iec[markers.iec$cluster_id == '13',]
de.14 <- markers.iec[markers.iec$cluster_id == '14',]
de.15 <- markers.iec[markers.iec$cluster_id == '15',]


p<-EnhancedVolcano(de.0, lab = de.0$gene, x = 'avg_log2FC', y = 'p_val_adj', title = 'iKO vs Control',  subtitle = 'Enterocyte Progenitor Cells (late) 1',
                   legendLabels=c('Not sig.','Log (base 2) FC','q-value', 'q-value & Log (base 2) FC'), drawConnectors = TRUE,
                   pCutoff = 0.05, FCcutoff = 0.05, pointSize = 3.0, axisLabSize = 12, legendLabSize = 10, legendIconSize = 3.0, labSize = 3.0)
ggsave(paste0('EPCL1_VolcanoPlot.png'),plot=p,dpi=300, scale = 2)

p<-EnhancedVolcano(de.1, lab = de.1$gene, x = 'avg_log2FC', y = 'p_val_adj', title = 'iKO vs Control',  subtitle = 'Mature Enterocyte Cells (proximal) 1',
                   legendLabels=c('Not sig.','Log (base 2) FC','q-value', 'q-value & Log (base 2) FC'), drawConnectors = TRUE,
                   pCutoff = 0.05, FCcutoff = 0.05, pointSize = 3.0, axisLabSize = 12, legendLabSize = 10, legendIconSize = 3.0, labSize = 3.0)
ggsave(paste0('MECP1_VolcanoPlot.png'),plot=p,dpi=300, scale = 2)

p<-EnhancedVolcano(de.2, lab = de.2$gene, x = 'avg_log2FC', y = 'p_val_adj', title = 'iKO vs Control',  subtitle = 'Mature Enterocyte Cells (proximal) 2',
                   legendLabels=c('Not sig.','Log (base 2) FC','q-value', 'q-value & Log (base 2) FC'), drawConnectors = TRUE,
                   pCutoff = 0.05, FCcutoff = 0.05, pointSize = 3.0, axisLabSize = 12, legendLabSize = 10, legendIconSize = 3.0, labSize = 3.0)
ggsave(paste0('MECP2_VolcanoPlot.png'),plot=p,dpi=300, scale = 2)

p<-EnhancedVolcano(de.3, lab = de.3$gene, x = 'avg_log2FC', y = 'p_val_adj', title = 'iKO vs Control',  subtitle = 'TA Cells',
                   legendLabels=c('Not sig.','Log (base 2) FC','q-value', 'q-value & Log (base 2) FC'), drawConnectors = TRUE,
                   pCutoff = 0.05, FCcutoff = 0.05, pointSize = 3.0, axisLabSize = 12, legendLabSize = 10, legendIconSize = 3.0, labSize = 3.0)
ggsave(paste0('TA_VolcanoPlot.png'),plot=p,dpi=300, scale = 2)

p<-EnhancedVolcano(de.4, lab = de.4$gene, x = 'avg_log2FC', y = 'p_val_adj', title = 'iKO vs Control',  subtitle = 'Mature Enterocyte Cells (distal)',
                   legendLabels=c('Not sig.','Log (base 2) FC','q-value', 'q-value & Log (base 2) FC'), drawConnectors = TRUE,
                   pCutoff = 0.05, FCcutoff = 0.05, pointSize = 3.0, axisLabSize = 12, legendLabSize = 10, legendIconSize = 3.0, labSize = 3.0)
ggsave(paste0('MECD_VolcanoPlot.png'),plot=p,dpi=300, scale = 2)

p<-EnhancedVolcano(de.5, lab = de.5$gene, x = 'avg_log2FC', y = 'p_val_adj', title = 'iKO vs Control',  subtitle = 'Immature Enterocyte Cells (proximal)',
                   legendLabels=c('Not sig.','Log (base 2) FC','q-value', 'q-value & Log (base 2) FC'), drawConnectors = TRUE,
                   pCutoff = 0.05, FCcutoff = 0.05, pointSize = 3.0, axisLabSize = 12, legendLabSize = 10, legendIconSize = 3.0, labSize = 3.0)
ggsave(paste0('IECP_VolcanoPlot.png'),plot=p,dpi=300, scale = 2)

p<-EnhancedVolcano(de.6, lab = de.6$gene, x = 'avg_log2FC', y = 'p_val_adj', title = 'iKO vs Control',  subtitle = 'Mature Enterocyte Cells (proximal) 3',
                   legendLabels=c('Not sig.','Log (base 2) FC','q-value', 'q-value & Log (base 2) FC'), drawConnectors = TRUE,
                   pCutoff = 0.05, FCcutoff = 0.05, pointSize = 3.0, axisLabSize = 12, legendLabSize = 10, legendIconSize = 3.0, labSize = 3.0)
ggsave(paste0('MECP3_VolcanoPlot.png'),plot=p,dpi=300, scale = 2)

p<-EnhancedVolcano(de.7, lab = de.7$gene, x = 'avg_log2FC', y = 'p_val_adj', title = 'iKO vs Control',  subtitle = 'Enterocyte Progenitor Cells (early) 1',
                   legendLabels=c('Not sig.','Log (base 2) FC','q-value', 'q-value & Log (base 2) FC'), drawConnectors = TRUE,
                   pCutoff = 0.05, FCcutoff = 0.05, pointSize = 3.0, axisLabSize = 12, legendLabSize = 10, legendIconSize = 3.0, labSize = 3.0)
ggsave(paste0('EPCE1_VolcanoPlot.png'),plot=p,dpi=300, scale = 2)

p<-EnhancedVolcano(de.8, lab = de.8$gene, x = 'avg_log2FC', y = 'p_val_adj', title = 'iKO vs Control',  subtitle = 'ISC',
                   legendLabels=c('Not sig.','Log (base 2) FC','q-value', 'q-value & Log (base 2) FC'), drawConnectors = TRUE,
                   pCutoff = 0.05, FCcutoff = 0.05, pointSize = 3.0, axisLabSize = 12, legendLabSize = 10, legendIconSize = 3.0, labSize = 3.0)
ggsave(paste0('VolcanoPlot_ISC.png'),plot=p,dpi=300, scale = 2)

p<-EnhancedVolcano(de.9, lab = de.9$gene, x = 'avg_log2FC', y = 'p_val_adj', title = 'iKO vs Control',  subtitle = 'Enterocyte Progenitor Cells (early) 2',
                   legendLabels=c('Not sig.','Log (base 2) FC','q-value', 'q-value & Log (base 2) FC'), drawConnectors = TRUE,
                   pCutoff = 0.05, FCcutoff = 0.05, pointSize = 3.0, axisLabSize = 12, legendLabSize = 10, legendIconSize = 3.0, labSize = 3.0)
ggsave(paste0('VolcanoPlot_EPCE2.png'),plot=p,dpi=300, scale = 2)

p<-EnhancedVolcano(de.10, lab = de.10$gene, x = 'avg_log2FC', y = 'p_val_adj', title = 'iKO vs Control',  subtitle = 'Enterocyte Progenitor Cells (late) 2',
                   legendLabels=c('Not sig.','Log (base 2) FC','q-value', 'q-value & Log (base 2) FC'), drawConnectors = TRUE,
                   pCutoff = 0.05, FCcutoff = 0.05, pointSize = 3.0, axisLabSize = 12, legendLabSize = 10, legendIconSize = 3.0, labSize = 3.0)
ggsave(paste0('VolcanoPlot_EPCL2.png'),plot=p,dpi=300, scale = 2)

p<-EnhancedVolcano(de.11, lab = de.11$gene, x = 'avg_log2FC', y = 'p_val_adj', title = 'iKO vs Control',  subtitle = 'Goblet Cells',
                   legendLabels=c('Not sig.','Log (base 2) FC','q-value', 'q-value & Log (base 2) FC'), drawConnectors = TRUE,
                   pCutoff = 0.05, FCcutoff = 0.05, pointSize = 3.0, axisLabSize = 12, legendLabSize = 10, legendIconSize = 3.0, labSize = 3.0)
ggsave(paste0('VolcanoPlot_Goblet.png'),plot=p,dpi=300, scale = 2)

p<-EnhancedVolcano(de.12, lab = de.12$gene, x = 'avg_log2FC', y = 'p_val_adj', title = 'iKO vs Control',  subtitle = 'Tuft Cells',
                   legendLabels=c('Not sig.','Log (base 2) FC','q-value', 'q-value & Log (base 2) FC'), drawConnectors = TRUE,
                   pCutoff = 0.05, FCcutoff = 0.05, pointSize = 3.0, axisLabSize = 12, legendLabSize = 10, legendIconSize = 3.0, labSize = 3.0)
ggsave(paste0('VolcanoPlot_Tuft.png'),plot=p,dpi=300, scale = 2)

p<-EnhancedVolcano(de.13, lab = de.13$gene, x = 'avg_log2FC', y = 'p_val_adj', title = 'iKO vs Control',  subtitle = 'Enteroendocrine Cells',
                   legendLabels=c('Not sig.','Log (base 2) FC','q-value', 'q-value & Log (base 2) FC'), drawConnectors = TRUE,
                   pCutoff = 0.05, FCcutoff = 0.05, pointSize = 3.0, axisLabSize = 12, legendLabSize = 10, legendIconSize = 3.0, labSize = 3.0)
ggsave(paste0('VolcanoPlot_Enteroendocrine.png'),plot=p,dpi=300, scale = 2)

p<-EnhancedVolcano(de.14, lab = de.14$gene, x = 'avg_log2FC', y = 'p_val_adj', title = 'iKO vs Control',  subtitle = 'Mature Enterocyte Cells (proximal) 4',
                   legendLabels=c('Not sig.','Log (base 2) FC','q-value', 'q-value & Log (base 2) FC'), drawConnectors = TRUE,
                   pCutoff = 0.05, FCcutoff = 0.05, pointSize = 3.0, axisLabSize = 12, legendLabSize = 10, legendIconSize = 3.0, labSize = 3.0)
ggsave(paste0('VolcanoPlot_MECP4.png'),plot=p,dpi=300, scale = 2)

p<-EnhancedVolcano(de.15, lab = de.15$gene, x = 'avg_log2FC', y = 'p_val_adj', title = 'iKO vs Control',  subtitle = 'Paneth Cells',
                   legendLabels=c('Not sig.','Log (base 2) FC','q-value', 'q-value & Log (base 2) FC'), drawConnectors = TRUE,
                   pCutoff = 0.05, FCcutoff = 0.05, pointSize = 3.0, axisLabSize = 12, legendLabSize = 10, legendIconSize = 3.0, labSize = 3.0)
ggsave(paste0('VolcanoPlot_Paneth.png'),plot=p,dpi=300, scale = 2)

## Summary Plot
summary.iec.up <- markers.iec %>% filter(avg_log2FC > 0 & p_val_adj < 0.05) %>% group_by(celltype) %>% summarise(Up = n())
summary.iec.down <- markers.iec %>% filter(avg_log2FC < 0 & p_val_adj < 0.05) %>% group_by(celltype) %>% summarise(Down = n())
summary.iec.full <-  as.data.frame(right_join(summary.iec.up,summary.iec.down,by='celltype'))
summary.iec.full[is.na(summary.iec.full)] <- 0 
m <-column_to_rownames(summary.iec.full, var = "celltype") %>% as.matrix()

library(reshape2)
dfm <- melt(summary.iec.full[,c('celltype','Up','Down')],id.vars = 1)
p.summary <- ggplot(dfm,aes(x = celltype,y = value)) + 
  geom_bar(aes(fill = variable),stat = "identity",position = "dodge") + theme(axis.text.x = element_text(size = 7, angle = 45, hjust = 1)) +
  xlab("") + ylab("Significantly expressed genes")

ggsave(paste0('Barplot_DE.png'),plot=p.summary,dpi=300, scale = 2)

## Heatmap
top10 <- markers.isc %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
p.heatmap <- DoHeatmap(isc.sct, features = top10$gene) + NoLegend()
SaveFigure(p.heatmap, "Heatmap_ISC_sct", width = 10, height = 10)

top10 <- markers.ta %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
p.heatmap <- DoHeatmap(ta.sct, features = top10$gene) + NoLegend()
SaveFigure(p.heatmap, "Heatmap_TA_sct", width = 10, height = 10)

top10 <- markers.tuft %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
p.heatmap <- DoHeatmap(tuft.sct, features = top10$gene) + NoLegend()
SaveFigure(p.heatmap, "Heatmap_Tuft_sct", width = 10, height = 10)

top10 <- markers.goblet %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
p.heatmap <- DoHeatmap(goblet.sct, features = top10$gene) + NoLegend()
SaveFigure(p.heatmap, "Heatmap_Goblet_sct", width = 10, height = 10)

top10 <- markers.eec %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
p.heatmap <- DoHeatmap(eec.sct, features = top10$gene) + NoLegend()
SaveFigure(p.heatmap, "Heatmap_EEC_sct", width = 10, height = 10)

top10 <- markers.paneth %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
p.heatmap <- DoHeatmap(paneth.sct, features = top10$gene) + NoLegend()
SaveFigure(p.heatmap, "Heatmap_Paneth_sct", width = 10, height = 10)


