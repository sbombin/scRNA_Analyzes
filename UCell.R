library(UCell)

set.seed(123)

pbmc <- AddModuleScore_UCell(pbmc, features = markers)

uscore <- AddModuleScore_UCell(pbmc, features = markers,storeRanks =TRUE)


signature.names <- paste0(names(markers), "_UCell")

p <- FeaturePlot(pbmc, reduction = "umap", features = signature.names, ncol = 3,
            order = T)

SaveFigure(p, "UCell_FeaturePlot", width = 15, height = 10)

p <- FeaturePlot(pbmc, features = signature.names, split.by = "group", max.cutoff = 3, cols = c("grey", "red"))
p
SaveFigure(p, "UCell_FeaturePlot_Group", width = 15, height = 20)


library(patchwork)
p <- VlnPlot(pbmc, features = signature.names, group.by = "seurat_clusters", split.by = "group", combine = FALSE)
p2 <- wrap_plots(plots = p, ncol = 2)
SaveFigure(p2, "UCell_VlnPlot_Group", width = 15, height = 10)


hep.up <- subpops %>% filter(Hep > 1 & symbol %in% genes)  
hep.up <- unique(hep.up$symbol)[1:600]
hep.up <- paste(hep.up, "+", sep="")
hep.dwn <- subpops %>% filter(Hep < -1 & symbol %in% genes)
hep.dwn <- unique(hep.dwn$symbol)[1:600]
hep.dwn <- paste(hep.dwn, "-", sep="")
hep.all <- c(hep.up,hep.dwn)

markers <- list(yap.top,notch.all,hep.all,pc.all,early.all,late.all)
names(markers) <- c('YAP','NOTCH','HEP','PC','TGFB_Early','TGFB_Late')


pbmc <- irGSEA.score(object = pbmc, assay = "RNA", slot = "data", msigdb = F, custom = TRUE, geneset = markers, method =  "UCell")
library(dplyr)
result.dge <- irGSEA.integrate(object = pbmc, 
                               group.by = "seurat_clusters",
                               metadata = NULL, col.name = NULL,
                               method = "UCell")

p1 <-irGSEA.heatmap.plot <- irGSEA.heatmap(object = result.dge, 
                                      method = "UCell",
                                      top = 50, 
                                      show.geneset = NULL)
SaveFigure(p1, "UCell_Heatmap", width = 15, height = 10)


