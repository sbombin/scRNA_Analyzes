##install.packages('BiocManager')
##BiocManager::install('multtest')
library(multtest)
library(metap)

immune.combined <- ReadObject("SeuratObj_Sham-BDL/seuratobj_intgrtn_CLR2")
DefaultAssay(immune.combined) <- "RNA"

nk.markers <- FindConservedMarkers(immune.combined, ident.1 = 6, grouping.var = "group", verbose = FALSE)

p1<- FeaturePlot(immune.combined, features = c("Mki67", "Pcna", "Ccna2", "Ccl2", "Cxcl1", "Cxcl2", "Sparc", "Tgfb2",
                                          "Cd14", "Fn1", "Cdk1", "Itgb6", "Ctgf", "Cyr61"), min.cutoff = "q9") 

SaveFigure((p1), "DE_Intgrtn_Sham-BDL", width = 20, height = 20)


plots <- VlnPlot(immune.combined, features = c("Mki67", "Pcna", "Ccna2", "Ccl2", "Cxcl1", "Cxcl2", "Sparc", "Tgfb2",
                                               "Cd14", "Fn1", "Cdk1", "Itgb6", "Ctgf", "Cyr61"),  
                 group.by = "sample", pt.size = 0, combine = FALSE)

wrap_plots(plots = plots, ncol = 1)
SaveFigure((plots), "VlnPlot2_Intgrtn_Sham-BDL", width = 20, height = 100)

plots <- VlnPlot(immune.combined, features = c("Mki67", "Pcna", "Ccna2", "Ccl2"), split.by = NULL, group.by = "sample", pt.size = 0, combine = FALSE)
wrap_plots(plots = plots, ncol = 1)
SaveFigure(wrap_plots(plots = plots, ncol = 1), "VlnPlots4_Intgrtn_Sham-BDL", width = 10, height = 15)

plots <- VlnPlot(immune.combined, features = c("Mki67", "Pcna", "Ccna2", "Ccl2", "Cxcl1", "Cxcl2", "Sparc", "Tgfb2",
                                               "Cd14", "Fn1", "Cdk1", "Itgb6", "Ctgf", "Cyr61"), split.by = NULL, group.by = "sample", pt.size = 0)
#wrap_plots(plots = plots, ncol = 1)
#SaveFigure(wrap_plots(plots = plots, ncol = 1), "VlnPlots4_Intgrtn_Sham-BDL", width = 10, height = 15)
SaveFigure((plots), "VlnPlot3_Intgrtn_Sham-BDL", width = 15, height = 15)

DefaultAssay(immune.combined) <- "integrated"

p2 <- DoHeatmap(immune.combined, features = c("Mki67", "Pcna", "Ccna2", "Ccl2", "Cxcl1", "Cxcl2", "Sparc", "Tgfb2",
                                              "Cd14", "Fn1", "Cdk1", "Itgb6", "Ctgf", "Cyr61"), group.by = "sample")

SaveFigure((p2), "DoHeatmap_Intgrtn_Sham-BDL", width = 15, height = 15)
                
