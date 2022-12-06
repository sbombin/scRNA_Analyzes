
library(enrichR)

pbmc <- ReadObject("seuratobj_BEC.scale_all")

get_goe<- function(cluster){
  DEenrichRPlot(object = pbmc, ident.1 = cluster,
                min.pct = 0.1, logfc.threshold = 0.25, test.use = "MAST", only.pos = FALSE,
                enrich.database = "KEGG_2019_Mouse", num.pathway = 50, max.genes =1000,
                return.gene.list = TRUE) %>% as.data.frame() %>% cbind(cluster_id = cluster, .)
}

# Iterate function across desired clusters
#goe <- map_dfr(c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15), get_goe) 
goe <- map_dfr(c(0,1,2,3,4,5,6,7,8,9,10,11), get_goe) 
write.csv(goe,"/Users/sbombin/Desktop/analysis_mm39.comb/Seurat_Analysis/BEC//GO_Enrichment_all_clusters.csv", row.names = FALSE)

## Plots
pgoe0 <- DEenrichRPlot(object = pbmc, ident.1 = 0, logfc.threshold = 0.8, test.use = "MAST", enrich.database = "KEGG_2019_Mouse", num.pathway = 50, max.genes =1000,
                       return.gene.list = FALSE, min.pct = 0.1)
pgoe1 <- DEenrichRPlot(object = pbmc, ident.1 = 1, logfc.threshold = 0.8, test.use = "MAST", enrich.database = "KEGG_2019_Mouse", num.pathway = 50, max.genes =1000,
                       return.gene.list = FALSE, min.pct = 0.1)
pgoe2 <- DEenrichRPlot(object = pbmc, ident.1 = 2, logfc.threshold = 0.8, test.use = "MAST", enrich.database = "KEGG_2019_Mouse", num.pathway = 50, max.genes =1000,
                       return.gene.list = FALSE, min.pct = 0.1)
pgoe3 <- DEenrichRPlot(object = pbmc, ident.1 = 3, logfc.threshold = 0.8, test.use = "MAST", enrich.database = "KEGG_2019_Mouse", num.pathway = 50, max.genes =1000,
                       return.gene.list = FALSE, min.pct = 0.1)
pgoe4 <- DEenrichRPlot(object = pbmc, ident.1 = 4, logfc.threshold = 0.8, test.use = "MAST", enrich.database = "KEGG_2019_Mouse", num.pathway = 50, max.genes =1000,
                       return.gene.list = FALSE, min.pct = 0.1)
pgoe5 <- DEenrichRPlot(object = pbmc, ident.1 = 5, logfc.threshold = 0.8, test.use = "MAST", enrich.database = "KEGG_2019_Mouse", num.pathway = 50, max.genes =1000,
                       return.gene.list = FALSE, min.pct = 0.1)
pgoe6 <- DEenrichRPlot(object = pbmc, ident.1 = 6, logfc.threshold = 0.8, test.use = "MAST", enrich.database = "KEGG_2019_Mouse", num.pathway = 50, max.genes =1000,
                       return.gene.list = FALSE, min.pct = 0.1)
pgoe7 <- DEenrichRPlot(object = pbmc, ident.1 = 7, logfc.threshold = 0.8, test.use = "MAST", enrich.database = "KEGG_2019_Mouse", num.pathway = 50, max.genes =1000,
                       return.gene.list = FALSE, min.pct = 0.1)
pgoe8 <- DEenrichRPlot(object = pbmc, ident.1 = 8, logfc.threshold = 0.8, test.use = "MAST", enrich.database = "KEGG_2019_Mouse", num.pathway = 50, max.genes =1000,
                       return.gene.list = FALSE, min.pct = 0.1)
pgoe9 <- DEenrichRPlot(object = pbmc, ident.1 = 9, logfc.threshold = 0.8, test.use = "MAST", enrich.database = "KEGG_2019_Mouse", num.pathway = 50, max.genes =1000,
                       return.gene.list = FALSE, min.pct = 0.1)
pgoe10 <- DEenrichRPlot(object = pbmc, ident.1 = 10, logfc.threshold = 0.8, test.use = "MAST", enrich.database = "KEGG_2019_Mouse", num.pathway = 50, max.genes =1000,
                        return.gene.list = FALSE, min.pct = 0.1)
pgoe11 <- DEenrichRPlot(object = pbmc, ident.1 = 11, logfc.threshold = 0.8, test.use = "MAST", enrich.database = "KEGG_2019_Mouse", num.pathway = 50, max.genes =1000,
                        return.gene.list = FALSE, min.pct = 0.1)
pgoe12 <- DEenrichRPlot(object = pbmc, ident.1 = 12, logfc.threshold = 0.8, test.use = "MAST", enrich.database = "KEGG_2019_Mouse", num.pathway = 50, max.genes =1000,
                        return.gene.list = FALSE, min.pct = 0.1)
pgoe13 <- DEenrichRPlot(object = pbmc, ident.1 = 13, logfc.threshold = 0.8, test.use = "MAST", enrich.database = "KEGG_2019_Mouse", num.pathway = 50, max.genes =1000,
                        return.gene.list = FALSE, min.pct = 0.1)
pgoe14 <- DEenrichRPlot(object = pbmc, ident.1 = 14, logfc.threshold = 0.8, test.use = "MAST", enrich.database = "KEGG_2019_Mouse", num.pathway = 50, max.genes =1000,
                        return.gene.list = FALSE, min.pct = 0.1)
pgoe15 <- DEenrichRPlot(object = pbmc, ident.1 = 15, logfc.threshold = 0.8, test.use = "MAST", enrich.database = "KEGG_2019_Mouse", num.pathway = 50, max.genes =1000,
                        return.gene.list = FALSE, min.pct = 0.1)

SaveFigure(pgoe0, "GOE_KEGG_cluster0", width = 20, height = 15)
SaveFigure(pgoe1, "GOE_KEGG_cluster1", width = 20, height = 15)
SaveFigure(pgoe2, "GOE_KEGG_cluster2", width = 20, height = 15)
SaveFigure(pgoe3, "GOE_KEGG_cluster3", width = 20, height = 15)
SaveFigure(pgoe4, "GOE_KEGG_cluster4", width = 20, height = 15)
SaveFigure(pgoe5, "GOE_KEGG_cluster5", width = 20, height = 15)
SaveFigure(pgoe6, "GOE_KEGG_cluster6", width = 20, height = 15)
SaveFigure(pgoe7, "GOE_KEGG_cluster7", width = 20, height = 15)
SaveFigure(pgoe8, "GOE_KEGG_cluster8", width = 20, height = 15)
SaveFigure(pgoe9, "GOE_KEGG_cluster9", width = 20, height = 15)
SaveFigure(pgoe10, "GOE_KEGG_cluster10", width = 20, height = 15)
SaveFigure(pgoe11, "GOE_KEGG_cluster11", width = 20, height = 15)
SaveFigure(pgoe12, "GOE_KEGG_cluster12", width = 20, height = 15)
SaveFigure(pgoe13, "GOE_KEGG_cluster13", width = 20, height = 15)
SaveFigure(pgoe14, "GOE_KEGG_cluster14", width = 20, height = 15)
SaveFigure(pgoe15, "GOE_KEGG_cluster15", width = 20, height = 15)



