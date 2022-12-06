library(tibble)

a <- DEenrichRPlot(object = pbmc, ident.1 = 0,  cells.1 = "Sham", cells.2 = "BDL",
                   min.pct = 0.1, logfc.threshold = 0.8, test.use = "MAST", only.pos = FALSE,
                   enrich.database = "KEGG_2019_Mouse", num.pathway = 50, max.genes =1000,return.gene.list = TRUE)
a1 <- as.data.frame(a)

de_markers1 <- FindMarkers(object = pbmc,ident.1 = "Sham", group.by = "group", min.pct = 0.1, logfc.threshold = 0.8,
                           test.use = "MAST", only.pos = FALSE,subset.ident ="0")
b <- rownames_to_column(de_markers1, var = "rowname")


p1 <- DEenrichRPlot(object = pbmc, ident.1 = 1,  cells.1 = "Sham", cells.2 = "BDL",
                    min.pct = 0.1, logfc.threshold = 0.8, test.use = "MAST", only.pos = FALSE,
                    enrich.database = "KEGG_2019_Mouse", num.pathway = 50, max.genes =1000,return.gene.list = FALSE)

p2 <- DEenrichRPlot(object = pbmc, ident.1 = "1",  cells.1 = "Sham", cells.2 = "BDL", subset.ident = "1",
                    min.pct = 0.1, logfc.threshold = 0.8, test.use = "MAST", only.pos = FALSE,
                    enrich.database = "KEGG_2019_Mouse", num.pathway = 50, max.genes =1000,return.gene.list = FALSE)

p3 <- DEenrichRPlot(object = pbmc, ident.1 = "1",  group.by = "group", subset.ident = "1",
                    min.pct = 0.1, logfc.threshold = 0.8, test.use = "MAST", only.pos = FALSE,
                    enrich.database = "KEGG_2019_Mouse", num.pathway = 50, max.genes =1000,return.gene.list = FALSE)


obj802 <- ReadObject("obj802.seurat")
obj877 <- ReadObject("obj877.seurat")

obj840 <- ReadObject("obj840.seurat")
obj805 <- ReadObject("obj805.seurat")

pbmc <- merge(obj802, y = c(obj840, obj877,  obj805), add.cell.ids = c("Sham", "Sham", "BDL", "BDL"), project = "SHAM-BDL")

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

pbmc <- FindNeighbors(pbmc, dims = 1:npcs)
pbmc <- FindClusters(pbmc, resolution = 0.6)
pbmc <- RunTSNE(pbmc, dims = 1:npcs)
pbmc <- RunUMAP(pbmc, dims = 1:npcs)