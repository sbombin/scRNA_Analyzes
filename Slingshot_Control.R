library(slingshot)
library(Seurat)
library(tradeSeq)

data_path <- "/Users/sbombin/Desktop/analysis_mm39.comb/Seurat_Analysis/IEC/"
fig_path <- "/Users/sbombin/Desktop/analysis_mm39.comb/Seurat_Analysis/IEC/Traj_Infer/Control/"
setwd("/Users/sbombin/Desktop/analysis_mm39.comb/Seurat_Analysis/IEC/")

object <- ReadObject("seuratobj.IEC.scale_all")

object <- RenameIdents(object, `0` = "EPL1", `1` = "MEP1", `2` = "MEP2",`3` = "TA", `4` = "MED", `5` = "IEP", `6` = "MEP3", `7` = "EPE1", `8` = "ISC",
                       `9` = "EPE2", `10` = "EPL2", `11` = "Goblet", `12` = "Tuft", `13` = "EEC", `14` = "MEP4", `15` = "Paneth")
#head(object@meta.data)

control <- subset(x = object, subset = group == "Tet1con" & percent.mt < 10) 

ctrl <- as.SingleCellExperiment(control, assay = "RNA")


ctrl1 <- slingshot(ctrl, reducedDim = 'UMAP', clusterLabels = colData(ctrl)$ident, start.clus = 'ISC', approx_points = FALSE)
ctrl2 <- slingshot(ctrl, reducedDim = 'UMAP', clusterLabels = colData(ctrl)$ident, approx_points = FALSE)
ctrl3 <- slingshot(ctrl, reducedDim = 'UMAP', clusterLabels = colData(ctrl)$ident, start.clus = 'ISC',
                   end.clus = c("Paneth","Goblet","EEC","Tuft","MED","MEP1","MEP2","MEP3","MEP4"), approx_points = FALSE, dist.method = "mnn")
#ctrl4 <- slingshot(ctrl, reducedDim = 'UMAP', clusterLabels = colData(ctrl)$ident, start.clus = 'ISC', omega = TRUE, approx_points = FALSE, dist.method = "mnn")

ctrl5 <- slingshot(ctrl, reducedDim = 'PCA', clusterLabels = colData(ctrl)$ident, start.clus = 'ISC',
                   end.clus = c("Paneth","Goblet","EEC","Tuft","MED","MEP1","MEP2","MEP3","MEP4"), approx_points = FALSE)
ctrl6 <- slingshot(ctrl, reducedDim = 'PCA', clusterLabels = colData(ctrl)$ident,start.clus = 'ISC', approx_points = FALSE)

ctrl1@colData@listData[["slingshot"]]@metadata[["lineages"]]
ctrl6@colData@listData[["slingshot"]]@metadata[["lineages"]]

# dev.off()
## UMAP_1a
plot.new() ## clean up device
par(mar=c(5, 5, 5, 7) + 0.1, xpd=TRUE)
plot(reducedDims(ctrl1)$UMAP, col = colors[ctrl1$ident], pch=16, cex = 0.5)
lines(SlingshotDataSet(ctrl1), lwd=2, type = 'lineages', col = 'black', show.constraints = TRUE)
legend("topright", pch = 16, col = colors, cex = 0.8, pt.cex=1.5, y.intersp=0.65, bty="n", xpd = TRUE,
       inset = c(-0.28, 0), legend = levels(factor(colData(ctrl1)$ident)))

p1 <- recordPlot()
SaveFigure(p1, "UMAP_Lineages_Cont_1c", width = 10, height = 8)

plot.new() ## clean up device
par(mar=c(5, 5, 5, 7) + 0.1, xpd=TRUE)
plot(reducedDims(ctrl1)$UMAP, col = colors[ctrl1$ident], pch=16, cex = 0.5)
lines(SlingshotDataSet(ctrl1), lwd=2, col = 'black', show.constraints = TRUE)
legend("topright", pch = 16, col = colors, cex = 0.8, pt.cex=1.5, y.intersp=0.65, bty="n", xpd = TRUE,
       inset = c(-0.28, 0), legend = levels(factor(colData(ctrl1)$ident)))

p2 <- recordPlot()
SaveFigure(p2, "UMAP_Curves_Cont_1c", width = 10, height = 8)

## UMAP_1b
plot.new() ## clean up device
par(mar=c(5, 5, 5, 7) + 0.1, xpd=TRUE)
plot(reducedDims(ctrl2)$UMAP, col = colors[ctrl2$ident], pch=16, cex = 0.5)
lines(SlingshotDataSet(ctrl2), lwd=2, type = 'lineages', col = 'black', show.constraints = TRUE)
legend("topright", pch = 16, col = colors, cex = 0.8, pt.cex=1.5, y.intersp=0.65, bty="n", xpd = TRUE,
       inset = c(-0.28, 0), legend = levels(factor(colData(ctrl2)$ident)))

p1 <- recordPlot()
SaveFigure(p1, "UMAP_Lineages_Cont_1d", width = 10, height = 8)

plot.new() ## clean up device
par(mar=c(5, 5, 5, 7) + 0.1, xpd=TRUE)
plot(reducedDims(ctrl2)$UMAP, col = colors[ctrl2$ident], pch=16, cex = 0.5)
lines(SlingshotDataSet(ctrl2), lwd=2, col = 'black', show.constraints = TRUE)
legend("topright", pch = 16, col = colors, cex = 0.8, pt.cex=1.5, y.intersp=0.65, bty="n", xpd = TRUE,
       inset = c(-0.28, 0), legend = levels(factor(colData(ctrl2)$ident)))

p2 <- recordPlot()
SaveFigure(p2, "UMAP_Curves_Cont_1d", width = 10, height = 8)

## UMAP_2a
plot.new() ## clean up device
par(mar=c(5, 5, 5, 7) + 0.1, xpd=TRUE)
plot(reducedDims(ctrl3)$UMAP, col = colors[ctrl3$ident], pch=16, cex = 0.5)
lines(SlingshotDataSet(ctrl3), lwd=2, type = 'lineages', col = 'black', show.constraints = TRUE)
legend("topright", pch = 16, col = colors, cex = 0.8, pt.cex=1.5, y.intersp=0.65, bty="n", xpd = TRUE,
       inset = c(-0.28, 0), legend = levels(factor(colData(ctrl3)$ident)))

p1 <- recordPlot()
SaveFigure(p1, "UMAP_Lineages_Cont_2c", width = 10, height = 8)

plot.new() ## clean up device
par(mar=c(5, 5, 5, 7) + 0.1, xpd=TRUE)
plot(reducedDims(ctrl3)$UMAP, col = colors[ctrl3$ident], pch=16, cex = 0.5)
lines(SlingshotDataSet(ctrl3), lwd=2, col = 'black', show.constraints = TRUE)
legend("topright", pch = 16, col = colors, cex = 0.8, pt.cex=1.5, y.intersp=0.65, bty="n", xpd = TRUE,
       inset = c(-0.28, 0), legend = levels(factor(colData(ctrl3)$ident)))

p2 <- recordPlot()
SaveFigure(p2, "UMAP_Curves_Cont_2c", width = 10, height = 8)

## UMAP_2b
plot.new() ## clean up device
par(mar=c(5, 5, 5, 7) + 0.1, xpd=TRUE)
plot(reducedDims(ctrl4)$UMAP, col = colors[ctrl4$ident], pch=16, cex = 0.5)
lines(SlingshotDataSet(ctrl4), lwd=2, type = 'lineages', col = 'black', show.constraints = TRUE)
legend("topright", pch = 16, col = colors, cex = 0.8, pt.cex=1.5, y.intersp=0.65, bty="n", xpd = TRUE,
       inset = c(-0.26, 0), legend = levels(factor(colData(ctrl4)$ident)))

p1 <- recordPlot()
SaveFigure(p1, "UMAP_Lineages_Cont_2b", width = 10, height = 8)

plot.new() ## clean up device
par(mar=c(5, 5, 5, 7) + 0.1, xpd=TRUE)
plot(reducedDims(ctrl4)$UMAP, col = colors[ctrl4$ident], pch=16, cex = 0.5)
lines(SlingshotDataSet(ctrl4), lwd=2, col = 'black', show.constraints = TRUE)
legend("topright", pch = 16, col = colors, cex = 0.8, pt.cex=1.5, y.intersp=0.65, bty="n", xpd = TRUE,
       inset = c(-0.26, 0), legend = levels(factor(colData(ctrl4)$ident)))

p2 <- recordPlot()
SaveFigure(p2, "UMAP_Curves_Cont_2b", width = 10, height = 8)

## PCA_1a
plot.new() ## clean up device
par(mar=c(5, 5, 5, 7) + 0.1, xpd=TRUE)
plot(reducedDims(ctrl5)$PCA, col = colors[ctrl5$ident], pch=16, cex = 0.5)
lines(SlingshotDataSet(ctrl5), lwd=2, type = 'lineages', col = 'black', show.constraints = TRUE)
legend("topright", pch = 16, col = colors, cex = 0.8, pt.cex=1.5, y.intersp=0.65, bty="n", xpd = TRUE,
       inset = c(-0.28, 0), legend = levels(factor(colData(ctrl5)$ident)))

p1 <- recordPlot()
SaveFigure(p1, "PCA_Lineages_Cont_1c", width = 10, height = 8)

plot.new() ## clean up device
par(mar=c(5, 5, 5, 7) + 0.1, xpd=TRUE)
plot(reducedDims(ctrl5)$PCA, col = colors[ctrl5$ident], pch=16, cex = 0.5)
lines(SlingshotDataSet(ctrl5), lwd=2, col = 'black', show.constraints = TRUE)
legend("topright", pch = 16, col = colors, cex = 0.8, pt.cex=1.5, y.intersp=0.65, bty="n", xpd = TRUE,
       inset = c(-0.28, 0), legend = levels(factor(colData(ctrl5)$ident)))

p2 <- recordPlot()
SaveFigure(p2, "PCA_Curves_Cont_1c", width = 10, height = 8)

## PCA_1d
plot.new() ## clean up device
par(mar=c(5, 5, 5, 7) + 0.1, xpd=TRUE)
plot(reducedDims(ctrl6)$PCA, col = colors[ctrl6$ident], pch=16, cex = 0.5)
lines(SlingshotDataSet(ctrl6), lwd=2, type = 'lineages', col = 'black', show.constraints = TRUE)
legend("topright", pch = 16, col = colors, cex = 0.8, pt.cex=1.5, y.intersp=0.65, bty="n", xpd = TRUE,
       inset = c(-0.28, 0), legend = levels(factor(colData(ctrl6)$ident)))

p1 <- recordPlot()
SaveFigure(p1, "PCA_Lineages_Cont_1d", width = 10, height = 8)

plot.new() ## clean up device
par(mar=c(5, 5, 5, 7) + 0.1, xpd=TRUE)
plot(reducedDims(ctrl6)$PCA, col = colors[ctrl6$ident], pch=16, cex = 0.5)
lines(SlingshotDataSet(ctrl6), lwd=2, col = 'black', show.constraints = TRUE)
legend("topright", pch = 16, col = colors, cex = 0.8, pt.cex=1.5, y.intersp=0.65, bty="n", xpd = TRUE,
       inset = c(-0.28, 0), legend = levels(factor(colData(ctrl6)$ident)))

p2 <- recordPlot()
SaveFigure(p2, "PCA_Curves_Cont_1d", width = 10, height = 8)

data3 <- SlingshotDataSet(ctrl3)
data5 <- SlingshotDataSet(ctrl5)

plot3d.SlingshotDataSet(data3, type = 'b')
plot3d.SlingshotDataSet(data5, type = 'curves', col = colors[ctrl5$ident], pch=16, cex = 0.5)


pt <- slingPseudotime(ctrl1)

