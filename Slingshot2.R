
plot.new() ## clean up device
par(mar=c(5, 5, 5, 7) + 0.1, xpd=TRUE)
plot(reducedDims(sce2)$UMAP, col = colors[sce2$ident], pch=16, cex = 0.5)
#lines(SlingshotDataSet(sce2), lwd=2, type = 'lineages', col = 'black', show.constraints = TRUE)
lines(SlingshotDataSet(sce2), lwd=2, col = 'black', show.constraints = TRUE)
legend("topright", pch = 16, col = colors, cex = 0.8, pt.cex=1.5, y.intersp=0.65, bty="n", xpd = TRUE,
       inset = c(-0.25, 0), legend = levels(factor(colData(sce2)$ident)))


dev.off()
plot.new() ## clean up device
par(mar=c(5, 5, 5, 7) + 0.1, xpd=TRUE)
plot(reducedDims(sce2)$UMAP, col = colors[sce2$ident], pch=16, cex = 0.5)
# lines(SlingshotDataSet(sce2), lwd=2, type = 'lineages', col = 'black', show.constraints = TRUE)
# lines(SlingshotDataSet(sce2), lwd=2, col = 'black', show.constraints = TRUE)
legend("topright", pch = 16, col = colors, cex = 0.8, pt.cex=1.5, y.intersp=0.65, bty="n", xpd = TRUE,
       inset = c(-0.25, 0), legend = levels(factor(colData(sce2)$ident)))

p <- recordPlot()
SaveFigure(p, "UMAP_Lineages_1b", width = 10, height = 8)
SaveFigure(p, "UMAP_Curves_1b", width = 10, height = 8)


n <- ncol(sce2)
L <- ncol(slingPseudotime(sce2))
plot(as.numeric(slingPseudotime(sce2)), jitter(rep(1:L, each=n)),pch=16, col = colors[sce2$ident])
legend("topright", pch = 16, col = colors, cex = 0.8, pt.cex=1.5, y.intersp=0.65, bty="n", xpd = TRUE,
       inset = c(-0.25, 0), legend = levels(factor(colData(sce2)$ident)))

p <- recordPlot()
SaveFigure(p, "Jitter_Cells_Lineages_1a", width = 12, height = 8)

plot(reducedDims(sce2)$UMAP, asp=0.5, col = 'grey75', pch = 16)
points(reducedDims(sce2)$UMAP, col = hcl.colors(100)[cut(slingPseudotime(sce2)[,5], 100)], pch=16)
legend('topleft', title = 'Lineage 1', col = hcl.colors(4), legend=c('0%','25%','50%','100%'), pch=16)


df2 <- (colData(sce2))
df2 <- df2[,c(1:19,21:28)] %>% as.data.frame()

df3 <- (colData(sce2))
df3 <- df3[,c(1:19,21:26)] %>% as.data.frame()

ggplot(df3, aes(x = slingPseudotime_6, y = ident, colour = ident)) +
  geom_quasirandom(groupOnX = FALSE) + theme_classic() +
  xlab("First Slingshot pseudotime") + ylab("cell type") +
  ggtitle("Cells ordered by Slingshot pseudotime") #+scale_colour_manual(values = pal)


## Print Clusters included in the lineage
#All lineages
sce5@colData@listData[["slingshot"]]@metadata[["lineages"]]
# One lineage
sce2@colData@listData[["slingshot"]]@metadata[["lineages"]][["Lineage1"]]

data3 <- SlingshotDataSet(sce2)
data4 <- SlingshotDataSet(sce2)
plot(data4, type = 'b')
slingBranchGraph(data4)
slingLineages(data4)
slingClusterLabels(data4)

plot3d.SlingshotDataSet(data3, type = 'curves')
pairs(data4, type = 'curves')

## Plot cells in lineages
nc <- 3
pt <- slingPseudotime(sce2)
nms <- colnames(pt)
nr <- ceiling(length(nms)/nc)
pal <- hcl.colors(100)

par(mfrow = c(nr, nc))
options(repr.plot.width=15, repr.plot.height=6)
for (i in nms) {
  colors2 <- pal[cut(pt[,i], breaks = 100)]
  plot(reducedDim(sce2, "UMAP"), col = colors2, pch = 16, main = i)
}

p <- recordPlot()
SaveFigure(p, "UMAP_Cells_Lineages_2a", width = 12, height = 10)



## Plot cells in lineages
dev.off()
plot.new()

nc <- 3
pt <- slingPseudotime(sce2)
nms <- colnames(pt)
nr <- ceiling(length(nms)/nc)
pal <- hcl.colors(100)

par(mfrow = c(nr, nc))
#par(mar=c(5, 5, 5, 7) + 0.1, xpd=TRUE)
options(repr.plot.width=15, repr.plot.height=6)
for (i in nms) {
  colors2 <- pal[cut(pt[,i], breaks = 100)]
  plot(reducedDims(sce2)$UMAP, col = 'grey75', pch = 16,main=i)
  points(reducedDim(sce2, "UMAP"), col = colors2, pch = 16)
}

p <- recordPlot()
SaveFigure(p, "UMAP_Cells_Lineages_1b", width = 12, height = 10)


## If stopped plotting
dev.off()
print(plot(1))



## Plot Only one lineage
# full plot
plot(umap, col = colvec)
# just lineage 1
plot(umap, col = scales::alpha(colvec, alpha = slingCurveWeights(sce)[,1]))
# with weights threshold
lin1 <- which(slingCurveWeights(sce)[,1] > .5)
plot(umap[lin1,], col = colvec[lin1])
