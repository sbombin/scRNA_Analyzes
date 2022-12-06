
data_path <- "/Users/sbombin/Desktop/analysis_mm39.comb/Seurat_Analysis/IEC/"
fig_path <- "/Users/sbombin/Desktop/analysis_mm39.comb/Seurat_Analysis/IEC/Traj_Infer/"
setwd("/Users/sbombin/Desktop/analysis_mm39.comb/Seurat_Analysis/IEC/")

object <- ReadObject("seuratobj.IEC.scale_all")

object <- RenameIdents(object, `0` = "EPL1", `1` = "MEP1", `2` = "MEP2",`3` = "TA", `4` = "MED", `5` = "IEP", `6` = "MEP3", `7` = "EPE1", `8` = "ISC",
                       `9` = "EPE2", `10` = "EPL2", `11` = "Goblet", `12` = "Tuft", `13` = "EEC", `14` = "MEP4", `15` = "Paneth")

full <- read.csv("~/Desktop/analysis_mm39.comb/Seurat_Analysis/IEC/Scanpy/scanpy_full/obs.csv", row.names=1)

### mtDNA check
object[["percent.mt"]] <- PercentageFeatureSet(object, pattern = "^mt-")
summary(object[["percent.mt"]])

#FeatureScatter(object, feature1 = "nCount_RNA", feature2 = "percent.mt", span = TRUE)
#VlnPlot(object, features = "percent.mt", group.by = "group")

a <- object[["percent.mt"]]
b <- a %>% filter(percent.mt < 10)

## Control
control <- subset(x = object, subset = group == "Tet1con" & percent.mt < 15) 
table(Idents(control))

mast <- FindAllMarkers(object, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25,test.use = "MAST")
roc <- FindAllMarkers(object, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25,test.use = "roc")

mast %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC) -> top10.mast
roc %>% group_by(cluster) %>% top_n(n = 10, wt = power) -> top10.roc

write.csv(top10.mast,"Top10_MAST_Control.csv", row.names = FALSE)
write.csv(top10.roc,"Top10_ROC_Control.csv", row.names = FALSE)

top10.mast <- read_csv("Traj_Infer/Exp_over_PsdT/Top10_MAST_Control.csv")
top10.roc <- read_csv("Traj_Infer/Exp_over_PsdT/Top10_ROC_Control.csv")

### Metagenes
my.data=GetAssayData(object, slot = "data") ## normalized counts data
clusters <- c("EPL1","MEP1","MEP2","TA","MED","IEP","MEP3","EPE1","ISC","EPE2","EPL2","Goblet", "Tuft","EEC","MEP4","Paneth")
#genes <- top10.mast$gene

#mast.mep1 <- dplyr::filter(top10.mast, cluster == "MEP1") 
#mast.mep1 <- mast.mep1$gene
#mast.mep1.mean <- colMeans(my.data[mast.mep1,],na.rm=T) 

### MAST Metagenes
mast.list <- lapply(clusters, FUN = function(x) {
  x <- dplyr::filter(top10.mast, cluster == x)
  x <- x$gene})
names(mast.list) <- clusters

mast.mean <- lapply(mast.list, FUN = function(x) {
  colMeans(my.data[x,],na.rm=T)})

## Combine all count tables from list to one dataframe
mast.df <- as.data.frame(mast.mean)
## Ad Pseudotipe to combined dataframe
mast.df <- merge(mast.df, full[,c(4,16,20)], by=0, all.x=TRUE)
## Plot Separately
ggplot(mast.df, aes(x=dpt_pseudotime, y=mast.df[,2],color = group)) + 
  geom_point()+ geom_smooth(method = 'loess',formula = y ~ x)+
  ggtitle("MAST Metagene") + ylab(paste0(names(mast.df[2]),'_Metagene_Expression'))
## Plot All
for(i in 2:17) {
  ggplot(mast.df, aes(x=dpt_pseudotime, y=mast.df[,i],color = group)) + 
    geom_point()+ geom_smooth(method = 'loess',formula = y ~ x)+
    ggtitle(paste0("MAST Metagene: ", names(mast.df[i]))) + ylab("Expression")
  ggsave(paste0('Scatter_MAST_',names(mast.df[i]),'.png'),dpi=300, scale =1.5)}

### ROC Metagenes
roc.list <- lapply(clusters, FUN = function(x) {
  x <- dplyr::filter(top10.roc, cluster == x)
  x <- x$gene})
names(roc.list) <- clusters

roc.mean <- lapply(roc.list, FUN = function(x) {
  colMeans(my.data[x,],na.rm=T)})

## Combine all count tables from list to one dataframe
roc.df <- as.data.frame(roc.mean)
## Ad Pseudotipe to combined dataframe
roc.df <- merge(roc.df, full[,c(4,16,20)], by=0, all.x=TRUE)

## Plot All
for(i in 2:17) {
  ggplot(roc.df, aes(x=dpt_pseudotime, y=roc.df[,i],color = group)) + 
    geom_point()+ geom_smooth(method = 'loess',formula = y ~ x)+
    ggtitle(paste0("ROC Metagene: ", names(roc.df[i]))) + ylab("Expression")
  ggsave(paste0('Scatter_ROC_',names(roc.df[i]),'.png'),dpi=300, scale =1.5)}

### Plot with smallest possible points
# MAST
for(i in 2:17) {
  ggplot(mast.df, aes(x=dpt_pseudotime, y=mast.df[,i],color = group)) + 
    geom_point(shape = ".")+ geom_smooth(method = 'loess',formula = y ~ x)+
    ggtitle(paste0("MAST Metagene: ", names(mast.df[i]))) + ylab("Expression")
  ggsave(paste0('Scatter_MAST_',names(mast.df[i]),'_dot001.png'),dpi=300, scale =1.2)}
# ROC
for(i in 2:17) {
  ggplot(roc.df, aes(x=dpt_pseudotime, y=roc.df[,i],color = group)) + 
    geom_point(shape = ".")+ geom_smooth(method = 'loess',formula = y ~ x)+
    ggtitle(paste0("ROC Metagene: ", names(roc.df[i]))) + ylab("Expression")
  ggsave(paste0('Scatter_ROC_',names(roc.df[i]),'_dot001.png'),dpi=300, scale =1.2)}

### Plot without dots
# MAST
for(i in 2:17) {
  ggplot(mast.df, aes(x=dpt_pseudotime, y=mast.df[,i],color = group)) + 
    geom_smooth(method = 'loess',formula = y ~ x)+
    ggtitle(paste0("MAST Metagene: ", names(mast.df[i]))) + ylab("Expression")
  ggsave(paste0('Scatter_MAST_',names(mast.df[i]),'_dot0.png'),dpi=300, scale =1.2)}
# ROC
for(i in 2:17) {
  ggplot(roc.df, aes(x=dpt_pseudotime, y=roc.df[,i],color = group)) + 
    geom_smooth(method = 'loess',formula = y ~ x) +
    ggtitle(paste0("ROC Metagene: ", names(roc.df[i]))) + ylab("Expression")
  ggsave(paste0('Scatter_ROC_',names(roc.df[i]),'_dot0.png'),dpi=300, scale =1.2)}

