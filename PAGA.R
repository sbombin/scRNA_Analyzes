library(Seurat)
library(SeuratWrappers)
library(SeuratDisk)

data_path <- "/Users/sbombin/Desktop/analysis_mm39.comb/Seurat_Analysis/IEC/"
fig_path <- "/Users/sbombin/Desktop/analysis_mm39.comb/Seurat_Analysis/IEC/Traj_Infer/"
setwd("/Users/sbombin/Desktop/analysis_mm39.comb/Seurat_Analysis/IEC/")

object <- ReadObject("seuratobj.IEC.scale_all")

object <- RenameIdents(object, `0` = "EPL1", `1` = "MEP1", `2` = "MEP2",`3` = "TA", `4` = "MED", `5` = "IEP", `6` = "MEP3", `7` = "EPE1", `8` = "ISC",
                       `9` = "EPE2", `10` = "EPL2", `11` = "Goblet", `12` = "Tuft", `13` = "EEC", `14` = "MEP4", `15` = "Paneth")
#head(object@meta.data)
#npcs <- fgetNPCs(object,MIN_CUM_SD)

## Control
control <- subset(x = object, subset = group == "Tet1con") 
control@meta.data$cells <- Idents(control)
## KO
ko <- subset(x = object, subset = group == "Tet1iKO") 
ko@meta.data$cells <- Idents(ko)
## Full
object@meta.data$cells <- Idents(object)

npcs2 <- fgetNPCs(control,MIN_CUM_SD)
control <- FindNeighbors(control, dims = 1:npcs)
## Rebuild Neighborhood Graph


### Convert to AnnData
#library(reticulate)
#ad <- import("anndata", convert = FALSE)
#control.ad <- Convert(from = control, to = "anndata", filename = "iec_controls.h5ad")
## Method 2
#library(SeuratDisk)
SaveH5Seurat(object, filename = "iec_full.h5Seurat")
Convert("iec_full.h5Seurat", dest = "h5ad")
########


install.packages("renv")
renv::init()
library(Seurat)
renv::use_python()


obs <- read.csv("~/Desktop/analysis_mm39.comb/Seurat_Analysis/IEC/Scanpy/filename/obs.csv", row.names=1)
p3 <- ggplot(obs, aes(dpt_pseudotime, reorder(cells, dpt_pseudotime, median), fill = cells)) +
  geom_boxplot()
p3

SaveFigure(p3, "BoxPlot_PsdT_Cont", width = 10, height = 8)


ctrl <- read.csv("~/Desktop/analysis_mm39.comb/Seurat_Analysis/IEC/Scanpy/scanpy_control/obs.csv", row.names=1)
p1 <- ggplot(ctrl, aes(dpt_pseudotime, reorder(cells, dpt_pseudotime, median), fill = cells)) +
  geom_boxplot()
p1

SaveFigure(p1, "BoxPlot_PsdT_Control", width = 10, height = 8)

ko <- read.csv("~/Desktop/analysis_mm39.comb/Seurat_Analysis/IEC/Scanpy/scanpy_ko/obs.csv", row.names=1)
p2 <- ggplot(ko, aes(dpt_pseudotime, reorder(cells, dpt_pseudotime, median), fill = cells)) +
  geom_boxplot()
p2

SaveFigure(p2, "BoxPlot_PsdT_KO", width = 10, height = 8)


p2 <- ggplot(ko, aes(dpt_pseudotime, y=cells,  fill = cells)) +
  geom_boxplot()
p2

#####
full <- read.csv("~/Desktop/analysis_mm39.comb/Seurat_Analysis/IEC/Scanpy/scanpy_full/obs.csv", row.names=1)

p3 <- ggplot(full, aes(dpt_pseudotime, reorder(cells, dpt_pseudotime, median), fill = cells)) +
  geom_boxplot()
p3

SaveFigure(p3, "BoxPlot_PsdT_Full", width = 10, height = 8)


## Add pseudotime to Seurat
object <- AddMetaData(object, metadata = full[,19:20])

p1 <- VlnPlot(object,features = "dpt_pseudotime", group.by = "group")
p1
SaveFigure(p1, "VlnPlot_PsdT_Group", width = 10, height = 8)

p2 <- VlnPlot(object,features = "dpt_pseudotime", group.by = "sample")
SaveFigure(p2, "VlnPlot_PsdT_Sample", width = 10, height = 8)

p <- FeatureScatter(object, feature1 = "dpt_pseudotime", feature2 = "Ascl2", group.by = "group", slot = 'data', span = NULL)
a <- as.data.frame(p$data)
zz <- ggplot(a, aes(x = dpt_pseudotime, y = Ascl2, color = colors)) + geom_point() + geom_smooth(method = 'loess',formula = y ~ x)
ggsave(paste0("Scatter",'_Ascl2.png'),plot=zz,dpi=300, scale = 1.5)      

p2 <- FeatureScatter(object, feature1 = "dpt_pseudotime", feature2 = "Sox9", group.by = "group", slot = 'data', span = NULL)

b <- as.data.frame(p2$data)
z <- p2$data

ggplot(b, aes(x = dpt_pseudotime, y = Sox9, color = colors)) + geom_point() + geom_smooth(method = 'loess',formula = y ~ x)
# ggsave(paste0(comparison,'_VolcanoPlot.png'),plot=p,dpi=300, scale = 2.0)

ggsave(paste0("Scatter",'_Sox9.png'),dpi=300, scale = 1.5)

x <- list(c("Sox9","Ascl2"))



genes <- c('Ascl2','Fabp1','Mafb','Atoh1','Fabp6','Mex3a','Cdx1','Gata4','Muc2','Cdx2','Hes1','Neurog3',
           'Chga','Hnf4a','Pou2f3','Maf','Lgr5','Sis','Dclk1','Lyz1','Sox4','Elf3','Lyz2','Sox9')
genes <- genes[genes %in% row.names(object@assays[["RNA"]]@counts)]

### Not Working
Expression <- function(x){
  p0 <- FeatureScatter(object, feature1 = "dpt_pseudotime", feature2 = x, group.by = "group", slot = 'data', span = NULL)
  data <- p0$data
  p <- ggplot(data, aes(x = dpt_pseudotime, y = x, color = colors)) + geom_point() + geom_smooth(method = 'loess',formula = y ~ x)
  ggsave(paste0("Scatter",x),dpi=300, scale = 1.5)
}
lapply(gene.list,Expression(x))
######


p0 <- lapply(genes, function(x) FeatureScatter(object, feature1 = "dpt_pseudotime", feature2 = x, group.by = "group", slot = 'data', span = NULL))
names(p0) <- genes
data <- lapply(p0, function(x) x$data)
colnames <- c("dpt_pseudotime","expression","colors")
data <- lapply(data, setNames, colnames)                                
                
                
#p <- lapply(seq_along(data2), function(z) ggplot(z, aes(x = dpt_pseudotime, y = expression, color = colors)) + geom_point() +
#              geom_smooth(method = 'loess',formula = y ~ x))


## Plot with title as list name
lapply(seq_along(data), function(i) {
  ggplot(data[[i]], aes(x=dpt_pseudotime, y=expression,color = colors)) + 
    geom_point(size = 0.001)+ geom_smooth(method = 'loess',formula = y ~ x)+
    ggtitle(names(data[i]))
  ggsave(paste0('Scatter_',names(data[i]),'_dot0.png'),dpi=300, scale =1.5)})

## Separately
ggplot(data[["Maf"]], aes(x=dpt_pseudotime, y=expression,color = colors)) + 
  geom_point()+ geom_smooth(method = 'loess',formula = y ~ x)+
  ggtitle(names(data)["Maf"])
ggsave(paste0('Scatter_Maf.png'),dpi=300, scale =1.5)


p <- ggplot(data[["Maf"]], aes(x=dpt_pseudotime, y=expression,color = colors)) + 
  geom_point(shape = ".") + geom_smooth(method = 'loess',formula = y ~ x)+
  ggtitle(names(data["Maf"]))
p
ggsave(paste0('Scatter_Maf_dot001.png'),plot = p,dpi=300, scale =1.5)


### Plot with smallest possible points
lapply(seq_along(data), function(i) {
  ggplot(data[[i]], aes(x=dpt_pseudotime, y=expression,color = colors)) + 
    geom_point(shape = ".")+ geom_smooth(method = 'loess',formula = y ~ x)+
    ggtitle(names(data[i]))
  ggsave(paste0('Scatter_',names(data[i]),'_dot001.png'),dpi=300, scale =1.2)})

### Plot without dots
lapply(seq_along(data), function(i) {
  ggplot(data[[i]], aes(x=dpt_pseudotime, y=expression,color = colors)) + 
    geom_smooth(method = 'loess',formula = y ~ x)+
    ggtitle(names(data[i]))
  ggsave(paste0('Scatter_',names(data[i]),'_dot0.png'),dpi=300, scale =1.2)})
