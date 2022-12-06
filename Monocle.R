
#library(monocle)
library(Seurat)
library(SeuratWrappers)
library(monocle3)


data_path <- "/Users/sbombin/Desktop/analysis_mm39.comb/Seurat_Analysis/IEC/"
fig_path <- "/Users/sbombin/Desktop/analysis_mm39.comb/Seurat_Analysis/IEC/Traj_Infer/Control/"
setwd("/Users/sbombin/Desktop/analysis_mm39.comb/Seurat_Analysis/IEC/")

object <- ReadObject("seuratobj.IEC.scale_all")
#object[["percent.mt"]] <- PercentageFeatureSet(object, pattern = "^mt-")
#summary(object[["percent.mt"]])

object <- RenameIdents(object, `0` = "EPL1", `1` = "MEP1", `2` = "MEP2",`3` = "TA", `4` = "MED", `5` = "IEP", `6` = "MEP3", `7` = "EPE1", `8` = "ISC",
                       `9` = "EPE2", `10` = "EPL2", `11` = "Goblet", `12` = "Tuft", `13` = "EEC", `14` = "MEP4", `15` = "Paneth")
#head(object@meta.data)
control <- subset(x = object, subset = group == "Tet1con") 
control@meta.data$clusters <- Idents(control)

### Convert to AnnData
library(reticulate)
ad <- import("anndata", convert = FALSE)
control.ad <- Convert(from = control, to = "anndata", filename = "iec_controls.h5ad")
## Method 2
library(SeuratDisk)
SaveH5Seurat(control, filename = "iec_controls.h5Seurat")
Convert("iec_controls.h5Seurat", dest = "h5ad")
########

## Native Monocle conversion function have not been released yet
#mon <- monocle::importCDS(control, import_all = TRUE)

cds <- as.cell_data_set(control)

# assign paritions
reacreate.partition <- c(rep(1,length(cds@colData@rownames)))
names(reacreate.partition) <- cds@colData@rownames
reacreate.partition <- as.factor(reacreate.partition)


cds@clusters$UMAP$partitions <- reacreate.partition

# Assign the cluster info 

list_cluster <- control@active.ident
cds@clusters$UMAP$clusters <- list_cluster


# Assign UMAP coordinate - cell embeddings
cds@int_colData@listData$reducedDims$UMAP <- control@reductions$umap@cell.embeddings

p1 <- plot_cells(cds, color_cells_by = "cluster", show_trajectory_graph = FALSE)
p2 <- plot_cells(cds, color_cells_by = "partition", show_trajectory_graph = FALSE)
wrap_plots(p1, p2)


cluster.before.trajectory <- plot_cells(cds, color_cells_by = 'cluster', show_trajectory_graph = FALSE,
                                        label_groups_by_cluster = FALSE, group_label_size = 5) + theme(legend.position = "right")


cds <- learn_graph(cds, use_partition = TRUE, verbose = FALSE)

p1 <- plot_cells(cds, color_cells_by = "cluster", label_groups_by_cluster=FALSE, label_leaves=TRUE,
           label_branch_points=TRUE,label_roots = TRUE,group_label_size = 4,show_trajectory_graph=TRUE)
p1
SaveFigure(p1, "UMAP_Lineage_Cont_2", width = 10, height = 8)


cds <- order_cells(cds, root_cells = colnames(cds[,clusters(cds) == "ISC"]))

p2 <- plot_cells(cds, color_cells_by = "pseudotime", group_cells_by = "cluster", label_cell_groups = FALSE, label_groups_by_cluster=FALSE,
           label_leaves=FALSE, label_branch_points=FALSE, label_roots = FALSE, trajectory_graph_color = "grey60")
p2

SaveFigure(p2, "UMAP_Pseudotime_Cont", width = 10, height = 8)

cds$monocle3_pseudotime <- pseudotime(cds)
data.pseudo <- as.data.frame(colData(cds))

p3 <- ggplot(data.pseudo, aes(monocle3_pseudotime, reorder(ident, monocle3_pseudotime, median), fill = ident)) +
  geom_boxplot()

p3

SaveFigure(p3, "BoxPlot_Pseudotime_Cont", width = 10, height = 8)


cds_sub <- choose_graph_segments(cds)
