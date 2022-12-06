## Slc25a1 = ENSMUSG00000003528.14
#my.data=GetAssayData(clean, slot = "data")
#a <- as.data.frame(my.data)
#a <- as.data.frame(rownames(my.data))

setwd('/Users/sbombin/Desktop/EICC/Projects/Kwong_Jennifer/')
data_path <- "/Users/sbombin/Desktop/EICC/Projects/Kwong_Jennifer/"
fig_path <- "/Users/sbombin/Desktop/EICC/Projects/Kwong_Jennifer/Cleaned/"

## Replace Ensemb IDs with Gene Symbols
library(Matrix)
cleaned <- ReadObject("gene_count_cleaned")
gene_annotate <- read_csv("gene_annotate.csv")
genes <- gene_annotate$gene_short_name
cleaned@Dimnames[[1]] <- genes

clean3 <- ReadObject("seurobj_cardiac_subtrj_clean")

clean3[["percent.mt"]] <- PercentageFeatureSet(clean3, pattern = "^MT-")
p <- VlnPlot(clean3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)
SaveFigure(p, "QC_Plot1", width = 10, height = 8)

p <- FeatureScatter(clean3, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
SaveFigure(p, "QC_Plot2", width = 10, height = 8)

clean<- subset(clean, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 5)
## Normalize 
clean <- NormalizeData(clean, normalization.method = "LogNormalize", scale.factor = 10000)
## Variable Features
clean <- FindVariableFeatures(clean, selection.method = "vst", nfeatures = 2000)
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(clean), 10)
# plot variable features with and without labels
p <- VariableFeaturePlot(clean)
p <- LabelPoints(plot = p, points = top10, repel = TRUE)
SaveFigure(p, "Variable_Features", width = 10, height = 8)
## Scale
all.genes <- rownames(clean)
clean <- ScaleData(clean, features = all.genes)
## PCA
clean <- RunPCA(clean, features = VariableFeatures(object = clean))
## Find optimal PCA dimensions
npcs <- fgetNPCs(clean,MIN_CUM_SD)
## Find Optimal number of PCs
p <- ElbowPlot(clean)
SaveFigure(p, "ElbowPlot", width = 10, height = 8)
clean <- JackStraw(clean, num.replicate = 100, dims = 50)
clean <- ScoreJackStraw(clean, dims = 1:50)
p <- JackStrawPlot(clean, dims = 1:50)
p
SaveFigure(p, "JackStrawPlot", width = 10, height = 8)

clean <- FindNeighbors(clean, dims = 1:42)
## Clustree to find resolution
library(clustree)
clustree <- FindClusters(clean, resolution = c(0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.2,1.4))
p <- clustree(clustree)
SaveFigure(p, "Clustree", width = 8, height = 10)
## Cluster 
clean <- FindClusters(clean, resolution = 0.8)
clean <- RunUMAP(clean, dims = 1:42)

#Idents(clean) = clean$RNA_snn_res.0.7

p <- DimPlot(clean, reduction = "umap", label = TRUE)
SaveFigure(p, "UMAP_Clusters", width = 10, height = 8)

p <- DimPlot(clean, label = FALSE, group.by = "development_stage", reduction = "umap", repel = FALSE)
p
SaveFigure(p, "UMAP_DevStage", width = 10, height = 8)

VlnPlot(clean, features = "Slc25a1", split.by = "development_stage")
VlnPlot(clean, features = "Slc25a1")


p1 <- FeaturePlot(clean, features = "Slc25a1", split.by = "development_stage", max.cutoff = 5, cols = c("grey", "red"))
SaveFigure(p1, "Slc25a1_DevelopStage_Cardiac", width = 20, height = 6)

p2 <- FeaturePlot(clean, features = "Slc25a1", split.by = "development_stage", max.cutoff = 5, cols = c("grey", "red"), combine = FALSE)
p2 <- wrap_plots(plots = p2, ncol = 2)
SaveFigure(p2, "Slc25a1_DevelopStage_Cardiac2", width = 10, height = 10)

p3 <- FeaturePlot(clean, features = "Slc25a1", cols = c("grey", "red"))
SaveFigure(p3, "UMAP_Slc25a1_Cardiac", width = 10, height = 8)

SaveObject(clean, "seurobj_cardiac_subtrj_clean_processed")

markers <- FindAllMarkers(clean, only.pos = FALSE, test.use = "MAST")
write.csv(markers,"DE-MAST_All_Markers_Cardiac.csv", row.names = TRUE)


### Part2 Labeling Clusters
###########################
my.data=GetAssayData(clean, slot = "data")
df <- as.data.frame(my.data)
genes <- unique(row.names(df)) %>% as.data.frame()

clean <- ReadObject("seurobj_cardiac_subtrj_clean_processed")

p <- DotPlot(clean, features = markers2) & coord_flip()
p
VlnPlot(clean, features = markers2)
VlnPlot(clean, features = 'ENSMUSG00000035458.15')


#######################################################
#### To use different resolutions for downstream analysis
Idents(clean) = clean$RNA_snn_res.0.8

DotPlot(clean, features = "Kdr") & coord_flip() 
# Pecam1 1.6.13

markers <- c('C1qa', 'Cav1','Cdh5','Dcn','Ehd3','Fabp4','Gsn','Myh11', 'Npr3', 'Pecam1', 'Tie1', 'Tnni3') # 'Hbb−b1'
fb <- c('Pdgfra', 'Pdgfrb','Col1a1', 'Tcf21','Ddr2','Lamc1','Dpep1','Pcsk6','Dkk3','Tbx20','Gstm5') # Fibroblast
end <- c('Cdh5', 'Emcn', 'Kdr','Cd31','Nfatc1',"Pecam1","Npr3") # Endothelial
cm <- c('Tnnt2', 'Tpma', 'Tnnc1', 'Tnni3', 'Actc1', 'Ryr2')  # Cardiomyocytes
ep <- c('Wt1','Sema3d','Tbx18','Scx','Tcf21', 'Aldh1a2','Upk3b','Upk1b','Tmem255a', 'Gpm6a', 'Dmkn')  # Epicardial
sm <- c('Acta2','Mastn1','Des')  ## Smooth muscle
genes <- c('Lum','Tnnt2','Isl1','Wt1','Epcam')
  
p2 <- VlnPlot(clean, features = end)
Idents(clean) = clean$RNA_snn_res.0.6
p2 <- VlnPlot(clean, features = cm)
Idents(clean) = clean$RNA_snn_res.0.7
p3 <- VlnPlot(clean, features = cm)

SaveFigure(p2, "VlnPlot_Endocardial_Markers", width = 12, height = 8)
SaveFigure(p2, "VlnPlot_Cardiomyocytes_Markers_clean06", width = 12, height = 8)
SaveFigure(p3, "VlnPlot_Cardiomyocytes_Markers_clean07", width = 12, height = 8)

p1 <- DotPlot(clean, features = end) & coord_flip() 
Idents(clean) = clean$RNA_snn_res.0.8
p3 <- VlnPlot(clean, features = genes) & coord_flip() 
Idents(clean) = clean$RNA_snn_res.0.7
p2 <- VlnPlot(clean, features = ) & coord_flip() 
p4 <- VlnPlot(full, features = genes) & coord_flip() 

p3 <- DotPlot(clean, features = "Slc25a1")# & coord_flip() 


SaveFigure(p1, "DotPlot_Endocardial_Markers", width = 12, height = 8)
SaveFigure(p2, "DotPlot_Cardiomyocytes_Markers_clean06", width = 12, height = 8)
SaveFigure(p3, "DotPlot_Slc25a1_res07", width = 8, height = 12)

#roc =object06, roc2 = clean07, roc3 =clean06, roc4 = clean09 
cluster <- dplyr::filter(roc2, cluster == '6')


#### Full raw cardiac
full <- ReadObject("seurobj_cardiac_processed")
VlnPlot(full, features = end)
DotPlot(full, features = ep) & coord_flip()
DimPlot(full, label = TRUE, group.by = "Sub_trajectory_name", reduction = "umap", repel = TRUE)



##############################
### Clean Resolution 0.8 final
##############################
clean <- ReadObject("seurobj_cardiac_subtrj_clean_processed")
Idents(clean) = clean$RNA_snn_res.0.8

clean2 <- RenameIdents(clean, `0` = "PC4", `1` = "EC3", `2` = "PC1",`3` = "CM4", `4` = "PC3", `5` = "CM1", `6` = "MF", `7` = "EC1", `8` = 'PC5',
                       `9` = "PC6", `10` = "PC2", `11` = "CM2", `12` = "CM3", `13` = "EC2")

p <- DimPlot(clean2, reduction = "umap", label = TRUE, repel = FALSE)
SaveFigure(p, "UMAP_Clusters_renm", width = 10, height = 8)

p1 <- DotPlot(clean2, features = ep) & coord_flip() 
p2 <- DotPlot(clean2, features = "Slc25a1", split.by = "development_stage", cols=c('blue','red','green','black','purple'))# & coord_flip() 
p3 <- DotPlot(clean2, features = "Slc25a1")# & coord_flip() 

SaveFigure(p1, "DotPlot_Cardiomyocytes_Markers_renm", width = 12, height = 8)
SaveFigure(p2, "DotPlot_Slc25a1_DevelopStage_renm", width = 8, height = 12)
SaveFigure(p3, "DotPlot_Slc25a1_renm", width = 8, height = 12)



meta <- clean2@meta.data
sum <- meta%>% group_by(RNA_snn_res.0.8,development_stage) %>% dplyr::count()

de <- read_csv("Cleaned/DE_Slc25a1_Full_Embryo_SubTrj.csv")
det <- dplyr::select(de,class,max.expr,status) %>% tibble::column_to_rownames(var ="status") %>% as.data.frame() %>% dplyr::rename(stage=class)

#p2 <- ggplot(det, aes(x=factor(stage), y=log(max.expr), fill = factor(stage))) + geom_boxplot(outlier.shape = NA) +
#  geom_jitter(shape=16, position=position_jitter(0.2)) + ggrepel::geom_text_repel(aes(label = row.names(det)))
#p2 

p2 <- ggplot(det, aes(x=factor(stage), y=log(max.expr), fill = factor(stage)))  +
  geom_dotplot(binaxis='y', stackdir='center',dotsize=1) + ggrepel::geom_text_repel(aes(label = row.names(det)))
p2 

ggsave(paste0('Slc25a1_MaxExpr_Embryo_SubTraj.png'),plot=p2,dpi=300, scale = 1.8)

de <- read_csv("Cleaned/DE_Slc25a1_Full_Embryo_MainTrj.csv")
det <- dplyr::select(de,class,max.expr,status) %>% tibble::column_to_rownames(var ="status") %>% as.data.frame() %>% dplyr::rename(stage=class)

p2 <- ggplot(det, aes(x=factor(stage), y=log(max.expr), fill = factor(stage)))  +
  geom_dotplot(binaxis='y', stackdir='center',dotsize=1) + ggrepel::geom_text_repel(aes(label = row.names(det)))
p2 

ggsave(paste0('Slc25a1_MaxExpr_Embryo_MainTraj.png'),plot=p2,dpi=300, scale = 1.2)


#+ aes(group=class)
#+ #+ scale_color_manual(values=c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2","#CC79A7"))

p1<-ggplot(det, aes(x=class, y=max.expr)) + geom_dotplot(binaxis='y', stackdir='center', dotsize=1) + aes(group=class)
p1


library(reshape2)
dfm <- melt(de[,c('class','max.expr','status')],id.vars = 1)
p.summary <- ggplot(dfm,aes(x = celltype,y = value)) + 
  geom_bar(aes(fill = variable),stat = "identity",position = "dodge") + theme(axis.text.x = element_text(size = 7, angle = 45, hjust = 1)) +
  xlab("") + ylab("Significantly expressed genes")

ggsave(paste0('Barplot_DE.png'),plot=p.summary,dpi=300, scale = 2)

#+ aes(group=class)

