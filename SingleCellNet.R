#source("Fxns.R")
## Downloading UMI count data
download.file("ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE92nnn/GSE92332/suppl/GSE92332_atlas_UMIcounts.txt.gz", destfile="GSE92332_atlas_UMIcounts.txt.gz")
## Reading UMI count data from file
atlas_umis = read.delim("GSE92332_atlas_UMIcounts.txt.gz")
info(sprintf("Data dimensions: %s" , paste(dim(atlas_umis), collapse = "x")))

#umi = read.delim("/Users/sbombin/Downloads//GSE92332_Org_RANKL_UMIcounts.txt")
library(singleCellNet)
pbmc <- ReadObject("seuratobj.EIC.scale_all")

seuratfile <- extractSeurat(pbmc, exp_slot_name = "counts")
sampTab = seuratfile$sampTab
expDat = seuratfile$expDat

#stTM <- utils_loadObject(fname = "/Users/sbombin/Downloads/sampTab_TM_053018.rda")
#expTm<- utils_loadObject(fname = "/Users/sbombin/Downloads/expTM_Raw_053018.rda")

#stList<-splitCommon(sampTab = stTM, ncells = 50, dLevel = "newAnn")

a <- as.data.frame(t(atlas_umis))
library(tidyverse)
library(tidyr)

b <- rownames_to_column(a, var = "rowname")
c <- b$rowname
d <- as.data.frame(c)
e <- separate(d, c, into = c("cell", "newAnn"), sep = "(?<=[A-Z])_")

f <- e %>% mutate(name = cell)
g <- column_to_rownames(f, var = "name")


#w <- gsub("(_[^_]+)_.*", "\\1", colnames(atlas_umis))
colnames(atlas_umis) <- stringr::str_extract(colnames(atlas_umis), "[^_]*_[^_]*")

#commonGenes<-intersect(rownames(atlas_umis), rownames(expTm))
#expTrain <- atlas_umis[commonGenes, ]
#expQuery <- expTm[commonGenes, ]

#stList<-splitCommon(sampTab = g, ncells = 50, dLevel = "newAnn")

#stTrain<-stList[[1]]
#expTrain <- expTrain[,stTrain$cell]

#stTest <- stList[[2]]
#expTest <- expTrainOrth[,stTest$cell]

stPark = sampTab
expPark = expDat
genesPark = rownames(expDat)

expTMraw = atlas_umis
stTM<-droplevels(g)

commonGenes = intersect(rownames(expTMraw), genesPark)
length(commonGenes)
expTMraw = expTMraw[commonGenes,]
## Split for training and assessment, and transform training data
set.seed(11) #can be any random seed number
stList = splitCommon(sampTab=stTM, ncells=100, dLevel="newAnn")
stTrain = stList[[1]]
expTrain = as.matrix(expTMraw[,rownames(stTrain)])
## Train the classifier
system.time(class_info<-scn_train(stTrain = stTrain, expTrain = expTrain, nTopGenes = 20, nRand = 70, nTrees = 1000,
                                  nTopGenePairs = 50, dLevel = "newAnn", colName_samp = "cell"))



#validate data
stTestList = splitCommon(sampTab=stList[[2]], ncells=100, dLevel="newAnn") #normalize validation data so that the assessment is as fair as possible
stTest = stTestList[[1]]
expTest = expTMraw[commonGenes,rownames(stTest)]

#predict
classRes_val_all = scn_predict(cnProc=class_info[['cnProc']], expDat=expTest, nrand = 50)

tm_heldoutassessment = assess_comm(ct_scores = classRes_val_all, stTrain = stTrain, stQuery = stTest,
                                   dLevelSID = "cell", classTrain = "newAnn", classQuery = "newAnn", nRand = 50)

p1 <- plot_PRs(tm_heldoutassessment)
p2 <- plot_metrics(tm_heldoutassessment)

#Create a name vector label used later in classification heatmap where the values are cell types/ clusters and names are the sample names
nrand = 50
sla = as.vector(stTest$newAnn)
names(sla) = as.vector(stTest$cell)
slaRand = rep("rand", nrand) 
names(slaRand) = paste("rand_", 1:nrand, sep='')
sla = append(sla, slaRand) #include in the random cells profile created

sc_hmClass(classMat = classRes_val_all,grps = sla, max=300, isBig=TRUE)
## Viusalize average top pairs genes expression for training data
gpTab = compareGenePairs(query_exp = expTest, training_exp = expTrain, training_st = stTrain, classCol = "newAnn", sampleCol = "cell",
                         RF_classifier = class_info$cnProc$classifier, numPairs = 20, trainingOnly= TRUE)

train = findAvgLabel(gpTab = gpTab, stTrain = stTrain, dLevel = "newAnn")

hm_gpa_sel(gpTab, genes = class_info$cnProc$xpairs, grps = train, maxPerGrp = 50)
## Apply to  query data
nqRand = 50
system.time(crParkall<-scn_predict(class_info[['cnProc']], expPark, nrand=nqRand))

sgrp = as.vector(stPark$description1)
names(sgrp) = as.vector(stPark$sample_name)
grpRand =rep("rand", nqRand)
names(grpRand) = paste("rand_", 1:nqRand, sep='')
sgrp = append(sgrp, grpRand)

# heatmap classification result
sc_hmClass(crParkall, sgrp, max=5000, isBig=TRUE, cCol=F, font=8)

# This classifies a cell with  the catgory with the highest classification score or higher than a classification score threshold 
# The annotation result can be found in a column named category in the query sample table.

z <- crPBMC[,1:15060]
stQuery <- assign_cate(classRes = z, sampTab = stPark, cThresh = 0.5)

final <- stQuery[c("orig.ident","sample", "group", "seurat_clusters", "category")]
write.csv(final,"/Users/sbombin/Desktop/analysis_mm39.comb/Seurat_Analysis/EIC//SCN_EIC_cell_assignment.csv", row.names = TRUE)

SaveFigure(p1, "SCN_Training", width = 12, height = 8)


df <- stQuery[c( "seurat_clusters", "category")]
sum <- df%>% group_by(seurat_clusters,category) %>% count()
#sum2 <- df%>% group_by(seurat_clusters,category) %>% summarise(count = n())
write.csv(sum,"/Users/sbombin/Desktop/analysis_mm39.comb/Seurat_Analysis/EIC//Cell-Cluster_Count.csv", row.names = FALSE)

library(viridis)
library(RColorBrewer)
library(randomcoloR)
nb.cols <- 16
#mycolors <- colorRampPalette(brewer.pal(12, "Set3"))(nb.cols)
palette <- distinctColorPalette(nb.cols)
manualcolors<-c('forestgreen', 'red2', 'orange', 'cornflowerblue', 
                 'darkolivegreen4', 'indianred1', 'tan4', 'darkblue', 
                'mediumorchid1','firebrick4',  'yellowgreen', 'lightsalmon', 'tan3',
                 'wheat4', '#DDAD4B',  
                'seagreen1', 'moccasin', 'mediumvioletred', 'seagreen','cadetblue1')

p2 <- ggplot(sum, aes(fill=category, y=n, x=seurat_clusters)) + 
  geom_bar(position="fill", stat="identity") + scale_fill_manual(values = mycolors)

#p3 <- ggplot(sum, aes(fill=category, y=n, x=seurat_clusters)) +  geom_bar(position="fill", stat="identity") + scale_color_viridis(discrete=TRUE, option="plasma")

SaveFigure(p2, "Cell-Cluster_Count2", width = 12, height = 8)

final <- final %>% mutate(cluster = seurat_clusters)

final$cluster <- gsub("^0$", "Enterocyte Progenitor (late) 1", final$cluster)
final$cluster <- gsub("^1$", "Mature Enterocyte (proximal) 1", final$cluster)
final$cluster <- gsub("^2$", "Mature Enterocyte (proximal) 2", final$cluster)
final$cluster <- gsub("^3$", "TA", final$cluster)
final$cluster <- gsub("^4$", "Mature Enterocyte (distal)", final$cluster)
final$cluster <- gsub("^5$", "Immature Enterocyte (proximal)", final$cluster)
final$cluster <- gsub("^6$", "Mature Enterocyte (proximal) 3", final$cluster)
final$cluster <- gsub("^7$", "Enterocyte Progenitor (early) 1", final$cluster)
final$cluster <- gsub("^8$", "ISC", final$cluster)
final$cluster <- gsub("^9$", "Enterocyte Progenitor (early) 2", final$cluster)
final$cluster <- gsub("^10$", "Enterocyte Progenitor (late) 2", final$cluster)
final$cluster <- gsub("^11$", "Goblet", final$cluster)
final$cluster <- gsub("^12$", "Tuft", final$cluster)
final$cluster <- gsub("^13$", "Enteroendocrine", final$cluster)
final$cluster <- gsub("^14$", "Mature Enterocyte (proximal) 4", final$cluster)
final$cluster <- gsub("^15$", "Paneth", final$cluster)

df <- final[c( "sample", "cluster")]
sum <- df%>% group_by(sample,cluster) %>% count()
write.csv(sum,"/Users/sbombin/Desktop/analysis_mm39.comb/Seurat_Analysis/IEC//Cell_Count_perMouse.csv", row.names = FALSE)
colors <- c("#CC6666", "#C7E2AE", "#7AD9DA", "#DC7850","#A5947E", "#D4DE5A", "#7B76D8", "#DFA8BE", "#006666",
            "#72DF9E","#80BDD7", "#D08DDB", "#8F4BE3", "#D8568A", "#D6B270", "#7FA081")
p3 <- ggplot(sum, aes(fill=cluster, y=n, x=sample)) + 
  geom_bar(position="fill", stat="identity") + scale_fill_manual(values = colors)
SaveFigure(p3, "Cell_Count_perMouse", width = 12, height = 8)

df <- final[c( "group", "cluster")]
sum <- df%>% group_by(group,cluster) %>% count()
write.csv(sum,"/Users/sbombin/Desktop/analysis_mm39.comb/Seurat_Analysis/IEC//Cell_Count_perGroup.csv", row.names = FALSE)
p4 <- ggplot(sum, aes(fill=cluster, y=n, x=group)) + 
  geom_bar(position="fill", stat="identity") + scale_fill_manual(values = colors)
SaveFigure(p4, "Cell_Count_perGroup", width = 8, height = 8)
