library(Seurat)
library(singleCellNet)

setwd("/Users/sbombin/Desktop/analysis_mm39.comb/Seurat_Analysis/BDL/")
data_path <- "/Users/sbombin/Desktop/analysis_mm39.comb/Seurat_Analysis/BDL/"
fig_path <- "/Users/sbombin/Desktop/analysis_mm39.comb/Seurat_Analysis/BDL/"

pbmc <- ReadObject("seuratobj.Sham-BDL.scale_all")
table(pbmc@active.ident)
## Separate Sham (0,3,4,7) and BDL (1,2,5,6,9,10)  Clusters
table(pbmc@meta.data$seurat_clusters) ## 693 cells

sham <- subset(x = pbmc, idents = c(0,3,4,7))
bdl <- subset(x = pbmc, idents = c(1,2,5,6,9,10))
table(bdl@active.ident)


seuratfile <- extractSeurat(sham, exp_slot_name = "counts")

## Training expression matrix
expTMraw = seuratfile$expDat
## Training Metadata file
stTM = seuratfile$sampTab
stTM$cell <- row.names(stTM)

stTM<-droplevels(stTM)


## Split for training and assessment, and transform training data
set.seed(100) #can be any random seed number
stList = splitCommon(sampTab=stTM, ncells=150, dLevel="seurat_clusters")
stTrain = stList[[1]]
expTrain = expTMraw[,rownames(stTrain)]
## Train
system.time(class_info<-scn_train(stTrain = stTrain, expTrain = expTrain, nTopGenes = 50, nRand = 70, nTrees = 5000,
                                  nTopGenePairs = 100, dLevel = "seurat_clusters", colName_samp = "cell"))


#validate data
stTestList = splitCommon(sampTab=stList[[2]], ncells=150, dLevel="seurat_clusters") #normalize validation data so that the assessment is as fair as possible
stTest = stTestList[[1]]
#expTest = expTMraw[commonGenes,rownames(stTest)]
expTest = expTMraw[,rownames(stTest)]

#predict
classRes_val_all = scn_predict(cnProc=class_info[['cnProc']], expDat=expTest, nrand = 50)

system.time(classRes_val_all <- scn_predict(class_info[['cnProc']], expTest, nrand = 50))

#test = classRes_val_all[1:4,]

#stTrain$cell <- row.names(stTrain)
#stTest$cell <- row.names(stTest)

  
tm_heldoutassessment = assess_comm(ct_scores = classRes_val_all, stTrain = stTrain, stQuery = stTest,
                                   dLevelSID = "cell", classTrain = "seurat_clusters", classQuery = "seurat_clusters", nRand = 50)

p <- plot_PRs(tm_heldoutassessment)
SaveFigure(p, "SCN_Train_Classifier", width = 12, height = 10)

### Classify Data

seuratfile2 <- extractSeurat(bdl, exp_slot_name = "counts")

## expression matrix
expQuery = seuratfile2$expDat
## Training Metadata file
stQuery = seuratfile2$sampTab
stQuery$cell <- row.names(stQuery)
## Classify dataset
system.time(crPBMC <- scn_predict(class_info[['cnProc']], expQuery, nrand = 50))

## Add annotation to metadata (method 1)
#z <- crPBMC[,1:5808]
#stQuery2 <- assign_cate(classRes = z, sampTab = stQuery, cThresh = 0.5) 
#test <- t(crPBMC) %>% as.data.frame()
#test$cell <-row.names(test)

## Add annotation to metadata (method 2)
stPark <- get_cate(classRes = crPBMC, sampTab = stQuery, dLevel = "seurat_clusters", sid = "cell", nrand = 50)

## Extract SCN results
result <- stPark[,18:21]

## Rename Clusters
result$seurat_clusters <- gsub("^1$", "BDL1", result$seurat_clusters)
result$seurat_clusters <- gsub("^2$", "BDL2", result$seurat_clusters)
result$seurat_clusters <- gsub("^5$", "BDL3", result$seurat_clusters)
result$seurat_clusters <- gsub("^6$", "BDL4", result$seurat_clusters)
result$seurat_clusters <- gsub("^9$", "BDL5", result$seurat_clusters)
result$seurat_clusters <- gsub("^10$", "BDL6", result$seurat_clusters)

result$category <- gsub("^0$", "Sham1", result$category)
result$category <- gsub("^3$", "Sham2", result$category)
result$category <- gsub("^4$", "Sham3", result$category)
result$category <- gsub("^7$", "Sham4", result$category)

write.csv(result,"SCN_Results_Table.csv", row.names = TRUE)

#p1 <- plot_attr(crPBMC, stQuery, nrand=50, dLevel="seurat_clusters", sid="cell")
#p1

library(dplyr)
sum <- result%>% group_by(seurat_clusters,category) %>% count()

write.csv(sum,"SCN_Summary_Table.csv", row.names = FALSE)

colors <- c("#CC6666", "#C7E2AE", "#7AD9DA", "#DC7850","#A5947E", "#D4DE5A", "#7B76D8", "#DFA8BE", "#006666",
            "#72DF9E","#80BDD7", "#D08DDB", "#8F4BE3", "#D8568A", "#D6B270", "#7FA081")
p3 <- ggplot(sum, aes(fill=category, y=n, x=seurat_clusters)) + 
  geom_bar(position="fill", stat="identity") + scale_fill_manual(values = colors)
p3
SaveFigure(p3, "SCN_Summary_BarPlot", width = 8, height = 6)
