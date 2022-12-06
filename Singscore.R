library(tidyverse)
library(singscore) 
library(ggbeeswarm)
library(tximport)
library(biomaRt)
library(ComplexHeatmap)
library(broom)

library(GSEABase)

setwd("/Users/sbombin/Desktop/analysis_mm39.comb/Seurat_Analysis/BEC/")
data_path <- "/Users/sbombin/Desktop/analysis_mm39.comb/Seurat_Analysis/BEC/"
fig_path <- "/Users/sbombin/Desktop/analysis_mm39.comb/Seurat_Analysis/BEC/"

pbmc <- ReadObject("seuratobj_BEC_cKO.scale_all")
#table(Idents(pbmc))
clst0 <- subset(x = pbmc, idents = 0)


library(irGSEA)

rankData <- rankGenes(pbmc)

my.data=GetAssayData(pbmc, slot = "data")
df <- as.data.frame(my.data)
genes <- unique(row.names(df)) #%>% as.data.frame()

library(readr)
yap <- read_csv("~/Desktop/EICC/Projects/Gracz/Singscore_Data/DongJ_YAPsignature_Cell_2007.csv")
yap <- dplyr::filter(yap, symbol %in% genes)
yap.top <- unique(yap$symbol)[1:600]
#yap.top <- paste(yap.top, "+", sep="")

notch <- read_csv("~/Desktop/EICC/Projects/Gracz/Singscore_Data/VillanuevaA_NotchTargets_Gastro_2012_S7.csv")
notch.up <- notch %>% dplyr::filter(Direction == 'Up'& symbol %in% genes) 
notch.up <- unique(notch.up$symbol)
#notch.up <- paste(notch.up, "+", sep="")

notch.dwn <- notch %>% dplyr::filter(Direction == 'Down' & symbol %in% genes) 
notch.dwn <- unique(notch.dwn$symbol)
#notch.dwn <- paste(notch.dwn, "-", sep="")

subpops <- read_csv('~/Desktop/EICC/Projects/Gracz/Singscore_Data/Schaub2018_S1.csv')
subpops <- subpops %>% dplyr::select(symbol, pC1, pC2, pC3, H1, H2, H3, HpC1, HpC2, HpC3, HpC4) 
subpops <- subpops %>%
  tidyr::gather(sample, score, -symbol)  %>%
  dplyr::mutate(group = case_when(sample %in% c('pC1', 'pC2', 'pC3') ~ 'PC', 
                           sample %in% c('H1', 'H2', 'H3') ~ 'Hep', 
                           sample %in% c('HpC1', 'HpC2', 'HpC3', 'HpC4') ~ 'HpC') ) %>%
  dplyr::group_by(group, symbol) %>%
  dplyr::summarise(u_score  = mean(score, na.rm = T) ) 
subpops <- subpops %>% pivot_wider(names_from = group, values_from = u_score) 

hep.up <- subpops %>% filter(Hep > 1 & symbol %in% genes)  
hep.up <- unique(hep.up$symbol)
#hep.up <- paste(hep.up, "+", sep="")
hep.dwn <- subpops %>% filter(Hep < -1 & symbol %in% genes)
hep.dwn <- unique(hep.dwn$symbol)
#hep.dwn <- paste(hep.dwn, "-", sep="")

pc.up   <- subpops %>% filter(PC > 1 & symbol %in% genes)
pc.up <- unique(pc.up$symbol)
#pc.up <- paste(pc.up, "+", sep="")
pc.dwn <- subpops %>% filter(PC < -1 & symbol %in% genes) 
pc.dwn <- unique(pc.dwn$symbol)
#pc.dwn <- paste(pc.dwn, "-", sep="")

tgfb <- read_csv('~/Desktop/EICC/Projects/Gracz/Singscore_Data/TGFB_GeneList.csv')
early.up <- dplyr::filter(tgfb, condition == 'early' & direction == 'up' &  symbol %in% genes)
early.dwn <- dplyr::filter(tgfb, condition == 'early' & direction == 'down' &  symbol %in% genes)
early.up <- unique(early.up$symbol)
#early.up <- paste(early.up, "+", sep="")
early.dwn <- unique(early.dwn$symbol)
#early.dwn <- paste(early.dwn, "-", sep="")
late.up <- dplyr::filter(tgfb, condition == 'late' & direction == 'up' &  symbol %in% genes)
late.dwn <- dplyr::filter(tgfb, condition == 'late' & direction == 'down' &  symbol %in% genes)
late.up <- unique(late.up$symbol)
#late.up <- paste(late.up, "+", sep="")
late.dwn <- unique(late.dwn$symbol)
#late.dwn <- paste(late.dwn, "-", sep="")

## Combine all markers in List and name list elements
notch.all <- c(notch.up,notch.dwn)
hep.all <- c(hep.up,hep.dwn)
pc.all <- c(pc.up,pc.dwn)
early.all <- c(early.up,early.dwn)
late.all <- c(late.up,late.dwn)

markers <- list(yap.top,notch.all,hep.all,pc.all,early.all,late.all)
names(markers) <- c('YAP','NOTCH','HEP','PC','TGFB_Early','TGFB_Late')

clst0 <- irGSEA.score(object = clst0, assay = "RNA", slot = "data", msigdb = F, custom = TRUE, geneset = markers, method =  "singscore") 

result.dge <- irGSEA.integrate(object = pbmc2, 
                               group.by = "sample",
                               metadata = NULL, col.name = NULL,
                               method = "singscore")

irGSEA.heatmap.plot <- irGSEA.heatmap(object = result.dge, 
                                      method = "RRA",
                                      top = 50, 
                                      show.geneset = NULL)
irGSEA.heatmap.plot

markers <- list(yap.top,notch.up,notch.dwn,hep.up,hep.dwn,pc.up,pc.dwn,early.up,early.dwn,late.up,late.dwn)
names(markers) <- c('yap.top','notch.up','notch.dwn','hep.up','hep.dwn','pc.up','pc.dwn','early.up','early.dwn','late.up','late.dwn')

markers = lapply(markers, unique)


## Convert Seurat to SingleCellExperiment
library(scater)

sce <- as.SingleCellExperiment(pbmc)


library(SeuratDisk)
#pbmc.loom <- as.loom(pbmc, filename = "BEC_cKO.loom", verbose = FALSE)
#pbmc.loom$close_all()
#se <- makeSummarizedExperimentFromLoom(BEC_cKO.loom, rownames_attr = NULL,colnames_attr = NULL)
