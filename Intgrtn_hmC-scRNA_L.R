library(ggpubr)
library(dplyr)
library(tibble)
library(readr)
library(ggplot2)
library(reshape2)

setwd('/Users/sbombin/Desktop/EICC/Projects/Gracz/Bulk-seq/')

results.de <- read_delim("/Users/sbombin/Library/CloudStorage/OneDrive-EmoryUniversity/Gracz_Adam/Seurat_out/IEC/Expression_plots_tables/DE-MAST_Tet1iKO_vs_Control_Full.csv",
                         delim = ",", escape_double = FALSE, trim_ws = TRUE) %>% as.data.frame()

counts.hmc <- read_delim("hmC-seq_2021/L_Mean_Normalized_Counts.csv", delim = ",", 
                         escape_double = FALSE, trim_ws = TRUE) %>% as.data.frame()
## Rename colums
colnames(counts.hmc) <- c('symbol','1055L','1069L','1262L', '1263L')

de.up <- dplyr::filter(results.de, p_val_adj < 0.05 & avg_log2FC >  0.5 & cluster_id == "8")
de.dwn <- dplyr::filter(results.de, p_val_adj < 0.05 & avg_log2FC <  -0.5 & cluster_id == "8")

body.up <- de.up[,"gene"]
body.dwn<- de.dwn[,"gene"]

### KO

## Counts Up
counts.up <- counts.hmc[c("symbol","1262L","1263L")] %>% subset(., subset = symbol %in% body.up) %>% na.omit() 
# Rearange table
counts.up <- melt(counts.up[,c('symbol','1262L','1263L')],id.vars = 1)
write.csv(counts.up,paste0('L_RNA-sign_up_hmC-GBody_KO.csv'),row.names = FALSE)

# Plot
#p <- ggplot(counts.up, aes(x=symbol, y=value)) + geom_boxplot() +  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5)
p1 <- ggplot(counts.up, aes(x=symbol, y=value)) + labs(y= "hmC", x = "gene") + geom_boxplot() + geom_jitter(shape=16, position=position_jitter(0.2))
p1 
ggsave(paste0('BoxPlot_L_scRNA-sign_up_hmC-GBody_KO', ".png"),plot=p1,dpi=300, scale = 2)

## Counts Down
counts.dwn <- counts.hmc[c("symbol","1262L","1263L")] %>% subset(., subset = symbol %in% body.dwn) %>% na.omit() 
# Rearange table
counts.dwn <- melt(counts.dwn[,c('symbol','1262L','1263L')],id.vars = 1)
write.csv(counts.dwn,paste0('L_RNA-sign_dwn_hmC-GBody_KO.csv'),row.names = FALSE)

# Plot
p2 <- ggplot(counts.dwn, aes(x=symbol, y=value)) + labs(y= "hmC", x = "gene") + geom_boxplot() + geom_jitter(shape=16, position=position_jitter(0.2))
p2 
ggsave(paste0('BoxPlot_L_scRNA-sign_dwn_hmC-GBody_KO', ".png"),plot=p2,dpi=300, scale = 2)

### Control

## Counts Up
counts.up <- counts.hmc[c("symbol","1055L","1069L")] %>% subset(., subset = symbol %in% body.up) %>% na.omit() 
# Rearange table
counts.up <- melt(counts.up[,c("symbol","1055L","1069L")],id.vars = 1)
write.csv(counts.up,paste0('L_RNA-sign_up_hmC-GBody_CTRL.csv'),row.names = FALSE)
## Plot
p1 <- ggplot(counts.up, aes(x=symbol, y=value)) + labs(y= "hmC", x = "gene") + geom_boxplot() + geom_jitter(shape=16, position=position_jitter(0.2))
p1 
ggsave(paste0('BoxPlot_L_scRNA-sign_up_hmC-GBody_CTRL', ".png"),plot=p1,dpi=300, scale = 2)

## Counts Down
counts.dwn <- counts.hmc[c("symbol","1055L","1069L")] %>% subset(., subset = symbol %in% body.dwn) %>% na.omit() 
# Rearange table
counts.dwn <- melt(counts.dwn[,c("symbol","1055L","1069L")],id.vars = 1)
write.csv(counts.dwn,paste0('L_RNA-sign_dwn_hmC-GBody_CTRL.csv'),row.names = FALSE)
# Plot
p2 <- ggplot(counts.dwn, aes(x=symbol, y=value)) + labs(y= "hmC", x = "gene") + geom_boxplot() + geom_jitter(shape=16, position=position_jitter(0.2))
p2 
ggsave(paste0('BoxPlot_L_scRNA-sign_dwn_hmC-GBody_CTRL', ".png"),plot=p2,dpi=300, scale = 2)


