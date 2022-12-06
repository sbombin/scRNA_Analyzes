
#counts = GetAssayData(pbmc, slot = "data")
#count <- as.data.frame(counts)
#metadata <- pbmc@meta.data
#metadata$cluster_id <- factor(pbmc@active.ident)

cts <- AggregateExpression(pbmc, group.by = c("seurat_clusters", "group"), assays = 'RNA', slot = "data", return.seurat = TRUE)

data =  GetAssayData(cts, slot = "data")

rank.data <- rankGenes(data)

score.notch <-  simpleScore(rank.data, upSet = notch.up, downSet = notch.dwn)
scr.notch <- score.notch %>% dplyr::select( TotalScore, TotalDispersion) %>% dplyr::rename(NOTCH_TotalScore = TotalScore, NOTCH_TotalDispersion = TotalDispersion) %>%
  tibble::rownames_to_column( var = "sample")

score.yap <-  simpleScore(rank.data, upSet = yap.top)
scr.yap <- score.yap %>% dplyr::select( TotalScore, TotalDispersion) %>% dplyr::rename(YAP_TotalScore = TotalScore, YAP_TotalDispersion = TotalDispersion) %>%
  tibble::rownames_to_column( var = "sample")

score.hep <-  simpleScore(rank.data, upSet = hep.up, downSet = hep.dwn)
scr.hep <- score.hep %>% dplyr::select( TotalScore, TotalDispersion) %>% dplyr::rename(HEP_TotalScore = TotalScore, HEP_TotalDispersion = TotalDispersion) %>%
  tibble::rownames_to_column( var = "sample")

score.pc <-  simpleScore(rank.data, upSet = pc.up, downSet = pc.dwn)
scr.pc <- score.pc %>% dplyr::select( TotalScore, TotalDispersion) %>% dplyr::rename(PC_TotalScore = TotalScore, PC_TotalDispersion = TotalDispersion) %>%
  tibble::rownames_to_column( var = "sample")

score.early <-  simpleScore(rank.data, upSet = early.up, downSet = early.dwn)
scr.early <- score.early %>% dplyr::select( TotalScore, TotalDispersion) %>% dplyr::rename(TGFB.ERL_TotalScore = TotalScore, TGFB.ERL_TotalDispersion = TotalDispersion) %>%
  tibble::rownames_to_column( var = "sample")

score.late <-  simpleScore(rank.data, upSet = late.up, downSet = late.dwn)
scr.late <- score.late %>% dplyr::select( TotalScore, TotalDispersion) %>% dplyr::rename(TGFB.Late_TotalScore = TotalScore, TGFB.Late_TotalDispersion = TotalDispersion) %>%
  tibble::rownames_to_column( var = "sample")

scr.all <- left_join(scr.notch,scr.yap, by = "sample", all=TRUE) %>% left_join(.,scr.hep, by = "sample", all=TRUE) %>% left_join(.,scr.pc, by = "sample", all=TRUE) %>%
   left_join(.,scr.early, by = "sample", all=TRUE) %>% left_join(.,scr.late, by = "sample", all=TRUE)
  
write.csv(scr.all,'SignScore_All.csv',row.names = FALSE)
write.csv(score.notch,'SignScore_Notch.csv',row.names = TRUE)
write.csv(score.yap,'SignScore_YAP.csv',row.names = TRUE)
write.csv(score.hep,'SignScore_HEP.csv',row.names = TRUE)
write.csv(score.pc,'SignScore_PC.csv',row.names = TRUE)
write.csv(score.early,'SignScore_TGFB_Early.csv',row.names = TRUE)
write.csv(score.late,'SignScore_TGFB_Late.csv',row.names = TRUE)

library(reshape2)
dfm <- melt(summary.iec.full[,c('celltype','Up','Down')],id.vars = 1)
p.summary <- ggplot(dfm,aes(x = celltype,y = value)) + 
  geom_bar(aes(fill = variable),stat = "identity",position = "dodge") + theme(axis.text.x = element_text(size = 7, angle = 45, hjust = 1)) +
  xlab("") + ylab("Significantly expressed genes")

ggsave(paste0('Barplot_DE.png'),plot=p.summary,dpi=300, scale = 2)



# Rearange table
counts.up <- melt(counts.up[,c('symbol','1262L','1263L')],id.vars = 1)
p1 <- ggplot(counts.up, aes(x=symbol, y=value)) + labs(y= "hmC", x = "gene") + geom_boxplot() + geom_jitter(shape=16, position=position_jitter(0.2))
p1 
