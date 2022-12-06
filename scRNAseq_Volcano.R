library(EnhancedVolcano)

RNAseq_Volcano <- function(results, pcut, fcut, title){
  
  # Remove NA q-values
  volres <- results[!is.na(results$p_val_adj),]
  
  p <- EnhancedVolcano(volres,
                       lab = rownames(volres),
                       x = 'avg_log2FC',
                       y = 'p_val_adj',
                       title = title,
                       legendLabels=c('Not sig.','Log (base 2) FC','q-value',
                                      'q-value & Log (base 2) FC'),
                       titleLabSize = 12,
                       subtitleLabSize = 9,
                       axisLabSize = 12,
                       legendLabSize = 10,
                       legendIconSize = 3.0,
                       drawConnectors = TRUE,
                       labSize = 3.0,
                       pCutoff = pcut,
                       FCcutoff = fcut,
                       xlim = c(min(volres$avg_log2FC,na.rm=TRUE)-1,
                                max(volres$avg_log2FC,na.rm=TRUE)+1),
                       ylim = c(0,max(-log10(volres$p_val_adj)+1, na.rm = TRUE) ))
  
  ggsave(paste0(title,'_VolcanoPlot.png'),plot=p,dpi=300, scale = 2)
}

RNAseq_Volcano(results=msc1.vs.msc2, pcut = 10e-32, fcut= 0.05, title="MSC1.vs.MSC2")



RNAseq_Volcano(results=msc2.vs.msc3, pcut = 0.05, fcut= 0.05, title="MSC2.vs.MSC3")


RNAseq_Volcano(results=msc4.vs.msc5, pcut = 0.05, fcut= 0.05, title="MSC4.vs.MSC5")

res =msc4.vs.msc5
p<-EnhancedVolcano(res,
                   lab = rownames(res),
                   x = 'avg_log2FC',
                   y = 'p_val_adj',
                   title = 'MSC4.vs.MSC5',
                   legendLabels=c('Not sig.','Log (base 2) FC','q-value',
                                  'q-value & Log (base 2) FC'),
                   drawConnectors = TRUE,
                   pCutoff = 10e-32,
                   FCcutoff = 0.05,
                   pointSize = 3.0,
                   axisLabSize = 12,
                   legendLabSize = 10,
                   legendIconSize = 3.0,
                   labSize = 3.0)
ggsave(paste0('MSC4.vs.MSC5_VolcanoPlot.png'),plot=p,dpi=300, scale = 2)

