###Volcano Plot
library(ggplot2)
library(ggrepel)
library(org.Mm.eg.db)
library(AnnotationDbi)
library(stringr)

setwd('/Users/indrikwijaya/Downloads/Subgroups 2/')
comp1 <- read.table('HGNC_complex I.txt',sep = '\t',header= T)
comp2 <- read.table('HGNC_complex II.txt',sep = '\t',header = T)
comp3 <- read.table('HGNC_complex III.txt',sep = '\t',header = T)
comp4 <- read.table('HGNC_complex IV.txt', sep = '\t',header = T)
comp5 <- read.table('HGNC_complex V.txt', sep = '\t',header = T)
comp_list <- list(comp1,comp2,comp3,comp4,comp5)

plot_volcano_complex <- function(res_tableDE, type, day,list_of_complexes,folder_label,xlims,ylims){
  res_tableDE <- data.frame(res_tableDE)
  res_tableDE$SYMBOL <- mapIds(org.Mm.eg.db,
                               keys = as.character(rownames(res_tableDE)),
                               column = 'SYMBOL',
                               keytype ='ACCNUM',
                               multiVals = 'first')
  threshold_DE <- res_tableDE$padj < 0.1
  
  res_tableDE$threshold <- threshold_DE
  #res_tableDE$threshold[res_tableDE$SYMBOL %in% list_of_genes] <- 'highlight'
  
  res_tableDE$genelabels <- ""
  res_tableDE$genelabels[res_tableDE$SYMBOL %in% list_of_complexes[[1]]$Approved.Symbol] <- 'CompI'
  res_tableDE$genelabels[res_tableDE$SYMBOL %in% list_of_complexes[[2]]$Approved.Symbol] <- 'CompII'
  res_tableDE$genelabels[res_tableDE$SYMBOL %in% list_of_complexes[[3]]$Approved.Symbol] <- 'CompIII'
  res_tableDE$genelabels[res_tableDE$SYMBOL %in% list_of_complexes[[4]]$Approved.Symbol] <- 'CompIV'
  res_tableDE$genelabels[res_tableDE$SYMBOL %in% list_of_complexes[[5]]$Approved.Symbol] <- 'CompV'
  comp_all <- rbind(list_of_complexes[[1]],list_of_complexes[[2]],list_of_complexes[[3]],
                        list_of_complexes[[4]],list_of_complexes[[5]])$Approved.Symbol
  
  ggplot(res_tableDE,aes(x = log2FoldChange, y = -log10(padj))) +
    
    # rna axes limit
    # ylim(-0.05,100) +
    # xlim(-5,5) +
    # 
    #te axes limit
    # ylim(-0.05,18) +
    # xlim(-3,3) +
    
    #rna enrich axes limit
    # ylim(-0.05,200) +
    # xlim(-3,6) +
    
    #rpf enrich axes limit
    # ylim(-0.05,300) +
    # xlim(-3,10) +
    
    #te enrich axes limit
    # ylim(-0.05,50) +
    # xlim(-6,6.5) +
    ylim(ylims[[1]],ylims[[2]]) +
    xlim(xlims[[1]],xlims[[2]]) +
    
    #threshold lines
    geom_vline(xintercept = -1, linetype = 'dashed', color = 'brown',size = 0.7) +
    geom_vline(xintercept = 1, linetype = 'dashed', color = 'brown', size = 0.7) +
    geom_hline(yintercept = 1, linetype = 'dashed', color = 'brown', size = 0.7) +
    
    #points
    geom_point(aes(colour = threshold),alpha = 0.5) +
    geom_point(data = subset(res_tableDE, SYMBOL %in% comp_all),
               aes(colour = genelabels,shape = genelabels),alpha = 1,size = 2)+
    # geom_text_repel(aes(x = log2FoldChange, y = -log10(padj), 
    #                     label = ifelse(genelabels==TRUE,res_tableDE$SYMBOL,""),size = 3),
    #                 box.padding = unit(0.2,'lines'),
    #                 point.padding = unit(0.2, 'lines'))+
    
    
    scale_colour_manual(breaks = c('CompI','CompII','CompIII','CompIV','CompV'),
                        labels = c('CompI','CompII','CompIII','CompIV','CompV'),
                      values = c("#0072B2", "#D55E00","#E69F00", "grey3", "coral4",'#F8766D','#00BFC4'))+
    labs(colour = 'test',shape ='test')+
    guides(colour = guide_legend(ncol = 2))+

    geom_text_repel(data = subset(res_tableDE, SYMBOL %in% comp_all), #& log2FoldChange > 1),
                    aes(label = SYMBOL), size = 3,
                    box.padding = unit(0.4, 'lines'),
                    point.padding = unit(0.4, 'lines'),
                    segment.size = 0.2, segment.colour = 'grey50')+
    
    theme_bw() + 
    
    ggtitle(paste(type,day,sep = '')) +
    xlab(bquote(~log[2]~ "FC")) +
    ylab(bquote(~-log[10]~italic(p-adj)))+
    
    theme(legend.position = c(0.15,0.92),
          legend.title = element_blank(),
          legend.background = element_rect(color = 'black',size = 0.5,linetype= 'solid'),
          legend.text = element_text(size = 12),
          plot.title = element_text(size = rel(1.5),hjust = 0.5),
          axis.title = element_text(size = rel(1.25)),
          axis.text.y = element_text(size = 12),
          axis.text.x = element_text(size=14))
  #save_dir <- '/Users/indrikwijaya/Desktop/de_analysis/enrichment/notnorm/no_er_genes/volcano_plots/'
  save_dir <- '/Users/indrikwijaya/Desktop/de_analysis/alpha_0.05/refseq_updated_shifted/'
  
  #check whether folder exists
  save_folder <- paste(save_dir, folder_label, sep = '') 
  ifelse(!dir.exists(save_folder), dir.create(save_folder), FALSE)
  
  filepath <- paste(save_folder, '/', type, toString(day), '.png',sep = '')
  ggsave(filepath,width = 9, height = 7)
}
