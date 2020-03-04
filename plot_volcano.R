###Volcano Plot
library(ggplot2)
library(ggrepel)
library(org.Mm.eg.db)
library(AnnotationDbi)

###res_tableDE: object from running results in DESeq
###type: rnabulk/rnacyt/..etc 
###list_of_genes: particular group of genes to be highlighted, can leave it as c() if none highlighted
###label: label of particular group of genes to be highlighted, can leave it as '' if none highlighted
###folder_label: name of folder to save plots
###xlims: c(x_min, x_max); ylims: c(y_min, y_max)
alpha <- 0.05
plot_volcano <- function(res_tableDE, type, day, 
                         list_of_genes, label, folder_label, 
                         xlims, ylims){
  
  res_tableDE <- data.frame(res_tableDE)
  res_tableDE$SYMBOL <- mapIds(org.Mm.eg.db,
                               keys = as.character(rownames(res_tableDE)),
                               column = 'SYMBOL',
                               keytype ='ACCNUM',
                               multiVals = 'first')

  threshold_DE <- res_tableDE$padj < alpha
  res_tableDE$threshold <- threshold_DE
  #res_tableDE$threshold[res_tableDE$SYMBOL %in% list_of_genes] <- 'highlight'
  
  res_tableDE$genelabels <- ""
  res_tableDE$genelabels[res_tableDE$SYMBOL %in% list_of_genes] <- label
  
  volcano_plot <- ggplot(res_tableDE,aes(x = log2FoldChange, y = -log10(padj))) +
    
    ylim(ylims[[1]],ylims[[2]]) +
    xlim(xlims[[1]], xlims[[2]]) +
    
    #threshold lines
    geom_vline(xintercept = -1, linetype = 'dashed', color = 'brown', size = 1.2) +
    geom_vline(xintercept = 1, linetype = 'dashed', color = 'brown', size = 1.2) +
    geom_hline(yintercept = -log10(0.05), linetype = 'dashed', color = 'brown', 
               size = 1.2) +
    
    #color points based on lfc and padj threshold
    geom_point(aes(colour = threshold), alpha = 0.5) +
    
    #color points in subgroup
    geom_point(data = subset(res_tableDE, SYMBOL %in% list_of_genes),
               aes(fill = genelabels), color = 'blue')+

    scale_fill_manual(name = '',
                       values = c(label = 'blue'))+
    
    #remove legend for colour
    guides(colour = F)
    
  
    #add gene labels for particular group from list_of_genes if no of genes <= 30
    subgroup_genes <- subset(res_tableDE, SYMBOL %in% list_of_genes)
    if(nrow(subgroup_genes) <= 30){
    volcano_plot <- volcano_plot + geom_text_repel(data = subgroup_genes,# & log2FoldChange > 1),
                    #data = subset(res_tableDE, SYMBOL %in% c('Htra2','Tomm7','Sqstm1')), #for rna autophagy
                    aes(label = SYMBOL), size = 3,
                    box.padding = unit(0.4, 'lines'),
                    point.padding = unit(0.4, 'lines'),
                    segment.size = 0.2, segment.colour = 'grey50')
      
    }
    
    volcano_plot + 
      theme_bw() + 
    
    ggtitle(paste(type,day, sep = '')) +
    xlab(bquote(~log[2]~ "FC")) +
    ylab(bquote(~-log[10]~italic(p-adj)))+
    
    theme(#legend.position = c(0.15,0.92),
          legend.position = 'none',
          #legend.title = element_blank(),
          #legend.background = element_rect(color = 'black', size = 0.5, linetype= 'solid'),
          #legend.text = element_text(size = 12),
          #plot.title = element_text(size = rel(1.5), hjust = 0.5),
          plot.title = element_blank(),
          #axis.title = element_text(size = 20),
          axis.title = element_blank(),
          #axis.text.y = element_text(size = 16),
          axis.text.y = element_blank(),
          axis.text.x = element_text(size = 42)
          )

  #save_dir <- '/Users/indrikwijaya/Desktop/de_analysis/alpha_0.05/enrichment_crubulk/notnorm/er_genes/volcano_plots/'
  save_dir <- '/Users/indrikwijaya/Desktop/de_analysis/alpha_0.05/volcano_plots_refseq_updated_shifted/'
  
  #check whether folder exists
  save_folder <- paste(save_dir, folder_label, sep='') 
  ifelse(!dir.exists(save_folder), dir.create(save_folder), FALSE)
  
  filepath <- paste(save_folder, '/', type, toString(day), '.png', sep='')
  ggsave(filepath, width=9, height=7)
}
