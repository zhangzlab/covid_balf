total_fig <- function(){
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(ggplot2)
  library(reshape2)
  setwd('/home/data/results/workspace/nCoV_rev/balf/Figure')
  nCoV.integrated = readRDS(file = "../nCoV.rds")
  
  dpi = 300
  png(file="1-umap_1.png", width = dpi*8, height = dpi*6, units = "px",res = dpi,type='cairo')
  pp = DimPlot(object = nCoV.integrated, reduction = 'umap',pt.size = 0.5,label = FALSE,label.size = 5)
  pp = pp + theme(axis.title = element_text(size = 15),axis.text =  element_text(size = 15,family = 'sans'),
                  legend.text = element_text(size = 15),axis.line = element_line(size = 0.8))
  print(pp)
  dev.off()
  
  ###umap of marker
  markers = c('TPPP3','KRT18','CD68','FCGR3B','CD1C','CLEC9A','LILRA4','TPSB2','CD3D','KLRD1','MS4A1','IGHG4')
  png(file="1-umap_marker_1.png", width = dpi*16, height = dpi*10, units = "px",res = dpi,type='cairo')
  pp_temp = FeaturePlot(object = nCoV.integrated, features = markers,cols = c("lightgrey","#ff0000"),combine = FALSE)
  plots <- lapply(X = pp_temp, FUN = function(p) p + labs(title='') + theme(axis.title = element_text(size = 22),axis.text = element_text(size = 22),
                                                                            plot.title = element_text(family = 'sans',face='italic',size=22),axis.line = element_line(size = 1.5),axis.ticks = element_line(size = 1.2),
                                                                            legend.text = element_text(size = 22),legend.key.height = unit(1.8,"line"),legend.key.width = unit(1.2,"line")))
  pp = CombinePlots(plots = plots,ncol = 4,legend = 'right')
  print(pp)
  dev.off()
  
  pdf(file="1-umap_marker.pdf", width = 16, height = 10)
  pp_temp = FeaturePlot(object = nCoV.integrated, features = markers,cols = c("lightgrey","#ff0000"),combine = FALSE)
  plots <- lapply(X = pp_temp, FUN = function(p) p + theme(axis.title = element_text(size = 18),axis.text = element_text(size = 18),
                                                           plot.title = element_text(family = 'sans',face='italic',size=20),
                                                           legend.text = element_text(size = 20),legend.key.height = unit(1.8,"line"),legend.key.width = unit(1.2,"line") ))
  pp = CombinePlots(plots = plots,ncol = 4,legend = 'right')
  print(pp)
  dev.off()
  
  pdf(file="1-marker_heatmap.pdf", width = 8, height = 7)
  pp = DotPlot(nCoV.integrated, features = rev(markers),cols = c('white','#F8766D'),dot.scale =5) + RotatedAxis()
  pp = pp + theme(axis.text.x = element_text(size = 12),axis.text.y = element_text(size = 12)) + labs(x='',y='') + 
    guides(color = guide_colorbar(title = 'Scale expression'),size = guide_legend(title = 'Percent expressed')) + 
    theme(axis.line = element_line(size = 0.6))
  print(pp)
  dev.off()
  
  ####绘制比例图
  library(ggplot2)
  nCoV.integrated[["cluster"]] <- Idents(object = nCoV.integrated)
  big.cluster = nCoV.integrated@meta.data
  organ.summary = table(big.cluster$sample_new,big.cluster$cluster)
  write.table(organ.summary,file = '1-nCoV-percentage-sample.txt',quote = FALSE,sep = '\t')
  
  library(ggplot2)
  library(dplyr)
  library(reshape2)
  setwd('/Users/yang/Downloads/work/nCoV/scrna_rev/balf/Figure/')
  organ.summary = read.delim2("1-nCoV-percentage-sample.txt",header = TRUE, stringsAsFactors = FALSE,check.names = FALSE, sep = "\t")
  organ.summary$group = rownames(organ.summary)
  organ.summary.dataframe = melt(organ.summary)
  colnames(organ.summary.dataframe) = c('group','cluster','cell')
  samples_name_new = c('HC1','HC2','HC3','HC4','O1','O2','O3','S1','C1','C2','C3','C4','C5')
  organ.summary.dataframe$group = factor(organ.summary.dataframe$group,labels = samples_name_new,levels = samples_name_new)
  organ.summary.dataframe$cell = as.numeric(organ.summary.dataframe$cell)
  
  ##Macrophage neutrophil plasma epithelial T NK mDC B pDC MAST
  new_order = c('0','1','2','3','4','5','7','8','10','11','12','18','21','22','23','26','15','19','13','16','25','28','31','6','9','14','17','20','27','29','30','24')
  organ.summary.dataframe$cluster = factor(organ.summary.dataframe$cluster,levels = new_order,labels = new_order)
  dpi = 300
  cols = c('#32b8ec','#60c3f0','#8ccdf1','#cae5f7','#92519c','#b878b0','#d7b1d2','#e7262a','#e94746','#eb666d','#ee838f','#f4abac','#fad9d9')
  png(file="1-nCoV-sample-percentage.png", width = dpi*8, height = dpi*4, units = "px",res = dpi,type='cairo')
  # Stacked barplot with multiple groups
  pp = ggplot(data=organ.summary.dataframe, aes(x=cluster, y=cell, fill=group)) + geom_bar(stat="identity",width = 0.6,position=position_fill(reverse = TRUE),size = 0.3,colour = '#222222') + labs(x='',y='Fraction of sample per cluster (%)') + 
    theme(axis.title.x = element_blank(),axis.text = element_text(size = 16),axis.title.y =element_text(size = 16), legend.text = element_text(size = 15),legend.key.height = unit(5,'mm'),
          legend.title = element_blank(),panel.grid.minor = element_blank()) + cowplot::theme_cowplot() + theme(axis.text.x = element_text(size = 10)) +
    scale_fill_manual(values = cols)
  print(pp)
  dev.off()
  
  
  nCoV.integrated1 <- RenameIdents(object = nCoV.integrated, 
                                   '13' = 'Epithelial','16' = 'Epithelial','25' = 'Epithelial','28' = 'Epithelial','31' = 'Epithelial',
                                   '0'='Macrophages','1'='Macrophages','2'='Macrophages','3'='Macrophages','4'='Macrophages','5'='Macrophages','7'='Macrophages','8'='Macrophages','10'='Macrophages','11'='Macrophages','12'='Macrophages','18'='Macrophages','21'='Macrophages','22'='Macrophages','23'='Macrophages','26'='Macrophages',
                                   '30'='Mast',
                                   '6'='T','9'='T','14'='T',
                                   '17'='NK',
                                   '19'='Plasma',
                                   '27'='B',
                                   '15'='Neutrophil',
                                   '20'='mDC',
                                   '29'='pDC',
                                   '24'='Doublets')
  
  nCoV.integrated1$celltype = Idents(nCoV.integrated1)
  nCoV.integrated1 = subset(nCoV.integrated1,subset = celltype != 'Doublets')
  nCoV_groups = c('Epithelial','Macrophages','Neutrophil','mDC','pDC','Mast','T','NK','B','Plasma')
  nCoV.integrated1$celltype = factor(nCoV.integrated1$celltype,levels = nCoV_groups,labels = nCoV_groups)
  Idents(nCoV.integrated1) = nCoV.integrated1$celltype
  
  ###HC/Mild/Severe
  png(file="1-nCoV-umap-split-group_1.png", width = dpi*11, height = dpi*4.3, units = "px",res = dpi,type='cairo')
  pp_temp = DimPlot(object = nCoV.integrated1, reduction = 'umap',label = FALSE, label.size = 6,split.by = 'group', ncol = 3,repel = TRUE,combine = TRUE)
  pp_temp = pp_temp + theme(axis.title = element_text(size = 17),axis.text = element_text(size = 17),
                            strip.text = element_text(family = 'arial',face='plain',size=17),
                            legend.text = element_text(size = 17),axis.line = element_line(size = 1),
                            axis.ticks = element_line(size = 0.8),
                            legend.key.height = unit(1.4,"line"))
  print(pp_temp)
  dev.off()
  
  ###sample
  png(file="1-nCoV-umap-split-sample.png", width = dpi*18, height = dpi*9, units = "px",res = dpi,type='cairo')
  pp_temp = DimPlot(object = nCoV.integrated1, reduction = 'umap',label = FALSE, split.by = 'sample_new', label.size = 7,ncol = 5,repel = TRUE, combine = TRUE,pt.size = 1.5)
  pp_temp = pp_temp + theme(axis.title = element_text(size = 22),axis.text = element_text(size = 22),
                            strip.text = element_text(family = 'sans',face='plain',size=22),
                            legend.text = element_text(size = 22),legend.key.height = unit(1.8,"line"),legend.key.width = unit(1.2,"line"),
                            axis.line = element_line(size = 1.2),axis.ticks = element_line(size = 1.2))
  print(pp_temp)
  dev.off()
  
  ####绘制比例图
  library(ggplot2)
  nCoV.integrated1[["cluster"]] <- Idents(object = nCoV.integrated1)
  big.cluster = nCoV.integrated1@meta.data
  organ.summary = table(big.cluster$group,big.cluster$cluster)
  write.table(organ.summary,file = '1-nCoV-percentage-group.txt',quote = FALSE,sep = '\t')
  
  library(ggplot2)
  library(dplyr)
  library(reshape2)
  organ.summary = read.delim2("1-nCoV-percentage-group.txt",header = TRUE, stringsAsFactors = FALSE,check.names = FALSE, sep = "\t")
  organ.summary$group = rownames(organ.summary)
  organ.summary.dataframe = melt(organ.summary)
  colnames(organ.summary.dataframe) = c('group','cluster','cell')
  organ.summary.dataframe = organ.summary.dataframe %>% filter(.,cluster %in% c('MP','Neutrophil','mDC','pDC','T','NK','B','Plasma'))
  organ.summary.dataframe$cell = as.numeric(organ.summary.dataframe$cell)
  organ.summary1 = organ.summary %>% select(c('MP','Neutrophil','mDC','pDC','T','NK','B','Plasma'))
  organ.summary2 <- round(organ.summary1 / rowSums(organ.summary1) * 100,2)
  write.table(organ.summary2,file = '1-nCoV-percentage-group-percent.txt',row.names = TRUE,quote = FALSE,sep='\t')
  
  dpi = 300
  png(file="1-nCoV-group-percentage.png", width = dpi*6, height = dpi*5, units = "px",res = dpi,type='cairo')
  # Stacked barplot with multiple groups
  pp = ggplot(data=organ.summary.dataframe, aes(x=group, y=cell, fill=cluster)) + geom_bar(stat="identity",width = 0.6,position=position_fill(reverse = TRUE),size = 0.3,colour = '#222222') + labs(x='',y='Fraction of cluster per group (%)') + 
    theme(axis.title.x = element_blank(),axis.text = element_text(size = 16),axis.title.y =element_text(size = 16), legend.text = element_text(size = 15),legend.key.height = unit(5,'mm'),
          legend.title = element_blank(),panel.grid.minor = element_blank()) + cowplot::theme_cowplot() + theme(axis.text.x = element_text(size = 10))
  
  print(pp)
  dev.off()
  
  pdf(file="1-nCoV-group-percentage.pdf", width = 6, height = 5)
  # Stacked barplot with multiple groups
  pp = ggplot(data=organ.summary.dataframe, aes(x=group, y=cell, fill=cluster)) + geom_bar(stat="identity",width = 0.6,position=position_fill(reverse = TRUE),size = 0.2,colour = '#222222') + labs(x='',y='Fraction of cluster per group (%)') + 
    theme(axis.title.x = element_blank(),axis.text = element_text(size = 16),axis.title.y =element_text(size = 16), legend.text = element_text(size = 15),legend.key.height = unit(5,'mm'),
          legend.title = element_blank(),panel.grid.minor = element_blank()) + cowplot::theme_cowplot() + theme(axis.text.x = element_text(size = 10))
  
  print(pp)
  dev.off()
  
  
  dpi = 300
  ###cell percentage for each sample and statistic analysis
  library(ggplot2)
  nCoV.integrated1[["cluster"]] <- Idents(object = nCoV.integrated1)
  big.cluster = nCoV.integrated1@meta.data
  organ.summary = as.data.frame.matrix(table(big.cluster$sample,big.cluster$cluster))
  write.table(organ.summary,file = '1-nCoV-percentage-sample.txt',quote = FALSE,sep = '\t')
  organ.summary1 = organ.summary %>% select(c('Macrophages','Neutrophil','mDC','pDC','T','NK','B','Plasma'))
  organ.summary2 <- round(organ.summary1 / rowSums(organ.summary1),3)
  write.table(organ.summary2,file = '1-nCoV-percentage-sample-percent.txt',row.names = TRUE,quote = FALSE,sep='\t')
  
  library(ggpubr)
  library(dplyr)
  library(ggplot2)
  library(cowplot)
  library(gridExtra)
  setwd('/Users/yang/Downloads/work/nCoV/scrna_rev/balf/Figure/')
  sample.percent = read.delim2("1-nCoV-percentage-sample-percent.txt",header = TRUE, stringsAsFactors = FALSE,check.names = FALSE, sep = "\t")
  nCoV_groups = c('Macrophages','Neutrophil','mDC','pDC','T','NK','B','Plasma')
  pplist = list()
  for(group_ in nCoV_groups){
    sample.percent_  = sample.percent %>% select(one_of(c('sample','group',group_)))
    colnames(sample.percent_) = c('sample','group','percent')
    sample.percent_$percent = as.numeric(sample.percent_$percent)
    sample.percent_ <- sample.percent_ %>% group_by(group) %>% mutate(upper =  quantile(percent, 0.75), 
                                                                      lower = quantile(percent, 0.25),
                                                                      mean = mean(percent),
                                                                      median = median(percent))
    print(group_)
    print(sample.percent_$median)
    
    pp1 = ggplot(sample.percent_,aes(x=group,y=percent)) + geom_jitter(shape = 21,aes(fill=group),width = 0.25) + stat_summary(fun.y=mean, geom="point", color="grey60") +
      theme_cowplot() +
      theme(axis.text = element_text(size = 6),axis.title = element_text(size = 6),legend.text = element_text(size = 6),
            legend.title = element_text(size = 6),plot.title = element_text(size = 6,face = 'plain'),legend.position = 'none') + labs(title = group_,y='Percentage') +
      geom_errorbar(aes(ymin = lower, ymax = upper),col = "grey60",width =  0.25)
    
    ###perform statistics analysis
    labely = max(sample.percent_$percent)
    compare_means(percent ~ group,  data = sample.percent_)
    my_comparisons <- list( c("HC", "M"), c("M", "S"), c("HC", "S") )
    pp1 = pp1 + stat_compare_means(comparisons = my_comparisons,size = 1,method = "t.test")
    pplist[[group_]] = pp1
    
  }
  pdf(file="sample-percentage.pdf", width = 6, height = 3)
  print(plot_grid(pplist[['Macrophages']],pplist[['Neutrophil']],
                  pplist[['mDC']],pplist[['pDC']],
                  pplist[['T']],pplist[['NK']],pplist[['B']],
                  pplist[['Plasma']],ncol = 4, nrow = 2))
  dev.off()
  
}