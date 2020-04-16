###TCR analysis
tcr_fig <- function(){
  ###NKT
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(ggplot2)
  library(reshape2)
  library(RColorBrewer)
  setwd('/home/data/results/workspace/nCoV_rev/balf/Figure')
  NKT.Integrated = readRDS("../NKT_T/6-NKT.rds")
  
  markers = c('TPPP3','KRT18','CD68','FCGR3B','CD1C','CLEC9A','LILRA4','TPSB2','CD3D','KLRD1','MS4A1','IGHG4')
  pdf(file="T-marker_heatmap.pdf", width = 8, height = 7)
  pp = DotPlot(NKT.Integrated, features = rev(markers),cols = c('white','#F8766D'),dot.scale =5) + RotatedAxis()
  pp = pp + theme(axis.text.x = element_text(size = 12),axis.text.y = element_text(size = 12)) + labs(x='',y='') + 
    guides(color = guide_colorbar(title = 'Scale expression'),size = guide_legend(title = 'Percent expressed')) + 
    theme(axis.line = element_line(size = 0.6))
  print(pp)
  dev.off()
  
  NKT.Integrated <- RenameIdents(object = NKT.Integrated, 
                                 '0' ='CCR7+ T','1'='ZNF683+ CD8 T','2' = 'GZMB+ CD8 T','3'='Proliferating T','4'='mixed T','5' ='Doublets',
                                 '6'='innate T','7'='PD1+ T','8'='Treg','9'='Doublets')
  NKT.Integrated$celltype = Idents(NKT.Integrated)
  nCoV_groups = c('CCR7+ T','ZNF683+ CD8 T','GZMB+ CD8 T','Proliferating T','innate T','PD1+ T','Treg','mixed T','Doublets')
  NKT.Integrated$celltype = factor(NKT.Integrated$celltype,levels = nCoV_groups,labels = nCoV_groups)
  Idents(NKT.Integrated) = NKT.Integrated$celltype
  NKT.Integrated1 = subset(NKT.Integrated,idents=c('CCR7+ T','ZNF683+ CD8 T','GZMB+ CD8 T','Proliferating T','innate T','PD1+ T','Treg','mixed T'))
  
  dpi = 300
  png(file="tcr-umap-nolabel.png", width = dpi*6, height = dpi*4, units = "px",res = dpi,type='cairo')
  pp = DimPlot(object = NKT.Integrated1, reduction = 'umap',label = FALSE,pt.size = 0.8,label.size = 4,repel = TRUE)
  pp = pp + theme(axis.title = element_text(size = 12),axis.text =  element_text(size = 12),
                  legend.text = element_text(size = 12),legend.key.height = unit(1.1,'line'),
                  axis.line = element_line(size = 1.0),axis.ticks = element_line(size = 1.0))
  print(pp)
  dev.off()
  
  ###heatmap
  nCoV_markers_rev = rev(c('CD3D','IL7R','CCR7','GZMA','CD8A','CD8B','NKG7','GZMB','CXCR3','GZMK','KLRD1','KLRC1','XCL1','MKI67','TYMS','KLRF1','EOMES','TRGV9','TRDV2','KLRB1','SLC4A10','CX3CR1',
                           'CD4','CXCR6','LAG3','MAF','PDCD1','HAVCR2','CTLA4','FOXP3','IL2RA'))
  nCoV_groups_rev = rev(c('CCR7+ T','ZNF683+ CD8 T','GZMB+ CD8 T','Proliferating T','innate T','PD1+ T','Treg','mixed T'))
  NKT.Integrated2 = NKT.Integrated1
  NKT.Integrated2 = subset(NKT.Integrated2,idents = nCoV_groups_rev)
  NKT.Integrated2$celltype = factor(NKT.Integrated2$celltype,levels = nCoV_groups_rev,labels = nCoV_groups_rev)
  Idents(NKT.Integrated2) = NKT.Integrated2$celltype
  pdf(file="tcr-heatmap-seurat.pdf", width = 9, height = 3.5)
  pp = DotPlot(NKT.Integrated2, features = nCoV_markers_rev,cols = c('white','#F8766D'),dot.scale=4.5) + RotatedAxis()
  pp = pp + theme(axis.text.x = element_text(size = 10),axis.text.y = element_text(size = 12)) + labs(x='',y='') + 
    guides(color = guide_colorbar(title = 'Scale expression'),size = guide_legend(title = 'Percent expressed')) + 
    theme(axis.line = element_line(size = 1))
  print(pp)
  dev.off()
  
  ##STEP0
  nkt.celltype = NKT.Integrated1@meta.data
  nkt.celltype$celltype = as.character(nkt.celltype$celltype)
  nkt.embedding = as.data.frame(Embeddings(NKT.Integrated1,reduction = 'umap'))
  nkt.embedding$ID = rownames(nkt.embedding)
  nkt.celltype = dplyr::left_join(nkt.celltype,nkt.embedding)
  nkt.celltype$celltype_big = nkt.celltype$celltype
  nkt.celltype[nkt.celltype$celltype == "CCR7+ CD4 T", "celltype_big"] <- "CD4 T"
  nkt.celltype[nkt.celltype$celltype == "ZNF683+ CD8 T", "celltype_big"] <- "CD8 T"
  nkt.celltype[nkt.celltype$celltype == "GZMB+ CD8 T", "celltype_big"] <- "CD8 T"
  write.table(nkt.celltype,file='7-tcr-nkt-meta.csv',row.names = FALSE,quote = FALSE,sep='\t')
  
  ##STEP2
  ####match meta and tcr,process with process.py first to add new barcode
  library(ggplot2)
  library(dplyr)
  setwd('/Users/yang/Downloads/work/nCoV/scrna_rev/tcr/')
  nkt.celltype = read.delim2("tcr-nkt-meta.csv",header = TRUE, stringsAsFactors = FALSE,check.names = FALSE, sep = "\t")
  nkt.celltype = nkt.celltype %>% dplyr::select(barcode_new,celltype,celltype_big,sample,sample_new,group,filter,UMAP_1,UMAP_2)
  colnames(nkt.celltype) = c('barcode','celltype','celltype_big','sample','sample_new','group','cellnum','UMAP_1','UMAP_2')
  tcr = read.delim2("TCR_both.tsv",header = TRUE, stringsAsFactors = FALSE,check.names = FALSE, sep = "\t")
  nkt.dataframe = dplyr::left_join(nkt.celltype,tcr)
  #nkt.dataframe = na.omit(nkt.dataframe)
  write.table(nkt.dataframe,file='nkt_tcr.csv',quote = FALSE,row.names = FALSE,sep='\t')
  
  ###STEP4
  ###preprocess with python first seperate to S/C and O group
  groups  = c('mild','severe')
  library(scales)
  for(group_ in groups){
    tcr.all = read.delim2(paste("tcr-",group_,'.csv',sep=''),header = TRUE, stringsAsFactors = FALSE,check.names = FALSE, sep = "\t")
    tcr.all$cloneCount = as.numeric(tcr.all$cloneCount)
    tcr.all$UMAP_1 = as.numeric(tcr.all$UMAP_1)
    tcr.all$UMAP_2 = as.numeric(tcr.all$UMAP_2)
    tcr.all = tcr.all[order(tcr.all$cloneCount),]
    
    maxcount = max(tcr.all$cloneCount)
    if(group_ == 'mild'){
      title_ = 'Moderate'
    }else{
      title_ = 'Severe'
    }
    
    dpi = 300
    png(file=paste("tcr-all-clone-",group_,".png",sep=""), width = dpi*6, height = dpi*4, units = "px",res = dpi,type='cairo')
    pp = ggplot(tcr.all,aes(x = UMAP_1,y=UMAP_2,colour=cloneCount)) + geom_point(size = 1.2)
    pp = pp + theme(axis.title = element_text(size = 12),axis.text =  element_text(size = 12),
                    legend.text = element_text(size = 12),legend.key.height = unit(1.1,'line'),
                    axis.line = element_line(size = 1.1),axis.ticks = element_line(size = 1.1)) + cowplot::theme_cowplot() + 
      scale_colour_gradientn('Clone count',colors=c("grey80","#FF9999","red"),values=rescale(c(0,10,maxcount)),limits=c(0,maxcount)) +
      ggtitle(title_)
    print(pp)
    dev.off()
    
    pdf(file=paste("tcr-all-clone-",group_,".pdf",sep=""), width = 6, height = 4)
    pp = ggplot(tcr.all,aes(x = UMAP_1,y=UMAP_2,colour=cloneCount)) + geom_point(size = 1.2)
    pp = pp + theme(axis.title = element_text(size = 12),axis.text =  element_text(size = 12),
                    legend.text = element_text(size = 12),legend.key.height = unit(1.1,'line'),
                    axis.line = element_line(size = 1.1),axis.ticks = element_line(size = 1.1)) + cowplot::theme_cowplot() + 
      scale_colour_gradientn('Clone count',colors=c("grey80","#FF9999","red"),values=rescale(c(0,10,maxcount)),limits=c(0,maxcount)) +
      ggtitle(title_)
    print(pp)
    dev.off()
    
  }
  
  ###calculate shannon index for all samples
  library(ggplot2)
  library(dplyr)
  setwd('/Users/yang/Downloads/work/nCoV/scrna_rev/tcr/')
  shannon_res = read.delim2('tcr-shannon.csv',header = TRUE, stringsAsFactors = FALSE,check.names = FALSE, sep = "\t")
  shannon_res$cnt = as.numeric(shannon_res$cnt)
  result_matrxi = matrix(0,length(unique(shannon_res$sample)),3)
  i = 1
  for(sample_ in unique(shannon_res$sample)){
    sample_i = shannon_res %>% filter(.,sample == sample_)
    sample_i$cloneFraction = sample_i$cnt / sum(sample_i$cnt)
    #shannon index
    shannon = -sum(sample_i$cloneFraction*log(sample_i$cloneFraction,base = exp(1)))
    
    #inversion simpson
    simpson = 1/sum(sample_i$cloneFraction*sample_i$cloneFraction)
    result_matrxi[i,1] = sample_
    result_matrxi[i,2] = shannon
    result_matrxi[i,3] = simpson
    
    i = i + 1
  }
  
  result_data = as.data.frame(result_matrxi)
  colnames(result_data) = c('sample','shannon','simpson')
  write.table(result_data,file = 'diversity.txt',row.names = FALSE,quote = FALSE)
  
  ###barplat for all cluster, count percentage
  library(ggplot2)
  library(dplyr)
  setwd('/Users/yang/Downloads/work/nCoV/scrna_rev/tcr/')
  nkt.tcr.all = read.delim2('nkt_tcr.csv',header = TRUE, stringsAsFactors = FALSE,check.names = FALSE, sep = "\t")
  nkt_groups = c('CD4 T','CD8 T','Proliferating T','gdT & MAIT','PD1+ T','Treg','mixed T')
  result.matrix = matrix(0,length(nkt_groups),3)
  i = 0
  for(nkt_group in nkt_groups){
    i = i + 1
    nkt.tcr.all.i = nkt.tcr.all %>% filter(.,celltype_big == nkt_group)
    nkt.tcr.all.i.tcr = nkt.tcr.all.i %>% filter(.,cdr3 != 'NA')
    result.matrix[i,1] = nkt_group
    result.matrix[i,2] = nrow(nkt.tcr.all.i)
    result.matrix[i,3] = nrow(nkt.tcr.all.i.tcr)
  }
  result.dataframe = as.data.frame(result.matrix)
  colnames(result.dataframe) = c('celltype','cellcount','tcrcount')
  write.table(result.dataframe,file = 'tcr-all-percentage.txt',row.names = FALSE,quote = FALSE,sep='\t')
  
  result.dataframe = read.delim2('tcr-all-percentage.txt',header = TRUE, stringsAsFactors = FALSE,check.names = FALSE, sep = "\t")
  result.dataframe$cellcount = as.numeric(result.dataframe$cellcount)
  result.dataframe$tcrcount = as.numeric(result.dataframe$tcrcount)
  result.dataframe$percentage = result.dataframe$tcrcount/result.dataframe$cellcount
  result.dataframe$celltype = factor(result.dataframe$celltype,levels = nkt_groups,labels = nkt_groups)
  pdf(file="tcr-all-percentage.pdf", width = 6, height = 4)
  pp = ggplot(result.dataframe,aes(x = celltype,y=percentage,fill =celltype )) + geom_bar(stat="identity",width = 0.6,size = 0.5,colour = '#222222')
  pp = pp + theme(axis.title = element_text(size = 15),axis.text = element_text(size = 15),
                  legend.text = element_text(size = 13)) + cowplot::theme_cowplot() +
    theme(axis.text.x = element_text(angle = 45,hjust = 1)) + guides(fill=FALSE) + ylab('Clonotype fraction') +
    xlab('T cell cluster') + theme(axis.line = element_line(size = 0.8))
  print(pp)
  dev.off()

  library(ggplot2)
  library(dplyr)
  setwd('/Users/yang/Downloads/work/nCoV/scrna_rev/tcr/')
  library(reshape2)
  library(gridExtra)
  library(cowplot)
  groups  = c('mild','severe')
  for(group_ in groups){
    tcr.count = read.delim2(paste("tcr-",group_,'-celltype-count.csv',sep=''),header = TRUE, stringsAsFactors = FALSE,check.names = FALSE, sep = "\t")
    tcr.clone = read.delim2(paste("tcr-",group_,'-celltype-clone.csv',sep=''),header = TRUE, stringsAsFactors = FALSE,check.names = FALSE, sep = "\t")
    tcr.count[,2:4] = sapply(tcr.count[,2:4],as.numeric)
    colnames(tcr.count) = c('celltype','1','2-5','>5')
    tcr.clone[,2:4] = sapply(tcr.clone[,2:4],as.numeric)
    colnames(tcr.clone) = c('celltype','1','2-5','>5')
    
    organ.summary = tcr.clone
    organ.summary.dataframe = melt(organ.summary)
    colnames(organ.summary.dataframe) = c('celltype','group','count')
    samples_name_new = rev(c('ZNF683+ CD8 T','GZMB+ CD8 T','Proliferating T','gdT & MAIT','CCR7+ CD4 T','PD1+ T','Treg','mixed T'))
    samples_name_new1 = rev(c('ZNF683+ CD8 T','GZMB+ CD8 T','Proliferating T','innate T','CCR7+ T','PD1+ T','Treg','mixed T'))
    organ.summary.dataframe$celltype = factor(organ.summary.dataframe$celltype,labels = samples_name_new1,levels = samples_name_new)
    organ.summary.dataframe$count = as.numeric(organ.summary.dataframe$count)
    
    organ.summary1 = tcr.count
    organ.summary.dataframe1 = melt(organ.summary1)
    colnames(organ.summary.dataframe1) = c('celltype','group','count')
    samples_name_new = rev(c('ZNF683+ CD8 T','GZMB+ CD8 T','Proliferating T','gdT & MAIT','CCR7+ CD4 T','PD1+ T','Treg','mixed T'))
    samples_name_new1 = rev(c('ZNF683+ CD8 T','GZMB+ CD8 T','Proliferating T','innate T','CCR7+ T','PD1+ T','Treg','mixed T'))
    organ.summary.dataframe1$celltype = factor(organ.summary.dataframe1$celltype,labels = samples_name_new1,levels = samples_name_new)
    organ.summary.dataframe1$count = as.numeric(organ.summary.dataframe1$count)
    
    pdf(file=paste("tcr-",group_,'-celltype.pdf',sep=''), width = 10, height = 2.6)
    pp1 = ggplot(data=organ.summary.dataframe, aes(x=celltype,y=count,fill = group)) + geom_bar(stat="identity",width = 0.6,position = position_stack(reverse = TRUE),size = 0.2,colour = '#222222') + labs(x='',y='Clonotype count') + 
      theme(axis.title.x = element_blank(),axis.text = element_text(size = 16),axis.title.y =element_text(size = 16), legend.text = element_text(size = 15),legend.key.height = unit(5,'mm'),
            legend.title = element_blank(),panel.grid.minor = element_blank()) + cowplot::theme_cowplot() + theme(axis.text.x = element_text(size = 10)) + 
      theme(legend.position = "left") + coord_flip() 
    
    pp2 = ggplot(data=organ.summary.dataframe1, aes(x=celltype,y=count,fill = group)) + geom_bar(stat="identity",width = 0.6,position = position_stack(reverse = TRUE),size = 0.2,colour = '#222222') + labs(x='',y='Cell count') + 
      theme(axis.title.x = element_blank(),axis.text.x = element_blank(),axis.text = element_text(size = 16),axis.title.y =element_text(size = 16), legend.text = element_text(size = 15),legend.key.height = unit(5,'mm'),
            legend.title = element_blank(),panel.grid.minor = element_blank()) + cowplot::theme_cowplot() + theme(axis.text.x = element_text(size = 10)) + 
       coord_flip() + scale_x_discrete(position = 'right')
    print(plot_grid(pp1, pp2,ncol = 2, nrow = 1))
    dev.off()
  }
  
  library(ggplot2)
  library(dplyr)
  setwd('/Users/yang/Downloads/work/nCoV/scrna_rev/tcr/')
  library(reshape2)
  library(gridExtra)
  library(cowplot)
  celltypes = c('CD4 T','CD8 T','Proliferating T')
  tcr.data = read.delim2('tcr-sample-celltype-result.csv',header = TRUE, stringsAsFactors = FALSE,check.names = FALSE, sep = "\t")
  for(celltype_ in celltypes){
    tcr.data.i = tcr.data %>% filter(.,celltype == celltype_)
    tcr.data.i$celltype = NULL
    pplist = list()
    for(type_ in unique(tcr.data.i$type)){
      tcr.data.j = tcr.data.i %>% filter(.,type==type_)
      tcr.data.j$type = NULL
      colnames(tcr.data.j) = c('sample','1','2-5','>5')
      tcr.data.j[,2:4] = sapply(tcr.data.j[,2:4],as.numeric)
      
      organ.summary = tcr.data.j
      organ.summary.dataframe = melt(organ.summary)
      colnames(organ.summary.dataframe) = c('sample','group','count')
      samples_name_new = rev(c('O1','O2','O3','S1','C1','C3','C4','C5'))
      organ.summary.dataframe$sample = factor(organ.summary.dataframe$sample,labels = samples_name_new,levels = samples_name_new)
      organ.summary.dataframe$count = as.numeric(organ.summary.dataframe$count)
      pp1 = ggplot(data=organ.summary.dataframe, aes(x=sample,y=count,fill = group)) + geom_bar(stat="identity",width = 0.6,position=position_fill(reverse = TRUE),size = 0.2,colour = '#222222') + 
        theme(axis.title.x = element_blank(),axis.text = element_text(size = 16),axis.title.y =element_text(size = 16), legend.text = element_text(size = 15),legend.key.height = unit(5,'mm'),
              legend.title = element_blank(),panel.grid.minor = element_blank()) + cowplot::theme_cowplot() + theme(axis.text.x = element_text(size = 10)) + 
        theme(legend.position = "left") + coord_flip() 
      if(type_ == 'clone'){
        pp1 = pp1 + labs(x='',y='Clonotype fraction') 
      }else{
        pp1 = pp1 + labs(x='',y='Cell fraction') 
      }
      print(type_)
      pplist[[type_]] = pp1
    }
    pdf(file=paste("tcr-",celltype_,'-sample.pdf',sep=''), width = 10, height = 3.4)
    print(plot_grid(pplist[['clone']], pplist[['count']],ncol = 2, nrow = 1))
    dev.off()
  }

  library(ggplot2)
  library(dplyr)
  setwd('/Users/yang/Downloads/work/nCoV/scrna_rev/tcr/')
  library(reshape2)
  library(gridExtra)
  library(cowplot)
  celltypes = c('CCR7+ CD4 T','ZNF683+ CD8 T','GZMB+ CD8 T','Proliferating T','gdT & MAIT','PD1+ T','Treg','mixed T')
  color_list = list('CCR7+ CD4 T'='#e3964e','ZNF683+ CD8 T'='#4fa65b','GZMB+ CD8 T'='#3f8dc1','Proliferating T'='#75659d','gdT & MAIT'='#c64e48','PD1+ T'='#97c9d3','Treg'='#df9db8','mixed T'='#ac7753')
  
  ###STEP1
  groups = c('mild','severe')
  for(group_ in groups){
    tcr.all = read.delim2(paste("tcr-",group_,'.csv',sep=''),header = TRUE, stringsAsFactors = FALSE,check.names = FALSE, sep = "\t")
    tcr.all$cloneCount = as.numeric(tcr.all$cloneCount)
    tcr.all = tcr.all[order(-tcr.all$cloneCount),]
    if(group_ == 'mild'){
      title_ = 'Ordinary'
    }else{
      title_ = 'Severe/Critical'
    }
    top20 = unique(tcr.all$cdr3)[1:20]
    tcr.all.i = tcr.all %>% filter(.,cdr3 %in% top20)
    samples = unique(tcr.all.i$sample_new)
    print(table(tcr.all.i$sample_new,tcr.all.i$cdr3))
    write.table(tcr.all.i,file=paste("trc-clonetype-top20-",group_,".txt",sep=""),row.names = FALSE,quote = FALSE,sep='\t')
  }
  
  ##STEP3
  groups = c('mild','severe')
  cols = c('#e3964e','#4fa65b','#3f8dc1','#75659d','#c64e48','#97c9d3','#df9db8','#ac7753' )
  for(group_ in groups){
    tcr.bar = read.delim2(paste("tcr-clonotype-bar-",group_,'.csv',sep=''),header = TRUE, stringsAsFactors = FALSE,check.names = FALSE, sep = "\t")
    ###bar plot
    organ.summary.dataframe = tcr.bar
    labels_temp = tcr.bar %>% select(clonoID,sample)
    labels_temp = unique(labels_temp)
    labels_temp = labels_temp %>% mutate(label=paste(clonoID,'\n',sample,sep=''))
    
    samples_name_new = c('CCR7+ CD4 T','ZNF683+ CD8 T','GZMB+ CD8 T','Proliferating T','gdT & MAIT','PD1+ T','Treg','mixed T')
    organ.summary.dataframe$celltype = factor(organ.summary.dataframe$celltype,labels = samples_name_new,levels = samples_name_new)
    organ.summary.dataframe$clonoID = factor(organ.summary.dataframe$clonoID,levels = c('C1','C2','C3','C4','C5','C6','C7','C8','C9','C10','C11','C12','C13','C14','C15','C16','C17','C18','C19','C20'),
                                             labels = labels_temp$label)
    organ.summary.dataframe$count = as.numeric(organ.summary.dataframe$count)
    pp1 = ggplot(data=organ.summary.dataframe, aes(x=clonoID,y=count,fill = celltype)) + geom_bar(stat="identity",width = 0.6,position = position_stack(reverse = TRUE),size = 0.2,colour = '#222222') + labs(x='',y='Clonotype count') + 
      theme(axis.title.x = element_blank(),axis.text = element_text(size = 16),axis.title.y =element_text(size = 16), legend.text = element_text(size = 15),legend.key.height = unit(5,'mm'),
            legend.title = element_blank(),panel.grid.minor = element_blank()) + cowplot::theme_cowplot() + theme(axis.text.x = element_text(size = 10)) +
      scale_fill_manual(values = cols) + scale_x_discrete(position = "top") + scale_y_continuous(trans = "reverse")
    pdf(file=paste("tcr-clonotype-bar-",group_,".pdf",sep=""), width = 12, height = 3)
    print(pp1)
    dev.off()
  }
    
    
  if(1==0){   ###pie plot
    pplist = list()
    blank_theme <- theme_minimal()+
      theme(
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.border = element_blank(),
        panel.grid=element_blank(),
        axis.ticks = element_blank(),
        plot.title=element_text(size=14, face="bold")
      )
    for(sample_ in samples){
      tcr.data.j = tcr.all.i %>% filter(.,sample_new==sample_)
      #tcr.data.j$sample_new = NULL
      tcr.data.j  = tcr.data.j %>% select(barcode,celltype)
      colnames(tcr.data.j) = c('barcode','celltype')
      tcr.data.j_cellcount = as.data.frame(table(tcr.data.j$celltype))
      colnames(tcr.data.j_cellcount) = c('celltype','count')
      tcr.data.j_cellcount$celltype = factor(tcr.data.j_cellcount$celltype,levels = celltypes,labels = celltypes)
      celltype_i = tcr.data.j_cellcount$celltype
      cols = c()
      for(item in celltypes){
        if(item %in% celltype_i){
          cols = c(cols,color_list[[item]])
        }
      }
    
      pp1 = ggplot(tcr.data.j_cellcount, aes(x="", y=count, fill=celltype)) + geom_bar(width = 1, stat = "identity") +
        coord_polar("y", start=0) + blank_theme + theme(axis.text.x=element_blank()) + ggtitle(paste(sample_," n = ",nrow(tcr.data.j),sep="")) +
        scale_fill_manual(values = cols)
      pplist[[sample_]] = pp1
    }
    
  }
  
  setwd('/Users/yang/Downloads/work/nCoV/scrna_rev/tcr/')
  tcr.top20 = read.delim2("nkt_tcr_full_plotdata.csv",header = TRUE, stringsAsFactors = FALSE,check.names = FALSE, sep = "\t")
  for(celltype_ in unique(tcr.top20$celltype)){
    print(celltype_)
    tcr.top20.1 = tcr.top20 %>% filter(.,celltype == celltype_)
    pplist = list()
    for(chain_ in unique(tcr.top20$chain)){
      mait = tcr.top20.1 %>% filter(.,chain == chain_)
      mait$count = as.numeric(mait$count)
      mait = mait[order(-mait$count),]
      mait$gene = factor(mait$gene,levels = mait$gene,labels = mait$gene)
      pp1 = ggplot(data=mait[1:40,], aes(x=gene,y=count)) + geom_bar(stat="identity",width = 0.6,position = position_stack(reverse = TRUE),size = 0.2,colour = '#222222') + labs(x='',y='Clonotype count') + 
        theme(axis.title.x = element_blank(),axis.text = element_text(size = 10),axis.title.y =element_text(size = 10), legend.text = element_text(size = 13),legend.key.height = unit(5,'mm'),
              legend.title = element_blank(),panel.grid.minor = element_blank()) + cowplot::theme_cowplot() + theme(axis.text.x = element_text(size = 10,angle = -45,hjust = 0),
                                                                                                                    plot.margin=unit(c(6,60,6,6),"points")) +
        ggtitle(paste(celltype_,' / ',chain_,sep=''))+ theme(plot.title = element_text(size = 12, face = "plain"))
      pplist[[chain_]] = pp1
    }
    
    pdf(file=paste("VJ-",celltype_,'.pdf',sep=''), width = 10, height = 6)
    print(plot_grid(pplist[['TRA']],pplist[['TRB']],ncol = 1, nrow = 2))
    dev.off()
  }
  
  ####match VDJ db
  library(dplyr)
  library(ggplot2)
  setwd('/Users/yang/Downloads/work/nCoV/scrna_rev/tcr/')
  nkt_tcr = read.delim2("nkt_tcr_full.csv",header = TRUE, stringsAsFactors = FALSE,check.names = FALSE, sep = "\t")
  nkt_tcr = nkt_tcr %>% filter(.,chain == 'TRB')
  VDJDB = read.delim2("VDJDB.tsv",header = TRUE, stringsAsFactors = FALSE,check.names = FALSE, sep = "\t")
  VDJDB = VDJDB %>% select(Gene,CDR3,V,J,`Epitope species`)
  colnames(VDJDB) = c('DB_gene','DB_cdr3','DB_V','DB_J','DB_species')
  VDJDB = VDJDB %>% mutate(id = paste(DB_gene,'_',DB_cdr3,sep=""))
  nkt_tcr = nkt_tcr %>% mutate(id=paste(chain,'_',cdr3,sep=''))
  nkt_tcr = dplyr::left_join(nkt_tcr,VDJDB)
  nkt_tcr = na.omit(nkt_tcr)
  nkt_tcr = nkt_tcr %>% select(barcode_new,celltype,sample_new,group,chain,v_gene,DB_V,j_gene,DB_J,DB_species,cdr3,DB_cdr3)
  write.table(nkt_tcr,file = 'nkt_tcr_full_VDJDB.csv',quote = FALSE,row.names = FALSE,sep='\t')
  
  
  ###gliph
  #CDR3b		TRBV	TRBJ	CDR3a		TRAV		TRAJ	PatientCounts
  #CAADTSSGANVLTF	TRBV30	TRBJ2-6	CALSDEDTGRRALTF	TRAV19		TRAJ5	09/02171
  library(dplyr)
  library(ggplot2)
  setwd('/Users/yang/Downloads/work/nCoV/scrna_rev/tcr/')
  nkt_tcr = read.delim2("nkt_tcr_full.csv",header = TRUE, stringsAsFactors = FALSE,check.names = FALSE, sep = "\t")
  nkt_tcr_trb = nkt_tcr %>% filter(.,chain == 'TRB')
  nkt_tcr_trb = nkt_tcr_trb %>% select(cdr3,v_gene,j_gene,sample_new,barcode_new)
  colnames(nkt_tcr_trb) = c('CDR3b','TRBV','TRBJ','PatientCounts','barcode_new')
  nkt_tcr_tra = nkt_tcr %>% filter(.,chain == 'TRA')
  nkt_tcr_tra = nkt_tcr_tra %>% select(cdr3,v_gene,j_gene,barcode_new)
  colnames(nkt_tcr_tra) = c('CDR3a','TRAV','TRAJ','barcode_new')
  nkt_tcr_both = dplyr::left_join(nkt_tcr_trb,nkt_tcr_tra)
  nkt_tcr_both = nkt_tcr_both %>% select('CDR3b','TRBV','TRBJ','CDR3a','TRAV','TRAJ','PatientCounts','barcode_new')
  nkt_tcr_both$EpitopeGene = ''
  nkt_tcr_both$EpitopeSpecies = ''
  
  VDJDB = read.delim2("VDJDB.tsv",header = TRUE, stringsAsFactors = FALSE,check.names = FALSE, sep = "\t")
  VDJDB_trb = VDJDB %>% filter(.,Gene == 'TRB')
  VDJDB_trb = VDJDB_trb %>% select(CDR3,V,J,`Epitope gene`,`Epitope species`)
  colnames(VDJDB_trb) = c('CDR3b','TRBV','TRBJ','EpitopeGene','EpitopeSpecies')
  VDJDB_trb$CDR3a = ''
  VDJDB_trb$TRAV = ''
  VDJDB_trb$TRAJ = ''
  VDJDB_trb$PatientCounts = ''
  VDJDB_trb$barcode_new = ''
  VDJDB_both = VDJDB_trb %>% select('CDR3b','TRBV','TRBJ','CDR3a','TRAV','TRAJ','PatientCounts','barcode_new','EpitopeGene','EpitopeSpecies')
  
  res = rbind(nkt_tcr_both,VDJDB_both)
  
  write.table(res,'gliph.csv',row.names = FALSE,quote = FALSE,sep='\t')
  
  
  ####compare clone >=2 and none clone = 2 in mild CD8 Tcell, process with python first
  library(ggplot2)
  library(dplyr)
  setwd('/Users/yang/Downloads/work/nCoV/scrna_rev/tcr/')
  nkt.tcr.all = read.delim2('tcr-mild-cd8T-clone.csv',header = TRUE, stringsAsFactors = FALSE,check.names = FALSE, sep = "\t")
  nkt.tcr.all = nkt.tcr.all %>% select(barcode,cloneGroup)
  colnames(nkt.tcr.all) = c('barcode_new','cloneGroup')
  ###
  tcr_nkt_meta = read.delim2('tcr-nkt-meta.csv',header = TRUE, stringsAsFactors = FALSE,check.names = FALSE, sep = "\t")
  tcr_nkt_meta = tcr_nkt_meta%>% select(ID,barcode_new)
  nkt.tcr.all = dplyr::left_join(nkt.tcr.all,tcr_nkt_meta)
  write.table(nkt.tcr.all,file='tcr-mild-cd8T-clone.meta',row.names = FALSE,quote = FALSE,sep='\t')
  
  
  ####load nkt data
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(ggplot2)
  library(reshape2)
  setwd('/home/data/results/workspace/nCoV_rev/balf/Figure')
  NKT.Integrated = readRDS("../NKT/6-NKT.rds")
  NKT.Integrated <- RenameIdents(object = NKT.Integrated, 
                                 '0' ='CCR7+ T','1'='CD8 T','2' = 'CD8 T','3'='Proliferating T','4'='NK','5' ='CD8 T',
                                 '6'='Treg','7'='Doublets','8'='NK','9'='innate T','10'='Proliferating T','11'='Uncertain','12'='Uncertain',
                                 '13'='Doublets')
  NKT.Integrated$celltype = Idents(NKT.Integrated)
  nCoV_groups = c('CCR7+ T','CD8 T','Proliferating T','NK','Treg','innate T','Uncertain','Doublets')
  NKT.Integrated$celltype = factor(NKT.Integrated$celltype,levels = nCoV_groups,labels = nCoV_groups)
  Idents(NKT.Integrated) = NKT.Integrated$celltype
  
  NKT.Integrated1 = NKT.Integrated
  tcr_nkt_meta = read.delim2('tcr-mild-cd8T-clone.meta',header = TRUE, stringsAsFactors = FALSE,check.names = FALSE, sep = "\t")
  NKT.Integrated1 = subset(NKT.Integrated1,cells = tcr_nkt_meta$ID)
  rownames(tcr_nkt_meta) = tcr_nkt_meta$ID
  tcr_nkt_meta$ID = NULL
  NKT.Integrated1 = AddMetaData(object = NKT.Integrated1, metadata = tcr_nkt_meta)
  
  DefaultAssay(NKT.Integrated1) = 'RNA'
  Idents(NKT.Integrated1) = 'cloneGroup'
  NKT.Integrated1@misc$markers <- FindAllMarkers(object = NKT.Integrated1, assay = 'RNA',only.pos = TRUE, test.use = 'MAST',logfc.threshold = 0)
  write.table(NKT.Integrated1@misc$markers,file=paste('4-mild-cd8T-clone-all.txt',sep=''),row.names = FALSE,quote = FALSE,sep = '\t')
  NKT.Integrated1@misc$markers <- FindAllMarkers(object = NKT.Integrated1, assay = 'RNA',only.pos = TRUE, test.use = 'MAST',logfc.threshold = 0)
  write.table(NKT.Integrated1@misc$markers,file=paste('4-mild-cd8T-clone.txt',sep=''),row.names = FALSE,quote = FALSE,sep = '\t')
  
  for(sample_ in c('O1','O2','O3')){
    print(sample_)
    temp_i = subset(NKT.Integrated1,subset = sample_new==sample_)
    print(table(temp_i@meta.data$cloneGroup))
    Idents(temp_i) = 'cloneGroup'
    print(table(temp_i@meta.data$cloneGroup))
    DefaultAssay(temp_i) = 'RNA'
    temp_i@misc$markers <- FindAllMarkers(object = temp_i, assay = 'RNA',only.pos = TRUE, test.use = 'MAST',logfc.threshold = 0)
    write.table(temp_i@misc$markers,file=paste('4-mild-cd8T-clone-',sample_,'.txt',sep=''),row.names = FALSE,quote = FALSE,sep = '\t')
  }
  
  library(ggrepel)
  library(ggplot2)
  library(dplyr)
  library(reshape2)
  setwd('/Users/yang/Downloads/work/nCoV/scrna_rev/balf/Figure/')
  dpi = 300
  nkt.markers <- read.delim2("4-mild-cd8T-clone-all.txt",header = TRUE, stringsAsFactors = FALSE,check.names = FALSE, sep = "\t")
  nkt.markers$p_val_adj = as.numeric(nkt.markers$p_val_adj)
  nkt.markers$avg_logFC = as.numeric(nkt.markers$avg_logFC)
  nkt.markers = nkt.markers %>% mutate(avg_logFC_new = ifelse(cluster == 'Y',avg_logFC,-avg_logFC))
  nkt.markers$color = '#BEBEBE'
  nkt.markers$group = 'Other'
  nkt.markers$label = nkt.markers$gene
  rownames(nkt.markers) = nkt.markers$gene
  group_7_up = nkt.markers %>% filter(., cluster == 'Y' & p_val_adj < 0.05 )
  group_7_down = nkt.markers %>% filter(.,cluster == 'N' & p_val_adj < 0.05 )
  print(dim(group_7_up))
  print(dim(group_7_down))
  
  ##包含在t vs macrophage差异上调的里面
  nkt.markers[group_7_up$gene,'color'] = '#d53d49'
  nkt.markers[group_7_up$gene,'group'] = 'Y'
  nkt.markers[group_7_down$gene,'color'] = '#418bbe'
  nkt.markers[group_7_down$gene,'group'] = 'N'
  
  nkt.markers$p_val_adj = -log10(nkt.markers$p_val_adj)
  top10 = nkt.markers %>% group_by(group) %>% top_n(n = 10, wt = p_val_adj)
  top10 = top10 %>% filter(.,group !='Other')
  my_gene_list = c('CD2','ZNF683','ITGA4','IKZF3','ITGB7','XCL1','XCL2','ID2','IRF1',
                   'CCR2','TOX','TNFSF14','CCL5','JAK1','GBP5','STAT1','LAT','HLA-E','CD3E',
                   'OASL','MX1','RPLP0','RPS15A','RPS26','RPS24','RPL7','RPL39','RPS13','EIF1',
                   'EIF1AY','DUT','TYMS','ATP5MC2','ATP5MG','UQCRB','LDHA',
                   'COX6C','NDUFS5','DNPH1','NDUFA12','HMGB2')
  #my_not_gene_list = c('MT−ATP6','MT−CO3','MT−ND4L','MT−CYB')
  my_gene_list = c('ITGAE')
  my_not_gene_list = c()
  nkt.markers$group = factor(nkt.markers$group,labels = c('Other','Y','N'),levels = c('Other','Y','N'))
  nkt.markers[!nkt.markers$gene %in% c(my_gene_list,top10$gene),'label'] = ''
  
  nkt.markers = with(nkt.markers, nkt.markers[order(group),])
  
  pdf(file = "4-mild-cd8T-clone-all.pdf",width = 3,height = 3)
  p <- ggplot(nkt.markers, aes(avg_logFC_new, p_val_adj,color = group,order = group)) + geom_point(size = 0.8)+
    labs(x='logFC',y='-Log10(p.adjust)')
  p <- p + theme(plot.title = element_text(hjust = 0.5)) + 
    theme(axis.text = element_text(size = 5),axis.title = element_text(size = 5),
          axis.title.x = element_text(margin = unit(c(4, 0, 0, 0), "mm")),
          axis.title.y = element_text(margin = unit(c(0, 4, 0, 0), "mm")),
          axis.line = element_line(colour = 'black',size = 0.6),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          panel.grid.major = element_line(colour = 'grey',size=0.4,linetype = 2),
          axis.ticks = element_line(colour = 'black',size = 0.6),
          legend.title = element_blank(),legend.text = element_text(size = 5,family = 'sans',face='plain'),
          legend.position = 'top',legend.key.size = unit(1.5, 'lines'))
  p = p + scale_color_manual(values = setNames(nkt.markers$color, nkt.markers$group))
  p = p + geom_text_repel(aes(label = nkt.markers$label),color= 'black',size = 1.8, box.padding = 0.25,point.padding = 0.5,segment.color = 'grey80')
  p = p + theme(legend.key=element_rect(fill=NA))
  print(p)
  dev.off()
  
  ###1. tissue resident score: cd8T mild group expanded vs non-expanded
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(ggplot2)
  library(reshape2)
  setwd('/home/data/results/workspace/nCoV_rev/balf/Figure/')
  library(sctransform)
  NKT.Integrated = readRDS("../NKT/6-NKT.rds")
  tcr_nkt_meta = read.delim2('tcr-mild-cd8T-clone.meta',header = TRUE, stringsAsFactors = FALSE,check.names = FALSE, sep = "\t")
  NKT.Integrated <- RenameIdents(object = NKT.Integrated, 
                                 '0' ='CCR7+ T','1'='CD8 T','2' = 'CD8 T','3'='Proliferating T','4'='NK','5' ='CD8 T',
                                 '6'='Treg','7'='Doublets','8'='NK','9'='innate T','10'='Proliferating T','11'='Uncertain','12'='Uncertain',
                                 '13'='Doublets')
  NKT.Integrated$celltype = Idents(NKT.Integrated)
  nCoV_groups = c('CCR7+ T','CD8 T','Proliferating T','NK','Treg','innate T','Uncertain','Doublets')
  NKT.Integrated$celltype = factor(NKT.Integrated$celltype,levels = nCoV_groups,labels = nCoV_groups)
  Idents(NKT.Integrated) = NKT.Integrated$celltype
  
  NKT.Integrated1 = NKT.Integrated
  NKT.Integrated1 = subset(NKT.Integrated1,cells = tcr_nkt_meta$ID)
  NKT.Integrated1@meta.data = dplyr::left_join(NKT.Integrated1@meta.data,tcr_nkt_meta)
  NKT.Integrated1 <- SCTransform(object = NKT.Integrated1, verbose = FALSE)
  write.table(t(NKT.Integrated1@assays$SCT@data), "4_T_sct.tsv",sep="\t", row.names = TRUE, quote = FALSE)
  
  ###now calculate score
  matrix.file = '4_T_sct.tsv'
  marker_p.file = 'marker_p.txt'
  marker_n.file = 'marker_n.txt'
  
  matrix = read.delim2(matrix.file,header = TRUE, stringsAsFactors = FALSE,check.names = FALSE, sep = "\t")
  matrix[]=lapply(matrix, as.numeric)
  marker_p = read.delim2(marker_p.file,header = FALSE, stringsAsFactors = FALSE,check.names = FALSE, sep = "\t")
  colnames(marker_p) = c('gene')
  marker_n = read.delim2(marker_n.file,header = FALSE, stringsAsFactors = FALSE,check.names = FALSE, sep = "\t")
  colnames(marker_n) = c('gene')
  
  ###now calculate the score
  totalUMI = as.data.frame(rowSums(matrix))
  colnames(totalUMI) = c('totalUMI')
  
  posi.matrix = matrix %>% select(one_of(marker_p$gene))
  positiveUMI = as.data.frame(rowSums(posi.matrix))
  colnames(positiveUMI) = c('posiUMI')
  nega.matrix = matrix %>% select(one_of(marker_n$gene))
  negativeUMI = as.data.frame(rowSums(nega.matrix))
  colnames(negativeUMI) = c('negaUMI')
  negativeUMI$negaUMI = -negativeUMI$negaUMI
  
  score.p = positiveUMI/totalUMI
  score.n = negativeUMI/totalUMI
  
  score.p.norm = (score.p - min(score.p) ) / (max(score.p) - min(score.p)) 
  score.n.norm = (score.n - min(score.n) ) / (max(score.n) - min(score.n)) 
  sig.score = score.p.norm*score.n.norm
  colnames(sig.score) = c('sigscore')
  sig.score$barcode = rownames(sig.score)
  sig.score = sig.score %>% select(one_of(c('barcode','sigscore')))
  write.table(sig.score,file = '4_T_sigscore.txt',quote = FALSE,row.names = FALSE,sep='\t')
  
  ###violin plot
  colnames(sig.score) = c('ID','sigscore')
  sig.plot.data = dplyr::left_join(sig.score,tcr_nkt_meta)
  write.table(sig.plot.data,file='4-mild-cd8T-clone-signature.txt',row.names = FALSE,quote = FALSE,sep='\t')
  
  setwd('/Users/yang/Downloads/work/nCoV/scrna_rev/balf/Figure/')
  sig.plot.data <- read.delim2("4-mild-cd8T-clone-signature.txt",header = TRUE, stringsAsFactors = FALSE,check.names = FALSE, sep = "\t")
  sig.plot.data$sigscore = as.numeric(sig.plot.data$sigscore)
  sig.plot.data$cloneGroup = factor(sig.plot.data$cloneGroup,levels = c('Y','N'),labels = c('Expanded','Non-expanded'))
  pdf(file = "4-mild-cd8T-clone-signature.pdf",width = 2.5,height = 1.5)
  p <- ggplot(sig.plot.data, aes(x=cloneGroup, y=sigscore,color=cloneGroup)) + geom_violin(trim=FALSE)
  p = p + theme_cowplot() + theme(axis.text = element_text(size = 5),axis.title = element_text(size = 5),
                legend.text = element_text(size = 5),legend.title = element_text(size = 5),legend.key.size = unit(3,'mm')) +
    scale_color_manual(values=c("#d53d49", "#418bbe"))  + stat_summary(fun.data=mean_sdl, mult=1, geom="pointrange", color="grey70") + 
    labs(x='',y='Tissue-resident core signatures score',color = "Group") 
  
  print(p)
  dev.off()
  
  ###t test
  library(lsr)
  expanded = sig.plot.data %>% filter(.,cloneGroup == 'Expanded')
  expanded = expanded$sigscore
  nonexpanded = sig.plot.data %>% filter(.,cloneGroup == 'Non-expanded')
  nonexpanded = nonexpanded$sigscore
  t.test(expanded,nonexpanded)
  cohensD(expanded,nonexpanded)
  ###p-value = 1.699e-14
  ###Cohen’s d 0.4704084
  
  
  ###
  ###2. tissue resident score: cd8T mild vs severe
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(ggplot2)
  library(reshape2)
  setwd('/home/data/results/workspace/nCoV_rev/balf/Figure/')
  library(sctransform)
  NKT.Integrated = readRDS("../NKT/6-NKT.rds")
  tcr_nkt_meta = read.delim2('tcr-mild-cd8T-clone.meta',header = TRUE, stringsAsFactors = FALSE,check.names = FALSE, sep = "\t")
  NKT.Integrated <- RenameIdents(object = NKT.Integrated, 
                                 '0' ='CCR7+ T','1'='CD8 T','2' = 'CD8 T','3'='Proliferating T','4'='NK','5' ='CD8 T',
                                 '6'='Treg','7'='Doublets','8'='NK','9'='innate T','10'='Proliferating T','11'='Uncertain','12'='Uncertain',
                                 '13'='Doublets')
  NKT.Integrated$celltype = Idents(NKT.Integrated)
  nCoV_groups = c('CCR7+ T','CD8 T','Proliferating T','NK','Treg','innate T','Uncertain','Doublets')
  NKT.Integrated$celltype = factor(NKT.Integrated$celltype,levels = nCoV_groups,labels = nCoV_groups)
  Idents(NKT.Integrated) = NKT.Integrated$celltype
  
  NKT.Integrated1 = NKT.Integrated
  NKT.Integrated1 = subset(NKT.Integrated1,subset = celltype=='CD8 T' & group != 'HC')
  NKT.Integrated1 <- SCTransform(object = NKT.Integrated1, verbose = FALSE)
  T_sct_mild_severe_meta = NKT.Integrated1@meta.data
  write.table(T_sct_mild_severe_meta, "4_T_sct_mild_severe_meta.tsv",sep="\t", row.names = TRUE, quote = FALSE)
  write.table(t(NKT.Integrated1@assays$SCT@data), "4_T_sct_mild_severe.tsv",sep="\t", row.names = TRUE, quote = FALSE)
  
  ###now calculate score
  matrix.file = '4_T_sct_mild_severe.tsv'
  marker_p.file = 'marker_p.txt'
  marker_n.file = 'marker_n.txt'
  
  matrix = read.delim2(matrix.file,header = TRUE, stringsAsFactors = FALSE,check.names = FALSE, sep = "\t")
  matrix[]=lapply(matrix, as.numeric)
  marker_p = read.delim2(marker_p.file,header = FALSE, stringsAsFactors = FALSE,check.names = FALSE, sep = "\t")
  colnames(marker_p) = c('gene')
  marker_n = read.delim2(marker_n.file,header = FALSE, stringsAsFactors = FALSE,check.names = FALSE, sep = "\t")
  colnames(marker_n) = c('gene')
  
  ###now calculate the score
  totalUMI = as.data.frame(rowSums(matrix))
  colnames(totalUMI) = c('totalUMI')
  
  posi.matrix = matrix %>% select(one_of(marker_p$gene))
  positiveUMI = as.data.frame(rowSums(posi.matrix))
  colnames(positiveUMI) = c('posiUMI')
  nega.matrix = matrix %>% select(one_of(marker_n$gene))
  negativeUMI = as.data.frame(rowSums(nega.matrix))
  colnames(negativeUMI) = c('negaUMI')
  negativeUMI$negaUMI = -negativeUMI$negaUMI
  
  score.p = positiveUMI/totalUMI
  score.n = negativeUMI/totalUMI

  score.p.norm = (score.p - min(score.p) ) / (max(score.p) - min(score.p)) 
  score.n.norm = (score.n - min(score.n) ) / (max(score.n) - min(score.n)) 
  sig.score = score.p.norm*score.n.norm
  colnames(sig.score) = c('sigscore')
  sig.score$barcode = rownames(sig.score)
  sig.score = sig.score %>% select(one_of(c('barcode','sigscore')))
  write.table(sig.score,file = '4_T_sct_mild_severe_sigscore.txt',quote = FALSE,row.names = FALSE,sep='\t')
  
  ###violin plot
  tcr_nkt_meta  = read.delim2("4_T_sct_mild_severe_meta.tsv",header = TRUE, stringsAsFactors = FALSE,check.names = FALSE, sep = "\t")
  tcr_nkt_meta = tcr_nkt_meta %>% select(ID,group)
  colnames(sig.score) = c('ID','sigscore')
  sig.plot.data = dplyr::left_join(sig.score,tcr_nkt_meta)
  write.table(sig.plot.data,file='4-cd8T_mild_severe-clone-signature.txt',row.names = FALSE,quote = FALSE,sep='\t')
  
  setwd('/Users/yang/Downloads/work/nCoV/scrna_rev/balf/Figure/')
  sig.plot.data <- read.delim2("4-cd8T_mild_severe-clone-signature.txt",header = TRUE, stringsAsFactors = FALSE,check.names = FALSE, sep = "\t")
  sig.plot.data$sigscore = as.numeric(sig.plot.data$sigscore)
  sig.plot.data$group = factor(sig.plot.data$group,levels = c('O','S/C'),labels = c('Moderate','Severe'))
  pdf(file = "4-cd8T-mild-severe-clone-signature.pdf",width = 2.5,height = 1.5)
  p <- ggplot(sig.plot.data, aes(x=group, y=sigscore,color=group)) + geom_violin(trim=FALSE)
  p = p + theme_cowplot() + theme(axis.text = element_text(size = 5),axis.title = element_text(size = 5),
                                  legend.text = element_text(size = 5),legend.title = element_text(size = 5),legend.key.size = unit(3,'mm')) +
    scale_color_manual(values=c("#d53d49", "#418bbe"))  + stat_summary(fun.data=mean_sdl, mult=1, geom="pointrange", color="grey70") + labs(x='',y='Tissue-resident core signatures score',color = "Group") 
  print(p)
  dev.off()
  
  ###t test
  expanded = sig.plot.data %>% filter(.,group == 'Moderate')
  expanded = expanded$sigscore
  nonexpanded = sig.plot.data %>% filter(.,group == 'Severe')
  nonexpanded = nonexpanded$sigscore
  t.test(expanded,nonexpanded)
  cohensD(expanded,nonexpanded)
  ###p-value = 4.567e-16
  ###Cohen’s d 0.3288864
  
  
  #####gliph
  library(dplyr)
  setwd('/Users/yang/Downloads/work/nCoV/scrna_rev/tcr/gliph/')
  vdjdb <- read.delim2("gliph.csv",header = TRUE, stringsAsFactors = FALSE,check.names = FALSE, sep = "\t")
  gliph <- read.delim2("gliph.csv-convergence-groups.txt",header = FALSE, stringsAsFactors = FALSE,check.names = FALSE, sep = "\t")
  gliph$V1 = as.numeric(gliph$V1)
  gliph = gliph %>% filter(.,V1 > 3)
  gliph = gliph[order(-gliph$V1),]
  
}





