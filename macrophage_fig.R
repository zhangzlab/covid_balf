macrophage_fig <- function(){
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(ggplot2)
  library(reshape2)
  setwd('/home/data/results/workspace/nCoV_rev/balf/Myeloid')
  Macrophage.Integrated = readRDS("2-Macrophage.rds")
  dpi = 300
  
  markers = c('AGER','SFTPC','SCGB3A2','TPPP3','KRT5',
              'CD68','FCN1','CD1C','TPSB2','CD14','MARCO','CXCR2',
              'CLEC9A','IL3RA',
              'CD3D','CD8A','KLRF1',
              'CD79A','IGHG4','MS4A1',
              'VWF','DCN',
              'FCGR3A','TREM2','KRT18')
  png(file="marker/violin_marker.png", width = dpi*30, height = dpi*24, units = "px",res = dpi,type='cairo')
  print(VlnPlot(object = Macrophage.Integrated, features = markers,pt.size = 0)+ theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5)))
  dev.off()
  png(file="marker/umap_marker.png", width = dpi*30, height = dpi*36, units = "px",res = dpi,type='cairo')
  print(FeaturePlot(object = Macrophage.Integrated, features = markers,cols = c("lightgrey","#ff0000")))
  dev.off()
  for(marker in markers){
    png(file=paste("marker/violin_",marker,".png",sep=''), width = dpi*8, height = dpi*3, units = "px",res = dpi,type='cairo')
    print(VlnPlot(object = Macrophage.Integrated, features = marker,pt.size = 0)+ theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5)))
    dev.off()
    png(file=paste("marker/umap_",marker,".png",sep=''), width = dpi*6, height = dpi*4, units = "px",res = dpi,type='cairo')
    print(FeaturePlot(object = Macrophage.Integrated, features = marker,cols = c("lightgrey","#ff0000")))
    dev.off()
  }
  
  
  library(ggplot2)
  pdf(file="marker_heatmap.pdf", width = 10, height = 8)
  pp = DotPlot(Macrophage.Integrated, features = rev(markers),cols = c('white','#F8766D'),dot.scale =5) + RotatedAxis()
  pp = pp + theme(axis.text.x = element_text(size = 12),axis.text.y = element_text(size = 12)) + labs(x='',y='') + 
    guides(color = guide_colorbar(title = 'Scale expression'),size = guide_legend(title = 'Percent expressed')) + 
    theme(axis.line = element_line(size = 0.6))
  print(pp)
  dev.off()
  
  ###
  Macrophage.Integrated1 = Macrophage.Integrated
  Macrophage.Integrated1$celltype = Idents(Macrophage.Integrated1)
  macrophage_groups = c('1','0','8','9','4','13','5','11','15','14','17','2','3','6','7','10','16')
  Macrophage.Integrated1 = subset(Macrophage.Integrated1,subset = celltype %in% macrophage_groups)
  Macrophage.Integrated1$celltype = factor(Macrophage.Integrated1$celltype,levels = macrophage_groups,labels = macrophage_groups)
  Idents(Macrophage.Integrated1) = Macrophage.Integrated1$celltype
  
  ###
  markers = c('FCN1','SPP1','FABP4')
  pdf(file="2-violin_marker1.pdf", width = 6, height = 9)
  pp_temp = VlnPlot(object = Macrophage.Integrated1, ncol = 1, features = markers,pt.size = 0,combine = FALSE)
  plots <- lapply(X = pp_temp, FUN = function(p) p + labs(x='') + theme(axis.text = element_text(size = 18),axis.text.x = element_text(size = 18,angle = 0,hjust = 0.5),axis.title = element_text(size = 18),plot.title = element_text(size = 18,family = 'sans',face = 'italic'),
                                                                        axis.line = element_line(size = 1),axis.ticks = element_line(size =1)))
  pp = CombinePlots(plots = plots,ncol = 1,legend = 'none')
  print(pp)  
  dev.off()
  
  dpi = 300
  markers = c('SPP1','FCN1','IL1B','MAFB','FABP4','MARCO','INHBA','TREM2')
  png(file="2-umap_marker.png", width = dpi*13, height = dpi*6, units = "px",res = dpi,type='cairo')
  pp_temp = FeaturePlot(object = Macrophage.Integrated1, ncol = 4, features = markers,cols = c("lightgrey","#ff0000"),combine = FALSE,pt.size = 0.8)
  plots <- lapply(X = pp_temp, FUN = function(p) p + theme(plot.title = element_text(family = 'sans',face='italic',size=16),legend.position = 'right') +
                    theme(axis.line = element_line(size = 0.8),axis.ticks = element_line(size = 0.8),
                          axis.text = element_text(size = 16),axis.title = element_text(size = 16),
                          legend.text = element_text(size =16)))
  pp = CombinePlots(plots = plots,ncol = 4,legend = 'right')
  print(pp)
  dev.off()
  
  ###
  DefaultAssay(disease_macrophage) = 'group'
  disease_macrophage = subset(Macrophage.Integrated,group != 'HC')
  Idents(disease_macrophage) = 'group'
  disease_macrophage@misc$markers <- FindAllMarkers(object = disease_macrophage, assay = 'RNA',only.pos = TRUE, test.use = 'MAST')
  write.table(disease_macrophage@misc$markers,file='2-macrophage_severe_vs_mild.txt',row.names = FALSE,quote = FALSE,sep = '\t')
  
  ###heatmap
  DefaultAssay(Macrophage.Integrated) <- "RNA"
  big.averageExpression1 = AverageExpression(object = Macrophage.Integrated,slot = "scale.data")
  big.averageExpression = big.averageExpression1$RNA
  big.averageExpression = big.averageExpression %>% select(one_of(nCoV_groups))
  big.averageExpression$gene = rownames(big.averageExpression)
  write.table(big.averageExpression,file = '2-nCoV-average-RNA.txt',row.names = FALSE,quote = FALSE,sep = '\t')
  
  rownames(big.averageExpression) = big.averageExpression$gene
  big.averageExpression$gene = NULL
  big.averageExpression = big.averageExpression %>% dplyr::select(one_of(nCoV_groups))
  big.averageExpression[,] = sapply(big.averageExpression[,], as.numeric)
  big.averageExpression$gene = rownames(big.averageExpression)
  nCoV_markers = c('FCN1','TREM2','SPP1','FABP4','INHBA','SIGLEC10','IL1B','MARCO','CD1C','LAMP3','MAFB')
  nCoV_marker_left = nCoV_markers
  big.averageExpression = big.averageExpression %>% filter(.,gene %in% nCoV_markers)
  big.averageExpression = na.omit(big.averageExpression[match(nCoV_markers, big.averageExpression$gene),])
  rownames(big.averageExpression) = big.averageExpression$gene
  
  ###sort FCN1/###sort TREM2
  big.averageExpression.tmp = big.averageExpression
  big.averageExpression.tmp$gene = NULL
  big.averageExpression.tmp = as.data.frame(t(big.averageExpression.tmp))
  big.averageExpression.tmp = big.averageExpression.tmp[order(-big.averageExpression.tmp$FCN1),]
  big.averageExpression = as.data.frame(t(big.averageExpression.tmp))
  big.averageExpression$gene = rownames(big.averageExpression)
  
  
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(ggplot2)
  library(reshape2)
  setwd('/home/data/results/workspace/nCoV_rev/balf/Figure')
  Macrophage.Integrated = readRDS("../Myeloid/2-Macrophage.rds")
  dpi = 300
  ###umap of macrophage
  png(file="2-umap-1.png", width = dpi*8, height = dpi*6, units = "px",res = dpi,type='cairo')
  pp = DimPlot(object = Macrophage.Integrated, reduction = 'umap',label = FALSE,pt.size = 0.8,label.size = 6,repel = TRUE)
  pp = pp + theme(axis.title = element_text(size = 16),axis.text =  element_text(size = 16),
                  legend.text = element_text(size = 16),axis.line = element_line(size = 1))
  print(pp)
  dev.off()
  
  ###group
  Macrophage.Integrated1 = Macrophage.Integrated
  #Group 1 includes C1 (FCN1hi only), 
  #Group 2 includes C13, C4, C8, C9, C0 (FCN1+SPP1+)
  #Group 3 includes C15, C5 and C11 (FCN1-SPP1+)
  #Group 4 includes C2, C3, C6, C7, C10, C16, C14, C17 (FABP4+)
  #Doublets C12 C18 C19
  Macrophage.Integrated1 <- RenameIdents(object = Macrophage.Integrated1, 
                                         '1' = 'Group1',   
                                         '13' = 'Group2','4' = 'Group2','8' = 'Group2','9' = 'Group2','0' = 'Group2',
                                         '15' = 'Group3','5' = 'Group3','11' = 'Group3',
                                         '2' = 'Group4','3' = 'Group4','6' = 'Group4','7' = 'Group4','10' = 'Group4',
                                         '16' = 'Group4','14' = 'Group4','17' = 'Group4',
                                         '12' = 'Doublets','18' = 'Doublets','19' = 'Doublets')
  Macrophage.Integrated1$celltype = Idents(Macrophage.Integrated1)
  Macrophage.Integrated1 = subset(Macrophage.Integrated1,subset = celltype!='Doublets')
  nCoV_groups = c('Group1','Group2','Group3','Group4')
  Macrophage.Integrated1$celltype = factor(Macrophage.Integrated1$celltype,levels = nCoV_groups,labels = nCoV_groups)
  Idents(Macrophage.Integrated1) = Macrophage.Integrated1$celltype
  
  ###calculate deg in four myeloid groups
  for(nCoV_group in nCoV_groups){
    Macrophage.Integrated1_i = subset(Macrophage.Integrated1,group != 'HC' & celltype == nCoV_group)
    DefaultAssay(Macrophage.Integrated1_i) = 'RNA'
    Idents(Macrophage.Integrated1_i) = 'group'
    Macrophage.Integrated1_i@misc$markers <- FindAllMarkers(object = Macrophage.Integrated1_i, assay = 'RNA',only.pos = TRUE, test.use = 'MAST')
    write.table(Macrophage.Integrated1_i@misc$markers,file=paste('2-macrophage_',nCoV_group,'_severe_vs_mild.txt',sep=''),row.names = FALSE,quote = FALSE,sep = '\t')
  }
  
  ###HC/Mild/Severe
  png(file="2-nCoV-umap-split-group.png", width = dpi*11, height = dpi*4.3, units = "px",res = dpi,type='cairo')
  pp_temp = DimPlot(object = Macrophage.Integrated1, reduction = 'umap',label = FALSE, pt.size = 0.8,label.size = 5,split.by = 'group', ncol = 3,repel = TRUE,combine = TRUE,
                    cols = c('#C680FC','#00bf7d','#00b0f6','#f8766d'))
  pp_temp = pp_temp + theme(axis.title = element_text(size = 20),axis.text = element_text(size = 20),
                            strip.text = element_text(family = 'sans',face='plain',size=20),
                            legend.text = element_text(size = 20),
                            legend.key.height = unit(1.5,'line'),
                            axis.line = element_line(size = 1.2),axis.ticks = element_line(size = 1.2))
  print(pp_temp)
  dev.off()
  
  FCN1_p = subset(Macrophage.Integrated1,subset = FCN1 > 0 & SPP1 == 0 & FABP4==0,slot = 'counts')
  SPP1_p = subset(Macrophage.Integrated1,subset = FCN1 == 0 & SPP1 > 0 & FABP4==0,slot = 'counts')
  FCN1_SPP1_p = subset(Macrophage.Integrated1,subset = FCN1 > 0 & SPP1 > 0 & FABP4==0,slot = 'counts')
  FABP4_p = subset(Macrophage.Integrated1,subset = FABP4 > 0 & FCN1==0 & SPP1==0,slot = 'counts')
  
  sample_info = as.data.frame(colnames(Macrophage.Integrated1))
  colnames(sample_info) = c('ID')
  rownames(sample_info) = sample_info$ID
  sample_info = sample_info%>%mutate(MARKER=ifelse(ID %in% colnames(FCN1_p),'High FCN1',
                                                   ifelse(ID %in% colnames(SPP1_p),'High SPP1',
                                                          ifelse(ID %in% colnames(FABP4_p),'High FABP4',
                                                                 ifelse(ID %in% colnames(FCN1_SPP1_p),'FCN1/SPP1 high both','Other')))))
  rownames(sample_info) = sample_info$ID
  mac_groups = c('High FCN1','FCN1/SPP1 high both','High SPP1','High FABP4','Other')
  sample_info$MARKER = factor(sample_info$MARKER,levels = mac_groups ,labels = mac_groups)
  Macrophage.Integrated1 = AddMetaData(object = Macrophage.Integrated1, metadata = sample_info)
  
  png(file="2-umap-FCN1_SPP1_FABP4.png", width = dpi*13, height = dpi*4.3, units = "px",res = dpi,type='cairo')
  pp = DimPlot(object = Macrophage.Integrated1, reduction = 'umap',label = FALSE, pt.size = 0.8,group.by = 'MARKER',cols = c('#C680FC','#00bf7d','#00b0f6','#f8766d','lightgrey'),
               split.by = 'group')
  pp= pp + theme(axis.title = element_text(size = 20),axis.text = element_text(size = 20),legend.text = element_text(size = 20),
                 legend.key.height = unit(1.5,'line'),strip.text = element_text(family = 'sans',face='plain',size=20),
                 axis.line = element_line(size = 1.2),axis.ticks = element_line(size = 1.2)
                 
  )
  print(pp)
  dev.off()
  
  
  ###violin plot
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(ggplot2)
  library(reshape2)
  setwd('/home/data/results/workspace/nCoV_rev/balf/Figure')
  Macrophage.Integrated = readRDS("../Myeloid/2-Macrophage.rds")
  Macrophage.Integrated2 = Macrophage.Integrated1
  markers = c('FCN1','SPP1','FABP4')
  pdf(file="2-violin_marker_FCN1_SPP1_FABP4.pdf", width = 6, height = 9)
  pp_temp = VlnPlot(object = Macrophage.Integrated2, ncol = 1, features = markers, pt.size = 0,
                    cols = c('#00b3f2','#00b3f2','#00b3f2','#00b3f2',
                             '#d277ff','#d277ff','#d277ff','#f8766d','#f8766d','#f8766d','#f8766d','#f8766d','#f8766d'),
                    combine = FALSE)
  plots <- lapply(X = pp_temp, FUN = function(p) p + labs(x='') + theme(plot.title = element_text(size = 20,family = 'sans',face = 'italic'),axis.text = element_text(size = 16)))
  pp = CombinePlots(plots = plots,ncol = 1,legend = 'none')
  print(pp)  
  dev.off()
  
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(ggplot2)
  library(reshape2)
  setwd('/home/data/results/workspace/nCoV_rev/balf/Figure')
  Macrophage.Integrated = readRDS("../Myeloid/2-Macrophage.rds")
  Macrophage.Integrated1 = Macrophage.Integrated
  #Group 1 includes C1 (FCN1hi only), 
  #Group 2 includes C13, C4, C8, C9, C0 (FCN1+SPP1+)
  #Group 3 includes C15, C5 and C11 (FCN1-SPP1+)
  #Group 4 includes C2, C3, C6, C7, C10, C16, C14, C17 (FABP4+)
  #Doublets C12 C18 C19
  Macrophage.Integrated1 <- RenameIdents(object = Macrophage.Integrated1, 
                                         '1' = 'Group1',   
                                         '13' = 'Group2','4' = 'Group2','8' = 'Group2','9' = 'Group2','0' = 'Group2',
                                         '15' = 'Group3','5' = 'Group3','11' = 'Group3',
                                         '2' = 'Group4','3' = 'Group4','6' = 'Group4','7' = 'Group4','10' = 'Group4',
                                         '16' = 'Group4','14' = 'Group4','17' = 'Group4',
                                         '12' = 'Doublets','18' = 'Doublets','19' = 'Doublets')
  Macrophage.Integrated1$celltype = Idents(Macrophage.Integrated1)
  Macrophage.Integrated1 = subset(Macrophage.Integrated1,subset = celltype!='Doublets')
  nCoV_groups = c('Group1','Group2','Group3','Group4')
  Macrophage.Integrated1$celltype = factor(Macrophage.Integrated1$celltype,levels = nCoV_groups,labels = nCoV_groups)
  Idents(Macrophage.Integrated1) = Macrophage.Integrated1$celltype
  
  DefaultAssay(Macrophage.Integrated1) <- "RNA"
  for(b_group in nCoV_groups){
    print(b_group)
    B.Integrated.subtype = Macrophage.Integrated1
    B.Integrated.subtype$celltype = as.character(Idents(B.Integrated.subtype))
    B.Integrated.subtype$celltype[B.Integrated.subtype$celltype != b_group] = 'other'
    B.Integrated.subtype$celltype = factor(B.Integrated.subtype$celltype,levels = c(b_group,'other'),labels = c(b_group,'other'))
    Idents(B.Integrated.subtype) = 'celltype'
    ###PART1 logfc
    B_gsea <- FindMarkers(object = B.Integrated.subtype,ident.1 = b_group,test.use = 'MAST',logfc.threshold = 0,min.pct = 0)
    write.table(B_gsea,file = paste('3-macrophage_gsea_',b_group,'_logfc.txt',sep = ''),sep='\t',quote = FALSE)
    ###PART2 signal to noise
    Expression = B.Integrated.subtype[['RNA']]@data
    metadata = B.Integrated.subtype@meta.data
    lc.cell = metadata%>% filter(.,celltype==b_group)
    hc.cell = metadata%>% filter(.,celltype=='other')
    lc.exp = Expression[,lc.cell$ID]
    hc.exp = Expression[,hc.cell$ID]
    u.exp = (rowMeans(lc.exp) - rowMeans(hc.exp))
    sd.exp = apply(lc.exp, 1, sd) + apply(hc.exp, 1, sd)  ##4844 0
    sintonoise = u.exp/sd.exp
    sintonoise = na.omit(sintonoise)
    sintonoise.data = as.data.frame(sintonoise)
    write.table(sintonoise.data,file = paste('3-macrophage_gsea_',b_group,'_sig2noise.txt',sep = ''),sep='\t',quote = FALSE)
  }
  
  library(dplyr)
  setwd('/Users/yang/Downloads/work/nCoV/scrna_rev/balf/Figure/GSEA/data')
  file_names<- list.files(".")
  for (i in 1:length(file_names)) {
    name<-gsub(".txt","",file_names[i])
    print(name)
    my_file = read.table(file_names[i],sep = "\t",header = TRUE,quote = "",stringsAsFactors = FALSE)
    name = gsub(" ","",name)
    name = gsub("-","_",name)
    if(dim(my_file)[2] == 6){  ###logfc
      my_file = my_file %>% dplyr::select(one_of(c('gene','avg_logFC')))
    }else{  ##signoise
    }
    colnames(my_file) = c('gene',name)
    my_file[,2] = as.numeric(my_file[,2])
    my_file = my_file[order(-my_file[,2]),]
    my_file = my_file[(my_file[,2] > 0.01 | my_file[,2] < -0.01),]
    write.table(my_file,file = paste('../rnk/',name,'.rnk',sep=''),sep = '\t',quote = FALSE,row.names = FALSE)
  }
  
  ###analysis with clusterprofier
  library(clusterProfiler)
  library(enrichplot)
  library(dplyr)
  setwd('/Users/yang/Downloads/work/nCoV/scrna_rev/balf/Figure/GSEA/')
  ###load gene set
  hallmark <- read.gmt('h.all.v7.0.symbols.gmt')
  file_names<- list.files("rnk")
  for (i in 1:length(file_names)) {
    name<-gsub(".rnk","",file_names[i])
    print(name)
    rnkList = read.table(paste('rnk/',file_names[i],sep=''),sep = "\t",header = TRUE,quote = "",stringsAsFactors = FALSE)
    colnames(rnkList) = c('gene','value')
    rnkList$value = as.numeric(rnkList$value)
    rnkList = arrange(rnkList,desc(value))
    glist <- rnkList[,2]
    names(glist) <- as.character(rnkList[,1])
    glist <- sort(glist,decreasing = T)
    gsea <- GSEA(geneList = glist, TERM2GENE=hallmark, verbose=FALSE, pvalueCutoff = 0.05,nPerm = 100000)
    write.table(gsea@result,file = paste('result/',name,'.txt',sep=''),sep='\t',quote = FALSE)
    if(length(gsea@result$ID) > 0 ){
      dpi = 300
      png(file = paste('result/',name,'.png',sep=''),width = dpi * 10,height = dpi * 6,units = "px",res = dpi,type = 'cairo')
      print(gseaplot2(gsea,geneSetID=gsea@result$ID,base_size = 10))
      dev.off()
    }
  }
  
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(ggplot2)
  library(reshape2)
  setwd('/home/data/results/workspace/nCoV_rev/balf/Figure')
  Macrophage.Integrated = readRDS("../Myeloid/2-Macrophage.rds")
  Macrophage.Integrated1 = Macrophage.Integrated
  #Group 1 includes C1 (FCN1hi only), 
  #Group 2 includes C13, C4, C8, C9, C0 (FCN1+SPP1+)
  #Group 3 includes C15, C5 and C11 (FCN1-SPP1+)
  #Group 4 includes C2, C3, C6, C7, C10, C16, C14, C17 (FABP4+)
  #Doublets C12 C18 C19
  Macrophage.Integrated1 <- RenameIdents(object = Macrophage.Integrated1, 
                                         '1' = 'Group1',   
                                         '13' = 'Group2','4' = 'Group2','8' = 'Group2','9' = 'Group2','0' = 'Group2',
                                         '15' = 'Group3','5' = 'Group3','11' = 'Group3',
                                         '2' = 'Group4','3' = 'Group4','6' = 'Group4','7' = 'Group4','10' = 'Group4',
                                         '16' = 'Group4','14' = 'Group4','17' = 'Group4',
                                         '12' = 'Doublets','18' = 'Doublets','19' = 'Doublets')
  Macrophage.Integrated1$celltype = Idents(Macrophage.Integrated1)
  Macrophage.Integrated1 = subset(Macrophage.Integrated1,subset = celltype!='Doublets')
  nCoV_groups = c('Group1','Group2','Group3','Group4')
  Macrophage.Integrated1$celltype = factor(Macrophage.Integrated1$celltype,levels = nCoV_groups,labels = nCoV_groups)
  Idents(Macrophage.Integrated1) = Macrophage.Integrated1$celltype
  DefaultAssay(Macrophage.Integrated1) <- "RNA"
  Macrophage.Integrated1@misc$markers <- FindAllMarkers(object = Macrophage.Integrated1, assay = 'RNA',only.pos = TRUE, test.use = 'MAST')
  write.table(Macrophage.Integrated1@misc$markers,file='3-markers_MAST.txt',row.names = FALSE,quote = FALSE,sep = '\t')
  hc.markers = read.delim2("3-markers_MAST.txt",header = TRUE, stringsAsFactors = FALSE,check.names = FALSE, sep = "\t")
  hc.markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_logFC) -> top50
  hvg = top50$gene
  var.genes = c(Macrophage.Integrated1@assays$RNA@var.features,top50$gene)
  Macrophage.Integrated1 <- ScaleData(Macrophage.Integrated1, verbose = FALSE, vars.to.regress = c("nCount_RNA", "percent.mito"),features = var.genes)
  tt1 = DoHeatmap(object = Macrophage.Integrated1, features = hvg,angle = 0,hjust = 0.5,size = 6) + NoLegend() 
  ggplot2::ggsave(file="3-feature2-1.pdf",plot = tt1,device = 'pdf',width = 12, height = 14, units = "in",dpi = dpi,limitsize = FALSE)
  
  ###now pathway analysis
  library(ggrepel)
  library(clusterProfiler)
  library(enrichplot)
  library(org.Hs.eg.db)
  setwd('/Users/yang/Downloads/work/nCoV/scrna_rev/balf/Figure/')
  hc.markers = read.delim2("3-markers_MAST.txt",header = TRUE, stringsAsFactors = FALSE,check.names = FALSE, sep = "\t")
  groups = unique(hc.markers$cluster)
  hc.markers$p_val_adj = as.numeric(hc.markers$p_val_adj )
  for(group_ in groups){
    hc.markers_ = hc.markers %>% filter(.,cluster==group_ & p_val_adj<0.05)
    gene_ = hc.markers_$gene
    B_fun.gene.1 <- bitr(gene_, fromType = "SYMBOL",toType = c("ENTREZID"),OrgDb = org.Hs.eg.db)
    B_fun.gene = B_fun.gene.1$ENTREZID
    B_fun_bp = enrichGO(gene = B_fun.gene,OrgDb = org.Hs.eg.db,ont = "BP",pAdjustMethod = "BH",pvalueCutoff  = 0.05,qvalueCutoff  = 0.05,readable= TRUE)
    write.table(B_fun_bp,file = paste('3-',group_,'_bp.txt',sep = ''),sep='\t',quote = FALSE)
    dpi = 300
    png(file = paste('3-',group_,'_bp.png',sep = ''),width = dpi * 10,height = dpi * 6,units = "px",res = dpi,type = 'cairo')
    print(barplot(B_fun_bp,showCategory=20,drop=T)+theme(axis.text.y = element_text(size = 16),legend.text = element_text(size = 16),legend.title = element_text(size = 16)))
    dev.off()
    B_fun_kegg <- enrichKEGG(gene = B_fun.gene,organism = 'hsa',pvalueCutoff = 0.05)
    write.table(B_fun_kegg,file = paste('3-',group_,'_kegg.txt',sep = ''),sep='\t',quote = FALSE)
    dpi = 300
    png(file = paste('3-',group_,'_kegg.png',sep = ''),width = dpi * 10,height = dpi * 6,units = "px",res = dpi,type = 'cairo')
    print(barplot(B_fun_kegg,showCategory=20,drop=T)+theme(axis.text.y = element_text(size = 16),legend.text = element_text(size = 16),legend.title = element_text(size = 16)))
    dev.off()
  }
  
  ###now union result
  for(group_ in groups){
    hc.markers_i = read.delim2(paste('3-',group_,'_bp_sy.txt',sep = ''),header = TRUE, stringsAsFactors = FALSE,check.names = FALSE, sep = "\t")
    hc.markers_i$p.adjust = as.numeric(hc.markers_i$p.adjust)
    hc.markers_i = hc.markers_i %>% filter(.,p.adjust < 0.05)
    hc.markers_i = hc.markers_i[1:10,]
    hc.markers_i = hc.markers_i %>% dplyr::select(ID,Description,p.adjust)
    hc.markers_i$p.adjust = -log10(hc.markers_i$p.adjust)
    hc.markers_i$Description = factor(hc.markers_i$Description,labels =hc.markers_i$Description,levels = hc.markers_i$Description)
    pdf(file = paste('3-',group_,'_bp_final.pdf',sep = ''),width = 5,height = 4)
    #png(file = paste('3-',group_,'_bp_final.png',sep = ''),width = dpi * 8,height = dpi * 4,units = "px",res = dpi,type = 'cairo')
    pp <- ggplot(data=hc.markers_i, aes(x=Description, y=p.adjust)) + geom_bar(stat="identity",width = 0.7,fill = "#FDA1A1") + coord_flip()
    pp = pp + theme_minimal() + labs(x='',y='-log10(p.adjust)') + scale_x_discrete(limits = rev(levels(hc.markers_i$Description)),position = 'right') + scale_y_continuous(expand = c(0,0))
    pp = pp + geom_text(data=hc.markers_i, aes(x=Description, y=max(hc.markers_i$p.adjust)/100),label=hc.markers_i$Description, size=5,hjust = 0)
    pp = pp + theme(axis.ticks.y = element_blank(),axis.text.y = element_blank(),
                    panel.grid = element_blank(),
                    axis.line.x = element_line(size = 0.5, linetype = "solid", colour = "darkgrey"),
                    axis.line.y = element_line(size = 0.5, linetype = "solid", colour = "black"),
                    axis.title.x = element_text(size = 16),axis.text.x = element_text(size = 14),
                    axis.ticks.x = element_line(linetype = "solid",colour = "darkgrey",size = 0.5),axis.ticks.length.x = unit(2,'mm'),
                    plot.title = element_text(hjust = 0.5,size = 16)) +
      ggtitle(group_)
    print(pp)
    dev.off()
  }
  
  ####GSEA analysis
  library(dplyr)
  library(reshape2)
  setwd('/Users/yang/Downloads/work/nCoV/scrna_rev/balf/Figure/GSEA/result')
  B_groups = c("Group1","Group2","Group3","Group4")
  hallmark_now = c()
  for(b_group in B_groups){
    filename = gsub(" ",'',paste('3_macrophage_gsea_',b_group,'_sig2noise.txt',sep = ''))
    b_i = read.table(filename,sep = "\t",header = TRUE,quote = "",stringsAsFactors = FALSE)
    b_i = b_i %>% dplyr::select(ID,NES)
    colnames(b_i) = c('ID',b_group)
    hallmark_now = union(hallmark_now,b_i$ID)
  }
  hallmark_now = as.data.frame(hallmark_now)
  colnames(hallmark_now) = c('ID')
  for(b_group in B_groups){
    filename = gsub(" ",'',paste('3_macrophage_gsea_',b_group,'_sig2noise.txt',sep = ''))
    print(filename)
    b_i = read.table(filename,sep = "\t",header = TRUE,quote = "",stringsAsFactors = FALSE)
    b_i = b_i %>% dplyr::select(ID,NES)
    colnames(b_i) = c('ID',b_group)
    hallmark_now = dplyr::left_join(hallmark_now,b_i)
  }
  library(gplots)
  library(reshape2)
  library(ggplot2)
  dpi = 300
  ###change hallmark_now ID to 
  hallmark_now$ID = gsub("HALLMARK_","",hallmark_now$ID)
  hallmark_now$ID = gsub("_"," ",hallmark_now$ID)
  write.table(hallmark_now,file='3-Macrophage_gsea.txt',row.names = FALSE,quote = FALSE,sep='\t')
  
  hallmark_now = read.delim2("3-Macrophage_gsea_small.txt",header = TRUE, stringsAsFactors = FALSE,check.names = FALSE, sep = "\t")
  hallmark_now$Group1 = as.numeric(hallmark_now$Group1)
  hallmark_now$Group2 = as.numeric(hallmark_now$Group2)
  hallmark_now$Group3 = as.numeric(hallmark_now$Group3)
  hallmark_now$Group4 = as.numeric(hallmark_now$Group4)
  my_data = melt(hallmark_now,id.vars = 'ID')
  
  
  my_data$ID = factor(my_data$ID,levels = hallmark_now$ID,labels = hallmark_now$ID)
  my_data$variable = factor(my_data$variable,levels = B_groups,labels = B_groups)
  my_data$value[is.na(my_data$value)] = 0
  my_data$value[my_data$value > 4] = 4
  my_data$value[my_data$value < -4] = -4
  pdf(file = '3-Macrophage_gsea.pdf',width = 8.5,height = 3.3)
  pp = ggplot(my_data,aes(x=ID,y=variable,fill=value)) +
    geom_tile(color = 'grey') + 
    labs(x="",y="")+
    scale_x_discrete(expand=c(0,0))+
    #scale_y_discrete(expand=c(0,0)) + 
    scale_y_discrete(expand=c(0,0),limits = (rev(levels(my_data$variable))),position = "right") + 
    scale_fill_gradientn(name="", limits=c(-4,4),
                         colours=colorRampPalette(c("#00a9ff","white", "#f8766d"))(n = 51))+
    theme(legend.position="top",legend.key.width=unit(0.5,"cm"),legend.key.height=unit(0.3,"cm"),
          axis.ticks=element_blank(),axis.text.x = element_text(angle = 315,size = 9,hjust = 0,vjust = 0.5),
          axis.text.y.right = element_text(size = 12)) + 
    theme(plot.margin=unit(c(1,3,0,0),'lines'))
  print(pp)
  dev.off()
  
  ###SCENIC using pyscenic docker to get AUC result
  library(SCENIC)
  library(AUCell)
  setwd('/home/data/results/workspace/nCoV/SCENIC')
  regulonAUC <- importAUCfromText(file.path("auc_mtx.csv"))
  regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]
  macro.cellinfo = read.delim2("/home/data/results/workspace/nCoV/SCENIC/3-macrophage-SCENIC-meta.csv",header = TRUE, stringsAsFactors = FALSE,check.names = FALSE, sep = "\t")
  rownames(macro.cellinfo) = macro.cellinfo$ID
  
  regulonActivity_byCellType <- sapply(split(rownames(macro.cellinfo), macro.cellinfo$celltype),
                                       function(cells) rowMeans(getAUC(regulonAUC)[,cells]))
  regulonActivity_byCellType_Scaled <- t(scale(t(regulonActivity_byCellType), center = T, scale=T))
  write.table(regulonActivity_byCellType_Scaled,file='3-SCENIC-heatmap.txt',row.names = TRUE,quote = FALSE,sep='\t')
  
  pdf(file="3-SCENIC-heatmap.pdf", width = 6, height = 13)
  pheatmap::pheatmap(regulonActivity_byCellType_Scaled, fontsize_row = 5, fontsize_col = 12, 
                     color=colorRampPalette(c("#00a9ff","white","#F8766D"))(100), breaks=seq(-3, 3, length.out = 100),
                     treeheight_row=10, treeheight_col=10, border_color=NA)
  dev.off()
  
  
}
