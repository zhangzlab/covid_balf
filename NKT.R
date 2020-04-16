NKT_fig <-function(){
  ##load macrophage
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(ggplot2)
  library(reshape2)
  setwd('/home/data/results/workspace/nCoV_rev/balf/Figure')
  Macrophage.Integrated = readRDS("../Myeloid/2-Macrophage.rds")
  Macrophage.Integrated1 = Macrophage.Integrated
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
  Macrophage.Integrated_severe = subset(Macrophage.Integrated1,subset = group =='S/C')
  macro_barcode_severe = colnames(Macrophage.Integrated_severe)
  Macrophage.Integrated_mild = subset(Macrophage.Integrated1,subset = group =='O')
  macro_barcode_mild = colnames(Macrophage.Integrated_mild)
  
  ##load NKT
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(ggplot2)
  library(reshape2)
  setwd('/home/data/results/workspace/nCoV_rev/balf/Figure')
  NKT.Integrated = readRDS("../NKT/6-NKT.rds")
  NKT.Integrated <- RenameIdents(object = NKT.Integrated, 
                                 '0' ='CD4 T','1'='CD8 T','2' = 'CD8 T','3'='Cycling T','4'='NK','5' ='CD8 T',
                                 '6'='Treg','7'='Doublets','8'='NK','9'='gdT & MAIT','10'='Cycling T','11'='Unidentified','12'='Unidentified',
                                 '13'='Doublets')
  NKT.Integrated$celltype = Idents(NKT.Integrated)
  nCoV_groups = c('CD4 T','CD8 T','Cycling T','NK','Treg','gdT & MAIT','Unidentified','Doublets')
  NKT.Integrated$celltype = factor(NKT.Integrated$celltype,levels = nCoV_groups,labels = nCoV_groups)
  Idents(NKT.Integrated) = NKT.Integrated$celltype
  NKT.Integrated = subset(NKT.Integrated,subset = celltype != 'Doublets' & celltype != 'Unidentified')
  NKT.Integrated_severe = subset(NKT.Integrated,subset = group == 'S/C')
  nkt_barcode_severe = colnames(NKT.Integrated_severe)
  NKT.Integrated_mild = subset(NKT.Integrated,subset = group == 'O')
  nkt_barcode_mild = colnames(NKT.Integrated_mild)
  
  ###load all cells
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(ggplot2)
  library(reshape2)
  setwd('/home/data/results/workspace/nCoV_rev/balf/Figure')
  nCoV.integrated = readRDS(file = "../nCoV.rds")
  nCoV.integrated1 <- RenameIdents(object = nCoV.integrated, 
                                   '13' = 'Epithelial','16' = 'Epithelial','25' = 'Epithelial','28' = 'Epithelial','31' = 'Epithelial',
                                   '0'='MP','1'='MP','2'='MP','3'='MP','4'='MP','5'='MP','7'='MP','8'='MP','10'='MP','11'='MP','12'='MP','18'='MP','21'='MP','22'='MP','23'='MP','26'='MP',
                                   '30'='Mast',
                                   '6'='T','9'='T','14'='T',
                                   '17'='NK',
                                   '19'='Plasma',
                                   '27'='B',
                                   '15'='Neutrophil',
                                   '20'='mDC',
                                   '29'='pDC',
                                   '24'='Doublets'
  )
  
  nCoV.integrated1$celltype = Idents(nCoV.integrated1)
  nCoV.integrated1 = subset(nCoV.integrated1,subset = celltype != 'Doublets')
  nCoV_groups = c('Epithelial','MP','Neutrophil','mDC','pDC','Mast','T','NK','B','Plasma')
  nCoV.integrated1$celltype = factor(nCoV.integrated1$celltype,levels = nCoV_groups,labels = nCoV_groups)
  Idents(nCoV.integrated1) = nCoV.integrated1$celltype
  
  ##
  nCoV.integrated_severe = subset(nCoV.integrated1,cells = c(macro_barcode_severe,nkt_barcode_severe))
  DefaultAssay(nCoV.integrated_severe) = 'RNA'
  Idents(nCoV.integrated_severe) = 'celltype'
  nCoV.integrated_severe = subset(nCoV.integrated_severe,celltype != 'NK')
  nCoV.integrated_severe@misc$markers <- FindAllMarkers(object = nCoV.integrated_severe, assay = 'RNA',only.pos = TRUE, test.use = 'MAST')
  write.table(nCoV.integrated_severe@misc$markers,file='4-macro_vs_nkt_severe.txt',row.names = FALSE,quote = FALSE,sep = '\t')
  
  nCoV.integrated_mild = subset(nCoV.integrated1,cells = c(macro_barcode_mild,nkt_barcode_mild))
  DefaultAssay(nCoV.integrated_mild) = 'RNA'
  Idents(nCoV.integrated_mild) = 'celltype'
  nCoV.integrated_mild = subset(nCoV.integrated_mild,celltype != 'NK')
  nCoV.integrated_mild@misc$markers <- FindAllMarkers(object = nCoV.integrated_mild, assay = 'RNA',only.pos = TRUE, test.use = 'MAST')
  write.table(nCoV.integrated_mild@misc$markers,file='4-macro_vs_nkt_mild.txt',row.names = FALSE,quote = FALSE,sep = '\t')
  
  ###NKT analysis
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(ggplot2)
  library(reshape2)
  setwd('/home/data/results/workspace/nCoV_rev/balf/Figure')
  NKT.Integrated = readRDS("../NKT/6-NKT.rds")
  
  markers = c('TPPP3','KRT18','CD68','FCGR3B','CD1C','CLEC9A','LILRA4','TPSB2','CD3D','KLRD1','MS4A1','IGHG4')
  pdf(file="NKT-marker_heatmap.pdf", width = 8, height = 7)
  pp = DotPlot(NKT.Integrated, features = rev(markers),cols = c('white','#F8766D'),dot.scale =5) + RotatedAxis()
  pp = pp + theme(axis.text.x = element_text(size = 12),axis.text.y = element_text(size = 12)) + labs(x='',y='') + 
    guides(color = guide_colorbar(title = 'Scale expression'),size = guide_legend(title = 'Percent expressed')) + 
    theme(axis.line = element_line(size = 0.6))
  print(pp)
  dev.off()
  
  NKT.Integrated <- RenameIdents(object = NKT.Integrated, 
                                 '0' ='CCR7+ T','1'='CD8 T','2' = 'CD8 T','3'='Proliferating T','4'='NK','5' ='CD8 T',
                                 '6'='Treg','7'='Doublets','8'='NK','9'='innate T','10'='Proliferating T','11'='Uncertain','12'='Uncertain',
                                 '13'='Doublets')
  NKT.Integrated$celltype = Idents(NKT.Integrated)
  nCoV_groups = c('CCR7+ T','CD8 T','Proliferating T','NK','Treg','innate T','Uncertain','Doublets')
  NKT.Integrated$celltype = factor(NKT.Integrated$celltype,levels = nCoV_groups,labels = nCoV_groups)
  Idents(NKT.Integrated) = NKT.Integrated$celltype
  
  NKT.Integrated11 = subset(NKT.Integrated,subset = celltype == 'NK' & group != 'HC')
  write.table(NKT.Integrated11@meta.data,file='nk-matadata.txt',row.names = FALSE,quote = FALSE,sep='\t')
  
  cytokines = c('CD33','CD14','CXCR6','CD163','IL3RA','CD27','CD19','CCR6','NCAM1','KLRB1','CD69','CD68','CXCR3','PDCD1','CCR7','ITGAE','CX3CR1','CD4','FCGR3A','IL7R','ITGAM','ITGAX','CD38','PTPRC','B3GAT1','CD3D','CD8A')
  pdf(file="4-heatmap-cytokine_cytof_sample.pdf", width = 10, height = 4)
  pp = DotPlot(NKT.Integrated, features = rev(cytokines),cols = c('white','#F8766D'),dot.scale =5) + RotatedAxis()
  pp = pp + theme(axis.text.x = element_text(size = 9),axis.text.y = element_text(size = 9)) + labs(x='',y='') + 
    guides(color = guide_colorbar(title = 'Scale expression'),size = guide_legend(title = 'Percent expressed')) + 
    theme(axis.line = element_line(size = 1))
  print(pp)
  dev.off()
  
  dpi = 300
  png(file="4-umap_1.png", width = dpi*6, height = dpi*4, units = "px",res = dpi,type='cairo')
  pp = DimPlot(object = NKT.Integrated, reduction = 'umap',label = FALSE,pt.size = 0.8,label.size = 4,repel = TRUE)
  pp = pp + theme(axis.title = element_text(size = 12),axis.text =  element_text(size = 12),
                  legend.text = element_text(size = 12),legend.key.height = unit(1.05,'line'),
                  axis.line = element_line(size = 0.8),axis.ticks = element_line(size = 0.8))
  print(pp)
  dev.off()
  
  ###HC/Mild/Severe
  NKT.Integrated1 = subset(NKT.Integrated,idents=c('CCR7+ T','CD8 T','Proliferating T','NK','Treg','innate T'))
  cols = c('#f8766d','#cd9600','#7cae00','#00be67','#00bfc4','#00a9ff')
  png(file="4-NKT-umap-split-group.png", width = dpi*11, height = dpi*4.3, units = "px",res = dpi,type='cairo')
  pp_temp = DimPlot(object = NKT.Integrated1, reduction = 'umap',label = FALSE, split.by = 'group', ncol = 3,repel = TRUE,combine = TRUE,pt.size = 1.2)
  pp_temp = pp_temp + theme(axis.title = element_text(size = 17),axis.text = element_text(size = 17),
                            strip.text = element_text(family = 'arial',face='bold',size=17),
                            legend.text = element_text(size = 17),legend.key.height = unit(1.5,'line'),
                            axis.line = element_line(size = 1.1),axis.ticks = element_line(size = 1.1)) + scale_colour_manual(values = cols)
  print(pp_temp)
  dev.off()
  
  
  NKT.Integrated1 = subset(NKT.Integrated,idents=c('CCR7+ T','CD8 T','Cycling T','NK','Treg','innate T'))
  DefaultAssay(NKT.Integrated1) <- "RNA"
  NKT.Integrated1@misc$markers <- FindAllMarkers(object = NKT.Integrated1, assay = 'RNA',only.pos = TRUE, test.use = 'MAST')
  write.table(NKT.Integrated1@misc$markers,file='4-markers_MAST.txt',row.names = FALSE,quote = FALSE,sep = '\t')
  
  ####percentage
  library(ggplot2)
  NKT.Integrated1[["cluster"]] <- Idents(object = NKT.Integrated1)
  big.cluster = NKT.Integrated1@meta.data
  organ.summary = table(big.cluster$sample_new,big.cluster$cluster)
  write.table(organ.summary,file = '4-NKT-percentage-sample.txt',quote = FALSE,sep = '\t')
  
  library(ggplot2)
  library(dplyr)
  library(reshape2)
  setwd('/Users/yang/Downloads/work/nCoV/scrna_rev/balf/Figure/')
  organ.summary = read.delim2("4-NKT-percentage-sample.txt",header = TRUE, stringsAsFactors = FALSE,check.names = FALSE, sep = "\t")
  organ.summary$group = rownames(organ.summary)
  organ.summary.dataframe = melt(organ.summary)
  colnames(organ.summary.dataframe) = c('group','cluster','cell')
  clusters = c('CCR7+ T','CD8 T','Proliferating T','NK','Treg','innate T')
  organ.summary.dataframe$cell = as.numeric(organ.summary.dataframe$cell)
  dpi = 300
  
  pdf(file="4-NKT-sample-percentage.pdf", width = 7, height = 3.6)
  # Stacked barplot with multiple groups
  pp = ggplot(data=organ.summary.dataframe, aes(x=group, y=cell, fill=cluster)) + geom_bar(stat="identity",width = 0.6,position=position_fill(reverse = FALSE),size = 0.5,colour = '#222222') + labs(x='',y='Fraction of cluster per sample (%)') + 
    theme(axis.title.x = element_blank(),axis.text = element_text(size = 14),axis.title.y =element_text(size = 14), legend.text = element_text(size = 14),legend.key.height = unit(5,'mm'),
          legend.title = element_blank(),panel.grid.minor = element_blank()) + cowplot::theme_cowplot() + theme(axis.text.x = element_text(angle = 45,size = 10,hjust = 0.5,vjust = 0.5))+
    scale_fill_manual(values = c('#f8766d','#cd9600','#7cae00','#00be67','#00bfc4','#00a9ff'))
  print(pp)
  dev.off()
  
  ####marker heatmap
  nCoV_markers_rev = rev(c('CD3D','IL7R','CCR7','GZMA','CD8A','CD8B','CXCR3','GZMK','MKI67','TYMS','NKG7','GZMB','KLRD1','KLRC1','XCL1','KLRF1','EOMES','CX3CR1',
                           'CD4','FOXP3','CTLA4','IL2RA','CXCR6','TRGV9','TRDV2','KLRB1','SLC4A10'))
  nCoV_groups_rev = rev(c('CCR7+ T','CD8 T','Cycling T','NK','Treg','innate T'))
  NKT.Integrated2 = NKT.Integrated
  NKT.Integrated2 = subset(NKT.Integrated2,idents = nCoV_groups_rev)
  NKT.Integrated2$celltype = factor(NKT.Integrated2$celltype,levels = nCoV_groups_rev,labels = nCoV_groups_rev)
  Idents(NKT.Integrated2) = NKT.Integrated2$celltype
  pdf(file="4-heatmap-seurat.pdf", width = 8.5, height = 3.5)
  pp = DotPlot(NKT.Integrated2, features = nCoV_markers_rev,cols = c('white','#F8766D'),dot.scale=5) + RotatedAxis()
  pp = pp + theme(axis.text.x = element_text(size = 10),axis.text.y = element_text(size = 12)) + labs(x='',y='') + 
    guides(color = guide_colorbar(title = 'Scale expression'),size = guide_legend(title = 'Percent expressed')) + 
    theme(axis.line = element_line(size = 1))
  print(pp)
  dev.off()
  
  
  ###
  NKT.Integrated1 = subset(NKT.Integrated,idents=c('CCR7+ T','CD8 T','Cycling T','NK','Treg','innate T'))
  hc.markers = read.delim2("4-markers_MAST.txt",header = TRUE, stringsAsFactors = FALSE,check.names = FALSE, sep = "\t")
  hc.markers %>% group_by(cluster) %>% top_n(n = 40, wt = avg_logFC) -> top50
  hvg = top50$gene
  var.genes = c(NKT.Integrated1@assays$RNA@var.features,top50$gene)
  NKT.Integrated1 <- ScaleData(NKT.Integrated1, verbose = FALSE, vars.to.regress = c("nCount_RNA", "percent.mito"),features = var.genes)
  dpi = 300
  tt1 = DoHeatmap(object = NKT.Integrated1, features = hvg,angle = 0,hjust = 0.5,size = 6) + NoLegend()
  ggplot2::ggsave(file="4-feature2-1.pdf",plot = tt1,device = 'pdf',width = 20, height = 24, units = "in",dpi = dpi,limitsize = FALSE)
  
  
  ###
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(ggplot2)
  library(reshape2)
  setwd('/home/data/results/workspace/nCoV_rev/balf/Figure')
  NKT.Integrated = readRDS("../NKT/6-NKT.rds")
  NKT.Integrated <- RenameIdents(object = NKT.Integrated, 
                                 '0' ='CCR7+ T','1'='CD8 T','2' = 'CD8 T','3'='Cycling T','4'='NK','5' ='CD8 T',
                                 '6'='Treg','7'='Doublets','8'='NK','9'='innate T','10'='Cycling T','11'='Uncertain','12'='Uncertain',
                                 '13'='Doublets')
  NKT.Integrated$celltype = Idents(NKT.Integrated)
  nCoV_groups = c('CCR7+ T','CD8 T','Cycling T','NK','Treg','innate T','Uncertain','Doublets')
  NKT.Integrated$celltype = factor(NKT.Integrated$celltype,levels = nCoV_groups,labels = nCoV_groups)
  Idents(NKT.Integrated) = NKT.Integrated$celltype
  NKT.Integrated1 = subset(NKT.Integrated,idents=c('CCR7+ T','CD8 T','Cycling T','NK','Treg','innate T'))
  NKT.Integrated1 = subset(NKT.Integrated1,subset = group != 'HC')
  NKT.Integrated.Mild = subset(NKT.Integrated,subset = group == 'O')
  NKT.Integrated.Severe = subset(NKT.Integrated,subset = group == 'C/S')
  NKT.Integrated.Mild <- FindVariableFeatures(object = NKT.Integrated.Mild, selection.method = "vst", nfeatures = 1000,verbose = FALSE)
  NKT.Integrated.Severe <- FindVariableFeatures(object = NKT.Integrated.Severe, selection.method = "vst", nfeatures = 1000,verbose = FALSE)
  hvg_Mild = NKT.Integrated.Mild@assays$RNA@var.features
  hvg_severe = NKT.Integrated.Severe@assays$RNA@var.features
  hvg = unique(c(hvg_Mild,hvg_severe))
  write.table(hvg,file = '4-hvg_1000.txt',row.names = FALSE,quote = FALSE)
  
  nkt_groups = c('CCR7+ T','CD8 T','Cycling T','NK','Treg','innate T')
  for(group_ in nkt_groups){
    NKT.Integrated.i = subset(NKT.Integrated1,subset = celltype == group_)
    DefaultAssay(NKT.Integrated.i) <- "RNA"
    Idents(NKT.Integrated.i) = 'group'
    NKT.Integrated.i@misc$markers <- FindAllMarkers(object = NKT.Integrated.i, assay = 'RNA',only.pos = TRUE, test.use = 'MAST',logfc.threshold = 0)
    write.table(NKT.Integrated.i@misc$markers,file=paste('4-NKT-all-',group_,'.txt',sep=''),row.names = FALSE,quote = FALSE,sep = '\t')
    NKT.Integrated.i@misc$markers <- FindAllMarkers(object = NKT.Integrated.i, assay = 'RNA',only.pos = TRUE, test.use = 'MAST')
    write.table(NKT.Integrated.i@misc$markers,file=paste('4-NKT-',group_,'.txt',sep=''),row.names = FALSE,quote = FALSE,sep = '\t')
  }
  
  ###add a latent.var
  nkt_groups = c('CD8 T')
  for(group_ in nkt_groups){
    NKT.Integrated.i = subset(NKT.Integrated1,subset = celltype == group_)
    DefaultAssay(NKT.Integrated.i) <- "RNA"
    Idents(NKT.Integrated.i) = 'group'
    NKT.Integrated.i@misc$markers <- FindAllMarkers(object = NKT.Integrated.i, assay = 'RNA',only.pos = TRUE, test.use = 'MAST',latent.vars = 'orig.ident')
    write.table(NKT.Integrated.i@misc$markers,file=paste('4-NKT-latent-',group_,'.txt',sep=''),row.names = FALSE,quote = FALSE,sep = '\t')
  }
  
  
  library(ggrepel)
  library(ggplot2)
  library(dplyr)
  library(reshape2)
  setwd('/Users/yang/Downloads/work/nCoV/scrna_rev/balf/Figure/')
  macro_nkt_deg_severe = read.delim2("4-macro_vs_nkt_severe.txt",header = TRUE, stringsAsFactors = FALSE,check.names = FALSE, sep = "\t")
  macro_nkt_deg_mild = read.delim2("4-macro_vs_nkt_mild.txt",header = TRUE, stringsAsFactors = FALSE,check.names = FALSE, sep = "\t")
  macro_nkt_deg_severe$p_val_adj = as.numeric(macro_nkt_deg_severe$p_val_adj)
  macro_nkt_deg_severe = macro_nkt_deg_severe %>% filter(.,p_val_adj < 0.05 & cluster == 'T')
  macro_nkt_deg_mild$p_val_adj = as.numeric(macro_nkt_deg_mild$p_val_adj)
  macro_nkt_deg_mild = macro_nkt_deg_mild %>% filter(.,p_val_adj < 0.05 & cluster == 'T')
  
  #hvg_used = hvg_500
  nkt_groups = c('CCR7+ T','CD8 T')
  dpi = 300
  for(group_ in nkt_groups){
    nkt.markers <- read.delim2(paste('4-NKT-all-',group_,'.txt',sep=''),header = TRUE, stringsAsFactors = FALSE,check.names = FALSE, sep = "\t")
    nkt.markers$p_val_adj = as.numeric(nkt.markers$p_val_adj)
    nkt.markers$avg_logFC = as.numeric(nkt.markers$avg_logFC)
    nkt.markers = nkt.markers %>% mutate(avg_logFC_new = ifelse(cluster == 'S/C',avg_logFC,-avg_logFC))
    nkt.markers$color = '#BEBEBE'
    nkt.markers$group = 'Other'
    nkt.markers$label = nkt.markers$gene
    rownames(nkt.markers) = nkt.markers$gene
    group_7_up = nkt.markers %>% filter(., cluster == 'S/C' & p_val_adj < 0.05 )
    group_7_down = nkt.markers %>% filter(.,cluster == 'O' & p_val_adj < 0.05 )
    print(dim(group_7_up))
    print(dim(group_7_down))
    nkt.markers[intersect(group_7_up$gene,macro_nkt_deg_severe$gene),'color'] = '#d53d49'
    nkt.markers[intersect(group_7_up$gene,macro_nkt_deg_severe$gene),'group'] = 'Severe'
    nkt.markers[intersect(group_7_down$gene,macro_nkt_deg_mild$gene),'color'] = '#418bbe'
    nkt.markers[intersect(group_7_down$gene,macro_nkt_deg_mild$gene),'group'] = 'Moderate'
    
    nkt.markers$p_val_adj = -log10(nkt.markers$p_val_adj)
    top10 = nkt.markers %>% group_by(group) %>% top_n(n = 5, wt = p_val_adj)
    top10 = top10 %>% filter(.,group !='Other')
    my_gene_list = c('CD2','ZNF683','ITGA4','IKZF3','ITGB7','XCL1','XCL2','ID2','IRF1',
                     'CCR2','TOX','TNFSF14','CCL5','JAK1','GBP5','STAT1','LAT','HLA-E','CD3E',
                     'OASL','MX1','RPLP0','RPS15A','RPS26','RPS24','RPL7','RPL39','RPS13','EIF1',
                     'EIF1AY','DUT','TYMS','ATP5MC2','ATP5MG','UQCRB','LDHA',
                     'COX6C','NDUFS5','DNPH1','NDUFA12','HMGB2')
    my_not_gene_list = c()
    nkt.markers$group = factor(nkt.markers$group,labels = c('Other','Severe','Moderate'),levels = c('Other','Severe','Moderate'))
    my_not_gene_list = c()
    nkt.markers[!nkt.markers$gene %in% c(my_gene_list,top10$gene),'label'] = ''
    nkt.markers = with(nkt.markers, nkt.markers[order(group),])
    
    pdf(file = paste('4-NKT-volcano-',group_,'.pdf',sep=''),width = 2.8,height = 2.5)
    p <- ggplot(nkt.markers, aes(avg_logFC_new, p_val_adj,color = group,order = group)) + geom_point(size = 0.8)+
      labs(x='logFC',y='-Log10(p.adjust)')
    p <- p + theme(plot.title = element_text(hjust = 0.5)) + 
      theme(axis.text = element_text(size = 5),axis.title = element_text(size = 5),
            axis.title.x = element_text(margin = unit(c(4, 0, 0, 0), "mm")),
            axis.title.y = element_text(margin = unit(c(0, 4, 0, 0), "mm")),
            axis.line = element_line(colour = 'black',size = 0.4),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            panel.grid.major = element_line(colour = 'grey',size=0.4,linetype = 2),
            axis.ticks = element_line(colour = 'black',size = 0.4),
            legend.title = element_blank(),legend.text = element_text(size = 5,family = 'sans',face='plain'),
            legend.position = 'top',legend.key.size = unit(1.5, 'lines'))
    p = p + scale_color_manual(values = setNames(nkt.markers$color, nkt.markers$group))
    p = p + geom_text_repel(aes(label = nkt.markers$label),color= 'black',size = 1.8, box.padding = 0.25,point.padding = 0.5,segment.color = 'grey80')
    p = p + theme(legend.key=element_rect(fill=NA))
    print(p)
    dev.off()
  }
  
  ###GO/KEGG analysis
  library(ggrepel)
  library(clusterProfiler)
  library(enrichplot)
  library(org.Hs.eg.db)
  cc = c()
  for(group_ in nkt_groups){
    print(group_)
    nkt.markers = read.delim2(paste('4-NKT-',group_,'.txt',sep=''),header = TRUE, stringsAsFactors = FALSE,check.names = FALSE, sep = "\t")
    print(dim(nkt.markers))
    nkt.markers$p_val_adj = as.numeric(nkt.markers$p_val_adj)
    nkt.markers$avg_logFC = as.numeric(nkt.markers$avg_logFC)
    up.markers_ = nkt.markers %>% filter(.,cluster=='S/C' & p_val_adj < 0.05)
    down.markers_ = nkt.markers %>% filter(.,cluster=='O' & p_val_adj < 0.05)
    
    up.gene_ = up.markers_$gene
    down.gene_ = down.markers_$gene
    up.gene_ = intersect(macro_nkt_deg_severe$gene,up.gene_)
    down.gene_ = intersect(macro_nkt_deg_mild$gene,down.gene_)
    print(length(up.gene_))
    print(length(down.gene_))
    write.table(up.gene_,file=paste('4-',group_,'severe_gene.txt',sep=''),row.names = FALSE,quote = FALSE,sep='\t')
    write.table(down.gene_,file=paste('4-',group_,'mild_gene.txt',sep=''),row.names = FALSE,quote = FALSE,sep='\t')
    up.B_fun.gene.1 <- bitr(up.gene_, fromType = "SYMBOL",toType = c("ENTREZID"),OrgDb = org.Hs.eg.db)
    up.B_fun.gene = up.B_fun.gene.1$ENTREZID
    down.B_fun.gene.1 <- bitr(down.gene_, fromType = "SYMBOL",toType = c("ENTREZID"),OrgDb = org.Hs.eg.db)
    down.B_fun.gene = down.B_fun.gene.1$ENTREZID
    
    up.B_fun_bp_cc = enrichGO(gene = up.B_fun.gene,OrgDb = org.Hs.eg.db,ont = "BP",pAdjustMethod = "BH",pvalueCutoff  = 0.05,qvalueCutoff  = 0.05,readable= TRUE)
    write.table(up.B_fun_bp_cc,file = paste('4-',group_,'_bp_severe.txt',sep = ''),sep='\t',quote = FALSE)
    dpi = 300
    png(file = paste('4-',group_,'_bp_severe.png',sep = ''),width = dpi * 10,height = dpi * 6,units = "px",res = dpi,type = 'cairo')
    print(barplot(up.B_fun_bp_cc,showCategory=20,drop=T)+theme(axis.text.y = element_text(size = 16),legend.text = element_text(size = 16),legend.title = element_text(size = 16)))
    dev.off()
    
    down.B_fun_bp_cc = enrichGO(gene = down.B_fun.gene,OrgDb = org.Hs.eg.db,ont = "BP",pAdjustMethod = "BH",pvalueCutoff  = 0.05,qvalueCutoff  = 0.05,readable= TRUE)
    write.table(down.B_fun_bp_cc,file = paste('4-',group_,'_bp_mild.txt',sep = ''),sep='\t',quote = FALSE)
    dpi = 300
    png(file = paste('4-',group_,'_bp_mild.png',sep = ''),width = dpi * 10,height = dpi * 6,units = "px",res = dpi,type = 'cairo')
    print(barplot(down.B_fun_bp_cc,showCategory=20,drop=T)+theme(axis.text.y = element_text(size = 16),legend.text = element_text(size = 16),legend.title = element_text(size = 16)))
    dev.off()
    
  }
  
  ###now union result
  for(group_ in nkt_groups){
    hc.markers_i = read.delim2(paste('4-',group_,'_bp.txt',sep = ''),header = TRUE, stringsAsFactors = FALSE,check.names = FALSE, sep = "\t")
    hc.markers_i$p.adjust = as.numeric(hc.markers_i$p.adjust)
    hc.markers_i = hc.markers_i %>% filter(.,p.adjust < 0.05)
    hc.markers_i = hc.markers_i[1:10,]
    hc.markers_i = hc.markers_i %>% dplyr::select(ID,Description,p.adjust)
    hc.markers_i$p.adjust = -log10(hc.markers_i$p.adjust)
    hc.markers_i$Description = factor(hc.markers_i$Description,labels =hc.markers_i$Description,levels = hc.markers_i$Description)
    png(file = paste('4-',group_,'_bp_final.png',sep = ''),width = dpi * 8,height = dpi * 6,units = "px",res = dpi,type = 'cairo')
    pp <- ggplot(data=hc.markers_i, aes(x=Description, y=p.adjust)) + geom_bar(stat="identity",width = 0.7,fill = "#FDA1A1") + coord_flip()
    pp = pp + theme_minimal() + labs(x='',y='-log10(p.adjust)') + scale_x_discrete(limits = rev(levels(hc.markers_i$Description)),position = 'right') + scale_y_continuous(expand = c(0,0))
    pp = pp + geom_text(data=hc.markers_i, aes(x=Description, y=max(hc.markers_i$p.adjust)/100),label=hc.markers_i$Description, size=5,hjust = 0)
    pp = pp + theme(axis.ticks.y = element_blank(),axis.text.y = element_blank(),
                    panel.grid = element_blank(),
                    axis.line.x = element_line(size = 0.5, linetype = "solid", colour = "darkgrey"),
                    axis.line.y = element_line(size = 0.5, linetype = "solid", colour = "black"),
                    axis.title.x = element_text(size = 16),axis.text.x = element_text(size = 14),
                    axis.ticks.x = element_line(linetype = "solid",colour = "darkgrey",size = 0.5),axis.ticks.length.x = unit(2,'mm'))
    print(pp)
    dev.off()
  }
  
  library(ggrepel)
  library(ggplot2)
  library(dplyr)
  library(reshape2)
  setwd('/Users/yang/Downloads/work/nCoV/scrna_rev/balf/Figure')
  nkt_groups = c('CD8 T')
  for(group_ in nkt_groups){
    ##Severe
    hc.markers_i = read.delim2(paste('4-',group_,'_bp_severe_sy.txt',sep = ''),header = TRUE, stringsAsFactors = FALSE,check.names = FALSE, sep = "\t")
    hc.markers_i$p.adjust = as.numeric(hc.markers_i$p.adjust)
    hc.markers_i = hc.markers_i %>% filter(.,p.adjust < 0.05)
    hc.markers_i = hc.markers_i[1:10,]
    hc.markers_i = hc.markers_i %>% dplyr::select(ID,Description,p.adjust)
    hc.markers_i$p.adjust = -log10(hc.markers_i$p.adjust)
    hc.markers_i$Description = factor(hc.markers_i$Description,labels =hc.markers_i$Description,levels = hc.markers_i$Description)
    pdf(file = paste('4-',group_,'_bp_severe.pdf',sep = ''),width = 2.4,height = 1.8)
    pp <- ggplot(data=hc.markers_i, aes(x=Description, y=p.adjust)) + geom_bar(stat="identity",width = 0.7,fill = "#FDA1A1") + coord_flip()
    pp = pp + theme_minimal() + labs(x='',y='-log10(p.adjust)') + scale_x_discrete(limits = rev(levels(hc.markers_i$Description)),position = 'right') + scale_y_continuous(expand = c(0,0))
    pp = pp + geom_text(data=hc.markers_i, aes(x=Description, y=max(hc.markers_i$p.adjust)/100),label=hc.markers_i$Description, size=1.6,hjust = 0)
    pp = pp + theme(axis.ticks.y = element_blank(),axis.text.y = element_blank(),
                    panel.grid = element_blank(),
                    axis.line.x = element_line(size = 0.5, linetype = "solid", colour = "darkgrey"),
                    axis.line.y = element_line(size = 0.5, linetype = "solid", colour = "black"),
                    axis.title.x = element_text(size = 5),axis.text.x = element_text(size = 5),
                    axis.ticks.x = element_line(linetype = "solid",colour = "darkgrey",size = 0.5),axis.ticks.length.x = unit(1,'mm'))
    print(pp)
    dev.off()
    
    ##Mild
    hc.markers_i = read.delim2(paste('4-',group_,'_bp_mild_sy.txt',sep = ''),header = TRUE, stringsAsFactors = FALSE,check.names = FALSE, sep = "\t")
    hc.markers_i$p.adjust = as.numeric(hc.markers_i$p.adjust)
    hc.markers_i = hc.markers_i %>% filter(.,p.adjust < 0.05)
    hc.markers_i = hc.markers_i[1:10,]
    hc.markers_i = hc.markers_i %>% dplyr::select(ID,Description,p.adjust)
    hc.markers_i$p.adjust = -log10(hc.markers_i$p.adjust)
    hc.markers_i$Description = factor(hc.markers_i$Description,labels =hc.markers_i$Description,levels = hc.markers_i$Description)
    pdf(file = paste('4-',group_,'_bp_mild.pdf',sep = ''),width = 2.4,height = 1.8)
    pp <- ggplot(data=hc.markers_i, aes(x=Description, y=p.adjust)) + geom_bar(stat="identity",width = 0.7,fill = "#B7B8FF") + coord_flip()
    pp = pp + theme_minimal() + labs(x='',y='-log10(p.adjust)') + scale_x_discrete(limits = rev(levels(hc.markers_i$Description)),position = 'right') + scale_y_continuous(expand = c(0,0))
    pp = pp + geom_text(data=hc.markers_i, aes(x=Description, y=max(hc.markers_i$p.adjust)/100),label=hc.markers_i$Description, size=1.6,hjust = 0)
    pp = pp + theme(axis.ticks.y = element_blank(),axis.text.y = element_blank(),
                    panel.grid = element_blank(),
                    axis.line.x = element_line(size = 0.5, linetype = "solid", colour = "darkgrey"),
                    axis.line.y = element_line(size = 0.5, linetype = "solid", colour = "black"),
                    axis.title.x = element_text(size = 5),axis.text.x = element_text(size = 5),
                    axis.ticks.x = element_line(linetype = "solid",colour = "darkgrey",size = 0.5),axis.ticks.length.x = unit(1,'mm'))
    print(pp)
    dev.off()
    
  }
  
}
