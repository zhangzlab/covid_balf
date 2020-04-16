cytokine_fig <- function(){
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(ggplot2)
  library(reshape2)
  library(RColorBrewer)
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
  
  cytokines = c('IL1B','IL4','IL5','IL6','IL10','IL12B','IFNA1',
                'IFNG','TNF','CXCL1','CXCL2','CXCL3','CXCL4','CXCL5','CXCL6','CXCL7','CXCL8','CXCL9','CXCL10','CXCL11','CXCL12','CXCL13','CXCL14','CXCL15','CXCL16','CXCL17',
                'XCL1','XCL2','CX3CL1','CCL1','CCL2','CCL3','CCL3L1','CCL4','CCL4L1','CCL5','CCL6','CCL7','CCL8','CCL9','CCL10','CCL12','CCL13','CCL14','CCL15','CCL18','CCL19','CCL20','CCL21','CCL22','CCL23','CCL24','CCL25','CCL26','CCL27','CCL28')
  
  DefaultAssay(Macrophage.Integrated1) = 'RNA'
  dpi = 600
  png(file="5-umap_cytokine.png", width = 14*600, height = 10*600,units = "px",res = dpi,type='cairo')
  pp_temp = FeaturePlot(object = Macrophage.Integrated1, features = cytokines,cols = c("lightgrey","#ff0000"),combine = FALSE)
  plots <- lapply(X = pp_temp, FUN = function(p) p + theme(axis.title = element_text(size = 18),axis.text = element_text(size = 18),
                                                           plot.title = element_text(family = 'sans',face='italic',size=18),
                                                           legend.text = element_text(size = 18),legend.key.height = unit(1.8,"line"),legend.key.width = unit(1.2,"line"),
                                                           axis.line = element_line(size = 1.25),axis.ticks = element_line(size = 1)))
  pp = CombinePlots(plots = plots,ncol = 4,legend = 'right')
  print(pp)
  dev.off()
  
  samples_name_rev = rev(c('HC1','HC2','HC3','HC4','O1','O2','O3','S1','C1','C2','C3','C4','C5'))
  Macrophage.Integrated1@meta.data$sample_new = factor(Macrophage.Integrated1@meta.data$sample_new,levels =samples_name_rev,labels = samples_name_rev)
  Idents(Macrophage.Integrated1) = 'sample_new'
  pdf(file="5-heatmap-cytokine_sample.pdf", width = 11, height = 4)
  pp = DotPlot(Macrophage.Integrated1, features = rev(cytokines),cols = c('white','#F8766D'),dot.scale =5) + RotatedAxis()
  pp = pp + theme(axis.text.x = element_text(size = 9),axis.text.y = element_text(size = 9)) + labs(x='',y='') + 
    guides(color = guide_colorbar(title = 'Scale expression'),size = guide_legend(title = 'Percent expressed')) + 
    theme(axis.line = element_line(size = 1))
  print(pp)
  dev.off()
  Idents(Macrophage.Integrated1) = 'celltype'
  
  ###perform de analysis of the cytokines between three groups
  cytokines_new = c('IL1B','IL4','IL5','IL6','IL10','IL12B','IFNA1',
                    'IFNG','TNF','CXCL1','CXCL2','CXCL3','CXCL5','CXCL6','CXCL8','CXCL9','CXCL10','CXCL11','CXCL12','CXCL13','CXCL14','CXCL16','CXCL17',
                    'XCL1','XCL2','CX3CL1','CCL1','CCL2','CCL3','CCL3L1','CCL4','CCL5','CCL7','CCL8','CCL13','CCL14','CCL15','CCL18','CCL19','CCL20','CCL22','CCL23','CCL24','CCL25','CCL26','CCL27','CCL28')
  cytokines_new_select = c('IL1B','IL6','TNF','CXCL1','CXCL2','CXCL3','CXCL8','CCL2','CCL3',
                           'CCL3L1','CCL4','CCL7','CCL8','CXCL9','CXCL10','CXCL11','CXCL16')
  Macrophage.Integrated1_violin = Macrophage.Integrated1
  Idents(Macrophage.Integrated1_violin) = 'group'
  pdf(file="5-violin_cytokine.pdf", width = 6, height = 3)
  pp_temp = VlnPlot(object = Macrophage.Integrated1_violin, ncol = 1, features = cytokines_new_select,pt.size = 0,combine = FALSE)
  plots <- lapply(X = pp_temp, FUN = function(p) p + labs(x='') + theme(axis.text = element_text(size = 5),axis.text.x = element_text(size = 5,angle = 0,hjust = 0.5),axis.title = element_text(size = 5),plot.title = element_text(size = 5,family = 'sans',face = 'italic'),
                                                                        axis.line = element_line(size = 0.4),axis.ticks = element_line(size = 0.4)))
  pp = CombinePlots(plots = plots,ncol = 6,legend = 'none')
  print(pp)  
  dev.off()
  
  Idents(Macrophage.Integrated1) = 'group'
  Macrophage.Integrated1_1 = subset(Macrophage.Integrated1,idents = c('O','S/C'))
  deg1 = FindAllMarkers(Macrophage.Integrated1_1,assay='RNA',features=cytokines_new,test.use = 'MAST',logfc.threshold = 0,min.pct = 0,only.pos = TRUE)
  
  Macrophage.Integrated1_2 = subset(Macrophage.Integrated1,idents = c('HC','O'))
  deg2 = FindAllMarkers(Macrophage.Integrated1_2,features=cytokines_new,assay='RNA',test.use = 'MAST',logfc.threshold = 0,min.pct = 0,only.pos = TRUE)
  
  Macrophage.Integrated1_3 = subset(Macrophage.Integrated1,idents = c('HC','S/C'))
  deg3 = FindAllMarkers(Macrophage.Integrated1_3,features=cytokines_new,assay='RNA',test.use = 'MAST',logfc.threshold = 0,min.pct = 0,only.pos = TRUE)
  write.table(deg1,file = '5-cytokine-deg-severe-mild.txt',quote = FALSE,sep='\t')
  write.table(deg2,file = '5-cytokine-deg-mild-hc.txt',quote = FALSE,sep='\t')
  write.table(deg3,file = '5-cytokine-deg-severe-hc.txt',quote = FALSE,sep='\t')
  
  deg1_o = deg1 %>% select(gene,cluster,p_val_adj) 
  colnames(deg1_o) = c('Gene','Severe vs Moderate Cluster','Severe vs Moderate adjust P')
  deg2_o = deg2 %>% select(gene,cluster,p_val_adj) 
  colnames(deg2_o) = c('Gene','Moderate vs HC Cluster','Moderate vs HC adjust P')
  deg3_o = deg3 %>% select(gene,cluster,p_val_adj) 
  colnames(deg3_o) = c('Gene','Severe vs HC Cluster','Severe vs HC adjust P')
  allgene = as.data.frame(unique(c(deg1_o$Gene,deg2_o$Gene,deg3_o$Gene)))
  colnames(allgene) = 'Gene'
  allgene = dplyr::left_join(allgene,deg1_o)
  allgene = dplyr::left_join(allgene,deg2_o)
  allgene = dplyr::left_join(allgene,deg3_o)
  write.table(allgene,file='5-cytokine.txt',sep='\t',row.names = FALSE,quote = FALSE)
  
  ###Macrophage: group1 .. group4
  ###NKT: CCR7+ T CD8 T NK Treg innate T
  ###neutrophill
  ###B
  ###plasma
  ###pDC
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(ggplot2)
  library(reshape2)
  library(RColorBrewer)
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
  
  ###NKT
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
  NKT.Integrated1 = subset(NKT.Integrated,idents=c('CCR7+ T','CD8 T','NK','Treg','innate T'))
  
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(ggplot2)
  library(reshape2)
  setwd('/home/data/results/workspace/nCoV_rev/balf/Figure')
  nCoV.integrated = readRDS(file = "../nCoV.rds")
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
  nCoV.integrated1 = subset(nCoV.integrated1,idents=c('Neutrophil','B','Plasma','pDC'))
  
  macro.celltype = Macrophage.Integrated1@meta.data
  nkt.celltype = NKT.Integrated1@meta.data
  other.celltype = nCoV.integrated1@meta.data
  macro.celltype = macro.celltype %>% select(ID,celltype)
  nkt.celltype = nkt.celltype %>% select(ID,celltype)
  other.celltype = other.celltype %>% select(ID,celltype)
  macro.nkt = rbind(macro.celltype,nkt.celltype)
  macro.nkt = rbind(macro.nkt,other.celltype)
  
  ###load nCoV
  nCoV.integrated = readRDS(file = "../nCoV.rds")
  nCoV.integrated2 = subset(nCoV.integrated, cells = macro.nkt$ID)
  nCoV.integrated2@meta.data = dplyr::left_join(nCoV.integrated2@meta.data,macro.nkt)
  Idents(nCoV.integrated2) = 'celltype'
  cytokines = c('CCR1','CCR2','CCR2B','CCR3','CCR4','CCR5','CCR6','CCR7','CCR8','CCR9','CCR10','CXCR1','CXCR2','CXCR3','CXCR3B','CXCR4','CXCR5','CXCR6','ACKR3','XCR1','CX3CR1')
  nCoV_groups = rev(c('Group1','Group2','Group3','Group4','CCR7+ T','CD8 T','NK','Treg','innate T','Neutrophil','B','Plasma','pDC'))
  nCoV.integrated2$celltype = factor(nCoV.integrated2$celltype,levels = nCoV_groups,labels = nCoV_groups)
  Idents(nCoV.integrated2) = nCoV.integrated2$celltype
  pdf(file="5-heatmap-cytokine_macro_nkt.pdf", width = 9, height = 4)
  pp = DotPlot(nCoV.integrated2, features = rev(cytokines),cols = c('white','#F8766D'),dot.scale =6.5) + RotatedAxis()
  pp = pp + theme(axis.text.x = element_text(size = 9),axis.text.y = element_text(size = 9)) + labs(x='',y='') + 
    guides(color = guide_colorbar(title = 'Scale expression'),size = guide_legend(title = 'Percent expressed')) + 
    theme(axis.line = element_line(size = 1))
  print(pp)
  dev.off()
}
