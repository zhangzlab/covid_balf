T_integration<- function(){
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(ggplot2)
  library(reshape2)
  dpi = 300
  setwd('/home/data/results/workspace/nCoV_rev/balf/NKT_T')
  NKT_T = readRDS(file = "../NKT/6-NKT.rds")
  NKT_T = subset(NKT_T,idents = c('0','1','2','3','5','6','9','10'))
  write.table(table(NKT_T@meta.data$sample),file='NKT_sample.txt',quote = FALSE,sep='\t',row.names = FALSE)
  write.table(table(NKT_T@meta.data$group),file='NKT_group.txt',quote = FALSE,sep='\t',row.names = FALSE)
  
  ##############T###############
  DefaultAssay(NKT_T) <- "RNA"
  SubseNKTs.list <- SplitObject(NKT_T, split.by = "sample")
  for (i in 1:length(SubseNKTs.list)) {
    SubseNKTs.list[[i]] <- NormalizeData(SubseNKTs.list[[i]], verbose = FALSE)
    SubseNKTs.list[[i]] <- FindVariableFeatures(SubseNKTs.list[[i]], selection.method = "vst", nfeatures = 2000,verbose = FALSE)
  }
  samples_name = c('C141','C142','C144','C143','C145','C148','C149','C152')
  reference.list <- SubseNKTs.list[samples_name]
  NKT.temp <- FindIntegrationAnchors(object.list = reference.list, dims = 1:50,k.filter = 140)
  NKT.Integrated_clean <- IntegrateData(anchorset = NKT.temp, dims = 1:50)
  
  ###first generate data and scaledata in RNA assay
  DefaultAssay(NKT.Integrated_clean) <- "RNA"
  NKT.Integrated_clean[['percent.mito']] <- PercentageFeatureSet(NKT.Integrated_clean, pattern = "^MT-")
  NKT.Integrated_clean <- NormalizeData(object = NKT.Integrated_clean, normalization.method = "LogNormalize", scale.factor = 1e4)
  NKT.Integrated_clean <- FindVariableFeatures(object = NKT.Integrated_clean, selection.method = "vst", nfeatures = 2000,verbose = FALSE)
  NKT.Integrated_clean <- ScaleData(NKT.Integrated_clean, verbose = FALSE, vars.to.regress = c("nCount_RNA", "percent.mito"))
  
  ##change to integrated assay
  DefaultAssay(NKT.Integrated_clean) <- "integrated"
  dpi = 300
  png(file="6-qc.png", width = dpi*16, height = dpi*8, units = "px",res = dpi,type='cairo')
  VlnPlot(object = NKT.Integrated_clean, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
  dev.off()
  
  png(file="6-umi-gene.png", width = dpi*6, height = dpi*5, units = "px",res = dpi,type='cairo')
  FeatureScatter(object = NKT.Integrated_clean, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  dev.off()
  
  # Run the standard workflow for visualization and clustering
  NKT.Integrated_clean <- ScaleData(NKT.Integrated_clean, verbose = FALSE, vars.to.regress = c("nCount_RNA", "percent.mito"))
  NKT.Integrated_clean <- RunPCA(NKT.Integrated_clean, verbose = FALSE)
  #visulaization pca result
  NKT.Integrated_clean <- ProjectDim(object = NKT.Integrated_clean)
  png(file="6-pca.png", width = dpi*10, height = dpi*6, units = "px",res = dpi,type='cairo')
  ElbowPlot(object = NKT.Integrated_clean,ndims = 50)
  dev.off()
  
  ###cluster
  NKT.Integrated_clean <- FindNeighbors(object = NKT.Integrated_clean, dims = 1:50)
  NKT.Integrated_clean <- FindClusters(object = NKT.Integrated_clean, resolution = 0.8)
  
  ###tsne
  NKT.Integrated_clean <- RunTSNE(object = NKT.Integrated_clean, dims = 1:50)
  NKT.Integrated_clean <- RunUMAP(NKT.Integrated_clean, reduction = "pca", dims = 1:50)
  png(file="6-NKT-tsne.png", width = dpi*10, height = dpi*8, units = "px",res = dpi,type='cairo')
  DimPlot(object = NKT.Integrated_clean, reduction = 'tsne',label = TRUE)
  dev.off()
  png(file="6-NKT-umap.png", width = dpi*10, height = dpi*8, units = "px",res = dpi,type='cairo')
  DimPlot(object = NKT.Integrated_clean, reduction = 'umap',label = TRUE)
  dev.off()
  saveRDS(NKT.Integrated_clean, file = "6-NKT.rds")
  dpi = 300
  png(file="6-NKT-umap-split-sample.png", width = dpi*20, height = dpi*18, units = "px",res = dpi,type='cairo')
  DimPlot(object = NKT.Integrated_clean, reduction = 'umap',label = TRUE, split.by = 'sample', ncol = 5)
  dev.off()
  
  png(file="6-NKT-umap-split-group.png", width = dpi*20, height = dpi*8, units = "px",res = dpi,type='cairo')
  DimPlot(object = NKT.Integrated_clean, reduction = 'umap',label = TRUE, split.by = 'group',ncol = 3)
  dev.off()
  
  png(file="6-NKT-umap-group-sample.png", width = dpi*10, height = dpi*8, units = "px",res = dpi,type='cairo')
  DimPlot(object = NKT.Integrated_clean, reduction = 'umap',label = TRUE, group.by = 'sample')
  dev.off()
  
  png(file="6-NKT-umap-group-group.png", width = dpi*10, height = dpi*8, units = "px",res = dpi,type='cairo')
  DimPlot(object = NKT.Integrated_clean, reduction = 'umap',label = TRUE, group.by = 'group')
  dev.off()
  
  DefaultAssay(NKT.Integrated_clean) <- "RNA"
  # find markers for every cluster compared to all remaining cells, report only the positive ones
  NKT.Integrated_clean@misc$markers <- FindAllMarkers(object = NKT.Integrated_clean, assay = 'RNA',only.pos = TRUE, test.use = 'MAST')
  write.table(NKT.Integrated_clean@misc$markers,file='6-NKT-markers_MAST.txt',row.names = FALSE,quote = FALSE,sep = '\t')
  ##add average expression information
  library(reshape2)
  NKT.Integrated_clean@misc$averageExpression = AverageExpression(object = NKT.Integrated_clean)
  write.table(NKT.Integrated_clean@misc$averageExpression$RNA,file='6-NKT-average_MAST.txt',row.names = TRUE,quote = FALSE,sep = '\t')
  
  dpi = 300
  png(file="6-NKT-feature.png", width = dpi*24, height = dpi*5, units = "px",res = dpi,type='cairo')
  VlnPlot(object = NKT.Integrated_clean, features = c("nFeature_RNA", "nCount_RNA"))
  dev.off()
  
  saveRDS(NKT.Integrated_clean, file = "6-NKT.rds")
  
  hc.markers = read.delim2("6-NKT-markers_MAST.txt",header = TRUE, stringsAsFactors = FALSE,check.names = FALSE, sep = "\t")
  hc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC) -> top10
  tt = DoHeatmap(object = NKT.Integrated_clean, features = top10$gene) + NoLegend()
  ggplot2::ggsave(file="6-NKT-feature2.pdf",plot = tt,device = 'pdf',width = 80, height = 40, units = "in",dpi = dpi,limitsize = FALSE)
  tt1 = DoHeatmap(object = subset(NKT.Integrated_clean, downsample = 500), features = top10$gene) + NoLegend()
  ggplot2::ggsave(file="6-NKT-feature2-1.pdf",plot = tt1,device = 'pdf',width = 20, height = 10, units = "in",dpi = dpi,limitsize = FALSE)
  
}
