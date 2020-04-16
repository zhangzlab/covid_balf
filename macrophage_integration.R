macrophage_integration <- function(){
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(ggplot2)
  library(reshape2)
  setwd('/home/data/results/workspace/nCoV_rev/balf/Myeloid')
  nCoV.integrated = readRDS(file = "../nCoV.rds")
  Macrophage = subset(nCoV.integrated,idents = c('0','1','2','3','4','5','7','8','10','11','12','18','21','22','23','26'))
  write.table(table(Macrophage@meta.data$sample_new),file='Macrophage_sample.txt',quote = FALSE,sep='\t',row.names = FALSE)
  write.table(table(Macrophage@meta.data$group),file='Macrophage_group.txt',quote = FALSE,sep='\t',row.names = FALSE)
  
  ##############Macrophage###############
  DefaultAssay(Macrophage) <- "RNA"
  SubseMacrophages.list <- SplitObject(Macrophage, split.by = "sample")
  for (i in 1:length(SubseMacrophages.list)) {
    SubseMacrophages.list[[i]] <- NormalizeData(SubseMacrophages.list[[i]], verbose = FALSE)
    SubseMacrophages.list[[i]] <- FindVariableFeatures(SubseMacrophages.list[[i]], selection.method = "vst", nfeatures = 2000,verbose = FALSE)
  }
  samples_name = c('C51','C52','C100','GSM3660650','C141','C142','C144','C143','C145','C146','C148','C149','C152')
  reference.list <- SubseMacrophages.list[samples_name]
  Macrophage.temp <- FindIntegrationAnchors(object.list = reference.list, dims = 1:50,k.filter = 115)
  Macrophage.Integrated <- IntegrateData(anchorset = Macrophage.temp, dims = 1:50)
  
  ###first generate data and scaledata in RNA assay
  DefaultAssay(Macrophage.Integrated) <- "RNA"
  Macrophage.Integrated[['percent.mito']] <- PercentageFeatureSet(Macrophage.Integrated, pattern = "^MT-")
  Macrophage.Integrated <- NormalizeData(object = Macrophage.Integrated, normalization.method = "LogNormalize", scale.factor = 1e4)
  Macrophage.Integrated <- FindVariableFeatures(object = Macrophage.Integrated, selection.method = "vst", nfeatures = 2000,verbose = FALSE)
  Macrophage.Integrated <- ScaleData(Macrophage.Integrated, verbose = FALSE, vars.to.regress = c("nCount_RNA", "percent.mito"))
  
  ##change to integrated assay
  DefaultAssay(Macrophage.Integrated) <- "integrated"
  dpi = 300
  png(file="2-qc.png", width = dpi*16, height = dpi*8, units = "px",res = dpi,type='cairo')
  VlnPlot(object = Macrophage.Integrated, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
  dev.off()
  
  png(file="2-umi-gene.png", width = dpi*6, height = dpi*5, units = "px",res = dpi,type='cairo')
  FeatureScatter(object = Macrophage.Integrated, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  dev.off()
  
  # Run the standard workflow for visualization and clustering
  Macrophage.Integrated <- ScaleData(Macrophage.Integrated, verbose = FALSE, vars.to.regress = c("nCount_RNA", "percent.mito"))
  Macrophage.Integrated <- RunPCA(Macrophage.Integrated, verbose = FALSE)
  #visulaization pca result
  Macrophage.Integrated <- ProjectDim(object = Macrophage.Integrated)
  png(file="2-pca.png", width = dpi*10, height = dpi*6, units = "px",res = dpi,type='cairo')
  ElbowPlot(object = Macrophage.Integrated,ndims = 50)
  dev.off()
  
  ###cluster
  Macrophage.Integrated <- FindNeighbors(object = Macrophage.Integrated, dims = 1:50)
  Macrophage.Integrated <- FindClusters(object = Macrophage.Integrated, resolution = 0.8)
  
  ###tsne and umap
  Macrophage.Integrated <- RunTSNE(object = Macrophage.Integrated, dims = 1:50)
  Macrophage.Integrated <- RunUMAP(Macrophage.Integrated, reduction = "pca", dims = 1:50)
  png(file="2-Macrophage-tsne.png", width = dpi*7, height = dpi*5, units = "px",res = dpi,type='cairo')
  DimPlot(object = Macrophage.Integrated, reduction = 'tsne',label = TRUE)
  dev.off()
  png(file="2-Macrophage-umap.png", width = dpi*7, height = dpi*5, units = "px",res = dpi,type='cairo')
  DimPlot(object = Macrophage.Integrated, reduction = 'umap',label = TRUE)
  dev.off()
  saveRDS(Macrophage.Integrated, file = "2-Macrophage.rds")
  dpi = 300
  png(file="2-Macrophage-umap-split-sample.png", width = dpi*16, height = dpi*14, units = "px",res = dpi,type='cairo')
  DimPlot(object = Macrophage.Integrated, reduction = 'umap',label = TRUE, split.by = 'sample', ncol = 4)
  dev.off()
  
  png(file="2-Macrophage-umap-split-group.png", width = dpi*16, height = dpi*6, units = "px",res = dpi,type='cairo')
  DimPlot(object = Macrophage.Integrated, reduction = 'umap',label = TRUE, split.by = 'group',ncol = 3)
  dev.off()
  
  png(file="2-Macrophage-umap-group-sample.png", width = dpi*7, height = dpi*5, units = "px",res = dpi,type='cairo')
  DimPlot(object = Macrophage.Integrated, reduction = 'umap',label = TRUE, group.by = 'sample')
  dev.off()
  
  png(file="2-Macrophage-umap-group-group.png", width = dpi*7, height = dpi*5, units = "px",res = dpi,type='cairo')
  DimPlot(object = Macrophage.Integrated, reduction = 'umap',label = TRUE, group.by = 'group')
  dev.off()
  
  DefaultAssay(Macrophage.Integrated) <- "RNA"
  # find markers for every cluster compared to all remaining cells, report only the positive ones
  Macrophage.Integrated@misc$markers <- FindAllMarkers(object = Macrophage.Integrated, assay = 'RNA',only.pos = TRUE, test.use = 'MAST')
  write.table(Macrophage.Integrated@misc$markers,file='2-Macrophage-markers_MAST.txt',row.names = FALSE,quote = FALSE,sep = '\t')
  ##add average expression information
  library(reshape2)
  Macrophage.Integrated@misc$averageExpression = AverageExpression(object = Macrophage.Integrated)
  write.table(Macrophage.Integrated@misc$averageExpression$RNA,file='2-Macrophage-average_MAST.txt',row.names = TRUE,quote = FALSE,sep = '\t')
  
  dpi = 300
  png(file="2-Macrophage-feature.png", width = dpi*24, height = dpi*5, units = "px",res = dpi,type='cairo')
  VlnPlot(object = Macrophage.Integrated, features = c("nFeature_RNA", "nCount_RNA"))
  dev.off()
  
  saveRDS(Macrophage.Integrated, file = "2-Macrophage.rds")
  
  hc.markers = read.delim2("2-Macrophage-markers_MAST.txt",header = TRUE, stringsAsFactors = FALSE,check.names = FALSE, sep = "\t")
  hc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC) -> top10
  tt1 = DoHeatmap(object = subset(Macrophage.Integrated, downsample = 500), features = top10$gene) + NoLegend()
  ggplot2::ggsave(file="2-Macrophage-feature2-1.pdf",plot = tt1,device = 'pdf',width = 20, height = 20, units = "in",dpi = dpi,limitsize = FALSE)
}
