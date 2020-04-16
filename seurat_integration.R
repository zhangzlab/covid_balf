##virus summary data
library(Seurat)
library(Matrix)
library(dplyr)
setwd('/home/data/results/workspace/nCoV_rev/balf')
samples = read.delim2("../balf_1.txt",header = TRUE, stringsAsFactors = FALSE,check.names = FALSE, sep = "\t")
nCoV.list = list()
result.matrix = matrix(0,dim(samples)[1],4)
i = 0
for(sample_ in samples$sample){
  print(sample_)
  i = i + 1
  print(i)
  result.matrix[i,1] = sample_
  samplesi = samples %>% filter(.,sample == sample_)
  if(samplesi$group == 'HC'){
    result.matrix[i,2] = 0
    result.matrix[i,3] = 0
    result.matrix[i,4] = 0
  }else{
    datadir = paste("/home/data/results/workspace/COVID_matrix/",sample_,"/outs/filtered_feature_bc_matrix/",sep="")
    tmp = Read10X(data.dir = datadir)
    tmp.gene = tmp['nCoV',]
    aa <- summary(tmp.gene)
    bb = sum(tmp.gene)
    result.matrix[i,2] = aa[[4]]
    result.matrix[i,3] = aa[[6]]
    result.matrix[i,4] = bb
  }
}
result.dataframe = as.data.frame(result.matrix)
colnames(result.dataframe) = c('sample','nCoV_mean','nCoV_max','nCoV_sum')
write.table(result.dataframe,file='statistics.txt',row.names = FALSE,quote = FALSE,sep='\t')

##data quality control
library(Seurat)
library(Matrix)
library(dplyr)
setwd('/home/data/results/workspace/nCoV_rev/balf')
samples = read.delim2("../balf_1.txt",header = TRUE, stringsAsFactors = FALSE,check.names = FALSE, sep = "\t")
result.matrix = matrix(0,dim(samples)[1],3)
i = 0
for(sample_s in samples$sample){
  i = i + 1
  print(sample_s)
  sample_i = samples %>% dplyr::filter(.,sample == sample_s)
  datadir = paste("/home/data/results/workspace/COVID_matrix/",sample_s,"/outs/filtered_feature_bc_matrix/",sep="")
  nCoV.data.i <- Read10X(data.dir = datadir)
  nCoV.seurat <- CreateSeuratObject(counts = nCoV.data.i, min.cells = 3, min.features = 200, project = sample_s)
  nCoV.seurat[['percent.mito']] <- PercentageFeatureSet(nCoV.seurat, pattern = "^MT-")
  print(dim(nCoV.seurat))
  
  dpi = 300
  png(file=paste('filter/',sample_s,"_qc.png",sep=''), width = dpi*16, height = dpi*8, units = "px",res = dpi,type='cairo')
  print(VlnPlot(object = nCoV.seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3))
  dev.off()
  
  png(file=paste('filter/',sample_s,"_umi-mito.png",sep=''), width = dpi*6, height = dpi*5, units = "px",res = dpi,type='cairo')
  print(FeatureScatter(object = nCoV.seurat, feature1 = "nCount_RNA", feature2 = "percent.mito"))
  dev.off()
  
  png(file=paste('filter/',sample_s,"_umi-gene.png",sep=''), width = dpi*6, height = dpi*5, units = "px",res = dpi,type='cairo')
  print(FeatureScatter(object = nCoV.seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA"))
  dev.off()
  
  result.matrix[i,1] = sample_s
  result.matrix[i,2] = dim(nCoV.seurat)[2]
  sample_i$nFeature_RNA_low = as.numeric(sample_i$nFeature_RNA_low)
  sample_i$nFeature_RNA_high = as.numeric(sample_i$nFeature_RNA_high)
  sample_i$nCount_RNA = as.numeric(sample_i$nCount_RNA)
  sample_i$percent.mito = as.numeric(sample_i$percent.mito)
  nCoV.seurat.filter <- subset(x = nCoV.seurat, subset = nFeature_RNA > sample_i$nFeature_RNA_low & nFeature_RNA < sample_i$nFeature_RNA_high 
                               & nCount_RNA > sample_i$nCount_RNA & percent.mito < sample_i$percent.mito)
  result.matrix[i,3] = dim(nCoV.seurat.filter)[2]
}
result.dataframe = as.data.frame(result.matrix)
colnames(result.dataframe) = c('sample','before','filter')
write.table(result.dataframe,file='filter/statistics_filter.txt',row.names = FALSE,quote = FALSE,sep='\t')

#umi statistics
library(Seurat)
library(Matrix)
library(dplyr)
setwd('/home/data/results/workspace/nCoV_rev/balf')
samples = read.delim2("../balf_1.txt",header = TRUE, stringsAsFactors = FALSE,check.names = FALSE, sep = "\t")
result.matrix = matrix(0,dim(samples)[1],3)
i = 0
for(sample_s in samples$sample){
  i = i + 1
  print(sample_s)
  sample_i = samples %>% dplyr::filter(.,sample == sample_s)
  datadir = paste("/home/data/results/workspace/COVID_matrix/",sample_s,"/outs/filtered_feature_bc_matrix/",sep="")
  nCoV.data.i <- Read10X(data.dir = datadir)
  nCoV.seurat <- CreateSeuratObject(counts = nCoV.data.i, min.cells = 3, min.features = 200, project = sample_s)
  nCoV.seurat[['percent.mito']] <- PercentageFeatureSet(nCoV.seurat, pattern = "^MT-")
  print(dim(nCoV.seurat))
  
  sample_i$nFeature_RNA_low = as.numeric(sample_i$nFeature_RNA_low)
  sample_i$nFeature_RNA_high = as.numeric(sample_i$nFeature_RNA_high)
  sample_i$nCount_RNA = as.numeric(sample_i$nCount_RNA)
  sample_i$percent.mito = as.numeric(sample_i$percent.mito)
  nCoV.seurat.filter <- subset(x = nCoV.seurat, subset = nFeature_RNA > sample_i$nFeature_RNA_low & nFeature_RNA < sample_i$nFeature_RNA_high 
                               & nCount_RNA > sample_i$nCount_RNA & percent.mito < sample_i$percent.mito)
  result.matrix[i,1] = sample_s
  result.matrix[i,2] = median(nCoV.seurat.filter@meta.data$nFeature_RNA)
  result.matrix[i,3] = median(nCoV.seurat.filter@meta.data$nCount_RNA)
}
result.matrix.dataframe = as.data.frame(result.matrix)
colnames(result.matrix.dataframe) = c('sample','median.gene','median.umi')
write.table(result.matrix.dataframe,file='statistics_tables1.txt',sep='\t',row.names = FALSE,quote = FALSE)

##integrate data with seurat v3
library(Seurat)
library(Matrix)
library(dplyr)
library(ggplot2)
setwd('/home/data/results/workspace/nCoV_rev/balf')
samples = read.delim2("../balf_1.txt",header = TRUE, stringsAsFactors = FALSE,check.names = FALSE, sep = "\t")
nCoV.list = list()
for(sample_s in samples$sample){
  print(sample_s)
  sample_i = samples %>% dplyr::filter(.,sample == sample_s)
  datadir = paste("/home/data/results/workspace/COVID_matrix/",sample_s,"/outs/filtered_feature_bc_matrix/",sep="")
  sample.tmp = Read10X(data.dir = datadir)
  sample.tmp.seurat <- CreateSeuratObject(counts = sample.tmp, min.cells = 3, min.features = 200,project = sample_s)
  sample.tmp.seurat[['percent.mito']] <- PercentageFeatureSet(sample.tmp.seurat, pattern = "^MT-")
  sample_i$nFeature_RNA_low = as.numeric(sample_i$nFeature_RNA_low)
  sample_i$nFeature_RNA_high = as.numeric(sample_i$nFeature_RNA_high)
  sample_i$nCount_RNA = as.numeric(sample_i$nCount_RNA)
  sample_i$percent.mito = as.numeric(sample_i$percent.mito)
  sample.tmp.seurat <- subset(x = sample.tmp.seurat, subset = nFeature_RNA > sample_i$nFeature_RNA_low & nFeature_RNA < sample_i$nFeature_RNA_high 
                                & nCount_RNA > sample_i$nCount_RNA & percent.mito < sample_i$percent.mito)
  sample.tmp.seurat <- NormalizeData(sample.tmp.seurat, verbose = FALSE)
  sample.tmp.seurat <- FindVariableFeatures(sample.tmp.seurat, selection.method = "vst", nfeatures = 2000,verbose = FALSE)
  nCoV.list[sample_s] = sample.tmp.seurat
}
nCoV <- FindIntegrationAnchors(object.list = nCoV.list, dims = 1:50)
nCoV.integrated <- IntegrateData(anchorset = nCoV, dims = 1:50,features.to.integrate = rownames(nCoV))

####add  sample info
sample_info = as.data.frame(colnames(nCoV.integrated))
colnames(sample_info) = c('ID')
rownames(sample_info) = sample_info$ID
sample_info$sample = nCoV.integrated@meta.data$orig.ident
sample_info = dplyr::left_join(sample_info,samples)
rownames(sample_info) = sample_info$ID
nCoV.integrated = AddMetaData(object = nCoV.integrated, metadata = sample_info)

###first generate data and scale data in RNA assay
DefaultAssay(nCoV.integrated) <- "RNA"
nCoV.integrated[['percent.mito']] <- PercentageFeatureSet(nCoV.integrated, pattern = "^MT-")
nCoV.integrated <- NormalizeData(object = nCoV.integrated, normalization.method = "LogNormalize", scale.factor = 1e4)
nCoV.integrated <- FindVariableFeatures(object = nCoV.integrated, selection.method = "vst", nfeatures = 2000,verbose = FALSE)
nCoV.integrated <- ScaleData(nCoV.integrated, verbose = FALSE, vars.to.regress = c("nCount_RNA", "percent.mito"))

##change to integrated assay
DefaultAssay(nCoV.integrated) <- "integrated"
dpi = 300
png(file="qc.png", width = dpi*16, height = dpi*8, units = "px",res = dpi,type='cairo')
VlnPlot(object = nCoV.integrated, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
dev.off()

png(file="umi-gene.png", width = dpi*6, height = dpi*5, units = "px",res = dpi,type='cairo')
FeatureScatter(object = nCoV.integrated, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
dev.off()

# Run the standard workflow for visualization and clustering
nCoV.integrated <- ScaleData(nCoV.integrated, verbose = FALSE, vars.to.regress = c("nCount_RNA", "percent.mito"))
nCoV.integrated <- RunPCA(nCoV.integrated, verbose = FALSE,npcs = 100)
nCoV.integrated <- ProjectDim(object = nCoV.integrated)
png(file="pca.png", width = dpi*10, height = dpi*6, units = "px",res = dpi,type='cairo')
ElbowPlot(object = nCoV.integrated,ndims = 100)
dev.off()

###cluster
nCoV.integrated <- FindNeighbors(object = nCoV.integrated, dims = 1:50)
nCoV.integrated <- FindClusters(object = nCoV.integrated, resolution = 1.2) 

###tsne and umap
nCoV.integrated <- RunTSNE(object = nCoV.integrated, dims = 1:50)
nCoV.integrated <- RunUMAP(nCoV.integrated, reduction = "pca", dims = 1:50)
png(file="tsne.png", width = dpi*8, height = dpi*6, units = "px",res = dpi,type='cairo')
DimPlot(object = nCoV.integrated, reduction = 'tsne',label = TRUE)
dev.off()
png(file="umap.png", width = dpi*8, height = dpi*6, units = "px",res = dpi,type='cairo')
DimPlot(object = nCoV.integrated, reduction = 'umap',label = TRUE)
dev.off()

DefaultAssay(nCoV.integrated) <- "RNA"
# find markers for every cluster compared to all remaining cells, report only the positive ones
nCoV.integrated@misc$markers <- FindAllMarkers(object = nCoV.integrated, assay = 'RNA',only.pos = TRUE, test.use = 'MAST')
write.table(nCoV.integrated@misc$markers,file='marker_MAST.txt',row.names = FALSE,quote = FALSE,sep = '\t')

dpi = 300
png(file="feature.png", width = dpi*24, height = dpi*5, units = "px",res = dpi,type='cairo')
VlnPlot(object = nCoV.integrated, features = c("nFeature_RNA", "nCount_RNA"))
dev.off()
saveRDS(nCoV.integrated, file = "nCoV.rds")

hc.markers = read.delim2("marker_MAST.txt",header = TRUE, stringsAsFactors = FALSE,check.names = FALSE, sep = "\t")
hc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC) -> top10
tt1 = DoHeatmap(object = subset(nCoV.integrated, downsample = 500), features = top10$gene) + NoLegend()
ggplot2::ggsave(file="marker_heatmap_MAST.pdf",plot = tt1,device = 'pdf',width = 20, height = 16, units = "in",dpi = dpi,limitsize = FALSE)

#draw heatmap
library(Seurat)
library(Matrix)
library(dplyr)
library(ggplot2)
library(reshape2)
setwd('/home/data/results/workspace/nCoV_rev/balf')
markers = c('AGER','SFTPC','SCGB3A2','TPPP3','KRT5',
            'CD68','FCN1','CD1C','TPSB2','CD14','MARCO','CXCR2',
            'CLEC9A','IL3RA',
            'CD3D','CD8A','KLRF1',
            'CD79A','IGHG4','MS4A1',
            'VWF','DCN',
            'FCGR3A','TREM2','KRT18')
hc.markers = read.delim2("marker_MAST.txt",header = TRUE, stringsAsFactors = FALSE,check.names = FALSE, sep = "\t")
hc.markers %>% group_by(cluster) %>% top_n(n = 30, wt = avg_logFC) -> top30
var.genes = c(nCoV.integrated@assays$RNA@var.features,top30$gene,markers)
nCoV.integrated <- ScaleData(nCoV.integrated, verbose = FALSE, vars.to.regress = c("nCount_RNA", "percent.mito"),features = var.genes)
saveRDS(nCoV.integrated, file = "nCoV.rds")

dpi = 300
png(file="tsne.png", width = dpi*8, height = dpi*6, units = "px",res = dpi,type='cairo')
DimPlot(object = nCoV.integrated, reduction = 'tsne',label = TRUE)
dev.off()
png(file="umap.png", width = dpi*8, height = dpi*6, units = "px",res = dpi,type='cairo')
DimPlot(object = nCoV.integrated, reduction = 'umap',label = TRUE)
dev.off()

dpi = 300
png(file="nCoV-umap-group-sample.png", width = dpi*8, height = dpi*6, units = "px",res = dpi,type='cairo')
DimPlot(object = nCoV.integrated, reduction = 'umap',label = FALSE, group.by = 'sample_new')
dev.off()

png(file="nCoV-umap-split-sample.png", width = dpi*16, height = dpi*16, units = "px",res = dpi,type='cairo')
DimPlot(object = nCoV.integrated, reduction = 'umap',label = TRUE, split.by = 'sample_new', ncol = 4)
dev.off()

png(file="nCoV-umap-group-group.png", width = dpi*8, height = dpi*6, units = "px",res = dpi,type='cairo')
DimPlot(object = nCoV.integrated, reduction = 'umap',label = FALSE, group.by = 'group')
dev.off()

png(file="nCoV-umap-split-group.png", width = dpi*12, height = dpi*4, units = "px",res = dpi,type='cairo')
DimPlot(object = nCoV.integrated, reduction = 'umap',label = TRUE, split.by = 'group', ncol = 3)
dev.off()

png(file="nCoV-umap-group-disease.png", width = dpi*8, height = dpi*6, units = "px",res = dpi,type='cairo')
DimPlot(object = nCoV.integrated, reduction = 'umap',label = FALSE, group.by = 'disease')
dev.off()

png(file="nCoV-umap-split-disease.png", width = dpi*10, height = dpi*4.5, units = "px",res = dpi,type='cairo')
DimPlot(object = nCoV.integrated, reduction = 'umap',label = TRUE, split.by = 'disease', ncol = 2)
dev.off()

####marker expression
dpi = 300
markers = c('AGER','SFTPC','SCGB3A2','TPPP3','KRT5',
            'CD68','FCN1','CD1C','TPSB2','CD14','MARCO','CXCR2',
            'CLEC9A','IL3RA',
            'CD3D','CD8A','KLRF1',
            'CD79A','IGHG4','MS4A1',
            'VWF','DCN',
            'FCGR3A','TREM2','KRT18','HBB')
#markers = c('HBB')
png(file="marker/violin_marker.png", width = dpi*30, height = dpi*24, units = "px",res = dpi,type='cairo')
print(VlnPlot(object = nCoV.integrated, features = markers,pt.size = 0)+ theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5)))
dev.off()
png(file="marker/umap_marker.png", width = dpi*30, height = dpi*36, units = "px",res = dpi,type='cairo')
print(FeaturePlot(object = nCoV.integrated, features = markers,cols = c("lightgrey","#ff0000")))
dev.off()
for(marker in markers){
  png(file=paste("marker/violin_",marker,".png",sep=''), width = dpi*8, height = dpi*3, units = "px",res = dpi,type='cairo')
  print(VlnPlot(object = nCoV.integrated, features = marker,pt.size = 0)+ theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5)))
  dev.off()
  png(file=paste("marker/umap_",marker,".png",sep=''), width = dpi*6, height = dpi*4, units = "px",res = dpi,type='cairo')
  print(FeaturePlot(object = nCoV.integrated, features = marker,cols = c("lightgrey","#ff0000")))
  dev.off()
}

library(ggplot2)
pdf(file="marker_heatmap.pdf", width = 10, height = 8)
pp = DotPlot(nCoV.integrated, features = rev(markers),cols = c('white','#F8766D'),dot.scale =5) + RotatedAxis()
pp = pp + theme(axis.text.x = element_text(size = 12),axis.text.y = element_text(size = 12)) + labs(x='',y='') + 
  guides(color = guide_colorbar(title = 'Scale expression'),size = guide_legend(title = 'Percent expressed')) + 
  theme(axis.line = element_line(size = 0.6))
print(pp)
dev.off()



