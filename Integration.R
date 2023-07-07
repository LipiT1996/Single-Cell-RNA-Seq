library(Seurat)
library(tidyverse)
library(ggplot2)
library(gridExtra)


directories <- list.dirs(path = "Integrate scRNA data/", recursive = F, full.names = F)

for (x in directories){
  name <- gsub("_filtered_feature_bc_matrix", "",x)
  
  cts <- ReadMtx(mtx = paste0('Integrate scRNA data/',x, '/matrix.mtx.gz'),
          features = paste0('Integrate scRNA data/', x, '/features.tsv.gz'), 
          cells = paste0('Integrate scRNA data/',x, '/barcodes.tsv.gz'))
  
  assign(name, CreateSeuratObject(counts = cts))
  
}

merged_seurat<- merge(HB17_background, y = c(HB17_PDX, 
                             HB17_tumor,
                             HB30_PDX, 
                             HB30_tumor, 
                             HB53_background, 
                             HB53_tumor), 
      add.cell.ids = ls()[3:9], 
      project = 'HB')
merged_seurat
View(merged_seurat@meta.data)


merged_seurat$sample <- row.names(merged_seurat@meta.data)


merged_seurat@meta.data<- separate(merged_seurat@meta.data, 
                                   col = 'sample', into = c('Patient', 'Type', 'Barcode'), 
                                   sep = '_')

unique(merged_seurat@meta.data$Type)

merged_seurat[["percentage.mt"]] <- PercentageFeatureSet(object = merged_seurat, pattern = "^MT-")

VlnPlot(merged_seurat, features = c('nCount_RNA', 'nFeature_RNA', 'percentage.mt'), ncol = 3)
FeatureScatter(merged_seurat, feature1 = 'nFeature_RNA', feature2 = 'nCount_RNA')+
  geom_smooth(method = 'lm')

merged_seurat_filtered <- subset(merged_seurat, subset = nFeature_RNA > 500 & 
         nCount_RNA > 800 &
         percentage.mt < 10)
View(merged_seurat_filtered)

#Check whether the data requires integration. Explore the data.

#Normalize the data
merged_seurat_filtered <- NormalizeData(object = merged_seurat_filtered) 

#identify Varible feature
merged_seurat_filtered <- FindVariableFeatures(object = merged_seurat_filtered)

#scalling of the data
merged_seurat_filtered <- ScaleData(object = merged_seurat_filtered)

#Linear Dimentionality Reducation
merged_seurat_filtered <- RunPCA(object = merged_seurat_filtered)
ElbowPlot(merged_seurat_filtered)
merged_seurat_filtered<- FindNeighbors(object = merged_seurat_filtered, dims = 1:20)

#clusteing
merged_seurat_filtered<- FindClusters(object = merged_seurat_filtered)
merged_seurat_filtered <- RunUMAP(object = merged_seurat_filtered, dims = 1:20)

#non-linear dimentionality reduction
p1 <- DimPlot(merged_seurat_filtered, reduction = 'umap', group.by = 'Patient')
p2 <- DimPlot(merged_seurat_filtered, reduction = 'umap', group.by = 'Type', 
        cols = c('red','green', 'blue'))

grid.arrange(p1, p2, ncol = 2, nrow = 2)

#integration
obj.list <- SplitObject(merged_seurat_filtered, split.by = 'Patient')

for(i in 1:length(obj.list)){
  obj.list[[i]] <- NormalizeData(object = obj.list[[i]])
  obj.list[[i]] <- FindVariableFeatures(object = obj.list[[i]])
}

# select integration features
featurs <- SelectIntegrationFeatures(object.list = obj.list)

# find integration anchors
anchors <- FindIntegrationAnchors(object.list = obj.list, 
                      anchor.features = featurs)

#integrate data
seurat.integrated <- IntegrateData(anchorset = anchors)

#scale data, runPCA and UMAP with integrated data
seurat.integrated<- ScaleData(object = seurat.integrated)
seurat.integrated <- RunPCA(object = seurat.integrated)
seurat.integrated <- RunUMAP(object = seurat.integrated, dims = 1:50)

p3<-DimPlot(object = seurat.integrated, reduction = 'umap', group.by = 'Patient')
p4<-DimPlot(object = seurat.integrated, reduction = 'umap', group.by = 'Type')

grid.arrange(p3, p4, ncol = 2)

grid.arrange(p1,p2,p3,p4, ncol = 2, nrow = 2)
