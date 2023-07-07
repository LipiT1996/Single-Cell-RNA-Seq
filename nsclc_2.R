#Single-Cell RNA-Seq data analysis

# Load Libraries ####
library(Seurat)
library(tidyverse)

# 1. Data #####
#The data used here is 10x Genomics data for non-small cell lung cancer(NSCLC)

nsclc <- Read10X_h5(filename = "../lipivikramthakker/Downloads/20k_NSCLC_DTC_3p_nextgem_intron_Multiplex_count_raw_feature_bc_matrix.h5")
#message: Genome matrix has multiple modalities, returning a list of matrices for this genome
str(nsclc)
#3 modalities: Gene expression, antibody capture and multiplexing capture
data <- nsclc$`Gene Expression`
data[1:10, 1:10]

# 2. Covert to Seurat objects ####
nsclc.seurat.obj <- CreateSeuratObject(count = data, project = 'NSCLC', min.cells = 2, min.features = 200)
str(nsclc.seurat.obj)
nsclc.seurat.obj
# 33623 features across 71880 samples within 1 assay 

# 3. Quality Control ####
# %MT READS- Mitochondrial genes %. higher quality cells contains higher MT%.
View(nsclc.seurat.obj@meta.data)
nsclc.seurat.obj[["Percentage.MT"]] <- PercentageFeatureSet(nsclc.seurat.obj, pattern = "^MT-")
#Violin Plot
VlnPlot(nsclc.seurat.obj, features = c("nFeature_RNA","nCount_RNA", "Percentage.MT"), ncol = 3)
#Feature Scatter plot
FeatureScatter(nsclc.seurat.obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")+
  geom_smooth(method = "lm")

# Filtering the data ####
nsclc.seurat.obj <- subset(nsclc.seurat.obj, subset = nFeature_RNA > 200 & 
                             nFeature_RNA < 2500 & 
                             Percentage.MT < 5)
nsclc.seurat.obj
# 33623 features across 54053 samples within 1 assay
FeatureScatter(nsclc.seurat.obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")+
  geom_smooth(method = "lm")

# Normalization ####
nsclc.seurat.obj <- NormalizeData(nsclc.seurat.obj, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features ####
nsclc.seurat.obj <- FindVariableFeatures(nsclc.seurat.obj, selection.method = "vst", nfeatures = 2000)

# Identify highly variable Gene ####
top10 <- head(VariableFeatures(nsclc.seurat.obj), 10)

# plot variable features with and without labels ####
plot1 <- VariableFeaturePlot(nsclc.seurat.obj)
LabelPoints(plot = plot1, points = top10, repel = TRUE)

# Scaling ####
all.genes <- rownames(nsclc.seurat.obj)
nsclc.seurat.obj <- ScaleData(nsclc.seurat.obj, features = all.genes)

str(nsclc.seurat.obj)

# Perform Linear dimensionality reduction ####
nsclc.seurat.obj <- RunPCA(nsclc.seurat.obj, features = VariableFeatures(object = nsclc.seurat.obj))

# visualize PCA results ####
print(nsclc.seurat.obj[["pca"]], dims = 1:5, nfeatures = 5)
DimHeatmap(nsclc.seurat.obj, dims = 1, cells = 500, balanced = TRUE)

# determine dimensionality of the data
ElbowPlot(nsclc.seurat.obj)

# 7. Clustering ------------
nsclc.seurat.obj <- FindNeighbors(nsclc.seurat.obj, dims = 1:15)

# understanding resolution
nsclc.seurat.obj <- FindClusters(nsclc.seurat.obj, resolution = c(0.1,0.3, 0.5, 0.7, 1))
View(nsclc.seurat.obj@meta.data)
DimPlot(nsclc.seurat.obj, group.by = "RNA_snn_res.0.1", label = TRUE)
DimPlot(nsclc.seurat.obj, group.by = "RNA_snn_res.0.3", label = TRUE)
DimPlot(nsclc.seurat.obj, group.by = "RNA_snn_res.0.5", label = TRUE)

# setting identity of clusters
Idents(nsclc.seurat.obj)
Idents(nsclc.seurat.obj) <- "RNA_snn_res.0.1"
Idents(nsclc.seurat.obj)

non-linear dimensionality reduction --------------
  # If you haven't installed UMAP, you can do so via reticulate::py_install(packages =
  # 'umap-learn')
  nsclc.seurat.obj <- RunUMAP(nsclc.seurat.obj, dims = 1:15)
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(nsclc.seurat.obj, reduction = "umap")

