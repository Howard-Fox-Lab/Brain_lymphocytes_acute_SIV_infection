library(Seurat)
library(SeuratDisk)
library(dplyr)
library(tidyr)
library(harmony)
library(ggplot2)

# Load data
files <- list.files(path = "./Mac_T/h5file/", 
                    recursive = F, full.names = F)
for (x in files){
  name <- gsub('_filtered_feature_bc_matrix.h5','',x)
  x <- Read10X_h5(paste0("./Mac_T/h5file/",x))
  assign(name,CreateSeuratObject(counts = x, min.cell = 10))
}

#Merge data
merged_seurat <- merge(Campto_rep1, 
                       y = c(Campto_rep2, CD95_rep1, CD95_rep2, NoApop_rep1,
                             NoApop_rep2, NoCoCul_rep1, NoCoCul_rep2),
                       add.cell.ids = ls()[1:8])

#Create a new column
merged_seurat$sample <- rownames(merged_seurat@meta.data)
merged_seurat$name <- paste0(merged_seurat$treatment, merged_seurat$replicate)
#Split columns
merged_seurat@meta.data <- separate(merged_seurat@meta.data, col = 'sample', into = c('treatment','replicate','bc'), sep = '_')

#Add MT percentage
merged_seurat[['percent.mt']] <- PercentageFeatureSet(merged_seurat, pattern = '^MT-')

#Plot qc matrics
plot(VlnPlot(merged_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
             ncol = 3, log = T, pt.size = 0.1))

#Save the merged_seurat in h5 file
SaveH5Seurat(merged_seurat, file = "./merged_seurat")

# Filter
#Before filter: 46655
merged_seurat_filtered <- subset(merged_seurat,subset = nCount_RNA > 1000
                                 & nFeature_RNA >500 & percent.mt< 15)

#After filter: 42390

#Standard workflow
merged_seurat_filtered <- NormalizeData(merged_seurat_filtered)
merged_seurat_filtered <- FindVariableFeatures(object = merged_seurat_filtered)
merged_seurat_filtered <- ScaleData(object = merged_seurat_filtered)
merged_seurat_filtered <- RunPCA(object = merged_seurat_filtered)
merged_seurat_filtered <- merged_seurat_filtered %>%
  RunHarmony(group.by.vars = 'name', plot_convergence = F)

#Elbowplot
#Plot qc matrics
plot(ElbowPlot(merged_seurat_filtered, ndims = 50))

merged_seurat_filtered <- merged_seurat_filtered %>%
  RunUMAP(reduction = 'harmony', dims = 1:30) %>%
  FindNeighbors(reduction = 'harmony', dims = 1:30) %>%
  FindClusters(res = 0.2)

SaveH5Seurat(merged_seurat_filtered, file = "./merged_seurat_filtered")

# Convert seurat h5 to h5ad for visualization in python
merged_seurat_filtered <- DietSeurat(object = merged_seurat_filtered, counts = T, data = T, scale.data = F, assays = "RNA", dimreducs =c('pca','tsne','umap','harmony'))
SaveH5Seurat(merged_seurat_filtered,"./seurat_allgenes", overwrite = T)
Convert('seurat_allgenes.h5seurat', dest = 'h5ad')


# Reclustering on just NoApop and NoCoCul groups:
# Subset the NoApop and NoCoCul
New <- subset(merged_seurat_filtered, subset = treatment == c("NoApop","NoCoCul"))

# Standard work flow
New <- FindVariableFeatures(object = New)
New <- ScaleData(object = New)
New <- RunPCA(object = New)
New <- New %>%
  RunHarmony(group.by.vars = 'name', plot_convergence = F)
New <- New %>%
  RunUMAP(reduction = 'harmony', dims = 1:30) %>%
  FindNeighbors(reduction = 'harmony', dims = 1:30) %>%
  FindClusters(res = 0.2)

markers <- FindAllMarkers(New, only.pos = TRUE,
                          min.pct = 0.25, logfc.threshold = 0.25)
markers <- markers%>% filter(avg_log2FC > 0.5 & p_val_adj < 0.05)

SaveH5Seurat(New, filename = "./New_clusters")
New <- DietSeurat(object = New, counts = T, data = T, scale.data = F, assays = "RNA", dimreducs =c('pca','tsne','umap','harmony'))
SaveH5Seurat(New, filename = "./New_clusters_2")
# Convert seurat h5 to h5ad for visualization in python
Convert('New_clusters_2.h5seurat', dest = 'h5ad')

