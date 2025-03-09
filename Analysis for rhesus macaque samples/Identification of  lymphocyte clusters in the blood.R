library(Seurat)
library(SeuratDisk)
library(harmony)
library(dplyr)
library(ggplot2)

# Combine the blood samples
files <- list.files(path = "./T_cells_acute", 
                    recursive = F, full.names = F)
files <- files[grep('blood_CD4', files)]
for (x in files){
  name <- gsub('_filtered_features_bc_matrix.h5','',x)
  x <- Read10X_h5(paste0("./",x))
  assign(name,CreateSeuratObject(counts = x))
}
# Merge data
blood <- merge(RM_92T_blood_CD4_H, 
                             y = c(RM_93T_blood_CD4_SIV,
                                   RM_94T_blood_CD4_SIV,
                                   RM_95T_blood_CD4_SIV),
                             add.cell.ids = ls()[1:4])
#Create a new column
blood$sample <- rownames(merged_seurat_blood@meta.data)
# Add meta data
blood@meta.data <- separate(merged_seurat_blood@meta.data, col = 'sample', into = c('species','name','tissue','sort','condition','barcode', sep = '_'))
# Add mitochondria percentage
mt <- read.csv('./mt.genes.txt') # mt.genes.txt contains rhesus macaque mitochondria genes
mt <- gsub("_","-",mt$mt)
blood$mtUMI <- Matrix::colSums(blood[which(rownames(blood) %in% mt),], na.rm = T)
blood$mitoPercent <- blood$mtUMI*100/blood$nCount_RNA

# Filter
blood <- subset(blood,subset = nCount_RNA > 400 & nCount_RNA < 20000
                    & nFeature_RNA >400 & nFeature_RNA < 10000 & mitoPercent< 15)
# General workflow
blood <- NormalizeData(object = blood)
blood <- FindVariableFeatures(object = blood)
blood <- ScaleData(object = blood)
gene.list <- read.csv('./gene.list.csv') # gene.list contains rhesus macaque curated genes
gene.list <- gene.list$Gene.Symbol
blood <- RunPCA(object = blood, features = gene.list)

blood <- blood %>%
  RunHarmony(group.by.vars = 'name', plot_convergence = F, kmeans_init_nstart=20, kmeans_init_iter_max=100)

# Determine PCA
ElbowPlot(blood, ndims = 50)

# Clustering
blood <- blood %>%
  RunUMAP(reduction = 'harmony', dims = 1:20) %>%
  FindNeighbors(reduction = 'harmony', dims = 1:20) %>%
  FindClusters(res = 0.2)

blood$condition <- gsub("H","Uninfected",blood$condition)
blood$condition <- gsub("SIV","Acute-Infected",blood$condition)

# Plot the embeddings
plot1 <- DimPlot(blood, reduction = 'umap', raster = F, label = T, label.size = 8)+
  theme(legend.text = element_text(size = 20))
plot2 <- DimPlot(blood, reduction = 'umap', group.by = 'condition', raster = F)+
  theme(legend.text = element_text(size = 20))
plot3 <- DimPlot(blood, reduction = 'umap', group.by = 'name', raster = F)+
  theme(legend.text = element_text(size = 20))



# Characterization
plot1 <- VlnPlot(blood, features = c("CD3D","CD3G","CD3E","IL7R","TCF7","NKG7"),raster = F, ncol = 3,
                 pt.size = 0)&theme(axis.text = element_text(size = 20), 
                                    axis.title.y = element_text(size = 20),
                                    axis.title.x = element_blank(),
                                    plot.title = element_text(size = 25))

plot2 <- FeaturePlot(blood, features = c("CD3D","CD3G","CD3E","IL7R","TCF7","NKG7"),raster = F, ncol = 3)&
  theme(axis.title = element_blank(), plot.title = element_text(size = 25), axis.text = element_text(size = 20))



plot3 <- VlnPlot(blood, features = c("CSF1R","C1QB","CD14","FCGR3","MAMU-DRA","MAMU-DRB1"),raster = F, ncol = 3,
                 pt.size = 0)&theme(axis.text = element_text(size = 20), 
                                    axis.title.y = element_text(size = 20),
                                    axis.title.x = element_blank(),
                                    plot.title = element_text(size = 25))

plot4 <- FeaturePlot(blood, features = c("CSF1R","C1QB","CD14","FCGR3","MAMU-DRA","MAMU-DRB1"),raster = F, ncol = 3)&
  theme(axis.title = element_blank(), plot.title = element_text(size = 25), axis.text = element_text(size = 20))

# Rename the cluster based on annotations
new.cluster.id <- c("Lymphocyte","Lymphocyte","Lymphocyte","Lymphocyte","Lymphocyte",
                    "Monocyte","Lymphocyte","Monocyte","Lymphocyte","Monocyte","Lymphocyte")
names(new.cluster.id) <- levels(blood_all)
blood_all <- RenameIdents(blood_all, new.cluster.id)

SaveH5Seurat(blood, filename = "./blood_all")

# Subset the blood lymphocyte
blood_T <- subset(blood, idents = c(5,7,9), invert = T)

# Rerun the workflow
blood_T <- FindVariableFeatures(object = blood_T)
blood_T <- ScaleData(object = blood_T)
gene.list <- read.csv('./gene.list.csv')
gene.list <- gene.list$Gene.Symbol
blood_T <- RunPCA(object = blood_T, features = gene.list)

blood_T <- blood_T %>%
  RunHarmony(group.by.vars = 'name', plot_convergence = F)

# Determine PC to use
ElbowPlot(blood_T, ndims = 50)

#cluster
blood_T <- blood_T %>%
  RunUMAP(reduction = 'harmony', dims = 1:20) %>%
  FindNeighbors(reduction = 'harmony', dims = 1:20) %>%
  FindClusters(res = 0.2)

SaveH5Seurat(blood_T, filename = "./blood_T")

# Plot the embeddings
plot1 <- DimPlot(blood_T, reduction = 'umap', raster = F, label = T, label.size = 8)+
  theme(axis.title = element_text(size = 20), axis.text = element_text(size = 20), legend.position = "none")
plot2 <- DimPlot(blood_T, reduction = 'umap', group.by = 'condition', raster = F)+
  theme(legend.text = element_text(size = 20))
plot3 <- DimPlot(blood_T, reduction = 'umap', group.by = 'name', raster = F)+
  theme(legend.text = element_text(size = 20))

# Find the markers
Idents(blood_T) <-blood_T$seurat_clusters
Markers <- FindAllMarkers(blood_T, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)


# Heatmap for markers
blood_T <- ScaleData(blood_T, features = rownames(blood_T))
Markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 0.5) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10
p1 <- DoHeatmap(blood_T, features = top10$gene, size = 8)+theme(axis.text.y = element_text(size = 20), legend.position = "bottom")


# Plot CD4 and CD8B
# CD4/CD8B expression ratio
df <- FetchData(blood_T, vars = c("CD4","CD8A","CD8B","seurat_clusters"))
df$CD4CD8B <- log2((df$CD4+1)/(df$CD8B+1))
blood_T$CD4CD8B <- df$CD4CD8B

library(RColorBrewer)
plot1 <- FeaturePlot(blood_T, features = "CD4")+
  theme(axis.title = element_blank(), plot.title = element_text(size = 25), axis.text = element_text(size = 20))
plot2 <- FeaturePlot(blood_T, features = "CD8B")+
  theme(axis.title = element_blank(), plot.title = element_text(size = 25), axis.text = element_text(size = 20))
plot3 <- FeaturePlot(blood_T, features = "CD8A")+
  theme(axis.title = element_blank(), plot.title = element_text(size = 25), axis.text = element_text(size = 20))
plot4 <- FeaturePlot(blood_T, features = "CD4CD8B")+scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))+
  theme(axis.title = element_blank(), plot.title = element_text(size = 25), axis.text = element_text(size = 20))

# Plot T cell markers
plot1 <- FeaturePlot(blood_T, features = c("GZMB","PRF1","LEF1","SELL","ITGB1","FOXP3","STMN1","MKI67","mitoPercent"), ncol = 3)&
  theme(axis.title = element_blank(), plot.title = element_text(size = 25), axis.text = element_text(size = 20))

