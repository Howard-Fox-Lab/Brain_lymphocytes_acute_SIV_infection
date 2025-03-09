library(Seurat)
library(SeuratDisk)
library(dplyr)
library(ggplot2)


files <- list.files(path = "./PLOS-data/", 
                    recursive = F, full.names = F)
for (x in files){
  name <- gsub('.h5','',x)
  x <- Read10X_h5(paste0("./PLOS-data/",x))
  assign(name,CreateSeuratObject(counts = x, min.cell = 10))
}

merged_seurat <- merge(RM_35488_brain_CD45_SIV_PLS, 
                       y = c(RM_37164_brain_CD45_SIV_PLS,RM_40742_brain_CD45_SIV_PLS,
                             RM_41812_brain_CD45_SIV_PLS,RM_C1_brain_CD45_H_PLS,
                             RM_C2_brain_CD45_H_PLS),
                       add.cell.ids = ls()[2:7])
merged_seurat$sample <- rownames(merged_seurat@meta.data)
merged_seurat@meta.data <- separate(merged_seurat@meta.data, col = 'sample', into = c('species','name','tissue','cells','condition','source','barcode'), sep = '_')

#Add Mt percentage
mt <- read.csv('./mt.genes.txt') #mt.genes contains mitochondria genes in rhesus macaques
mt <- gsub("_","-",mt$mt)
merged_seurat$mtUMI <- Matrix::colSums(merged_seurat[which(rownames(merged_seurat) %in% mt),], na.rm = T)
merged_seurat$mitoPercent <- merged_seurat$mtUMI*100/merged_seurat$nCount_RNA

#Filtering
merged_seurat <- subset(merged_seurat,subset = nCount_RNA > 400 & nCount_RNA < 20000
              & nFeature_RNA >400 & nFeature_RNA < 10000 & mitoPercent< 15)


# Performed workflow
merged_seurat <- NormalizeData(object = merged_seurat)
merged_seurat <- FindVariableFeatures(object = merged_seurat)
merged_seurat <- ScaleData(object = merged_seurat)
gene.list <- read.csv('./gene.list.csv')
gene.list <- gene.list$Gene.Symbol
merged_seurat <- RunPCA(object = merged_seurat, features = gene.list)

merged_seurat <- merged_seurat %>%
  RunHarmony(group.by.vars = 'name', plot_convergence = F, kmeans_init_nstart=20, kmeans_init_iter_max=100)

# Perform clustering
merged_seurat <- merged_seurat %>%
  RunUMAP(reduction = 'harmony', dims = 1:20) %>%
  FindNeighbors(reduction = 'harmony', dims = 1:20) %>%
  FindClusters(res = 0.2)

plot1 <- DimPlot(merged_seurat, reduction = 'umap', raster = F, label = T, label.size = 8)+
  theme(axis.title = element_text(size = 20), axis.text = element_text(size = 20), legend.position = "none")
plot2 <- DimPlot(merged_seurat, reduction = 'umap', group.by = 'name', raster = F)+
  theme(legend.text = element_text(size = 20))

plot3 <- VlnPlot(merged_seurat, features = c("CD3D","CD3G","CD3E","IL7R","TCF7","NKG7"),raster = F, ncol = 3,
                 pt.size = 0)&theme(axis.text = element_text(size = 20), 
                                    axis.title.y = element_text(size = 20),
                                    axis.title.x = element_blank(),
                                    plot.title = element_text(size = 25))

plot4 <- VlnPlot(merged_seurat, features = c("GPR34","P2RY12","CX3CR1","C1QB","MAMU-DRA","MAMU-DRB1"), raster = F, ncol = 3,
                 pt.size = 0)&theme(axis.text = element_text(size = 20), 
                                    axis.title.y = element_text(size = 20),
                                    axis.title.x = element_blank(),
                                    plot.title = element_text(size = 25))

new.cluster.id <- c("Lymphocyte","Lymphocyte","CAM","Lymphocyte","Lymphocyte","CAM","Lymphocyte",
                    "CAM","Lymphocyte","Microglia","CAM","CAM","others")
names(new.cluster.id) <- levels(merged_seurat)
merged_seurat <- RenameIdents(merged_seurat, new.cluster.id)
merged_seurat$Label <- Idents(merged_seurat)

# Only kept lymphocyte populations
merged_seurat <- subset(merged_seurat, idents = "Lymphocyte")
SaveH5Seurat(merged_seurat, "./PLS_merged.h5seurat")
