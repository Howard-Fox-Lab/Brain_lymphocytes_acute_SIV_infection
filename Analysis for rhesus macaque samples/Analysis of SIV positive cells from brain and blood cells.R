library(Seurat)
library(SeuratDisk)
library(dplyr)
library(ggplot2)


#Read in blood and brain T data
blood_T <- LoadH5Seurat("./blood_T.h5seurat")
brain_T <- LoadH5Seurat("./brain_T.h5seurat")
PLS_T <- LoadH5Seurat("./PLS_merged.h5seurat")

#SIV+ cells in the blood
SIV <- FetchData(blood_T, vars = c("SIV","ident","condition"))
SIV.pos <- filter(SIV, SIV>0)
cell <- rownames(SIV.pos)
SIV.blood <- subset(blood_T, cells =cell)
SIV.blood$Tissue <- rep("Blood",707)
SIV.blood$sample <- paste(SIV.blood$name,SIV.blood$sort, sep = "_")

#SIV+ cells in the brain
SIV <- FetchData(brain_T, vars = c("SIV","ident","condition"))
SIV.pos <- filter(SIV, SIV>0)
cell <- rownames(SIV.pos)
SIV.brain <- subset(brain_T, cells =cell)
SIV.brain$Tissue <- rep("Brain",148)

#Merge the siv in blood and brain
SIV <- merge(SIV.blood, SIV.brain)

#General clustering workflow
SIV <- NormalizeData(object = SIV)
SIV <- FindVariableFeatures(object = SIV)
SIV <- ScaleData(object = SIV)

gene.list <- read.csv('./gene.list.csv')
gene.list <- gene.list$Gene.Symbol
SIV <- RunPCA(object = SIV, features = gene.list)

SIV <- SIV %>%
  RunHarmony(group.by.vars = 'sample', plot_convergence = F)

# Determine PCA
ElbowPlot(SIV, ndims = 50)

# Perform clustering
SIV <- SIV %>%
  RunUMAP(reduction = 'harmony', dims = 1:10) %>%
  FindNeighbors(reduction = 'harmony', dims = 1:10) %>%
  FindClusters(res = 0.1)

# Characterize each cluster
plot2 <- FeaturePlot(SIV, features = c("CD4","CD8B","GZMB","SELL","MKI67","P2RY12"),raster = F, ncol = 3)&
  theme(axis.title = element_blank(), plot.title = element_text(size = 25), axis.text = element_text(size = 20))

# Rename clusters
new.cluster.id <- c("CD4+ CTL","CD4+ Naive","CD4+ Proliferating","Myeloid T")
names(new.cluster.id) <- levels(SIV)
SIV <- RenameIdents(SIV, new.cluster.id)
SIV$Label <- Idents(SIV)


plot1 <- DimPlot(SIV, reduction = 'umap', raster = F, label = F, split.by = "Tissue")+
  theme(axis.title = element_text(size = 20), axis.text = element_text(size = 20), legend.position = "none")

SaveH5Seurat(SIV, filename = "/work/foxlab/xiaoke/seurat/T_cells_acute/SIV")
