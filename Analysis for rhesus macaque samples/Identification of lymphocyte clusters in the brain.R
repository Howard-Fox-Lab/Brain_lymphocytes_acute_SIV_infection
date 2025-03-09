library(Seurat)
library(SeuratDisk)
library(dplyr)
library(harmony)
library(ggplot2)

files <- list.files(path = "./T_cells_acute", 
                    recursive = F, full.names = F)
files <- c(files[grep('brain_CD45', files)], files[grep('brain_CD45CD11B', files)])

for (x in files){
  name <- gsub('_filtered_features_bc_matrix.h5','',x)
  x <- Read10X_h5(paste0("./T_cells_acute",x))
  assign(name,CreateSeuratObject(counts = x, min.cell = 10))
}

#Merge data
merged_seurat <- merge(RM_104T_brain_CD45CD11B_H_1, 
                       y = c(RM_111T_brain_CD45CD11B_H_1,RM_92T_brain_CD45CD11B_H_1,RM_93T_brain_CD45CD11B_SIV_1,
                             RM_94T_brain_CD45CD11B_SIV_1,RM_95T_brain_CD45CD11B_SIV_1,RM_92T_brain_CD45_H_1,
                             RM_92T_brain_CD45_H_2,RM_93T_brain_CD45_SIV_1,RM_93T_brain_CD45_SIV_2,
                             RM_94T_brain_CD45_SIV_1,RM_94T_brain_CD45_SIV_2,RM_95T_brain_CD45_SIV_1,
                             RM_95T_brain_CD45_SIV_2),
                       add.cell.ids = ls()[1:14])

#Create a new column
merged_seurat$sample <- rownames(merged_seurat@meta.data)
# Add meta data
merged_seurat@meta.data <- separate(merged_seurat@meta.data, col = 'sample', into = c('species','name','tissue','sort','condition',"number",'barcode', sep = '_'))
# Add mitochondria percentage
mt <- read.csv('./mt.genes.txt') # mt.genes.txt contains rhesus macaque mitochondria genes
mt <- gsub("_","-",mt$mt)
merged_seurat$mtUMI <- Matrix::colSums(merged_seurat[which(rownames(merged_seurat) %in% mt),], na.rm = T)
merged_seurat$mitoPercent <- merged_seurat$mtUMI*100/merged_seurat$nCount_RNA

merged_seurat$sample <- paset(paste(merged_seurat$name,Brain_all$sort,sep = "_"),Brain_all$number, sep = "_")

table(merged_seurat$sample)

# Filtering
Brain_all <- subset(merged_seurat,subset = nCount_RNA > 400 & nCount_RNA < 20000
                       & nFeature_RNA >400 & nFeature_RNA < 10000 & mitoPercent< 15)
table(Brain_all$sample)

# Standard workflow
Brain_all <- NormalizeData(object = Brain_all)
Brain_all <- FindVariableFeatures(object = Brain_all)
Brain_all <- ScaleData(object = Brain_all)
gene.list <- read.csv('/work/foxlab/xiaoke/seurat/gene.list.csv') #gene.list contains rhesus macaque curated genes
gene.list <- gene.list$Gene.Symbol
Brain_all <- RunPCA(object = Brain_all, features = gene.list)

Brain_all <- Brain_all %>%
  RunHarmony(group.by.vars = 'sample', plot_convergence = F, kmeans_init_nstart=20, kmeans_init_iter_max=100)

# Determine PCA
ElbowPlot(Brain_all, ndims = 50)

# Perform clustering
Brain_all <- Brain_all %>%
  RunUMAP(reduction = 'harmony', dims = 1:20) %>%
  FindNeighbors(reduction = 'harmony', dims = 1:20) %>%
  FindClusters(res = 0.2)

# Save the file
Brain_all$condition <- gsub("H","Uninfected",Brain_all$condition)
Brain_all$condition <- gsub("SIV","Acute-Infected",Brain_all$condition)

SaveH5Seurat(Brain_all, filename = "./brain_all_filtered")

# Plot the embeddings
plot1 <- DimPlot(Brain_all, reduction = 'umap', raster = F, label = T, label.size = 8)+
  theme(axis.title = element_text(size = 20), axis.text = element_text(size = 20), legend.position = "none")
plot2 <- DimPlot(Brain_all, reduction = 'umap', group.by = 'condition', raster = F, label = F)+
  theme(axis.title = element_text(size = 20), axis.text = element_text(size =20), legend.position = c(0.5,0))
plot3 <- DimPlot(Brain_all, reduction = 'umap', group.by = 'sample', raster = F)+
  theme(legend.text = element_text(size = 20))

# Plot Markers
# Lmphocyte
plot1 <- VlnPlot(Brain_all, features = c("CD3D","CD3G","CD3E","IL7R","TCF7","NKG7"),raster = F, ncol = 3,
                 pt.size = 0)&theme(axis.text = element_text(size = 20), 
                                    axis.title.y = element_text(size = 20),
                                    axis.title.x = element_blank(),
                                    plot.title = element_text(size = 25))
plot2 <- FeaturePlot(Brain_all, features = c("CD3D","CD3G","CD3E","IL7R","TCF7","NKG7"),raster = F, ncol = 3)&
  theme(axis.title = element_blank(), plot.title = element_text(size = 25), axis.text = element_text(size = 20))

#Myeloid cells
plot1 <- VlnPlot(Brain_all, features = c("GPR34","P2RY12","CX3CR1","C1QB","MAMU-DRA","MAMU-DRB1"), raster = F, ncol = 3,
                 pt.size = 0)&theme(axis.text = element_text(size = 20), 
                                    axis.title.y = element_text(size = 20),
                                    axis.title.x = element_blank(),
                                    plot.title = element_text(size = 25))

plot2 <- FeaturePlot(Brain_all, features = c("GPR34","P2RY12","CX3CR1","C1QB","MAMU-DRA","MAMU-DRB1"),raster = F, ncol = 3)&
  theme(axis.title = element_blank(), plot.title = element_text(size = 25), axis.text = element_text(size = 20))

#B cells
plot1 <- VlnPlot(Brain_all, features = c("CD79A", "CD19", "PAX5",'MS4A1'), raster = F,ncol = 2,
                 pt.size = 0)&theme(axis.text = element_text(size = 20), 
                                    axis.title.y = element_text(size = 20))

plot2 <- FeaturePlot(Brain_all, features = c("CD79A", "CD19", "PAX5",'MS4A1'),raster = F, ncol = 2)

# Find the markers
Markers <- FindAllMarkers(Brain_all, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# Rename the clusters based on annotations
new.cluster.id <- c("Microglia","Microglia","Lymphocyte","Lymphocyte","CAM","Lymphocyte",
                    "Microglia","CAM","Cluster8","Microglia","Others","Others")
names(new.cluster.id) <- levels(Brain_all)
Brain_all <- RenameIdents(Brain_all, new.cluster.id)
Brain_all$Label <- Idents(Brain_all)

SaveH5Seurat(Brain_all, filename = "./brain_all_filtered", overwrite = T)


# Plot the UMI count and gene count
Idents(Brain_all) <- "Label"
plot1 <- VlnPlot(Brain_all, features = c("nCount_RNA","nFeature_RNA"), raster = F, ncol = 1, log = T,
                 pt.size = 0)&theme(axis.text = element_text(size = 20), 
                                    axis.title.y = element_text(size = 20),
                                    axis.title.x = element_blank(),
                                    plot.title = element_text(size = 25))

#Subset the brain lymphocyte
brain_T <- subset(Brain_all, idents = "Lymphocyte")

# PCA
brain_T <- FindVariableFeatures(object = brain_T)
brain_T <- ScaleData(object = brain_T)
gene.list <- read.csv('./gene.list.csv')
gene.list <- gene.list$Gene.Symbol
brain_T <- RunPCA(object = brain_T, features = gene.list)

# Determine PC to use
ElbowPlot(brain_T, ndims = 50)

brain_T <- brain_T %>%
  RunHarmony(group.by.vars = 'sample', plot_convergence = F)

#cluster
brain_T <- brain_T %>%
  RunUMAP(reduction = 'harmony', dims = 1:20) %>%
  FindNeighbors(reduction = 'harmony', dims = 1:20) %>%
  FindClusters(res = 0.2)

SaveH5Seurat(brain_T, filename = "./brain_T")

# Plot the embeddings
plot1 <- DimPlot(brain_T, reduction = 'umap', raster = F, label = T, label.size = 8)+
  theme(axis.title = element_text(size = 20), axis.text = element_text(size = 20), legend.position = "none")
plot2 <- DimPlot(brain_T, reduction = 'umap', group.by = 'condition', raster = F)+
  theme(legend.text = element_text(size = 20))
plot3 <- DimPlot(brain_T, reduction = 'umap', group.by = 'name', raster = F)+
  theme(legend.text = element_text(size = 20))

#Find the markers
Markers <- FindAllMarkers(brain_T, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)


#Heatmap
brain_T <- ScaleData(brain_T, features = rownames(brain_T))
Markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 0.5) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10
p1 <- DoHeatmap(brain_T, features = top10$gene, size = 8)+theme(axis.text.y = element_text(size = 20))

# CD4/CD8B expression ratio
df <- FetchData(brain_T, vars = c("CD4","CD8A","CD8B","seurat_clusters"))
df$CD4CD8B <- log2((df$CD4+1)/(df$CD8B+1))
brain_T$CD4CD8B <- df$CD4CD8B

library(RColorBrewer)
plot1 <- FeaturePlot(brain_T, features = "CD4")+
  theme(axis.title = element_blank(), plot.title = element_text(size = 25), axis.text = element_text(size = 20))
plot2 <- FeaturePlot(brain_T, features = "CD8B")+
  theme(axis.title = element_blank(), plot.title = element_text(size = 25), axis.text = element_text(size = 20))
plot3 <- FeaturePlot(brain_T, features = "CD8A")+
  theme(axis.title = element_blank(), plot.title = element_text(size = 25), axis.text = element_text(size = 20))
plot4 <- FeaturePlot(brain_T, features = "CD4CD8B")+scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))+
  theme(axis.title = element_blank(), plot.title = element_text(size = 25), axis.text = element_text(size = 20))

# Plot markers
p1 <- FeaturePlot(brain_T, features = c( "GZMK","LTB","ITGA1",
                                         "LOC720538","LOC711031","IL23R",
                                         "GZMB","PRF1","CCL5",
                                         "KLRC3","FCGR3","GNLY",
                                        "STMN1","MKI67","TOP2A"), ncol = 3,raster = F)&
  theme(axis.title = element_blank(), plot.title = element_text(size = 25), axis.text = element_text(size = 20))

