library(Seurat)
library(SeuratDisk)
library(dplyr)
library(ggplot2)

PLS <- LoadH5Seurat("./PLS_merged.h5seurat")
PLS$condition[which(PLS$condition %in% "SIV")] <- "Acute-Infected"
PLS$condition[which(PLS$condition %in% "H")] <- "Uninfected"
PLS$duration <- rep("21-days",nrow(PLS_comb@meta.data))


# All cells were lymphocyte and no cells need to be removed
# Merge data from brain
brain_T <- LoadH5Seurat("./brain_T.h5seurat")
brain_T$duration <- rep("12-days",nrow(Comb@meta.data))
Comb_brain <- merge(brain_T, PLS)
Comb_brain$duration[which(Comb$name %in% c("C1","C2"))] <- "ctl_21days"
Comb_brain$duration[which(Comb$name %in% c("92T","104T","111T"))] <- "ctl_12days"
SaveH5Seurat(Comb, filename = "./PLSComb_brain")

#---------------Clustering-----------------
Comb_brain <- NormalizeData(object = Comb_brain)
Comb_brain <- FindVariableFeatures(object = Comb_brain)
Comb_brain <- ScaleData(object = Comb_brain)

gene.list <- read.csv('./gene.list.csv')
gene.list <- gene.list$Gene.Symbol
Comb_brain <- RunPCA(object = Comb_brain, features = gene.list)

library(harmony)
Comb_brain <- Comb_brain %>%
  RunHarmony(group.by.vars = 'sample', plot_convergence = F, kmeans_init_nstart=20, kmeans_init_iter_max=100)

# Determine PCA
ElbowPlot(Comb_brain, ndims = 50)

# Perform clustering
Comb_brain <- Comb_brain %>%
  RunUMAP(reduction = 'harmony', dims = 1:20) %>%
  FindNeighbors(reduction = 'harmony', dims = 1:20) %>%
  FindClusters(res = 0.2)

# Plot the embeddings
plot1 <- DimPlot(Comb_brain, reduction = 'umap', raster = F, label = T, label.size = 8)+
  theme(legend.text = element_text(size = 20))
plot2 <- DimPlot(Comb_brain, reduction = 'umap', group.by = 'duration', raster = F)+
  theme(legend.text = element_text(size = 20), plot.title = element_blank())
plot3 <- DimPlot(Comb_brain, reduction = 'umap', group.by = 'name', raster = F)+
  theme(legend.text = element_text(size = 20), plot.title = element_blank())

# Findmarkers
Markers <- FindAllMarkers(Comb_brain, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# Assign annotation for new clusters
library(pheatmap)
df <- table(Comb_brain$Label, Comb_brain$seurat_clusters)
df <- scale(df)
plot4 <- pheatmap(log10(df+10),cluster_rows = T,cluster_cols = F,
                  color = colorRampPalette(c('white','blue'))(10))+
  theme(legend.position = "bottom")


#Rename the clusters
new.cluster.id <- c("TRM","CD4+ CTL","CD8+ Effector","EM","γδ T cells","Proliferating","NK cell","Myeloid T cell")
names(new.cluster.id) <- levels(Comb_brain)
Comb_brain <- RenameIdents(Comb_brain, new.cluster.id)
Comb_brain$Label <- Idents(Comb_brain)

plot1 <- DimPlot(Comb_brain, reduction = 'umap', raster = F, label = F)+
  theme(legend.text = element_text(size = 20))

SaveH5Seurat(Comb_brain, filename = "./PLSComb_brain", overwrite = T)


# Population change
pt <- table(Idents(Comb_brain), Comb_brain$duration)
pt <- as.data.frame(pt)
pt <- pt[which(pt$Var2 %in% c("21-days","ctl_21days")),]
pt <- pt[which(pt$Var2 %in% c("12-days","ctl_12days")),]

plot1 <- ggplot(pt, aes(x = Var2, y = Freq, fill = Var1)) +
  theme_bw(base_size = 15) +
  geom_col(position = "fill", width = 0.5) +
  theme(legend.title = element_blank())+
  theme(axis.text = element_text(size = 15, face = 'bold'),
        axis.title = element_blank(),
        legend.key.size = unit(1, 'cm'),
        legend.text = element_text(size = 18))

########### Pseudobulk analyses-------------------
library(dplyr)
library(ggrepel)
library(tidyr)
Comb_brain <- subset(Comb_brain, condition == "Acute-Infected")
pseudo <- AggregateExpression(Comb_brain, assays = "RNA", return.seurat = T,
                                  group.by = c("duration","name"))
pseudo$celltype <- colnames(pseudo)
pseudo@meta.data <- separate(pseudo@meta.data, col = 'celltype', into = c('duration','name','celltype'), sep = '_')
pseudo$sample <- paste0(pseudo$duration, pseudo$celltype)
Idents(pseudo) <- 'name'

df <- pseudo[["RNA"]]@counts
df <- round(df)
pseudo[["RNA"]]@counts <- df
bulk.de <- FindMarkers(object = pseudo, 
                       ident.1 = "21-days", 
                       ident.2 = "12-days",
                       test.use = "DESeq2")
Markers <- bulk.de %>%
  filter(p_val_adj < 0.05 & abs(avg_log2FC) > 2)


bulk.de$DE <- rep("NO", nrow(bulk.de))
bulk.de$DE[which(rownames(bulk.de) %in% rownames(Markers))] <- "DN"
bulk.de$DE[which(bulk.de$avg_log2FC>2 & bulk.de$p_val_adj < 0.05)] <- 'UP'

LB <- c("KEG06-t18","KEG06-p09","KEG06-p05","KEG06-t13","KEG06-t06","KEG06-p01",
        "CCL20","CXCL10","ITGB6","CXCR5","CCR7","CCR6","IL23A","FOS","NFKB1","FOSB",
        "IFI6")
bulk.de$Label <- rep(NA, nrow(bulk.de))
bulk.de[which(rownames(bulk.de) %in% LB),"Label"] <- rownames(bulk.de)[which(rownames(bulk.de) %in% LB)]

p1 <- ggplot(data=bulk.de, aes(x=avg_log2FC, y=-log10(p_val_adj), col=DE, label = Label)) +
  geom_point() + 
  theme_minimal() +
  scale_color_manual(values=c("blue", "grey", "red")) +
  geom_text_repel(max.overlaps = Inf) +
  geom_vline(xintercept=c(-2, 2), col="black") +
  geom_hline(yintercept=-log10(0.05), col="black")+
  xlim(-max(abs(bulk.de$avg_log2FC)),max(abs(bulk.de$avg_log2FC)))+
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5, size = 15))+
  ggtitle("21-day SIV infection")+
  annotate(geom = "text", x = -7.5, y = 250, label = "Downregulated", size = 14/.pt)+
  annotate(geom = "text", x = 7.5, y = 250, label = "Upregulated",size = 14/.pt)

#Combine SIV positive cells from brain (PLS, FOX) and blood
#Merged with blood_T 
blood_T <- LoadH5Seurat("./blood_T.h5seurat")
blood_T$duration <- rep("12-days",nrow(Comb@meta.data))
comb_PLS <- merge(Comb_brain, blood_T)
SaveH5Seurat(comb_PLS, "./merged_PLS_FOX")

#Recluster the SIV+ Cells
#Remove the myeloid T cell and just kept the lymphocyte clusters
Comb_PLS <- subset(Cob_PLS, Label == "Myeloid T cell", invert = T)
#Identify SIV+ cells
SIV <- FetchData(Comb_PLS, vars = c('SIV','duration',"tissue"))
SIV <- filter(SIV,SIV>0) #total 985 cells
siv.cells <- rownames(SIV)
SIV_all <- subset(Cob_PLS, cells= siv.cells)

#Reclustering on SIV+ cells
SIV_all <- NormalizeData(object = SIV_all)
SIV_all <- FindVariableFeatures(object = SIV_all)
SIV_all <- ScaleData(object = SIV_all)

gene.list <- read.csv('./gene.list.csv')
gene.list <- gene.list$Gene.Symbol
SIV_all <- RunPCA(object = SIV_all, features = gene.list)

SIV_all <- SIV_all %>%
  RunHarmony(group.by.vars = 'sample', plot_convergence = F)

# Determine PCA
ElbowPlot(SIV_all, ndims = 50)

# Perform clustering
SIV_all <- SIV_all %>%
  RunUMAP(reduction = 'harmony', dims = 1:10) %>%
  FindNeighbors(reduction = 'harmony', dims = 1:10) %>%
  FindClusters(res = 0.1)

Markers <- FindAllMarkers(SIV_all, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# Rename the clusters
new.cluster.id <- c("Cytotoic","EM/TRM","Naive","Proliferating")
names(new.cluster.id) <- levels(SIV_all)
SIV_all <- RenameIdents(SIV_all, new.cluster.id)

SaveH5Seurat(SIV_all, "./SIV_all")

# Plot marker genes
p1 <- FeaturePlot(SIV_all, features = c('LEF1','GZMB','MKI67','ICOS'), ncol = 2)

# Plot UMAP
plot1 <- DimPlot(SIV_all, split.by = 'tissue',reduction = 'umap', label = F)+
  theme(axis.title = element_blank(), axis.text = element_blank(), 
        axis.line = element_blank(),axis.ticks = element_blank())

# Plot percentage of SIV+ cells in different tissue
pt <- table(Idents(SIV_all), SIV_all$tissue)
pt <- as.data.frame(pt)
plot1 <- ggplot(pt, aes(x = Var2, y = Freq, fill = Var1)) +
  theme_bw(base_size = 15) +
  geom_col(position = "fill", width = 0.5) +
  theme(legend.title = element_blank())+
  theme(axis.text.x = element_text(size = 15, face = 'bold'),
        axis.title = element_blank(), plot.title = element_text(face = 'bold', size = 20, hjust = 0.5))

# Subset for only brain infected cells
SIV_brain <- subset(SIV_all, tissue == 'brain')
# Plot UMAP
plot2 <- DimPlot(SIV_brain, split.by = 'duration',reduction = 'umap', label = F)+
  theme(axis.title = element_blank(), axis.text = element_blank(), 
        axis.line = element_blank(),axis.ticks = element_blank())
# Plot percentage of SIV+ cells in the brain between 12-days and 21-days infection
pt <- table(Idents(SIV_brain), SIV_brain$duration)
pt <- as.data.frame(pt)
plot1 <- ggplot(pt, aes(x = Var2, y = Freq, fill = Var1)) +
  theme_bw(base_size = 15) +
  geom_col(position = "fill", width = 0.5) +
  theme(legend.title = element_blank())+
  theme(axis.text.x = element_text(size = 15, face = 'bold'),
        axis.title = element_blank(), plot.title = element_text(face = 'bold', size = 20, hjust = 0.5))

