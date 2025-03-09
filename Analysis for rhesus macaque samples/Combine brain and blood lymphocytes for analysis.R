library(Seurat)
library(SeuratDisk)
library(harmony)
library(dplyr)
library(ggplot2)

blood_T <- LoadH5Seurat("./blood_T.h5seurat")
brain_T <- LoadH5Seurat("./brain_T.h5seurat")

Comb <- merge(blood_T, brain_T)
Comb@meta.data <- Comb@meta.data[,-c(14,15)]
Comb$sample <- paste0(Comb$name, Comb$sort)

# Reperform the workflow
Comb <- FindVariableFeatures(object = Comb)
Comb <- ScaleData(object = Comb)
gene.list <- read.csv('./gene.list.csv') # The gene.list for curated rhesus macaque genes
gene.list <- gene.list$Gene.Symbol
Comb <- RunPCA(object = Comb, features = gene.list)
Comb <- Comb %>%
  RunHarmony(group.by.vars = 'sample', plot_convergence = F, kmeans_init_nstart=20, kmeans_init_iter_max=100)

# Determine PCA
ElbowPlot(Comb, ndims = 50)

# Perform clustering
Comb <- Comb %>%
  RunUMAP(reduction = 'harmony', dims = 1:20) %>%
  FindNeighbors(reduction = 'harmony', dims = 1:20) %>%
  FindClusters(res = 0.2)

plot1 <- DimPlot(Comb, reduction = 'umap', raster = F, label = T, label.size = 8)+
  theme(axis.title = element_text(size = 20), axis.text = element_text(size = 20), legend.position = "none")
plot2 <- DimPlot(Comb, reduction = 'umap', raster = F, label = F, group.by = "sample")+
  theme(axis.title = element_text(size = 20), axis.text = element_text(size = 20))
plot3 <- DimPlot(Comb, reduction = 'umap', raster = F, label = T, label.size = 8, group.by = "tissue")+
  theme(axis.title = element_text(size = 20), axis.text = element_text(size = 20), legend.position = "none")+
  ggtitle("")


SaveH5Seurat(Comb, filename = "/work/foxlab/xiaoke/seurat/T_cells_acute/Comb")

# Find the markers
Markers <- FindAllMarkers(Comb, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# IN Sankey Diagram
library(networkD3)
library(dplyr)
library(htmlwidgets)
BL <- subset(Comb, subset = tissue == "blood")
df.bl <- as.data.frame(table(BL$seurat_clusters, BL$Label))
BR <- subset(Comb, subset = tissue == "brain")
df.br <- as.data.frame(table(BR$seurat_clusters, BR$Label))

df.bl$target <- df.bl$Var1
df.bl$source <- df.bl$Var2
df.bl <- df.bl[,-c(1,2)]

df.br$target <- df.br$Var2
df.br$source <- df.br$Var1
df.br <- df.br[,-c(1,2)]

df <- rbind(df.bl,df.br)

nodes <- data.frame(name = c(as.character(df$targe),as.character(df$source)) %>% unique())
nodes$group <- as.factor(c(rep("new",10), rep("brain",8), rep("blood",8)))
my_color <- 'd3.scaleOrdinal() .domain(["blood", "brain","new"]) .range(["darkred", "darkblue","darkgreen"])'

df$IDsource <- match(df$source, nodes$name)-1 
df$IDtarget <- match(df$target, nodes$name)-1

p <- sankeyNetwork(Links = df, Nodes = nodes,
                   Source = "IDsource", Target = "IDtarget",
                   Value = "Freq", NodeID = "name", colourScale = my_color, NodeGroup = "group",
                   fontSize = 20, nodeWidth = 30)

saveWidget(p, file = "./sankey.html")

#Trajectory Analyses
library(monocle3)
library(SeuratWrappers)
Idents(Comb) <- Comb$Label
#Remove unrelated clusters
sub <- subset(Comb, idents = c("MT","Myeloid T cell","CD8+ cells"), invert = T)
Idents(sub) <- sub$seurat_clusters

#Convert seurat object to cell data set
cds <- as.cell_data_set(sub)
fData(cds)$gene_sn <- rownames(fData(cds))
#Cluster cells by using clustering info in seurat
#Assign partitiion 
partition <- c(rep(1, length(cds@colData@rownames)))
names(partition) <- cds@colData@rownames
partition <- as.factor(partition)
cds@clusters$UMAP$partitions <- partition

#Assign cluster info
list.cluster <- sub@active.ident
cds@clusters$UMAP$clusters <- list.cluster

list.cluster <- sub$Label
cds@clusters$UMAP$clusters <- list.cluster 

#Assign UMAP coordinate -cell embeddings
cds@int_colData@listData$reducedDims$UMAP <- sub@reductions$umap@cell.embeddings

#Trajectory analysis
cds <- learn_graph(cds, use_partition = FALSE)

#pseudotime
cds <- order_cells(cds, reduction_method = 'UMAP', root_cells = colnames(cds[,clusters(cds)=='0']))
plot1 <- plot_cells(cds, color_cells_by = 'cluster',
                    label_cell_group = F,
                    label_branch_points = F,
                    label_roots = F,
                    label_leaves = T,
                    graph_label_size = 5)

#Plot the pseudotime in boxplot
cds$mono3_pseudotime <- pseudotime (cds)
data.pseudo <- as.data.frame(colData(cds))
plot1 <- ggplot(data.pseudo,aes(mono3_pseudotime, reorder(Label, mono3_pseudotime, median), fill= tissue))+
  geom_boxplot()+ylab('cluster')+xlab('pseudotime')+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text = element_text(size = 20), axis.title.x = element_text(size =20), legend.title = element_blank(),
        legend.key.size = unit(2,"cm"), axis.title.y = element_blank(), legend.text = element_text(size = 20))


#Find the markers between brain and blood cells

#For CD4 CTLs
CD4CTL <- subset(Comb, idents = "1")
df <- FetchData(CD4CTL, vars = c("CD4","CD8A","CD8B","seurat_clusters"))
df$CD4CD8B <- log2((df$CD4+1)/(df$CD8B+1))
CD4CTL$CD4CD8B <- df$CD4CD8B
CD4CTL <- subset(CD4CTL, subset = CD4CD8B>0)
#Brain:877 cells Blood:2252 cells
Idents(CD4CTL) <- CD4CTL$tissue
Markers <- FindAllMarkers(CD4CTL, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

Markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 0.5) %>%
  slice_head(n = 20) %>%
  ungroup() -> top20

CD4CTL <- ScaleData(object = CD4CTL,features = rownames(CD4CTL))
p1 <- DoHeatmap(CD4CTL, features = top20$gene)+
  theme(axis.text.y = element_text(size = 15))+
  scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n =4, name = "RdBu")))

#For CD4 EM
EM <- subset(Comb, idents = "3")
df <- FetchData(EM, vars = c("CD4","CD8A","CD8B","seurat_clusters"))
df$CD4CD8B <- log2((df$CD4+1)/(df$CD8B+1))
EM$CD4CD8B <- df$CD4CD8B
EM <- subset(EM, subset = CD4CD8B>0)
Idents(EM) <- EM$tissue

Markers <- FindAllMarkers(EM, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

Markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 0.5) %>%
  slice_head(n = 20) %>%
  ungroup() -> top20

EM <- ScaleData(object = EM,features = rownames(EM))

p1 <- DoHeatmap(EM, features = top20$gene)+
  theme(axis.text.y = element_text(size = 15), legend.position = "bottom", 
        legend.text = element_text(size =20), legend.title = element_text(size = 20),
        legend.size = unit(1,"cm"))+
  scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n =4, name = "RdBu")))
