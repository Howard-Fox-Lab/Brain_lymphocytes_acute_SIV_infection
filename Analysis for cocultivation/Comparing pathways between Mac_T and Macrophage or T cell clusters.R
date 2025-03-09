library(Seurat)
library(SeuratDisk)
library(clusterProfiler)
library(ggplot2)
library(dplyr)
library(tidyr)
library(enrichplot)
library("org.Hs.eg.db")

merged_seurat_filtered <- LoadH5Seurat("./merged_seurat_filtered.h5seurat")
Idents(merged_seurat_filtered) <- merged_seurat_filtered$seurat_clusters

# Pathway analysis for macrophage clusters plus Mac_T cluster
# Subset for macrophage cluters plus Mac_T cluster
subset <- subset(merged_seurat_filtered, idents = c(1,3,5))
# Perform DEGs
Idents(subset) <- subset$Label
markers <- FindAllMarkers(subset, only.pos = TRUE,
                          min.pct = 0.25, logfc.threshold = 0.25)
markersList <- markers %>% 
  filter(avg_log2FC > 0.5 & p_val_adj < 0.05) %>%
  group_by(cluster)

gene.cluster <- markersList[,6:7]
sample <- split(gene.cluster$gene, gene.cluster$cluster)
sample$`SPP1_Mac` <- bitr(sample$`SPP1_Mac`, fromType = "SYMBOL", toType = 'ENTREZID', OrgDb = "org.Hs.eg.db")
sample$`HLAII_Mac` <- bitr(sample$`HLAII_Mac`, fromType = "SYMBOL", toType = 'ENTREZID', OrgDb = "org.Hs.eg.db")
sample$`Mac_T` <- bitr(sample$`Mac_T`, fromType = "SYMBOL", toType = 'ENTREZID', OrgDb = "org.Hs.eg.db")

genelist <- list('SPP1_Mac'= sample$`SPP1_Mac`$ENTREZID,
                 'HLAII_Mac'= sample$`HLAII_Mac`$ENTREZID,
                 'Mac_T'= sample$`Mac_T`$ENTREZID)
                 
Gocluster2 <- compareCluster(geneCluster = genelist, fun = 'enrichGO', OrgDb = "org.Hs.eg.db", ont = "BP")

plot2 <- dotplot(Gocluster2, showCategory=8, label_format = function(x) stringr::str_wrap(x, width = 50))+
  theme(axis.title = element_blank(),axis.text.x = element_text(angle = 45, color = "black"))+
  scale_fill_gradient(high = "blue",low = "red")

# Pathway analysis for T cell clusters plus Mac_T cluster
# Subset for T cell cluters plus Mac_T cluster
subset <- subset(merged_seurat_filtered, idents = c(0,2,5,7))

# Perform DEGs
Idents(subset) <- subset$Label
markers <- FindAllMarkers(subset, only.pos = TRUE,
                          min.pct = 0.25, logfc.threshold = 0.25)
markersList <- markers %>% 
  filter(avg_log2FC > 0.5 & p_val_adj < 0.05) %>%
  group_by(cluster)

gene.cluster <- markersList[,6:7]
sample <- split(gene.cluster$gene, gene.cluster$cluster)
sample$`Naive_T` <- bitr(sample$`Naive_T`, fromType = "SYMBOL", toType = 'ENTREZID', OrgDb = "org.Hs.eg.db")
sample$`Proli_T` <- bitr(sample$`Proli_T`, fromType = "SYMBOL", toType = 'ENTREZID', OrgDb = "org.Hs.eg.db")
sample$`Mac_T` <- bitr(sample$`Mac_T`, fromType = "SYMBOL", toType = 'ENTREZID', OrgDb = "org.Hs.eg.db")
sample$`Activated_T` <- bitr(sample$`Activated_T`, fromType = "SYMBOL", toType = 'ENTREZID', OrgDb = "org.Hs.eg.db")

genelist <- list('Naive_T'= sample$`Naive_T`$ENTREZID,
                 'Proli_T'= sample$`Proli_T`$ENTREZID,
                 'Mac_T'= sample$`Mac_T`$ENTREZID,
                 'Activated_T' = sample$`Activated_T`$ENTREZID)

Gocluster2 <- compareCluster(geneCluster = genelist, fun = 'enrichGO', OrgDb = "org.Hs.eg.db", ont = "BP")

plot2 <- dotplot(Gocluster2, showCategory=8, label_format = function(x) stringr::str_wrap(x, width = 50))+
  theme(axis.title = element_blank(),axis.text.x = element_text(angle = 45, color = "black"))+
  scale_fill_gradient(high = "blue",low = "red")
