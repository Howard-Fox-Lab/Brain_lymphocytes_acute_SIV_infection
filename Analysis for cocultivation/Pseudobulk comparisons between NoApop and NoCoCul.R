library(Seurat)
library(SeuratDisk)
library(muscat)
library(limma)
library(SingleCellExperiment)
library(ggplot2)


merged_seurat_filtered <- LoadH5Seurat("./merged_seurat_filtered.h5seurat")

# Subset for Mac_T cluster
Idents(merged_seurat_filtered) <- merged_seurat_filtered$seurat_clusters
Mac_T <- subset(merged_seurat_filtered, idents = '5')
Idents(Mac_T) <- Mac_T$Label
# Muscat comparison
# Covert seuratobject to sce
sce <- as.SingleCellExperiment(Mac_T)
# Data preparation
sce <- prepSCE(sce,
               kid = "Label",
               gid = "treatment",
               sid = "name",
               drop = T)

nk <- length(kids <- levels(sce$cluster_id))
ns <- length(sids <- levels(sce$sample_id))
names(kids) <- kids
names(sids) <- sids

# Generate pseudobulk data
pb <- aggregateData(sce,
                    assay = "logcounts", fun = "sum",
                    by = c("cluster_id", "sample_id"))

# Pseudobulk-level MDS plot 
# MDS: multi-dimensional scaling
trace(pbMDS, edit = T) # change the code to color the plot by group_id
pb_mds <- pbMDS(pb)

# Prepare for running DS analysis
ei <- metadata(sce)$experiment_info
mm <- model.matrix(~ 0 + ei$group_id)
dimnames(mm) <- list(ei$sample_id, levels(ei$group_id))
contrast <- makeContrasts("NoApop-NoCoCul", levels = mm)

# Run DS analysis
res <- pbDS(pb, design = mm, contrast = contrast)

# Access results table for 1st comparison
tbl <- res$table[[1]]
names(tbl)

# Filtering the results
tbl_fil <- lapply(tbl, function(u) {
  u <- dplyr::filter(u, p_adj.loc < 0.05, abs(logFC) > 0.25)
  dplyr::arrange(u, p_adj.loc)
})

tbl_fil <- as.data.frame(tbl_fil)
# nb. of DS genes & % of total by cluster
n_de <- vapply(tbl_fil, nrow, numeric(1))
p_de <- format(n_de / nrow(sce) * 100, digits = 3)
data.frame("#DS" = n_de, "%DS" = p_de, check.names = FALSE)

# Plot DS gene (top 20)
trace(pbHeatmap, edit = T) # change the code to color the plot by group_id
p1 <- pbHeatmap(sce, res, k = "Mac_T")
p2 <- ComplexHeatmap::draw(p1, heatmap_legend_side = "bottom", annotation_legend_side = "bottom")

# FeaturePlot of some of the genes
geneList <- c("IER3","CXCL3","CCL4","CCL3","CCL4L2","CXCL2")
p1 <- FeaturePlot(Mac_T, features = geneList, split.by = "treatment", ncol =5)&theme(plot.margin = unit(c(0,0,0,0),"cm"))

# Visualize the expression of CCL in brain_T
brain_all <- LoadH5Seurat("./brain_all_filtered.h5seurat")
Idents(brain_all) <- brain_all$Label

Cluster8 <- subset(brain_all, idents = "Cluster8")
gene1_expr <- FetchData(Cluster8, vars = "CCL3")
gene2_expr <- FetchData(Cluster8, vars = "CCL4L1")

# Define categories
Cluster8$GeneExpressionGroup <- "None"
Cluster8$GeneExpressionGroup[gene1_expr > 0 & gene2_expr == 0] <- "CCL3+"
Cluster8$GeneExpressionGroup[gene1_expr == 0 & gene2_expr > 0] <- "CCL4+"
Cluster8$GeneExpressionGroup[gene1_expr > 0 & gene2_expr > 0] <- "Both+"

# Calculate proportions
group_counts <- table(Cluster8$GeneExpressionGroup)
total_cells <- sum(group_counts)
group_percent <- round((group_counts / total_cells) * 100, 2)

plot_data <- FetchData(Cluster8, vars = c("CCL3", "CCL4L1", "GeneExpressionGroup"), slot = "data")

p <- ggplot(plot_data, aes(x = CCL3, y = CCL4L1, color = GeneExpressionGroup)) +
  geom_point(alpha = 0.6) +
  scale_color_manual(values = c("None" = "gray", "CCL3+" = "blue", "CCL4+" = "red", "Both+" = "purple")) +
  theme_bw() +
  labs(title = "",
       subtitle = paste("None:", group_percent["None"], "% |",
                        "CCL3+:", group_percent["CCL3+"], "% |",
                        "CCL4+:", group_percent["CCL4+"], "% |",
                        "Both+:", group_percent["Both+"], "%"),
       x = "CCL3 Expression",
       y = "CCL4L1 Expression") +
  theme(legend.position = "none", 
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 12),
        panel.border = element_rect(colour = "black", linewidth = 1.5))

# In Mac_T cluster
gene1_expr <- FetchData(Mac_T, vars = "CCL3")
gene2_expr <- FetchData(Mac_T, vars = "CCL4")

# Define categories
Mac_T$GeneExpressionGroup <- "None"
Mac_T$GeneExpressionGroup[gene1_expr > 0 & gene2_expr == 0] <- "CCL3+"
Mac_T$GeneExpressionGroup[gene1_expr == 0 & gene2_expr > 0] <- "CCL4+"
Mac_T$GeneExpressionGroup[gene1_expr > 0 & gene2_expr > 0] <- "Both+"

# Calculate proportions
group_counts <- table(Mac_T$GeneExpressionGroup)
total_cells <- sum(group_counts)
group_percent <- round((group_counts / total_cells) * 100, 2)

plot_data <- FetchData(Mac_T, vars = c("CCL3", "CCL4", "GeneExpressionGroup"), slot = "data")

p <- ggplot(plot_data, aes(x = CCL3, y = CCL4, color = GeneExpressionGroup)) +
  geom_point(alpha = 0.6) +
  scale_color_manual(values = c("None" = "gray", "CCL3+" = "blue", "CCL4+" = "red", "Both+" = "purple")) +
  theme_bw() +
  labs(title = "",
       subtitle = paste("None:", group_percent["None"], "% |",
                        "CCL3+:", group_percent["CCL3+"], "% |",
                        "CCL4+:", group_percent["CCL4+"], "% |",
                        "Both+:", group_percent["Both+"], "%"),
       x = "CCL3 Expression",
       y = "CCL4 Expression") +
  theme(legend.position = "none", 
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 12),
        panel.border = element_rect(colour = "black", linewidth = 1.5))

# Plot the percentage of each treatment with none expression of CCL3 and CCL4
data <- data.frame(Category = c("Campto","CD95","NoApop","NoCoCul"),
                   Count = c(24,4,4,138))
data$Percentage <- round((data$Count / sum(data$Count)) * 100, 1)
data$Label <- paste0(data$Category, "\n", data$Percentage, "%")

p <- ggplot(data, aes(x = "", y = Count, fill = Category)) +
  geom_bar(stat = "identity", width = 1) +  # Create a bar chart
  coord_polar(theta = "y") +  # Convert to pie chart
  theme_void() +
  theme(legend.position = "none")
