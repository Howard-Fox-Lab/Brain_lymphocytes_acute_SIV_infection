library(Seurat)
library(SeuratDisk)
library(muscat)
library(limma)
library(SingleCellExperiment)
library(dplyr)

# For brain_T
brain_T <- LoadH5Seurat("./brain_T.h5seurat")
brain_T$condition[which(brain_T$condition == "Acute-Infected")] <- "AcuteInfected"

### DS analysis for brain_T
# Covert seuratobject to sce
sce <- as.SingleCellExperiment(brain_T)
# Data preparation
sce <- prepSCE(sce,
               kid = "Label",
               gid = "condition",
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
pb_mds <- pbMDS(pb)

# Prepare for running DS analysis
ei <- metadata(sce)$experiment_info
mm <- model.matrix(~ 0 + ei$group_id)
dimnames(mm) <- list(ei$sample_id, levels(ei$group_id))
contrast <- makeContrasts("AcuteInfected-Uninfected", levels = mm)

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

# nb. of DS genes & % of total by cluster
n_de <- vapply(tbl_fil, nrow, numeric(1))
p_de <- format(n_de / nrow(sce) * 100, digits = 3)
data.frame("#DS" = n_de, "%DS" = p_de, check.names = FALSE)

df <- bind_rows(tbl_fil)

# Filtering based on expression frequencies
frq <- calcExprFreqs(sce, assay = "counts", th = 0)
gids <- levels(sce$group_id)
frq10 <- vapply(as.list(assays(frq)), 
                function(u) apply(u[, gids] > 0.1, 1, any), 
                logical(nrow(sce)))
tbl_fil2 <- lapply(kids, function(k)
  dplyr::filter(tbl_fil[[k]], 
                gene %in% names(which(frq10[, k]))))
df <- bind_rows(tbl_fil2)

# Visualize the results
library(UpSetR)
library(purrr)
library(ggplot2)
de_gs_by_k <- map(tbl_fil2, "gene")
p1 <- upset(fromList(de_gs_by_k), text.scale = 2)

# Heatmap
p2 <- pbHeatmap(sce, res, top_n = 10)

# top-20 DS genes for EM
trace(pbHeatmap, edit = T)
p3 <- pbHeatmap(sce, res, k = "EM")
p3 <- ComplexHeatmap::draw(p3, heatmap_legend_side = "bottom", annotation_legend_side = "bottom")

# For blood_T
blood_T <- LoadH5Seurat("./blood_T.h5seurat")
blood_T$condition[which(blood_T$condition == "Acute-Infected")] <- "AcuteInfected"
# Pseudobulk-level MDS plot 
# MDS: multi-dimensional scaling
pb_mds <- pbMDS(pb)
# Since blood_T only has one sample in uninfected condition, the comparison cannot be ran