library(Seurat)
library(SeuratDisk)
library(dplyr)
library(ggplot2)
library(DoubletFinder)

# Modify the functions to accomdate the used seurat version
trace(paramSweep, edit = T)
trace(doubletFinder, edit = T)
# Change seu@assays$RNA$counts to seu@assays$RNA@counts

# Load data
Brain_all <- LoadH5Seurat("./brain_all_filtered.h5seurat")
# Create a new column in metadata for storing the doublet info
Brain_all$Doublet <- rep("NA",nrow(Brain_all@meta.data))

# Since the doublet finder need to be ran on each sample separately, we have to split the run
# This code just demonstrate for one sample, for other samples, we changed "sample" in subset
sample1 <- subset(Brain_all, sample == "95T_CD45CD11B_1")
sample1 <- NormalizeData(sample1)
sample1 <- FindVariableFeatures(sample1, selection.method = "vst", nfeatures = 2000)
sample1 <- ScaleData(sample1)
sample1 <- RunPCA(sample1)
sample1 <- RunUMAP(sample1, dims = 1:20)

#pK identification (no ground-truth)
sweep.res.list <- paramSweep(sample1, PCs = 1:20, sct = F)
sweep.stats <- summarizeSweep(sweep.res.list, GT = F)
bcmvn <- find.pK(sweep.stats)
pK <- bcmvn %>%
  filter(BCmetric == max(BCmetric)) %>%
  select(pK)
pK <- as.numeric(as.character(pK[[1]]))

#Homotypic Doublet proportion estimate
annotations <- sample1$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.039*nrow(sample1@meta.data))
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

#Run DoubletFinder 
sample1 <- doubletFinder(sample1,
                         PCs = 1:20,
                         pN = 0.25,
                         pK = pK,
                         nExp = nExp_poi.adj,
                         reuse.pANN =F,
                         sct = F)

# The doublet information will be embedded in metadata
Brain_all$Doublet[which(rownames(Brain_all@meta.data) %in% rownames(sample1@meta.data))] <- sample1@meta.data[,ncol(sample1@meta.data)]
SaveH5Seurat(Brain_all, "/work/foxlab/xiaoke/seurat/T_cells_acute/Merge/brain_all_filtered", overwrite = T)

# Plot the doublet
plot1 <- DimPlot(Brain_all, group.by = "Doublet", raster = F, cols = c("orange","black"))+
  theme(plot.title = element_blank(), legend.text = element_text(size = 18))
# Plot percentage of doublet in each phenotype for all brain cells
pt <- table(Idents(Brain_all), Brain_all$Doublet)
pt <- as.data.frame(pt)
plot1 <- ggplot(pt, aes(x = Var1, y = Freq, fill = Var2)) +
  theme_bw(base_size = 15) +
  geom_col(position = "fill", width = 0.5) +
  theme(legend.title = element_blank())+
  theme(axis.text = element_text(size = 15, face = 'bold', angle = 45),
        axis.title = element_blank(),
        legend.key.size = unit(1, 'cm'),
        legend.text = element_text(size = 18),
        legend.position = "bottom")+
  scale_fill_manual(values = c("orange","black"))

# Plot percentage of doublet in each phenotype for all brain cells
brain_T <- LoadH5Seurat("./brain_T.h5seurat")
pt <- table(brain_T$Label, brain_T$Doublet)
pt <- as.data.frame(pt)
plot2 <- ggplot(pt, aes(x = Var1, y = Freq, fill = Var2)) +
  theme_bw(base_size = 15) +
  geom_col(position = "fill", width = 0.5) +
  theme(legend.title = element_blank())+
  theme(axis.text = element_text(size = 15, face = 'bold', angle = 45),
        axis.title = element_blank(),
        legend.key.size = unit(1, 'cm'),
        legend.text = element_text(size = 18),
        legend.position = "bottom")+
  scale_fill_manual(values = c("orange","black"))


