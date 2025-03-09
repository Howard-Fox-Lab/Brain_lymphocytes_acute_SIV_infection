library(Seurat)
library(SeuratDisk)
library(dplyr)
library(ggplot2)

brain_T <- LoadH5Seurat("./brain_T.h5seurat")

# Rename the clusters based on annotations
new.cluster.id <- c("TRM","CD4+ CTL","Proliferating","CD8+ Effector","EM","γδ T cells","Myeloid T cell","NK cell")
names(new.cluster.id) <- levels(brain_T)
brain_T <- RenameIdents(brain_T, new.cluster.id)
brain_T$Label <- Idents(brain_T)

# Calculate percetage of each cluster in infected vs uninfected conditions
pt <- table(Idents(brain_T), brain_T$condition)
pt <- as.data.frame(pt)
plot1 <- ggplot(pt, aes(x = Var2, y = Freq, fill = Var1)) +
  theme_bw(base_size = 15) +
  geom_col(position = "fill", width = 0.5) +
  theme(legend.title = element_blank())+
  theme(axis.text = element_text(size = 20, face = 'bold'),
        axis.title = element_blank(),
        legend.key.size = unit(1, 'cm'),
        legend.text = element_text(size = 18))

# Compare the CD8 Effector and EM
CD8 <- subset(brain_T, idents = c("CD8+ Effector","EM"))
Markers <- FindAllMarkers(CD8, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
plot1 <- VlnPlot(CD8, features = c("KLRC3","GNLY","CCL5","FCGR3",
                                   "GZMB","GZMH","GZMA","GZMK",
                                   "PAG1","THEMIS","BCL11B","ICOS"),raster = F, ncol = 4,pt.size = 0)&
  theme(axis.text = element_text(size = 20),axis.title = element_blank(), plot.title = element_text(size = 25))

# Compare the CD8 Effector and CD4+CTL
cyto <- subset(brain_T, idents = c("CD8+ Effector","CD4+ CTL"))
Markers <- FindAllMarkers(cyto, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
plot1 <- VlnPlot(cyto, features = c("IL7R","CD6","RORA","ITM2A","KLRC3","GNLY","GZMA","FCGR3"),raster = F, ncol = 4,pt.size = 0)&
  theme(axis.text = element_text(size = 20),axis.title = element_blank(),plot.title = element_text(size = 25))

# Compare the TCR expression between each lymphocyte cluster
plot1 <- RidgePlot(brain_T, features = c("LOC710951","LOC703029","CD3D","CD3E","GZMB","IL7R"),ncol = 3)&
  theme(axis.title = element_blank(),legend.position = "none", plot.title = element_text(size = 25))

# Compare the UMI count and gene count between each lymphocyte cluster
plot1 <- VlnPlot(brain_T,features = c("nCount_RNA","nFeature_RNA"),raster = F, ncol = 1,pt.size = 0, log = T)&
  theme(axis.text = element_text(size = 20),axis.title = element_blank(),plot.title = element_text(size =25))
png(filename ="/work/foxlab/xiaoke/seurat/plots/Figure 1C.png", width = 1500, height = 3000, units = 'px',res = 300)
plot1
dev.off()

# Compare the Infected vs Uninfected conditions
library(cowplot)
theme_set(theme_cowplot())
options(ggrepel.max.overlaps = Inf)
#-------------For TRM------------------
TRM <- subset(brain_T, idents = "TRM")
Idents(TRM) <- TRM$condition
Markers <- FindAllMarkers(TRM, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

df <- FetchData(TRM, vars = c("condition","CD4CD8B"))
df.SIV <- filter(df, condition == "Acute-Infected")
df.H <- filter(df, condition == "Uninfected")

avg <- as.data.frame(log1p(AverageExpression(TRM, verbose = F)$RNA))
avg$gene <- rownames(TRM)

gene.to.label <- c("GAMB","GZMA","GNLY","IFI27","IFI6","STMN1","TOP2A","MKI67","MAMU-DRA","CD74")
colnames(avg) <- c("H","SIV","gene")

p1 <- ggplot(avg, aes(H, SIV)) +
  geom_point()+
  geom_point(data = subset(avg, gene%in%gene.to.label), color = 'red') + 
  ggtitle("TRM (Brain)")+
  xlim (0,5)+
  ylim(0,5)+
  xlab("Uninfected")+
  ylab("Acute-Infected")+
  theme(axis.text = element_text(size = 18))

p1 <- LabelPoints(plot = p1,col= 'red',points = gene.to.label, 
                  repel = TRUE, xnudge = 0, ynudge = 0, fontface = "bold", size = 5)

#--------------For CD4+ CTL------------------
CD4 <- subset(brain_T, idents = "CD4+ CTL")
Idents(CD4) <- CD4$condition
Markers <- FindAllMarkers(CD4, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

gene.to.label <- c("ADORA2B","MAMU-DRA","MAMU-DQA1","RPL37","RPS20","IFI6","RPL35","MAMU-DRB1","CD74","IFI27")
avg <- as.data.frame(log1p(AverageExpression(CD4, verbose = F)$RNA))
avg$gene <- rownames(CD4)
colnames(avg) <- c("H","SIV","gene")

p1 <- ggplot(avg, aes(H, SIV)) +
  geom_point()+
  geom_point(data = subset(avg, gene%in%gene.to.label), color = 'red') + 
  ggtitle("CD4+ CTL (Brain)")+
  xlim (0,5)+
  ylim(0,5)+
  xlab("Uninfected")+
  ylab("Acute-Infected")+
  theme(axis.text = element_text(size = 18))

p1 <- LabelPoints(plot = p1,col= 'red',points = gene.to.label, 
                  repel = TRUE, xnudge = 0, ynudge = 0, fontface = "bold", size = 5)

plot1 <- VlnPlot(CD4, features = c("GZMB","PRF1","NKG7"),raster = F, ncol = 3,
                 pt.size = 0)&theme(axis.text = element_text(size = 20), 
                                    axis.title.y = element_text(size = 20))

#-------------For CD8+ effector--------
CD8E <- subset(brain_T, idents = "CD8+ Effector")
Idents(CD8E) <- CD8E$condition
Markers <- FindAllMarkers(CD8E, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

gene.to.label <- c("GAMB","GZMA","GNLY","IFI27","IFI6","STMN1","TOP2A","MKI67","MAMU-DRA","CD74","KLRC1","KLRC2","KLRC3","FCER1G")
avg <- as.data.frame(log1p(AverageExpression(CD8E, verbose = F)$RNA))
avg$gene <- rownames(CD8E)
colnames(avg) <- c("H","SIV","gene")

p1 <- ggplot(avg, aes(H, SIV)) +
  geom_point()+
  geom_point(data = subset(avg, gene%in%gene.to.label), color = 'red') + 
  ggtitle("CD8+ Effector (Brain)")+
  xlim (0,5)+
  ylim(0,5)+
  xlab("Uninfected")+
  ylab("Acute-Infected")+
  theme(axis.text = element_text(size = 18))

p1 <- LabelPoints(plot = p1,col= 'red',points = gene.to.label, 
                  repel = TRUE, xnudge = 0, ynudge = 0, fontface = "bold", size = 5)

#--------------For EM------------------------
EM <- subset(brain_T, idents = "EM")
Idents(EM) <- EM$condition

plot1 <- DimPlot(EM, reduction = 'umap', raster = F, label = T, label.size = 8, split.by = "condition")+
  theme(axis.title = element_text(size = 20), axis.text = element_text(size = 20), legend.position = "none")

plot1 <- FeaturePlot(EM, features = "CD8B")+
  theme(axis.title = element_blank(), plot.title = element_text(size = 25), axis.text = element_text(size = 20))
plot2 <- FeaturePlot(EM, features = "CD4CD8B")+scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))+
  theme(axis.title = element_blank(), plot.title = element_text(size = 25), axis.text = element_text(size = 20))

Markers <- FindAllMarkers(EM, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

gene.to.label <- c("GAMB","GZMA","GNLY","IFI27","IFI6","STMN1","TOP2A","MKI67","MAMU-DRA","CD74","CD8A","CD8B")
avg <- as.data.frame(log1p(AverageExpression(EM, verbose = F)$RNA))
avg$gene <- rownames(EM)
colnames(avg) <- c("H","SIV","gene")

p1 <- ggplot(avg, aes(H, SIV)) +
  geom_point()+
  geom_point(data = subset(avg, gene%in%gene.to.label), color = 'red') + 
  ggtitle("EM (Brain)")+
  xlim (0,5)+
  ylim(0,5)+
  xlab("Uninfected")+
  ylab("Acute-Infected")+
  theme(axis.text = element_text(size = 18))

p1 <- LabelPoints(plot = p1,col= 'red',points = gene.to.label, 
                  repel = TRUE, xnudge = 0, ynudge = 0, fontface = "bold", size = 5)

#--------------For Myeloid T cell------------------
Myeloid <- subset(brain_T, idents = "Myeloid T cell")
Idents(Myeloid) <- Myeloid$condition
Markers <- FindAllMarkers(Myeloid, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

gene.to.label <- c("GZMB","APOE","IFI27","CTSC","MKI67","IFI6","CD5L","STMN1")
avg <- as.data.frame(log1p(AverageExpression(Myeloid, verbose = F)$RNA))
avg$gene <- rownames(Myeloid)
colnames(avg) <- c("H","SIV","gene")

p1 <- ggplot(avg, aes(H, SIV)) +
  geom_point()+
  geom_point(data = subset(avg, gene%in%gene.to.label), color = 'red') + 
  ggtitle("Myeloid T cell (Brain)")+
  xlim (0,5)+
  ylim(0,5)+
  xlab("Uninfected")+
  ylab("Acute-Infected")+
  theme(axis.text = element_text(size = 18))

p1 <- LabelPoints(plot = p1,col= 'red',points = gene.to.label, 
                  repel = TRUE, xnudge = 0, ynudge = 0, fontface = "bold", size = 5)

#-----------For Proliferating cell------------
pro <- subset(brain_T, idents = "Proliferating")
Idents(pro) <- pro$condition
Markers <- FindAllMarkers(pro, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

gene.to.label <- c("IFI27","GZMA","GZMK","CD38","GZMB","MAMU-DRA","IL12RB2","KLRC3","CD74")
avg <- as.data.frame(log1p(AverageExpression(pro, verbose = F)$RNA))
avg$gene <- rownames(pro)
colnames(avg) <- c("H","SIV","gene")

p1 <- ggplot(avg, aes(H, SIV)) +
  geom_point()+
  geom_point(data = subset(avg, gene%in%gene.to.label), color = 'red') + 
  ggtitle("Proliferating Cell (Brain)")+
  xlim (0,5)+
  ylim(0,5)+
  xlab("Uninfected")+
  ylab("Acute-Infected")+
  theme(axis.text = element_text(size = 18))

p1 <- LabelPoints(plot = p1,col= 'red',points = gene.to.label, 
                  repel = TRUE, xnudge = 0, ynudge = 0, fontface = "bold", size = 5)

#----------------For NK cell---------------------------
NK <- subset(brain_T, idents = "NK cell")
Idents(NK) <- NK$condition
Markers <- FindAllMarkers(NK, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

gene.to.label <- c("GZMB","IFI27","CD74","MAMU-DRA","GNLY","CX3CR1","GZMH","PRF1","KLRC1","KLRC2","KLRC3","FCER1G")
avg <- as.data.frame(log1p(AverageExpression(NK, verbose = F)$RNA))
avg$gene <- rownames(NK)
colnames(avg) <- c("H","SIV","gene")

p1 <- ggplot(avg, aes(H, SIV)) +
  geom_point()+
  geom_point(data = subset(avg, gene%in%gene.to.label), color = 'red') + 
  ggtitle("NK cell (Brain)")+
  xlim (0,5)+
  ylim(0,5)+
  xlab("Uninfected")+
  ylab("Acute-Infected")+
  theme(axis.text = element_text(size = 18))

p1 <- LabelPoints(plot = p1,col= 'red',points = gene.to.label, 
                  repel = TRUE, xnudge = 0, ynudge = 0, fontface = "bold", size = 5)

#----------------For gamma/delta T -----------------
gd <- subset(brain_T, idents = "γδ T cells")
Idents(gd) <- gd$condition
write.csv(Markers, file = '/work/foxlab/xiaoke/seurat/plots/gd_marker.csv')

gene.to.label <- c("GZMB","IFI27","GNLY","KLRB1","IFI27L2","CCL5","CD52","IFI6","NKG7")
avg <- as.data.frame(log1p(AverageExpression(gd, verbose = F)$RNA))
avg$gene <- rownames(gd)
colnames(avg) <- c("H","SIV","gene")

p1 <- ggplot(avg, aes(H, SIV)) +
  geom_point()+
  geom_point(data = subset(avg, gene%in%gene.to.label), color = 'red') + 
  ggtitle("γδ  cell (Brain)")+
  xlim (0,5)+
  ylim(0,5)+
  xlab("Uninfected")+
  ylab("Acute-Infected")+
  theme(axis.text = element_text(size = 18))

p1 <- LabelPoints(plot = p1,col= 'red',points = gene.to.label, 
                  repel = TRUE, xnudge = 0, ynudge = 0, fontface = "bold", size = 5)
#--------------------------------------------
# Comparing the the expression of cytotoxic genes
Uninfect <- subset(brain_T, condition == "Uninfected")
Infect <- subset(brain_T, condition == "Acute-Infected")

genes <- c("GZMA","GZMB","GZMK","GZMH","GNLY","MAMU-DRA","CD74","IFI27","IFI27L2","IFI6")
p1 <- DotPlot(Infect, features = genes, cols = c("blue","red"),scale = F)+
  theme(axis.title = element_blank(), axis.text.x = element_text(angle = 90), 
        plot.margin = unit(c(0.2,0.2,0.2,0.2),"cm"), plot.title = element_text(hjust = 0.5))+
  ggtitle("Acute-Infected")
png(filename ="/work/foxlab/xiaoke/seurat/plots/infect.png", width = 2000, height = 1500, units = 'px',res = 300)
p1
dev.off()

p2 <- DotPlot(Uninfect, features = genes, cols = c("blue","red"), scale = F)+
  theme(axis.title = element_blank(), axis.text.x = element_text(angle = 90), 
        plot.margin = unit(c(0.2,0.2,0.2,0.2),"cm"), plot.title = element_text(hjust = 0.5))+
  ggtitle("Uninfected")

# Analysis on SIV+ cells
SIV <- FetchData(brain_T, vars = c("SIV","ident","condition"))
SIV <- filter(SIV, condition == "Acute-Infected")
SIV.pos <- filter(SIV, SIV>0)
cell <- rownames(SIV.pos)
plot1 <- DimPlot(brain_T, cells.highlight = cell, sizes.highlight = 2)

# Heatmap for CD4 gene
avg <- AverageExpression(brain_T, return.seurat = T)
p1 <- DoHeatmap(avg, features = c("CD4","CCR5","ITGAL","ITGB2","ITGA4","ITGB1"), draw.lines = F)+
  theme(axis.text.y = element_text(size = 20), legend.position = "bottom")+
  scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n =4, name = "RdBu")))
png(filename ="/work/foxlab/xiaoke/seurat/plots/heatmap_CD4.png", width = 4500, height = 4000, units = 'px',res = 300)
p1
dev.off()


