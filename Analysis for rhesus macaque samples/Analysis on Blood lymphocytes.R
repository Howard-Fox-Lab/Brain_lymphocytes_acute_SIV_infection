library(Seurat)
library(SeuratDisk)
library(dplyr)
library(ggplot2)

blood_T <- LoadH5Seurat("./blood_T.h5seurat")
# Rename the clusters based on annotation
new.cluster.id <- c("CD4+ Naive","CD4+ EM","CD4+ CTL(C2)","CD4+ Treg","CD4+ Proliferating","CD8+ cells","MT", "CD4+ CTL(C7)")
names(new.cluster.id) <- levels(blood_T)
blood_T <- RenameIdents(blood_T, new.cluster.id)
blood_T$Label <- Idents(blood_T)
SaveH5Seurat(blood_T, filename = "./blood_T", overwrite = T)

# Calculate the percentage of each cluster between uninfected and acute-infected conditions
pt <- table(Idents(blood_T), blood_T$condition)
pt <- as.data.frame(pt)
plot1 <- ggplot(pt, aes(x = Var2, y = Freq, fill = Var1)) +
  theme_bw(base_size = 15) +
  geom_col(position = "fill", width = 0.5) +
  theme(legend.title = element_blank())+
  theme(axis.text = element_text(size = 15, face = 'bold'),
        axis.title = element_blank(),
        legend.key.size = unit(1, 'cm'),
        legend.text = element_text(size = 18))

#Find the markers between C2 and C7
CTL <- subset(blood_T, idents = c(2,7))
Markers <- FindAllMarkers(CTL, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

plot1 <- VlnPlot(CTL, features = c("GZMB","NKG7","CCL5","LEF1","LTB","SELL"),raster = F, ncol = 3,pt.size = 0)&
  theme(axis.text = element_text(size = 20),axis.title = element_blank(), plot.title = element_text(size = 25))

# Compare the Infected vs Uninfected conditions
library(cowplot)
theme_set(theme_cowplot())
options(ggrepel.max.overlaps = Inf)
#-------------For CD4+ Naive-------------------------
Naive <- subset(blood_T, idents = "CD4+ Naive")
Idents(Naive) <- Naive$condition
Markers <- FindAllMarkers(Naive, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

gene.to.label <- c("LOC693820","LOC114677644","LOC694372","MAMU-A","RPL37","RPL41","RPL37A","CD3D","CD7")
avg <- as.data.frame(log1p(AverageExpression(Naive, verbose = F)$RNA))
avg$gene <- rownames(Naive)
colnames(avg) <- c("H","SIV","gene")

p1 <- ggplot(avg, aes(H, SIV)) +
  geom_point()+
  geom_point(data = subset(avg, gene%in%gene.to.label), color = 'red') + 
  ggtitle("CD4+ Naive (Blood)")+
  xlim (0,5)+
  ylim(0,5)+
  xlab("Uninfected")+
  ylab("Acute-Infected")+
  theme(axis.text = element_text(size = 18))

p1 <- LabelPoints(plot = p1,col= 'red',points = gene.to.label, 
                  repel = TRUE, xnudge = 0, ynudge = 0, fontface = "bold", size = 5)

#-------------For CD4+ EM-------------------------
EM <- subset(blood_T, idents = "CD4+ EM")
Idents(EM) <- EM$condition
Markers <- FindAllMarkers(EM, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

gene.to.label <- c("LOC694372",'LOC719250',"MAMU-A","RPL37","RPS9","RPL41","RPL34","RPL36AL","RPS28","CD3D","CD7")
avg <- as.data.frame(log1p(AverageExpression(EM, verbose = F)$RNA))
avg$gene <- rownames(EM)
colnames(avg) <- c("H","SIV","gene")

p1 <- ggplot(avg, aes(H, SIV)) +
  geom_point()+
  geom_point(data = subset(avg, gene%in%gene.to.label), color = 'red') + 
  ggtitle("CD4+ EM (Blood)")+
  xlim (0,5)+
  ylim(0,5)+
  xlab("Uninfected")+
  ylab("Acute-Infected")+
  theme(axis.text = element_text(size = 18))

p1 <- LabelPoints(plot = p1,col= 'red',points = gene.to.label, 
                  repel = TRUE, xnudge = 0, ynudge = 0, fontface = "bold", size = 5)

#-------------For CD4+ Treg-------------------------
Treg <- subset(blood_T, idents = "CD4+ Treg")
Idents(Treg) <- Treg$condition
Markers <- FindAllMarkers(Treg, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

gene.to.label <- c("LOC114677644","LOC694372",'LOC719250',"MAMU-A","RPL37","RPL32","RPL36AL","CD7","CD3D")
avg <- as.data.frame(log1p(AverageExpression(Treg, verbose = F)$RNA))
avg$gene <- rownames(Treg)
colnames(avg) <- c("H","SIV","gene")

p1 <- ggplot(avg, aes(H, SIV)) +
  geom_point()+
  geom_point(data = subset(avg, gene%in%gene.to.label), color = 'red') + 
  ggtitle("CD4+ Treg (Blood)")+
  xlim (0,5)+
  ylim(0,5)+
  xlab("Uninfected")+
  ylab("Acute-Infected")+
  theme(axis.text = element_text(size = 18))

p1 <- LabelPoints(plot = p1,col= 'red',points = gene.to.label, 
                  repel = TRUE, xnudge = 0, ynudge = 0, fontface = "bold", size = 5)

#-------------For CD4+ CTL (C2/7)-------------------------
CTL <- subset(blood_T, idents = "CD4+ CTL(C7)")
Idents(CTL) <- CTL$condition
Markers <- FindAllMarkers(CTL, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

gene.to.label <- c("CCL5","RPL32","RPL37","GZMB","NKG7","GZMM","MAMU-DRA","MAMU-DQA1","IFI27L2","IFI6")
avg <- as.data.frame(log1p(AverageExpression(CTL, verbose = F)$RNA))
avg$gene <- rownames(CTL)
colnames(avg) <- c("H","SIV","gene")

p1 <- ggplot(avg, aes(H, SIV)) +
  geom_point()+
  geom_point(data = subset(avg, gene%in%gene.to.label), color = 'red') + 
  ggtitle("CD4+ CTL_C7 (Blood)")+
  xlim (0,5)+
  ylim(0,5)+
  xlab("Uninfected")+
  ylab("Acute-Infected")+
  theme(axis.text = element_text(size = 18))

p1 <- LabelPoints(plot = p1,col= 'red',points = gene.to.label, 
                  repel = TRUE, xnudge = 0, ynudge = 0, fontface = "bold", size = 5)
#----------------------------------------
# SIV+ cells in the blood
SIV <- FetchData(blood_T, vars = c("SIV","ident","condition"))
SIV <- filter(SIV, condition == "Acute-Infected")
SIV.pos <- filter(SIV, SIV>0)
cell <- rownames(SIV.pos)
plot1 <- DimPlot(blood_T, cells.highlight = cell, sizes.highlight = 2)

# Compare SIV+ vs SIV- cells
Infect <- subset(blood_T, subset = condition == "Acute-Infected")

SIV <- FetchData(Infect, vars = c("SIV"))
siv.cell <- rownames(filter(SIV, SIV>0))
SIV$Infection <- rep("SIV-",18308)
id <- which(rownames(SIV)%in%siv.cell)
SIV[id,2]<-"SIV+"
Infect$Infection <- SIV$Infection
Idents(Infect) <- Infect$Infection

Markers <- FindAllMarkers(Infect, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

Markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 0.5) %>%
  slice_head(n = 20) %>%
  ungroup() -> top20

Infect <- ScaleData(object = Infect,features = rownames(Infect))

p1 <- DoHeatmap(Infect, features = top20$gene)+
  theme(axis.text.y = element_text(size = 15))+
  scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n =4, name = "RdBu")))

# Heatmap for CD4,CCR5 and integrin genes
library(pheatmap)
avg <- AverageExpression(blood_T, return.seurat = F)
gene <- c("CD4","CCR5","ITGAL","ITGB2","ITGA4","ITGB1")
df <- data.frame("phenotype" = c("non-CTL","non-CTL","CTL","non-CTL","CTL","CTL","non-CTL","CTL"))
rownames(df) <- colnames(avg$RNA)
df2 <- avg$RNA[which(rownames(avg$RNA) %in% gene),]
df2 <- df2[match(gene,rownames(df2)),] 

p1 <- pheatmap(df2,cluster_rows = F, cluster_cols = T, annotation_names_col = F,
               annotation_col = df, scale = "row", show_colnames = T)

# Comparing the the expression of cytotoxic genes
Uninfect <- subset(blood_T, condition == "Uninfected")
Infect <- subset(blood_T, condition == "Acute-Infected")

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
png(filename ="/work/foxlab/xiaoke/seurat/plots/uninfect.png", width = 2000, height = 1500, units = 'px',res = 300)
p2
dev.off()

# use auto annotation for blood T cells
Idents(blood_T) <- blood_T$seurat_clusters
seurat.markers <- FindAllMarkers(blood_T, method = "MAST")
library(scMayoMap)
mayomap.obj <- scMayoMap(data = seurat.markers, database = scMayoMapDatabase, tissue = "blood")
p1 <- scMayoMap.plot(mayomap.obj)
