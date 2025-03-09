library(clusterProfiler)
library(enrichplot)
library(org.Mmu.eg.db)
library(ggplot2)
library(tidyverse)

# Load the markers for brain compared to blood
MemoryT_brain <- read.csv("./EMmarkers_brain.csv")
CD4CTL_brain <- read.csv("./CTLmarkers_brain.csv")

ids_MT <- bitr(MemoryT_brain[,"X"], fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Mmu.eg.db")
ids_CTL <- bitr(CD4CTL_brain[,"X"], fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Mmu.eg.db")

MemoryT_brain_filter <- MemoryT_brain[match(ids_MT$SYMBOL, MemoryT_brain$X),]
MemoryT_brain_filter$ENTREZID <- ids_MT[match(MemoryT_brain_filter$X, ids_MT$SYMBOL),"ENTREZID"]


CD4CTL_brain_filter <- CD4CTL_brain[match(ids_CTL$SYMBOL, CD4CTL_brain$X),]
CD4CTL_brain_filter$ENTREZID <- ids_CTL[match(CD4CTL_brain_filter$X, ids_CTL$SYMBOL),"ENTREZID"]

geneList_MT <- MemoryT_brain_filter[,"avg_log2FC"]
names(geneList_MT) <- as.character(MemoryT_brain_filter[,"ENTREZID"])
geneList_CTL <- CD4CTL_brain_filter[,"avg_log2FC"]
names(geneList_CTL) <- as.character(CD4CTL_brain_filter[,"ENTREZID"])

geneList_MT <- sort(geneList_MT, decreasing = T)
geneList_CTL <- sort(geneList_CTL, decreasing = T)

# GO gsea 
gsea_MT <- gseGO(geneList = geneList_MT,
                 OrgDb = org.Mmu.eg.db,
                 ont = "BP",
                 minGSSize = 100,
                 maxGSSize = 500,
                 pvalueCutoff = 0.05,
                 verbose = F,
                 eps = 0)
ridgeplot(gsea_MT,showCategory = 10)+scale_fill_gradient(low = "red", high = "blue")+
  theme(axis.text.y = element_text(size = 17))
edo_MT <- setReadable(gsea_MT, OrgDb = "org.Mmu.eg.db")

gsea_CTL <- gseGO(geneList = geneList_CTL,
                 OrgDb = org.Mmu.eg.db,
                 ont = "BP",
                 minGSSize = 100,
                 maxGSSize = 500,
                 pvalueCutoff = 0.05,
                 verbose = F,
                 eps = 0)
edo_CTL <- setReadable(gsea_CTL, OrgDb = "org.Mmu.eg.db")

# Plot NES
gsea_MT <- read.csv("./gsea_memoryT_brain.csv")
gsea_CTL <- read.csv("./gsea_CTL_brain.csv")

gsea_MT %>%
  arrange(desc(abs(NES))) %>%
  slice_head(n =20)%>%
  mutate(Description = str_wrap(Description, width = 45))%>%
  ungroup() -> top_MT

gsea_CTL %>%
  arrange(desc(abs(NES))) %>%
  slice_head(n =20)%>%
  mutate(Description = str_wrap(Description, width = 45))%>%
  ungroup() -> top_CTL

ggplot(top_MT)+
  theme_bw()+
  geom_col()+
  aes(x = NES, y = reorder(Description,NES), fill = p.adjust)+
  scale_fill_gradient(low = "red",high = "blue")+
  theme(axis.text = element_text(size = 18, color = "black"),
        axis.title.x = element_text(size = 18, color = "black"),
        axis.title.y = element_blank(),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 15),
        plot.margin = unit(c(0.5,0.2,0.2,0.5),"cm"))

ggplot(top_CTL)+
  theme_bw()+
  geom_col()+
  aes(x = NES, y = reorder(Description,NES), fill = p.adjust)+
  scale_fill_gradient(low = "red",high = "blue")+
  theme(axis.text = element_text(size = 18, color = "black"),
        axis.title.x = element_text(size = 18, color = "black"),
        axis.title.y = element_blank(),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 15),
        plot.margin = unit(c(0.5,0.2,0.2,0.5),"cm"))



