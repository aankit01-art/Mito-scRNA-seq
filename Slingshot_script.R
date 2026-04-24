


################## Script for using slingshot


library(Seurat)
library(ggplot2)
library(cowplot)
library(dplyr)
library(readxl)
library(stringi)
library(purrr)
library(reshape2)
library(slingshot)
library(RColorBrewer)
library(grDevices)


MTGplana.integrated <- readRDS(file = "MTGplana.integrated.RDS")
MTGplana.integratedpiwi_p4hb_clusters <- subset(MTGplana.integrated, ident=c(0,1,2,7,10,12))
olig = (MTGplana.integratedpiwi_p4hb_clusters)
sce <- slingshot(Embeddings(olig, "umap")[,1:2], clusterLabels=Idents(olig), start.clus='0')
sce <- SlingshotDataSet(sce)
seu <- SCTransform(olig)
seu <- RunPCA(seu, npcs = 6, verbose = FALSE)
seu <- RunUMAP(seu, dims = 1:6, seed.use = 4867)


pdf(paste("plot_tracjectory_cluster.pdf"), width=10, height=10)
DimPlot(seu, reduction = "umap",
group.by = "seurat_clusters", pt.size = 0.5, label = TRUE, repel = TRUE) +
scale_color_brewer(type = "qual", palette = "Set1")
dev.off()


dimred <- seu@reductions$umap@cell.embeddings
clustering <- Idents(seu)
lineages <- getLineages(data = dimred, clusterLabels = clustering)
curves <- getCurves(lineages, approx_points = 300, thresh = 0.01, stretch = 0.8, allow.breaks = TRUE, shrink = 0.99)
pdf(paste("plot_tracjectory_clusters.pdf"), width=10, height=10)
plot(dimred, col = pals[clustering], asp = 1, pch = 16)
lines(curves, lwd = 2, col = "black")
dev.off()


