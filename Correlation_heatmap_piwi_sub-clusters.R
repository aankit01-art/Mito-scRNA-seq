
#Ankit Arora
#Planaria project

######## Script to perform the Correlation coefficient Heatmap for interested 11 genes obtained from mito-scRNA-seq data.


library(Seurat)
library(ggplot2)
library(cowplot)
library(dplyr)
library(readxl)
library(stringi)
library(purrr)
library(reshape2)


setwd("")
MTGplana.integrated <- readRDS(file = "scRNAseq_Planaria_iPSCsDopa.integrated.RDS")
pathto.outPlots <- ""
Idents(MTGplana.integrated) <- "seurat_clusters"
MTGplana.integratedpiwi_clusters <- subset(MTGplana.integrated, ident=c(1,7))

Idents(MTGplana.integratedpiwi_new_clusters) <- "seurat_clusters"
DefaultAssay(MTGplana.integratedpiwi_new_clusters) <- "integrated"
MTGplana.integratedpiwi_new_clusters <- RunUMAP(MTGplana.integratedpiwi_new_clusters, dims = 1:30)
MTGplana.integratedpiwi_new_clusters <- FindNeighbors(MTGplana.integratedpiwi_new_clusters, reduction = "pca", dims = 1:30, nn.eps = 0.5)
MTGplana.integratedpiwi_new_clusters <- FindClusters(MTGplana.integratedpiwi_new_clusters, resolution = 0.2, n.start = 10)


DefaultAssay(iPSCsDopa.integratedpiwi_new_clusters) <- "RNA"
neoblastmarkers=read.delim("Selected_genes_DE_SMESG-ID.txt", stringsAsFactors=F)

av.exp <- AverageExpression(MTGplana.integratedpiwi_new_clusters,features = unique(neoblastmarkers[,2]))
av.exp.RNA <- av.exp$RNA

df2 <- data.matrix(av.exp.RNA)
final_df <- as.data.frame(t(df2))
colnames(final_df)[1] <- "PYR1"
colnames(final_df)[2] <- "TRAP1"
colnames(final_df)[3] <- "HSP7C"
colnames(final_df)[4] <- "PESC"
colnames(final_df)[5] <- "RL37a"
colnames(final_df)[6] <- "H1B"
colnames(final_df)[7] <- "DDX43"
colnames(final_df)[8] <- "Unidentified-1"
colnames(final_df)[9] <- "Unidentified-2"
colnames(final_df)[10] <- "CEBPZ"         
colnames(final_df)[11] <- "E2F5"
corr_mat <- round(cor(final_df),2)
melted <- melt(corr_mat)

ggplot(data = melted, aes(x = Var1, y = Var2, fill = value)) +
geom_tile() + 
geom_text(aes(Var2, Var1, label = round(value, 2)))



pdf(paste("Heatmap_correlation_coefficient_genes.pdf"))
ggplot(data = melted, aes(x = Var1, y = Var2, fill = value)) +
geom_tile() +
scale_fill_gradient2(name="Correlation")+
labs(x = "X2 Enriched Markers", y = "X2 Enriched Markers")+
theme(axis.text.x = element_text(angle = 45,hjust=1))+
geom_text(aes(label = round(value, 2)))
dev.off()


