


######################## Script for Co-expression or Gene regulatory network (GRN) analysis


library(WGCNA)
library(hdWGCNA)
library(Seurat)
library(cowplot)
library(patchwork)

MTGplana.integrated <- readRDS(file = "MTGplana.integrated.RDS")


MTGplana.integrated <- SetupForWGCNA(
MTGplana.integrated,
gene_select = "fraction",
fraction = 0.05,
wgcna_name = "tutorial")


MTGplana.integrated <- MetacellsByGroups(
MTGplana.integrated,
group.by = c("seurat_clusters"),
reduction = "pca",
k = 25,
max_shared = 10,
ident.group = "seurat_clusters")


MTGplana.integrated <- NormalizeMetacells(MTGplana.integrated)


MTGplana.integrated <- SetDatExpr(
MTGplana.integrated,
group_name = '7',
group.by='seurat_clusters',
assay = 'RNA',
slot = 'data')


MTGplana.integrated <- TestSoftPowers(
MTGplana.integrated,
networkType = 'signed')


MTGplana.integrated <- ConstructNetwork(
MTGplana.integrated,soft_power=9,
setDatExpr=FALSE,
tom_name = 'INH')

pdf(paste("INH_dendo.pdf"), width=10, height=10)
PlotDendrogram(MTGplana.integrated, main='INH hdWGCNA Dendrogram')
dev.off()

MTGplana.integrated <- ScaleData(MTGplana.integrated, features=VariableFeatures(scRNAseq_Planaria_iPSCsDopa.integrated))


MTGplana.integrated <- ModuleEigengenes(
MTGplana.integrated)


hMEs <- GetMEs(scRNAseq_Planaria_iPSCsDopa.integrated)
MEs <- GetMEs(scRNAseq_Planaria_iPSCsDopa.integrated, harmonized=FALSE)

MTGplana.integrated <- ModuleConnectivity(
MTGplana.integrated)



 <- ResetModuleNames(
  ,
  new_name = "INH-M"
)


MEs <- GetMEs(MTGplana.integrated, harmonized=TRUE)
mods <- colnames(MEs); mods <- mods[mods != 'grey']
MTGplana.integrated@meta.data <- cbind(MTGplana.integrated@meta.data, MEs)



pdf(paste("INH_feature-plot.pdf"), width=10, height=10)
wrap_plots(plot_list, ncol=6)
dev.off()


pdf(paste(pathto.outPlots,"INH_dot-plot.pdf"), width=10, height=10)
DotPlot(MTGplana.integrated, features=mods, group.by = 'seurat_clusters')
dev.off()
