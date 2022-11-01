library(Matrix)
library(Seurat)
library(glmGamPoi)
library(ggplot2)
library(ggrepel)
library(dplyr)


data <- readRDS(snakemake@input$seurat)
hvfinfo <- HVFInfo(data)
topGenes <- min(nrow(hvfinfo), 10000)
variableGenes <- hvfinfo[order(-hvfinfo$residual_variance), ][1:topGenes, ]
VariableFeatures(data) <- rownames(variableGenes)

top20_variable_genes <- head(VariableFeatures(data), 20)
vfPlot <- VariableFeaturePlot(data) %>%
  LabelPoints(points = top20_variable_genes, repel = TRUE) +
  scale_y_log10()

ggsave(snakemake@output$variable_feature_plot, plot=vfPlot, width=6, height=4)

npcs <- min(ncol(data) - 1, 50)
data <- RunPCA(data,
               features=VariableFeatures(data),
               assay = "SCT",
               verbose = FALSE,
               npcs = npcs,
               rev.pca = TRUE,
               reduction.name = "pca.rev",
               reduction.key="PCR_")

E <- data@reductions$pca.rev@feature.loadings
saveRDS(data, file = snakemake@output$rev_pca)