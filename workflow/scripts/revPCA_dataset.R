library(Matrix)
library(Seurat)
library(glmGamPoi)
library(ggplot2)
library(ggrepel)
library(dplyr)

data <- readRDS(snakemake@input$seurat)
DefaultAssay(data) <- "RNA"

data <- SCTransform(
  data,
  new.assay.name = "SCT.full",
  variable.features.n = min(nrow(data), 10000),
  method = "glmGamPoi",
  ncells=min(100000, ncol(data)),
  vars.to.regress = c("percent.mt"),
  verbose = T,
  conserve.memory = T
)

top20_variable_genes <- head(VariableFeatures(data), 20)
vfPlot <- VariableFeaturePlot(data) %>%
  LabelPoints(points = top20_variable_genes, repel = TRUE) +
  scale_y_log10()

ggsave(snakemake@output$variable_feature_plot, plot=vfPlot, width=6, height=4)

npcs <- min(ncol(data) - 1, 50)
data <- RunPCA(data,
               assay = "SCT.full",
               verbose = FALSE,
               npcs = npcs,
               rev.pca = TRUE,
               reduction.name = "pca.rev",
               reduction.key="PCR_")

E <- data@reductions$pca.rev@feature.loadings
saveRDS(data, file = snakemake@output$rev_pca)