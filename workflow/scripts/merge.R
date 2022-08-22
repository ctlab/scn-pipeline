library(Seurat)
library(ggplot2)

object.list <- sapply(snakemake@input$objects, readRDS)
features <- SelectIntegrationFeatures(object.list = object.list, nfeatures = 3000)
object.list <- PrepSCTIntegration(object.list = object.list, anchor.features = features)

gc()

anchors <- FindIntegrationAnchors(object.list = object.list, 
                                  normalization.method = "SCT",
                                  anchor.features = features)

objects.combined <- IntegrateData(anchorset = anchors, normalization.method = "SCT")
objects.combined <- RunPCA(objects.combined, verbose = FALSE)

rm(anchors, object.list)
gc()

ElbowPlot(objects.combined, ndims = 50)
ggsave(snakemake@output$elbow_plot, width=6, height=4)

objects.combined <- RunTSNE(objects.combined,
                            reduction = "pca",
                            dims = 1:30, 
                            tsne.method = "FIt-SNE",
                            nthreads = 4, 
                            max_iter = 2000)

objects.combined <- RunUMAP(objects.combined, 
                            reduction = "pca", 
                            dims = 1:30)

## CLUSTERING

resolutions <- snakemake@params$resolutions
defaultResolution <- snakemake@params$default_resolution

objects.combined <- FindNeighbors(object = objects.combined, reduction="pca", dims = 1:30)
objects.combined <- FindClusters(object = objects.combined, resolution = resolutions)


## Default resolution to be 0.6

allIdents <- paste0('integrated_snn_res.', resolutions)
defaultIdent <- paste0('integrated_snn_res.', defaultResolution)

Idents(objects.combined) <- objects.combined[[defaultIdent]]

DimPlot(objects.combined, split.by = 'orig.ident', reduction = "tsne") + theme(aspect.ratio = 1)
ggsave(snakemake@output$tsne_plot, width=8, height=6)

DimPlot(objects.combined, split.by = 'orig.ident', reduction = "umap") + theme(aspect.ratio = 1)
ggsave(snakemake@output$umap_plot, width=8, height=6)

## SAVING: DATASET

saveRDS(objects.combined, file = snakemake@output$seurat)