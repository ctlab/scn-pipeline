library(ggplot2)
library(flexmix)
library(glmGamPoi)
library(Seurat)
library(SeuratWrappers)
library(jsonlite)

theme_set(theme_bw(base_size = 8))
set.seed(1)

counts <- readRDS(snakemake@input$filtered_counts)
data <- CreateSeuratObject(counts, project=snakemake@params$sample, min.cells = 3, min.features = 10)

data <- PercentageFeatureSet(data, pattern = "^Mt\\.|^MT\\.|^mt\\.|^Mt-|^MT-|^mt-", col.name = "percent.mt")
data[['percent.mt_log10']] <- log10(data[['percent.mt']] + 1)
data[['nCount_RNA_log10']] <- log10(data[['nCount_RNA']] + 1)
data[['nFeature_RNA_log10']] <- log10(data[['nFeature_RNA']] + 1)


## QC Plots
vlnFeaturePlots <- VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size=0.1)
ggsave(snakemake@output$vln_feature_plots, plot=vlnFeaturePlots, width=6, height=4)
  
umiMtPlot <- FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "percent.mt") + theme(aspect.ratio = 1)
ggsave(snakemake@output$umi_mt_plot, plot=umiMtPlot, width=6, height=4)

umiFeaturesLog10Plot <- FeatureScatter(data, feature1 = "nCount_RNA_log10", feature2 = "nFeature_RNA_log10") +
  theme(aspect.ratio = 1)
ggsave(snakemake@output$umi_features_log10_plot, plot=umiFeaturesLog10Plot, width=6, height=4)

umiFeaturesPlot <- FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + theme(aspect.ratio = 1)
ggsave(snakemake@output$umi_features_plot, plot=umiFeaturesPlot, width=6, height=4)

## Number of cells before

cells.before <- dim(data)[2]

data <- RunMiQC(data, 
                percent.mt = "percent.mt", 
                nFeature_RNA = "nFeature_RNA", 
                posterior.cutoff = 0.75, 
                model.slot = "flexmix_model")

tryCatch({
  miqcPlotProb <- PlotMiQC(data, color.by = "miQC.probability") +
    scale_color_gradient(low = "grey", high = "purple")
  ggsave(snakemake@output$miqc_plot_prob, plot=miqcPlotProb, width=6, height=4)
}, error = function(e) {
  message(e)
  pdf(snakemake@output$miqc_plot_prob, width=6, height=4)
  dev.off()
})

miqcPlotKeep <- FeatureScatter(data, feature1 = "nFeature_RNA", feature2 = "percent.mt", group.by = "miQC.keep") +
  theme(aspect.ratio = 1)
ggsave(snakemake@output$miqc_plot_keep, plot=miqcPlotKeep, width=6, height=4)

## FILTER MT CONTENT

data <- subset(data, miQC.keep == "keep")

## Number of cells after

seurat_stats <- list()

cells.after <- dim(data)[2]
print(paste0("cells.before:",cells.before))
print(paste0("cells.after:",cells.after))
print(paste0("cell.diff:", cells.before-cells.after))

seurat_stats$cells_before_mt_filtering <- cells.before
seurat_stats$cells_after_mt_filtering <- cells.after
seurat_stats$cells_filtered_mt <- cells.before-cells.after

## plots after filtering

## QC Plots

vlnFeaturePLotsAfterLog <- VlnPlot(data, features = c("nFeature_RNA_log10","nCount_RNA_log10","percent.mt"), pt.size=0.1)
ggsave(snakemake@output$vln_feature_plots_after_log, plot=vlnFeaturePLotsAfterLog, width=6, height=4)

vlnFeaturePlotsAfter <- VlnPlot(data, features = c("nFeature_RNA","nCount_RNA","percent.mt"), pt.size=0.1)
ggsave(snakemake@output$vln_feature_plots_after, plot=vlnFeaturePlotsAfter, width=6, height=4)

umiMtPlotAfter <- FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "percent.mt") + theme(aspect.ratio = 1)
ggsave(snakemake@output$umi_mt_plot_after, plot=umiMtPlotAfter, width=6, height=4)

umiFeaturesLog10PlotAfter <- FeatureScatter(data, feature1 = "nCount_RNA_log10", feature2 = "nFeature_RNA_log10") +
  theme(aspect.ratio = 1)
ggsave(snakemake@output$umi_features_log10_plot_after, plot=umiFeaturesLog10PlotAfter, width=6, height=4)

umiFeaturesPlotAfter <- FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
  theme(aspect.ratio = 1)
ggsave(snakemake@output$umi_features_plot_after, plot=umiFeaturesPlotAfter, width=6, height=4)


## NORMALIZATION

data <- SCTransform(
  data,
  method = "glmGamPoi",
  ncells=min(100000, ncol(data)),
  vars.to.regress = c("percent.mt"),
  verbose = T,
  conserve.memory = T
)

## PCA
data <- RunPCA(object = data,
               features = VariableFeatures(object = data), 
               npcs=50)

elbowPlot <- ElbowPlot(data, ndims = 50)
ggsave(snakemake@output$elbow_plot, plot=elbowPlot, width=6, height=4)

## TSNE

data <- RunTSNE(data, 
                dims = 1:20, 
                tsne.method = "FIt-SNE",
                nthreads = 4, 
                max_iter = 2000)


## UMAP

data <- RunUMAP(data, dims = 1:20)


## CLUSTERING

resolutions <- snakemake@params$resolutions
defaultResolution <- snakemake@params$default_resolution

data <- FindNeighbors(object = data, dims = 1:20)
data <- FindClusters(object = data, resolution = resolutions)


## Default resolution to be 0.6

allIdents <- paste0('SCT_snn_res.', resolutions)
defaultIdent <- paste0('SCT_snn_res.', defaultResolution)

seurat_stats$clustering <- list()
for (ident in allIdents) {
  seurat_stats$clustering[[ident]] <- length(levels(data@meta.data[[ident]]))
}

Idents(data) <- data[[defaultIdent]]

tsnePlot <- DimPlot(data, reduction = "tsne") + theme(aspect.ratio = 1)
ggsave(snakemake@output$tsne_plot, plot=tsnePlot, width=6, height=4)

umapPlot <- DimPlot(data, reduction = "umap") + theme(aspect.ratio = 1)
ggsave(snakemake@output$umap_plot, plot=umapPlot, width=6, height=4)


## SAVING: DATASET

write(toJSON(seurat_stats, auto_unbox = T, pretty = T),
      snakemake@output$seurat_stats)

saveRDS(data, file = snakemake@output$seurat)