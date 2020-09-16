

draw_plots <- function(path, data) {
  VlnPlot(
    data,
    features = "nFeature_RNA",
    pt.size = 0.1
    ) + theme(aspect.ratio = 1)
  ggsave(paste0(path, 'vln_plot_features.pdf'))
  VlnPlot(
    data,
    features = "nCount_RNA",
    pt.size = 0.1
  ) + theme(aspect.ratio = 1)
  ggsave(paste0(path, 'vln_plot_umi.pdf'))
  VlnPlot(
    data,
    features = "percent.mito",
    pt.size = 0.1
  ) + theme(aspect.ratio = 1)
  ggsave(paste0(path, 'vln_plot_mt.pdf'))
  
  
  
  FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "percent.mito") + theme(aspect.ratio = 1)
  ggsave(paste0(path, 'umi_mt_plot.pdf'))
  
  FeatureScatter(data, feature1 = "nCount_RNA_log10", feature2 = "nFeature_RNA_log10") + theme(aspect.ratio = 1)
  ggsave(paste0(path, 'umi_features_log10_plot.pdf'))
  
  FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + theme(aspect.ratio = 1)
  ggsave(paste0(path, 'umi_features_plot.pdf'))
  
  ElbowPlot(data, ndims = 50) + theme(aspect.ratio = 1)
  ggsave(paste0(path, 'elbow_plot.pdf'))
}