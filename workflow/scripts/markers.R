library(ggplot2)
library(Seurat)

data <- readRDS(snakemake@input$seurat)
outputMarkerDir <- snakemake@output$marker_dir
dir.create(outputMarkerDir, recursive = T)

for (resolution in snakemake@params$resolutions) {
  
  identName <- paste0('SCT_snn_res.', resolution)
  Idents(data) <- data[[identName]]
  cluster.averages <- AverageExpression(object = data, assays = 'SCT', slot = 'data')[[1]]
  write.table(cluster.averages, 
              file=file.path(outputMarkerDir, paste0("clusters_", resolution, "_average.tsv")))
  
  deResults <- FindAllMarkers(object = data,
                              assay='SCT',
                              slot='data',
                              
                              only.pos = FALSE,
                              
                              min.pct = 0,
                              logfc.threshold = 0,
                              
                              max.cells.per.ident = 20,
                              return.thresh=1.01,
                              
                              test.use = 'wilcox')
  
  write.table(deResults, 
              file.path(outputMarkerDir, paste0("markers_", resolution, ".tsv")),
              sep="\t", 
              quote=F,
              row.names=F)
  
}