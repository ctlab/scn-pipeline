library(ggplot2)
library(Seurat)

data <- readRDS(snakemake@input$seurat)

DefaultAssay(data) <- "RNA"
data <- NormalizeData(data)

averagePCT <- function(seurat,
                       ident, # columns from @meta.data
                   assay="RNA",
                       slot="counts",
                       threshold=0) {
  seurat <- SetIdent(seurat, value=ident)
  allLevels <- levels(Idents(seurat))
  data <- GetAssayData(seurat, assay=assay, slot=slot)

  results <- matrix(nrow=nrow(data), ncol=length(allLevels),
                    dimnames = list(rownames(data), allLevels))

  for (i in 1:length(allLevels)) {
    ident <- allLevels[i]
    cell.ids <- which(Idents(seurat) == ident)
    results[, i] <- round(rowSums(data[, cell.ids, drop=F] > threshold) / length(cell.ids), digits=3)
  }
  return(results)
}

averageExpression <- function(seurat,
                              ident,
                              assay="RNA",
                              slot="data") {
  cluster.averages <- AverageExpression(object = seurat, group.by=ident, assays = assay, slot = slot)
  return(cluster.averages[[1]])
}

allMarkers <- function(seurat,
                       ident,
                       assay="RNA",
                       slot="data") {

  Idents(seurat) <- seurat[[ident]]
  whole.markers <- FindAllMarkers(object = seurat,
                                  assay=assay,
                                  slot=slot,
                                  only.pos = TRUE,
                                  min.pct = 0.10,
                                  test.use = 'wilcox',
                                  max.cells.per.ident = 3e3,
                                  random.seed = 42)
  return(whole.markers)
}



for (index in 1:length(snakemake@params$resolutions)) {

  resolution <- snakemake@params$resolutions[index]
  identName <- paste0('SCT_snn_res.', resolution)

  clusterAverages <- averageExpression(data, identName)
  write.table(clusterAverages, file=snakemake@output$clusters_avg[index])

  clusterPCTs <- averagePCT(data, identName)
  write.table(clusterPCTs, file=snakemake@output$clusters_pct[index])

  deResults <- allMarkers(data, identName)
  write.table(deResults,
              snakemake@output$markers[index],
              sep="\t",
              quote=F,
              row.names=F)
  
}