

filter_umi <- function(dataset, path){
  peaks <- peakfinder(dataset[['nCount_RNA_log10']][[1]])
  umi_dist <- as.data.frame(dataset[['nCount_RNA_log10']][[1]])
  colnames(umi_dist) <- 'nCount_RNA_log10'
  peaks <- peakfinder(dataset[['nCount_RNA_log10']][[1]])
  umi_dist <- as.data.frame(dataset[['nCount_RNA_log10']][[1]])
  colnames(umi_dist) <- 'nCount_RNA_log10'
  ggplot(umi_dist, aes(x = nCount_RNA_log10)) +
    geom_histogram(color = "black", fill = "white", bins = 30) +
    geom_vline(xintercept=peaks, colour = "red") +
    ggtitle('nCount_RNA_log10 distribution before filtration')
  ggsave(paste0(paste0(path, unique(dataset$sample)), '_numi_hist_before_filtration.pdf'))
  
  if (length(peaks) == 1) {
    print("There is only one peak in nUMI distribution, so dataset have not been filtered")
    quit(save="no")
  }
  
  if (length(peaks) > 1) {
    expr <- FetchData(object = dataset, vars = 'nCount_RNA_log10')
    dataset <- dataset[, which(x = expr > peaks[1])]
    umi_dist <- as.data.frame(dataset[['nCount_RNA_log10']][[1]])
    colnames(umi_dist) <- 'nCount_RNA_log10'
    ggplot(umi_dist, aes(x = nCount_RNA_log10)) +
      geom_histogram(color = "black", fill = "white", bins = 30) +
      ggtitle('nCount_RNA_log10 distribution after filtration')
    ggsave(paste0(paste0(path, unique(dataset$sample)), '_numi_hist_after_filtration.pdf'))
    dataset
  }
}
