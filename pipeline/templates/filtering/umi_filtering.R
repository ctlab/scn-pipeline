
{% if AnalysisType == "single" %}
peaks <- peakfinder(whole[['nCount_RNA_log10']][[1]])
umi_dist <- as.data.frame(whole[['nCount_RNA_log10']][[1]])
colnames(umi_dist) <- 'nCount_RNA_log10'
peaks <- peakfinder(whole[['nCount_RNA_log10']][[1]])
umi_dist <- as.data.frame(whole[['nCount_RNA_log10']][[1]])
colnames(umi_dist) <- 'nCount_RNA_log10'
ggplot(umi_dist, aes(x = nCount_RNA_log10)) +
  geom_histogram(color = "black", fill = "white", bins = 30) +
  geom_vline(xintercept=peaks, colour = "red") +
  theme(aspect.ratio = 1) +
  ggtitle('nCount_RNA_log10 distribution before filtration')
ggsave(paste0(path, 'numi_hist_before_filtration.pdf'))

if (length(peaks) == 1) {
  print("There is one peak in nUMI distribution, so dataset have not been filtered")
  quit(save="no")
}

if (length(peaks) > 1) {
  whole <-
    subset(x = whole, subset = nCount_RNA_log10 > peaks[1])
  umi_dist <- as.data.frame(whole[['nCount_RNA_log10']][[1]])
  colnames(umi_dist) <- 'nCount_RNA_log10'
  ggplot(umi_dist, aes(x = nCount_RNA_log10)) +
    geom_histogram(color = "black", fill = "white", bins = 30) +
    theme(aspect.ratio = 1) +
    ggtitle('nCount_RNA_log10 distribution after filtration')
  ggsave(paste0(path, 'numi_hist_after_filtration.pdf'))
}
{% elif AnalysisType == "many" %}
test <- sapply(whole, function(x) filter_umi(x, path))
{% endif %}