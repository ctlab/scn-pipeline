
{% if AnalysisType == "single" %}
mt_dist <- as.data.frame(whole[['scaled_mito']][[1]])
colnames(mt_dist) <- 'scaled_mito'
ggplot(mt_dist, aes(x = scaled_mito)) +
  geom_histogram(color = "black", fill = "white", bins = 30) +
  geom_vline(xintercept=get_conf_interval(whole, 'scaled_mito')[2], colour = "red") +
  ggtitle('scaled_mito distribution before filtration')
ggsave(paste0(path, 'mt_content_hist_before_filtration.pdf'))
whole <-
  subset(
    x = whole,
    subset = scaled_mito < get_conf_interval(whole, 'scaled_mito')[2]
  )
{% elif AnalysisType == "many" %}
whole <- sapply(whole, function(x) filter_mito(x, path))
{% endif %}
