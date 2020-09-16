

filter_mito <- function(dataset, path){
  mt_dist <- as.data.frame(dataset[['scaled_mito']][[1]])
  mt_pers <- as.data.frame(dataset[['percent.mito']][[1]])
  scaled_mito_percentage <- scale(dataset[['percent.mito']][[1]])
  colnames(mt_dist) <- 'scaled_mito'
  colnames(mt_pers) <- 'percent.mito'
  filtration_coord <- get_conf_interval(dataset, 'scaled_mito')[2] * 
    attr(scaled_mito_percentage, 'scaled:scale') + 
    attr(scaled_mito_percentage, 'scaled:center')
  ggplot(mt_pers, aes(x = percent.mito)) +
    geom_histogram(color = "black", fill = "white", bins = 30) +
    geom_vline(xintercept=filtration_coord, colour = "red") +
    theme(aspect.ratio = 1) +
    ggtitle('percent.mito distribution before filtration')
  ggsave(paste0(paste0(path, unique(dataset$sample)), '_mt_content_hist_before_filtration.pdf'))
  expr <- FetchData(object = dataset, vars = 'scaled_mito')
  dataset <- dataset[, which(x = expr < get_conf_interval(dataset, 'scaled_mito')[2])]
  dataset
}
