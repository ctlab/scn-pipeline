

options(future.globals.maxSize = 10000 * 1024^2)

get_df <- function(path) {
  data <- get(load(path))
  data <- subset(x = data, features = (rowSums(as.matrix(GetAssayData(object = data, slot = "counts"))) > round(ncol(data) * 0.001, 0)))
  data <- add_metadata(data)
  data$sample <- gsub('.*/', '', gsub('/counts.RData', '', path))
  data
}

get_whole_obj <- function(pathes) {
  objects <- lapply(pathes, get_df)
  names(objects) <- sapply(objects, function(x) unique(x$sample))
  objects
}

whole <- get_whole_obj(c("{{ Objects|join('", "') }}"))