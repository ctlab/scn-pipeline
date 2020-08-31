

options(future.globals.maxSize = 10000 * 1024^2)

whole <- c()
names <- c()

{% for Object in Objects %}
data <- get(load("{{ Object }}"))
data <- subset(x = data, features = (rowSums(as.matrix(GetAssayData(object = data, slot = "counts"))) > round(ncol(data) * 0.001, 0)))
data <- add_metadata(data)
data$sample <- "{{ SampleIds[loop.index0] }}"
assign(paste0('data_', "{{ SampleIds[loop.index0] }}"), data)
whole <- c(get(paste0('data_', "{{ SampleIds[loop.index0] }}")), whole)
names <- c("{{ SampleIds[loop.index0] }}", names) # extract SRS from /path/SRA*_SRS*.sparse.RData
{% endfor %}
whole <- as.list(whole)
names(whole) <- names
rm(list=ls(pattern='data'))

if (length(whole) > 20) {
  ref <- 10
} else {
  ref <- NULL
}
