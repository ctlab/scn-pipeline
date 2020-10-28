
options(future.globals.maxSize = 8000 * 1024^2)

{% if Objects and not test_mode and not panglao %}
whole <- get(load("{{ Object }}"))
whole <- subset(x = whole, features = (rowSums(as.matrix(GetAssayData(object = whole, slot = "counts"))) > round(ncol(whole) * 0.001, 0)))
whole <- add_metadata(whole)

{% elif Objects and and panglao %}
whole <- get(load("{{ Object }}"))
whole@Dimnames[[1]] <-
  make.names(gsub("_ENS.*", "", rownames(whole)), unique = T)
whole <- CreateSeuratObject(whole, min.cells = 3, project = "{{ SampleId }}")
whole <- subset(x = whole, features = (rowSums(as.matrix(GetAssayData(object = whole, slot = "counts"))) > round(ncol(whole) * 0.001, 0)))
whole <- add_metadata(whole)


{% elif Objects and test_mode and not panglao %}

whole <- get(load("{{ Object }}"))
whole <- subset(x = whole, features = (rowSums(as.matrix(GetAssayData(object = whole, slot = "counts"))) > round(ncol(whole) * 0.001, 0)))
{% endif %}

{% if PathsToCounts and not test_mode and not panglao %}

fdata <- Read10X(data.dir = "{{ PathToCounts }}/{{ SampleId }}")
whole <- CreateSeuratObject(
  counts = fdata,
  min.cells = 2,
  min.features = 200,
  project = "{{ RunName }}"
)
whole <- add_metadata(whole)
{% elif PathsToCounts and test_mode %}

fdata <- Read10X(data.dir = "{{ PathToCounts }}/{{ SampleId }}")
whole <- CreateSeuratObject(
  counts = fdata,
  min.cells = 2,
  min.features = 200,
  project = "{{ RunName }}"
)
{% endif%}