
{% if AnalysisType == "single" %}
whole <- FindNeighbors(object = whole, dims = 1:{{ Clustering.GraphBased.PcaComponents }})
whole <- FindClusters(object = whole, resolution = c({{ Clustering.GraphBased.Resolution }}))
{% elif AnalysisType == "many" %}
whole.integrated <- FindNeighbors(object = whole.integrated, dims = 1:{{ Clustering.GraphBased.PcaComponents }})
whole.integrated <- FindClusters(object = whole.integrated, resolution = c({{ Clustering.GraphBased.Resolution }}))
{% endif %}