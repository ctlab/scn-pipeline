
{% if AnalysisType == "single" %}
cluster.averages <- AverageExpression(object = whole, assays = 'SCT', slot = 'data')
sapply(names(cluster.averages), 
       function(x) write.table(cluster.averages[[x]], file=paste0(x, "_clusters.tsv")))
{% elif AnalysisType == "many" %}
cluster.averages <- AverageExpression(object = whole.integrated, assays = 'SCT', slot = 'data')
sapply(names(cluster.averages), 
       function(x) write.table(cluster.averages[[x]], file=paste0(x, "_clusters.tsv")))

{% endif %}