
{% if AnalysisType == "single" %}
cells.before <- dim(GetAssayData(object = whole, slot = "counts"))[2]
{% elif AnalysisType == "many" %}
cells.before <- sapply(whole, function(x) dim(GetAssayData(object = x, slot = "counts"))[2])
{% endif %}
