
{% if AnalysisType == "single" %}
cells.after <- length(colnames(x = whole))
print(paste0("cells.before:",cells.before))
print(paste0("cells.after:",cells.after))
print(paste0("cell.diff:", cells.before-cells.after))

{% elif AnalysisType == "many" %}
cells.after <- sapply(whole, function(x) length(colnames(x = x)))
cells.diff <- cells.before-cells.after
rbind(cells.before, cells.after, cells.diff)
{% endif %}