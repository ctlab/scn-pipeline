

{% if AnalysisType == "single" %}
{% if db == "GEO" %}
file_out <- paste0(whole@project.name, '.RData')
{% elif db == "MTAB" %}
file_out <- paste0(whole@project.name, ".RData")
{% endif %}
save('whole', file = file_out)
{% elif AnalysisType == "many" %}
file_out <- paste0("{{ RunName }}", '.RData')
save(list = c('whole.integrated', 'whole.features', 'whole.anchors'), file = file_out)
{% endif %}