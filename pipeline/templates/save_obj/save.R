

{% if AnalysisType == "single" %}
save('whole', file = "{{ SampleId }}.RData")
{% elif AnalysisType == "many" %}
save(list = c('whole.integrated', 'whole.features', 'whole.anchors'), file = "{{ RunName }}.RData")
{% endif %}