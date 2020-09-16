
{% if AnalysisType == "single" %}
whole <- RunUMAP(whole, dims = 1:{{ PcaComponentsUmap }})
{% elif AnalysisType == "many" %}
whole.integrated <- RunUMAP(whole.integrated, dims = 1:{{ PcaComponentsUmap }})
{% endif %}