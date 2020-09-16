gc()

{% if AnalysisType == "single" %}
whole <- RunPCA(object = whole, features = VariableFeatures(object = whole), npcs={{ PcaComponentsTotal }})
{% elif AnalysisType == "many" %}
whole.integrated <- RunPCA(whole.integrated, verbose = FALSE)
{% endif %}