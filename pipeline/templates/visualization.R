
{% if AnalysisType == "single" and not test_mode %}
draw_plots(path, whole)
{% elif AnalysisType == "many" and not test_mode %}
draw_plots(path, whole.integrated)
{% endif %}