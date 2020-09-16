{% if AnalysisType == "single" and not test_mode %}

whole <-
  SCTransform(
    whole,
    ncells=min(100000, ncol(whole)),
    vars.to.regress = c("percent.mito"),
    verbose = T,
    conserve.memory = T
  )

{% elif AnalysisType == "single" and test_mode %}

whole <-
  SCTransform(
    whole,
    ncells=min(100000, ncol(whole)),
    verbose = T,
    conserve.memory = T
  )

{% elif AnalysisType == "many" %}
whole <- sapply(whole, function(x) SCTransform(
  x,
  ncells=min(100000, ncol(x)),
  vars.to.regress = c("percent.mito"),
  verbose = T,
  conserve.memory = T
))
{% endif %}

