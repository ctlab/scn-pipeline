
{% if AnalysisType == "single" and not test_mode %}
whole <-
  RunTSNE(whole, dims = 1:{{ PcaComponentsTsne }}, tsne.method = "FIt-SNE",
          fast_tsne_path = "{{ PathToFastTsne }}", nthreads = {{ Threads }}, max_iter = 2000)
{% elif AnalysisType == "many" and not test_mode %}
whole.integrated <- RunTSNE(whole.integrated, dims = 1:{{ PcaComponentsTsne }}, tsne.method = "FIt-SNE",
                            fast_tsne_path = "{{ PathToFastTsne }}", nthreads = {{ Threads }}, max_iter = 2000)
{% endif %}