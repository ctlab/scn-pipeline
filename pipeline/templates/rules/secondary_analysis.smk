rule run_analysis:
    """
    Run Seurat processing using count matrix from the get_count_matrix rule.
    """
{% if AnalysisType == "single" %}
    input: rules.get_count_matrix.output.count_mat
{% else %}
    input: {{ Objects }}
{% endif %}

    output: h5ad=temp("{{ RunName }}.h5ad"), rda="{{ RunName }}.RData"
    benchmark: "benchmarks/analysis.txt"
    log: "logs/seurat.log"
    conda: "{{ PathToCondaYml }}"
    singularity: "{{ PathToSinImage }}"
    shell: "Rscript scripts/analysis.R 2> {log}"