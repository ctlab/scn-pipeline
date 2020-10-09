rule add_uns_to_h5ad:
    """
    Add uns information to h5ad object after Seurat processing
    """
    input: rules.run_analysis.output.h5ad
    output: "{{ RunName }}_with_uns.h5ad"
    conda: "{{ PathToCondaYml }}"
{% if AnalysisType == "single" %}
    params: kallisto="kallisto.sh", s_d="sample_description.csv", summary="{{ SummaryFile }}"
{% else %}
    params: kallisto="{{ FirstSample }}/kallisto.sh", s_d="{{ FirstSample }}/sample_description.csv", summary="{{ SummaryFile }}"
{% endif %}

    shell:
         "python scripts/add_uns.py --h5 {input} --h5_out {output} --kallisto_script {params.kallisto} \
         --s_d {params.s_d} --summary_file {params.summary}"