rule add_uns_to_h5ad:
    """
    Add uns information to h5ad object after Seurat processing
    """
    input: rules.run_analysis.output.h5ad
    output: f"{sample_id}_with_uns.h5ad"
    conda: "{{ PathToCondaYml }}"
    params: kallisto="kallisto.sh", s_d="sample_description.csv"
    shell:
         "python scripts/add_uns.py --h5 {input} --h5_out {output} --kallisto_script {params.kallisto} --s_d {params.s_d}"