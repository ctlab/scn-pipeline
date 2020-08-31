rule run_analysis:
    """
    Run Seurat processing using count matrix from the get_count_matrix rule.
    """
    input: rules.get_count_matrix.output.count_mat
    benchmark: "benchmarks/analysis.txt"
    log: "logs/seurat.log"
    output: h5ad=temp(f"{sample_id}.h5ad"), rda=f"{sample_id}.RData"
    conda: "{{ PathToCondaYml }}"
    singularity: "{{ PathToSinImage }}"
    shell: "Rscript scripts/analysis.R 2> {log}"