rule get_count_matrix:
    """
    Extract count matrix which will be used in Seurat analysis using DropletUtils package.
    """
    input: rules.kallisto.output.out
    benchmark: "benchmarks/get_count_matrix.txt"
    log: "logs/get_count_matrix.log"
    output: count_mat="counts.RData"
    singularity: "{{ PathToSinImage }}"
    conda: "{{ PathToCondaYml }}"
    shell: "Rscript scripts/get_count_matrix.R 2> {log}"