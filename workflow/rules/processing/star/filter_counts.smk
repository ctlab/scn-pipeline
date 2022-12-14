include: "star_run.smk"

rule filter_counts_star:
    input:
        mtx=rules.run_star.output.solo_filtered_matrix,
        genes=rules.run_star.output.solo_filtered_features,
        barcodes=rules.run_star.output.solo_filtered_barcodes
    output:
        filtered_counts="data/samples/{dataset}/{sample}/star/filtered_counts.rds",
    threads: 4
    priority: 3
    resources:
        mem_mb=16000
    log: "logs/{dataset}/{sample}/star/filter_counts.log"
    benchmark: "logs/{dataset}/{sample}/star/filter_counts.benchmark"
    conda: "../../../envs/filtering_counts.yaml"
    script: "filter_counts.R"
