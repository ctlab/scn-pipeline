include: "../analysis/analysis.smk"

rule rev_pca_sample:
    input:
        seurat = rules.seurat_analysis.output.seurat
    output:
        rev_pca="data/samples/{dataset}/{sample}/rev_pca/rev_pca.rds",
        variable_feature_plot="data/samples/{dataset}/{sample}/rev_pca/variable_feature_plot.pdf"
    threads: 4
    resources: mem_mb=32000
    log: "logs/{dataset}/{sample}/rev_pca_sample.log"
    benchmark: "logs/{dataset}/{sample}/rev_pca_sample.benchmark"
    conda: "../../envs/seurat_analysis.yaml"
    script: "../../scripts/revPCA_sample.R"

rule rev_pca_dataset:
    input:
        seurat = rules.merge_samples.output.seurat
    output:
        rev_pca="data/datasets/{dataset}/rev_pca/rev_pca.rds",
        variable_feature_plot="data/datasets/{dataset}/rev_pca/variable_feature_plot.pdf"
    threads: 4
    resources: mem_mb=128000
    log: "logs/{dataset}/rev_pca_sample.log"
    benchmark: "logs/{dataset}/rev_pca_sample.benchmark"
    conda: "../../envs/seurat_analysis.yaml"
    script: "../../scripts/revPCA_dataset.R"