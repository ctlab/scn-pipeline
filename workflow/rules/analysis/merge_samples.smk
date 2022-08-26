from workflow.scripts.DependencyDispatcher import DependencyDispatcher
import os

RESOLUTIONS = [0.2, 0.4, 0.6, 0.8, 1.0]
DEFAULT_RESOLUTION = 0.6

def merge_samples_input(wildcards):
    dataset = wildcards.dataset
    samples = DependencyDispatcher(config).get_sample_names(wildcards)

    return [
        config["out_dir"] + f"/data/samples/{dataset}/{sample}/seurat.rds" for sample in samples
    ]

def merge_samples_forced_input(wildcards):
    input_files = merge_samples_input(wildcards)
    ret_val = [
        file for file in input_files if os.path.exists(file)
    ]
    return ret_val

rule merge_samples:
    input:
        objects = merge_samples_input
    output:
        seurat = config["out_dir"] + "/data/datasets/{dataset}/seurat.rds",
        elbow_plot= report(config["out_dir"] + "/data/datasets/{dataset}/plots/elbow_plot.pdf"),
        tsne_plot= report(config["out_dir"] + "/data/datasets/{dataset}/plots/tsne_plot.pdf"),
        umap_plot= report(config["out_dir"] + "/data/datasets/{dataset}/plots/umap_plot.pdf")
    params:
        resolutions = RESOLUTIONS,
        default_resolution = DEFAULT_RESOLUTION,
        dataset="{dataset}"
    threads: 4
    resources:
        mem_mb=64000
    log: config['logs_dir'] + "/{dataset}/merge_samples.log"
    benchmark: config['logs_dir'] + "/{dataset}/merge_samples.benchmark"
    conda: "../../envs/seurat_analysis.yaml"
    script: "../../scripts/merge.R"

use rule merge_samples as merge_samples_forced with:
    input:
        objects=merge_samples_forced_input
    output:
        seurat=config["out_dir"] + "/data/datasets/{dataset}/forced/seurat.rds",
        elbow_plot=report(config["out_dir"] + "/data/datasets/{dataset}/forced/plots/elbow_plot.pdf"),
        tsne_plot=report(config["out_dir"] + "/data/datasets/{dataset}/forced/plots/tsne_plot.pdf"),
        umap_plot=report(config["out_dir"] + "/data/datasets/{dataset}/forced/plots/umap_plot.pdf")

