from workflow.scripts.DependencyDispatcher import DependencyDispatcher
import os

RESOLUTIONS = [0.2, 0.4, 0.6, 0.8, 1.0]
DEFAULT_RESOLUTION = 0.6
dispatcher = DependencyDispatcher(config)

rule merge_samples:
    input:
        objects = dispatcher.get_seurat_objects
    output:
        seurat = "data/datasets/{dataset}/seurat.rds",
        elbow_plot= report("data/datasets/{dataset}/plots/elbow_plot.pdf"),
        tsne_plot= report("data/datasets/{dataset}/plots/tsne_plot.pdf"),
        umap_plot= report("data/datasets/{dataset}/plots/umap_plot.pdf")
    params:
        resolutions = RESOLUTIONS,
        default_resolution = DEFAULT_RESOLUTION,
        dataset="{dataset}"
    threads: 4
    priority: 5
    resources:
        mem_mb=64000
    log: "logs/{dataset}/merge_samples.log"
    benchmark: "logs/{dataset}/merge_samples.benchmark"
    conda: "../../envs/seurat_analysis.yaml"
    script: "../../scripts/merge.R"

use rule merge_samples as merge_samples_forced with:
    input:
        objects=dispatcher.get_seurat_objects_forced
    output:
        seurat="data/datasets/{dataset}/forced/seurat.rds",
        elbow_plot=report("data/datasets/{dataset}/forced/plots/elbow_plot.pdf"),
        tsne_plot=report("data/datasets/{dataset}/forced/plots/tsne_plot.pdf"),
        umap_plot=report("data/datasets/{dataset}/forced/plots/umap_plot.pdf")

