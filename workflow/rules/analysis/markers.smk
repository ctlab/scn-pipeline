from pathlib import Path

include: "seurat_analysis.smk"
include: "merge_samples.smk"

RESOLUTIONS = [0.2, 0.4, 0.6, 0.8, 1.0]
DEFAULT_RESOLUTION = 0.6

rule markers:
    input:
        seurat = rules.seurat_analysis.output.seurat
    output:
        marker_dir = directory("data/samples/{dataset}/{sample}/markers"),
        markers = [
            f"data/samples/{{dataset}}/{{sample}}/markers/markers_{resolution}.tsv"
            for resolution in RESOLUTIONS
        ],
        clusters_avg = [
            f"data/samples/{{dataset}}/{{sample}}/markers/clusters_{resolution}_average.tsv"
            for resolution in RESOLUTIONS
        ],
        clusters_pct = [
            f"data/samples/{{dataset}}/{{sample}}/markers/clusters_{resolution}_pct.tsv"
            for resolution in RESOLUTIONS
        ]
    params:
        resolutions = RESOLUTIONS
    threads: 4
    resources:
        mem_mb=16000
    log: "logs/{dataset}/{sample}/markers.log"
    benchmark: "logs/{dataset}/{sample}/markers.benchmark"
    conda: "../../envs/seurat_analysis.yaml"
    script: "../../scripts/markers.R"

use rule markers as markers_merged with:
    input:
        seurat = rules.merge_samples.output.seurat
    output:
        marker_dir = directory("data/datasets/{dataset}/markers"),
        markers=[
            f"data/datasets/{{dataset}}/markers/markers_{resolution}.tsv"
            for resolution in RESOLUTIONS
        ],
        clusters_avg=[
            f"data/datasets/{{dataset}}/markers/clusters_{resolution}_average.tsv"
            for resolution in RESOLUTIONS
        ],
        clusters_pct=[
            f"data/datasets/{{dataset}}/markers/clusters_{resolution}_pct.tsv"
            for resolution in RESOLUTIONS
        ]
    log: "logs/{dataset}/markers.log"
    benchmark: "logs/{dataset}/markers.benchmark"