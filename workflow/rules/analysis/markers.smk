include: "seurat_analysis.smk"
include: "merge_samples.smk"

RESOLUTIONS = [0.2, 0.4, 0.6, 0.8, 1.0]
DEFAULT_RESOLUTION = 0.6

rule markers:
    input:
        seurat = rules.seurat_analysis.output.seurat
    output:
        marker_dir = directory(config["out_dir"] + "/data/samples/{dataset}/{sample}/markers"),
        markers = [
            config["out_dir"] + "/data/samples/{dataset}/{sample}/markers/markers_" + f"{resolution}.tsv"
            for resolution in RESOLUTIONS
        ],
        clusters_avg = [
            config["out_dir"] + "/data/samples/{dataset}/{sample}/markers/clusters_" + f"{resolution}_average.tsv"
            for resolution in RESOLUTIONS
        ],
        clusters_pct = [
            config["out_dir"] + "/data/samples/{dataset}/{sample}/markers/clusters_" + f"{resolution}_pct.tsv"
            for resolution in RESOLUTIONS
        ]
    params:
        resolutions = RESOLUTIONS
    threads: 4
    resources:
        mem_mb=16000
    log: config['logs_dir'] + "/{dataset}/{sample}/markers.log"
    benchmark: config['logs_dir'] + "/{dataset}/{sample}/markers.benchmark"
    conda: "../../envs/seurat_analysis.yaml"
    script: "../../scripts/markers.R"

use rule markers as markers_merged with:
    input:
        seurat = rules.merge_samples.output.seurat
    output:
        marker_dir = directory(config["out_dir"] + "/data/datasets/{dataset}/markers"),
        markers=[
            config["out_dir"] + "/data/datasets/{dataset}/markers/markers_" + f"{resolution}.tsv"
            for resolution in RESOLUTIONS
        ],
        clusters_avg=[
            config["out_dir"] + "/data/datasets/{dataset}/markers/clusters_" + f"{resolution}_average.tsv"
            for resolution in RESOLUTIONS
        ],
        clusters_pct=[
            config["out_dir"] + "/data/datasets/{dataset}/markers/clusters_" + f"{resolution}_pct.tsv"
            for resolution in RESOLUTIONS
        ]
    log: config['logs_dir'] + "/{dataset}/markers.log"
    benchmark: config['logs_dir'] + "/{dataset}/markers.benchmark"