include: "seurat_analysis.smk"
include: "merge_samples.smk"

RESOLUTIONS = [0.2, 0.4, 0.6, 0.8, 1.0]
DEFAULT_RESOLUTION = 0.6

rule markers:
    input:
        seurat = rules.seurat_analysis.output.seurat
    output:
        marker_dir = directory(config["out_dir"] + "/data/samples/{dataset}/{sample}/markers")
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
        marker_dir = directory(config["out_dir"] + "/data/datasets/{dataset}/markers")
    log: config['logs_dir'] + "/{dataset}/markers.log"
    benchmark: config['logs_dir'] + "/{dataset}/markers.benchmark"


rule markers_default:
    input:
        marker_dir = rules.markers.output.marker_dir
    output:
        markers = config["out_dir"] + "/data/samples/{dataset}/{sample}/markers.tsv",
        clusters = config["out_dir"] + "/data/samples/{dataset}/{sample}/clusters.tsv",
        pct = config["out_dir"] + "/data/samples/{dataset}/{sample}/pct.tsv",
    params:
        default_resolution = DEFAULT_RESOLUTION
    log: config['logs_dir'] + "/{dataset}/{sample}/markers_default.log"
    benchmark: config['logs_dir'] + "/{dataset}/{sample}/markers_default.benchmark"
    shell: """
    cp "{input.marker_dir}/clusters_{params.default_resolution}_pct.tsv" {output.pct}
    cp "{input.marker_dir}/clusters_{params.default_resolution}_average.tsv" {output.clusters}
    cp "{input.marker_dir}/markers_{params.default_resolution}.tsv" {output.markers}
    echo `date`> {log}
    echo "cp successfull" >> {log}
    """

use rule markers_default as markers_default_merged with:
    input:
        marker_dir = rules.markers_merged.output.marker_dir
    output:
        markers = config["out_dir"] + "/data/datasets/{dataset}/markers.tsv",
        clusters = config["out_dir"] + "/data/datasets/{dataset}/clusters.tsv",
        pct = config["out_dir"] + "/data/datasets/{dataset}/pct.tsv",
    params:
        default_resolution = DEFAULT_RESOLUTION
    log: config['logs_dir'] + "/{dataset}/markers_default.log"
    benchmark: config['logs_dir'] + "/{dataset}/markers_default.benchmark"