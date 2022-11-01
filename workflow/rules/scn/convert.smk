from workflow.scripts.DependencyDispatcher import DependencyDispatcher
import os

include: "geseca.smk"
include: "../analysis/analysis.smk"
include: "../preparation/find_single_cell.smk"

rule convert_to_scn:
    input:
        study_meta=rules.parse_geo_meta.output.file_study,
        sample_meta=rules.parse_geo_meta.output.file_sample,
        seurat=rules.seurat_analysis.output.seurat,
        markers=rules.markers.output.markers[2]
    output:
        descriptor="data/samples/{dataset}/{sample}/scn/dataset.json",
        plot_data="data/samples/{dataset}/{sample}/scn/plot_data.json",
        exp_data="data/samples/{dataset}/{sample}/scn/exp_data.json",
        data_h5="data/samples/{dataset}/{sample}/scn/data.h5",
        markers="data/samples/{dataset}/{sample}/scn/markers.json"
    params:
        species=DependencyDispatcher(config).get_species,
        level='sample',
        dataset='{dataset}',
        sample='{sample}',
        token='{sample}'
    log: "logs/{dataset}/{sample}/convert_to_scn.log"
    benchmark: "logs/{dataset}/{sample}/convert_to_scn.benchmark"
    threads: 1
    resources:
        mem_mb=64000
    conda: "../../envs/scn.yaml"
    script: "../../scripts/convert.R"


use rule convert_to_scn as convert_to_scn_merged with:
    input:
        study_meta=rules.parse_geo_meta.output.file_study,
        sample_meta=rules.parse_geo_meta.output.file_sample,
        seurat=rules.merge_samples.output.seurat,
        markers=rules.markers_merged.output.markers[2]
    output:
        descriptor="data/datasets/{dataset}/scn/dataset.json",
        plot_data="data/datasets/{dataset}/scn/plot_data.json",
        exp_data="data/datasets/{dataset}/scn/exp_data.json",
        data_h5="data/datasets/{dataset}/scn/data.h5",
        markers="data/datasets/{dataset}/scn/markers.json"
    params:
        species=DependencyDispatcher(config).get_common_species,
        level='dataset',
        dataset='{dataset}',
        token='{dataset}'
    log: "logs/{dataset}/convert_to_scn.log"
    benchmark: "logs/{dataset}/convert_to_scn.benchmark"

rule upload_to_box:
    input:
        descriptor=rules.convert_to_scn.output.descriptor,
        seurat=rules.seurat_analysis.output.seurat,
        star_summary = rules.run_star.output.solo_summary,
        rev_pca = rules.rev_pca_sample.output.rev_pca,
        star_raw_counts=[
            rules.run_star.output.solo_raw_barcodes,
            rules.run_star.output.solo_raw_features,
            rules.run_star.output.solo_raw_matrix
        ],
        star_filtered_counts=[
            rules.run_star.output.solo_filtered_barcodes,
            rules.run_star.output.solo_filtered_features,
            rules.run_star.output.solo_filtered_matrix
        ]
    output:
        box_receipt = "data/samples/{dataset}/{sample}/box_receipt.txt"
    params:
        to_upload=lambda wildcards, output: os.path.split(output["box_receipt"])[0],
        box_path="dump_2022/data/{dataset}_{sample}"
    threads: 4
    conda: "../../envs/git.yaml"
    log: "logs/{dataset}/{sample}/upload_to_box.log"
    benchmark: "logs/{dataset}/{sample}/upload_to_box.benchmark"
    shell: """
    rclone copy {params.to_upload} remote:{params.box_path} -vv --transfers=16 --log-file={log}
    grep INFO {log} > {output.box_receipt}
    """

use rule upload_to_box as upload_to_box_merged with:
    input:
        descriptor=rules.convert_to_scn_merged.output.descriptor,
        seurat=rules.merge_samples.output.seurat,
        rev_pca=rules.rev_pca_dataset.output.rev_pca
    output:
        box_receipt = "data/datasets/{dataset}/dropbox_receipt.txt"
    params:
        to_upload=lambda wildcards, output: os.path.split(output["box_receipt"])[0],
        box_path="dump_2022/data/{dataset}"
    conda: "../../envs/git.yaml"
    log: "logs/{dataset}/upload_to_box.log"
    benchmark: "logs/{dataset}/upload_to_box.benchmark"

