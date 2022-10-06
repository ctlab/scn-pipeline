from workflow.scripts.DependencyDispatcher import DependencyDispatcher

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
    conda: "../../envs/seurat_analysis.yaml"
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
        species=DependencyDispatcher(config).get_species,
        level='dataset',
        dataset='{dataset}',
        token='{dataset}'
    log: "logs/{dataset}/convert_to_scn.log"
    benchmark: "logs/{dataset}/convert_to_scn.benchmark"

