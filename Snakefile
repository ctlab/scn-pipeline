from snakemake.utils import min_version
from workflow.scripts.DependencyDispatcher import DependencyDispatcher
import os

min_version("7.0")
configfile: "configs/config.yaml"

config["datasets"] = os.path.join(os.getcwd(), "configs", "datasets.json")
config["templates"] = os.path.join(os.getcwd(), "workflow", "templates")
config["samples"] = os.path.join(config["out_dir"], "sample_description.json")
config["resources"] = os.path.join(config["out_dir"], "resources")
config["logs_dir"] = os.path.join(config["out_dir"], "logs")


include: "workflow/rules/resources/get_all_resources.smk"
include: "workflow/rules/preparation/get_all_meta.smk"
include: "workflow/rules/processing/get_file.smk"
include: "workflow/rules/scn/convert.smk"

wildcard_constraints:
    run=r"SRR\d+",
    filename=r"[\\.|\w]+"

dispatcher = DependencyDispatcher(config)
datasets = dispatcher.get_datasets()

datasets_full = []
samples_full = []

for dataset in datasets.keys():
    for sample in datasets[dataset].samples.keys():
        datasets_full.append(dataset)
        samples_full.append(sample)

DEFAULT_RESOLUTION = 2

workdir: config["out_dir"]

rule process_all:
    input:
        seurat=expand(rules.seurat_analysis.output.seurat, zip, dataset=datasets_full, sample=samples_full),
        markers=expand(rules.markers.output.markers[DEFAULT_RESOLUTION], zip, dataset=datasets_full, sample=samples_full),
        pct=expand(rules.markers.output.clusters_pct[DEFAULT_RESOLUTION], zip, dataset=datasets_full, sample=samples_full),
        average=expand(rules.markers.output.clusters_avg[DEFAULT_RESOLUTION],zip,dataset=datasets_full,sample=samples_full),
        merged=expand(rules.merge_samples.output.seurat, dataset=datasets_full),
        merged_markers=expand(rules.markers_merged.output.markers[DEFAULT_RESOLUTION], dataset=datasets_full),
        merged_pct=expand(rules.markers_merged.output.clusters_pct[DEFAULT_RESOLUTION], dataset=datasets_full),
        merged_average=expand(rules.markers_merged.output.clusters_avg[DEFAULT_RESOLUTION], dataset=datasets_full)


rule scn:
    input:
        plot_data=expand(rules.convert_to_scn.output.plot_data, zip, dataset=datasets_full, sample=samples_full),
        plot_data_merged=expand(rules.convert_to_scn_merged.output.plot_data, dataset=datasets_full)

localrules: process_all, scn