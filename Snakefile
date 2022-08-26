from snakemake.utils import min_version
from workflow.scripts.DependencyDispatcher import DependencyDispatcher
import os


min_version("7.0")
configfile: "configs/config.yaml"

config["datasets"] = config.get("datasets", os.path.join(os.getcwd(), "configs", "datasets.json"))
config["samples"] = config.get("samples", os.path.join(config["out_dir"], "meta", "sample_description.json"))
config["templates"] = config.get("templates", os.path.join(os.getcwd(), "workflow", "templates"))
config["resources"] = config.get("resources", os.path.join(config["out_dir"], "resources"))
config["logs_dir"] = config.get("logs_dir", os.path.join(config["out_dir"], "logs"))

wildcard_constraints:
    run="SRR\d+",
    filename="\w+(\\.\w+)*"


include: "workflow/rules/resources/get_all_resources.smk"
include: "workflow/rules/preparation/get_all_meta.smk"
include: "workflow/rules/processing/get_file.smk"
include: "workflow/rules/analysis/analysis.smk"
include: "workflow/rules/utils/utils.smk"


dispatcher = DependencyDispatcher(config)
datasets = dispatcher.get_datasets()

datasets_full = []
samples_full = []

for dataset in datasets.keys():
    for sample in datasets[dataset].samples.keys():
        datasets_full.append(dataset)
        samples_full.append(sample)

rule process_all:
    input:
        seurat=expand(rules.seurat_analysis.output.seurat, zip, dataset=datasets_full, sample=samples_full),
        markers=expand(rules.markers_default.output.markers, zip, dataset=datasets_full, sample=samples_full),
        pct=expand(rules.markers_default.output.pct, zip, dataset=datasets_full, sample=samples_full),
        merged=expand(rules.merge_samples.output.seurat, dataset=datasets_full),
        merged_markers=expand(rules.markers_default_merged.output.markers, zip, dataset=datasets_full),
        merged_pct=expand(rules.markers_default_merged.output.pct, zip, dataset=datasets_full)

localrules: process_all, print_sample_names