from snakemake.utils import min_version
from workflow.scripts.DependencyDispatcher import DependencyDispatcher
import os


min_version("7.0")
configfile: "configs/config.yaml"

config["datasets"] = config.get("templates", os.path.join(os.getcwd(), "config", "datasets.json"))
config["samples"] = config.get("samples", os.path.join(config["out_dir"], "meta", "sample_description.json"))
config["templates"] = config.get("templates", os.path.join(os.getcwd(), "workflow", "templates"))
config["resources"] = config.get("resources", os.path.join(config["out_dir"], "resources"))
config["logs_dir"] = config.get("logs_dir", os.path.join(config["out_dir"], "logs"))

wildcard_constraints:
    run="SRR\d+",
    filename="\w+(\\.\w+)*"


include: "workflow/rules/resources/get_all_resources.smk"
include: "workflow/rules/preparation/get_all_meta.smk"
include: "workflow/rules/processing/star/star.smk"
include: "workflow/rules/processing/get_file.smk"
include: "workflow/rules/analysis/analysis.smk"


dispatcher = DependencyDispatcher(config)
dataset_samples = []
datasets = dispatcher.get_datasets()

for dataset in datasets.keys():
    for sample in datasets[dataset].samples.keys():
        dataset_samples.append((dataset, sample))

rule process_all:
    input: [config["out_dir"] + f"/data/samples/{dataset}/{sample}/star/seurat.rds" for dataset, sample in dataset_samples]