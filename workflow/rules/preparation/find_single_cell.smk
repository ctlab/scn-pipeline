import json

SPECIES = ['mm', 'hs']

checkpoint find_all_datasets:
    output:
        datasets = "meta/dump/{start_date}_{end_date}.geo_datasets.json"
    log: "logs/dump/{start_date}_{end_date}.find_all_datasets.log"
    benchmark: "logs/dump/{start_date}_{end_date}.find_all_datasets.benchmark"
    threads: 1
    params:
        species=SPECIES,
        start_date="{start_date}",
        end_date="{end_date}"
    conda: "../../envs/define_technology.yaml"
    script: "../../scripts/FindAllDatasets.py"


rule parse_geo_meta:
    output:
        file_study="meta/dump/datasets/{dataset}/study.tsv",
        file_sample="meta/dump/datasets/{dataset}/sample.tsv"
    log: "logs/dump/datasets/{dataset}/parse_geo_meta.log"
    benchmark: "logs/dump/datasets/{dataset}/.parse_geo_meta.benchmark"
    threads: 1
    params:
        dataset='{dataset}'
    conda: "../../envs/define_technology.yaml"
    script: "../../scripts/ParseGeoMeta.py"


def gather_all_meta_input_study(wildcards):
    datasets = []
    with checkpoints.find_all_datasets.get(start_date=wildcards.start_date, end_date=wildcards.end_date).output.datasets.open() as f:
        datasets = json.load(f)
    return expand(rules.parse_geo_meta.output.file_study, dataset=datasets)

def gather_all_meta_input_sample(wildcards):
    datasets = []
    with checkpoints.find_all_datasets.get(start_date=wildcards.start_date, end_date=wildcards.end_date).output.datasets.open() as f:
        datasets = json.load(f)
    return expand(rules.parse_geo_meta.output.file_sample, dataset=datasets)


rule gather_all_meta:
    input:
        study=gather_all_meta_input_study,
        sample=gather_all_meta_input_sample
    output:
        study="meta/dump/{start_date}_{end_date}.study.tsv",
        sample="meta/dump/{start_date}_{end_date}.sample.tsv",
    threads: 1
    conda: "../../envs/define_technology.yaml"
    wildcard_constraints:
        start_date='(\d\d\d\d)(\d\d)(\d\d)',
        end_date='(\d\d\d\d)(\d\d)(\d\d)'
    script: "../../scripts/ConcatTables.py"


rule find_all_single_cell:
    input:
        rules.gather_all_meta.output.study,
        rules.gather_all_meta.output.sample