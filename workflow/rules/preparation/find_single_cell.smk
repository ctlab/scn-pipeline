import json

SPECIES = ['mm', 'hs']

checkpoint find_all_datasets:
    output:
        datasets = config['out_dir'] + "/meta/dump/datasets.json"
    log: config['logs_dir'] + "/dump/find_all_datasets.log"
    benchmark: config['logs_dir'] + "/dump/find_all_datasets.benchmark"
    threads: 1
    params:
        species=SPECIES,
        start_date="2022/09/01",
        end_date="2022/09/5"
    conda: "../../envs/define_technology.yaml"
    script: "../../scripts/FindAllDatasets.py"


rule parse_geo_meta:
    output:
        file_study=config['out_dir'] + "/meta/dump/datasets/{dataset}/study.tsv",
        file_sample=config['out_dir'] + "/meta/dump/datasets/{dataset}/sample.tsv"
    threads: 1
    params:
        dataset='{dataset}'
    conda: "../../envs/define_technology.yaml"
    script: "../../scripts/FindSingleCell.py"


def gather_all_meta_input_study(wildcards):
    datasets = []
    with checkpoints.find_all_datasets.get().output.datasets.open() as f:
        datasets = json.load(f)
    return expand(rules.parse_geo_meta.output.file_study, dataset=datasets)

def gather_all_meta_input_sample(wildcards):
    datasets = []
    with checkpoints.find_all_datasets.get().output.datasets.open() as f:
        datasets = json.load(f)
    return expand(rules.parse_geo_meta.output.file_sample, dataset=datasets)


rule gather_all_meta:
    input:
        study=gather_all_meta_input_study,
        sample=gather_all_meta_input_sample
    output:
        study=config['out_dir'] + "/meta/dump/study.tsv",
        sample=config['out_dir'] + "/meta/dump/sample.tsv",
    threads: 1
    conda: "../../envs/define_technology.yaml"
    shell: """
    cat {input.sample} > {output.sample}
    cat {input.study} > {output.study}
    """


rule find_all_single_cell:
    input:
        rules.gather_all_meta.output.study,
        rules.gather_all_meta.output.sample