include: 'get_dataset_meta.smk'

import json
import os
from workflow.scripts.Classes import parse_sample_descriptions


def get_all_runs(wildcards):
    all_runs = []
    json_file = checkpoints.get_meta_single.get(dataset=wildcards.dataset).output.ffq_json
    with open(json_file) as f:
        data = json.load(f)
        dataset = parse_sample_descriptions(data)[wildcards.dataset]
        for sample_id, sample in dataset.samples.items():
            for run in sample.get_all_runs():
                all_runs.append(run.accession)

    return [os.path.join(config['ncbi_dir'], "sra", f"{run}.sra") for run in all_runs]

rule define_tech:
    input:
        get_all_runs,
        ffq_json=rules.get_meta_single.output.ffq_json,
        whitelist_10x_v1=rules.get_whitelists.output.whitelist_10x_v1,
        whitelist_10x_v2=rules.get_whitelists.output.whitelist_10x_v2,
        whitelist_10x_v3=rules.get_whitelists.output.whitelist_10x_v3,
    output: ffq_full=config['out_dir'] + "/meta/{dataset}/ffq.json"
    params:
        dataset="{dataset}",
        header_dir=config['out_dir'] + "/meta/{dataset}/headers/",
        ncbi_dir=config['ncbi_dir']
    # log: config['logs_dir'] + "/{dataset}/define_technology.log"
    benchmark: config['logs_dir'] + "/{dataset}/define_technology.benchmark"
    threads: 4
    resources:
        mem_mb=8000
    conda: "../../envs/define_technology.yaml"
    script: "../../scripts/DefineTechnology.py"
