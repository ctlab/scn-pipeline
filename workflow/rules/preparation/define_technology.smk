include: 'get_dataset_meta.smk'

rule define_tech:
    input:
        ffq_json=rules.get_meta_single.output.ffq_json,
        whitelist_10x_v1=rules.get_whitelists.output.whitelist_10x_v1,
        whitelist_10x_v2=rules.get_whitelists.output.whitelist_10x_v2,
        whitelist_10x_v3=rules.get_whitelists.output.whitelist_10x_v3,
    output: ffq_full="./meta/{dataset}/ffq.json"
    params:
        dataset="{dataset}",
        header_dir="./meta/{dataset}/headers/",
        ncbi_dir=config['ncbi_dir']
    log: ".logs/{dataset}/define_technology.log"
    benchmark: ".logs/{dataset}/define_technology.benchmark"
    threads: 4
    resources:
        mem_mb=8000
    conda: "../../envs/define_technology.yaml"
    script: "../../scripts/DefineTechnology.py"
