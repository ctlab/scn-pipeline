
rule get_meta_single:
    params: dataset="{dataset}"
    output: ffq_json=config['out_dir'] +"/meta/{dataset}/ffq_raw.json",
    log: config['logs_dir'] + "/{dataset}/ffq.log"
    benchmark: config['logs_dir'] + "/{dataset}/ffq.benchmark"
    resources:
        mem_mb=4000
    conda: "../../envs/ffq.yaml"
    shell: """
    ffq -o {output.ffq_json} {params.dataset} 2> {log}
    if grep "error_msg" {output.ffq_json}; then exit 1; fi
    """
