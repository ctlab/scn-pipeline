import json

DATASET_FILE = config['datasets']
DATASETS = json.load(open(DATASET_FILE, "r"))

include: "define_technology.smk"

rule get_all_meta:
  input: expand(rules.define_tech.output.ffq_full, dataset=DATASETS)
  output: config['samples']
  log: "./logs/concatenate_sample_descriptions.log"
  conda: "../../envs/jq.yaml"
  resources:
    mem_mb=2000
  shell: """
  jq -s add {input} > {output} 2> {log}
  """