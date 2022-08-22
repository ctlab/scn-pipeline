rule split_sample_descriptions:
  input: config["out_dir"] +  "/meta/{dataset}/sample_description.csv"
  output: config["out_dir"] + "/data/samples/{dataset}/{sample}/sample_description.csv"
  params: sample="{sample}", dataset="{dataset}"
  wildcard_constraints:
    dataset="GSE\d+",
    sample="GSM\d+"
  log: config['logs_dir'] + "/{dataset}/{sample}/split_sample_descriptions.log"
  benchmark: config['logs_dir'] + "/{dataset}/{sample}/split_sample_descriptions.benchmark"
  conda: "../../envs/isolate_sample.yaml"
  script: "../../scripts/isolate_sample.R"