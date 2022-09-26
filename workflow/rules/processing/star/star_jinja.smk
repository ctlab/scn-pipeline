from workflow.scripts.DependencyDispatcher import DependencyDispatcher
from pathlib import Path
dispatcher = DependencyDispatcher(config)

rule render_star_script:
  input: Path(config["templates"], "star.bash")
  params:
    idx=dispatcher.star_index,
    barcode=dispatcher.get_barcode_reads,
    cdna=dispatcher.get_cdna_reads,
    index=dispatcher.get_index_reads,
    bam=dispatcher.get_bam,
    technology=dispatcher.get_technology,
    version=dispatcher.get_version,
    processing_mode=dispatcher.get_processing_mode,
    whitelist=dispatcher.get_whitelist,
    threads=4,
    out_dir="data/samples/{dataset}/{sample}/star/solo/"
  output: "data/samples/{dataset}/{sample}/star/star.bash"
  threads: 1
  resources:
    mem_mb=4000
  benchmark: "logs/{dataset}/{sample}/kallisto/render_star_script.benchmark"
  log: "logs/{dataset}/{sample}/kallisto/render_star_script.log"
  conda: "jinja"
  script: "../../../scripts/RenderJinja.py"