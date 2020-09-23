rule kallisto:
    input: fq_dir=rules.get_fastq.output.fq_dir
    output: tsv=protected("out/kallisto/{srr}/abundance.tsv")
    log: "logs/kallisto/{srr}.log"
    conda: "{{ PathToCondaYml }}"
    params: out_dir = lambda wildcards: f'out/kallisto/{wildcards.srr}', index=config['index']
    shell:
         """
         bash scripts/kallisto.sh {wildcards.srr} {input.fq_dir} {params.index} {params.out_dir} 2> {log}
         """