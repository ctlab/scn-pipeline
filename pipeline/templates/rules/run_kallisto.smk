{% if not cell_ranger %}

rule kallisto:
    input: rules.prepare_kallisto.output
    conda: "{{ PathToCondaYml }}"
    output: out=temp("bus_out/output.bus")
    log: "logs/kallisto.log"
    benchmark: "benchmarks/kallisto.txt"
    shell: "bash kallisto.sh 2> {log}"

{% endif %}

{% if cell_ranger %}

rule kallisto:
{% if not fq_dump %}
    input: rules.define_version.output.kallisto_script
    output: temp("R1.gz"), temp("R2.gz"), out=temp("bus_out/correct_output.bus"), uncor_out=temp("bus_out/output.bus")
{% elif fq_dump %}
    input: expand("{accession}_1.fastq.gz", accession=runs), expand("{accession}_2.fastq.gz", accession=runs)
    output: temp("R1.gz"), temp("R2.gz"), out=temp("bus_out/correct_output.bus"), uncor_out=temp("bus_out/output.bus")
{% endif %}
    log: "logs/kallisto.log"
    benchmark: "benchmarks/kallisto.txt"
    conda: "{{ PathToCondaYml }}"
    singularity: "{{ PathToSinImage }}"
    shell: "bash kallisto.sh 2> {log}"
{% endif %}