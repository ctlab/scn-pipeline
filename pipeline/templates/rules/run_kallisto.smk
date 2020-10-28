{% if cell_ranger and bam %}

def is_merged():
    if len(bam_files.keys()) > 1:
        return temp('merged.bam')
    else:
        return []

def get_tmp_bams():
    return [temp(f'file{idx}.bam') for idx in range(0, len(bam_files.keys()))]

{% endif %}


{% if not cell_ranger %}

rule kallisto:
    input: rules.prepare_kallisto.output
    conda: "{{ PathToCondaYml }}"
    output: out=temp("bus_out/output.bus")
    log: "logs/kallisto.log"
    benchmark: "benchmarks/kallisto.txt"
    threads: {{ KallistoThreads }}
    shell: "bash kallisto.sh 2> {log}"

{% endif %}

{% if cell_ranger %}

rule kallisto:
{% if not fq_dump and not bam %}
    input: rules.define_version.output.kallisto_script
    output: temp("R1.gz"), temp("R2.gz"), out=temp("bus_out/correct_output.bus"), uncor_out=temp("bus_out/output.bus")
{% elif bam %}
    input: rules.define_version.output.kallisto_script
    output: is_merged(), get_tmp_bams(),
            out=temp("bus_out/correct_output.bus"), uncor_out=temp("bus_out/output.bus")
{% elif fq_dump %}
    input: expand("{accession}_1.fastq.gz", accession=runs), expand("{accession}_2.fastq.gz", accession=runs)
    output: temp("R1.gz"), temp("R2.gz"), out=temp("bus_out/correct_output.bus"), uncor_out=temp("bus_out/output.bus")
{% endif %}
    log: "logs/kallisto.log"
    benchmark: "benchmarks/kallisto.txt"
    conda: "{{ PathToCondaYml }}"
    singularity: "{{ PathToSinImage }}"
    threads: {{ KallistoThreads }}
    shell: "bash kallisto.sh 2> {log}"
{% endif %}