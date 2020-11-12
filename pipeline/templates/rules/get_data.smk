{% if cell_ranger and bam and db == 'GEO' and not fq_dump %}


rule get_bam_header:
    output: temp("{accession}_tmp.bam")
    conda: "{{ PathToCondaYml }}"
    singularity: "{{ PathToSinImage }}"
    log: "logs/get_bam_{accession}.log"
    benchmark: "benchmarks/get_bam_{accession}.txt"
    params: link=lambda wildcards, output: bam_files[f'{wildcards.accession}']
    shell: "wget -q -O - {params.link} | head -1000 > {output} || true"

{% elif cell_ranger and not bam and not fq_dump %}

rule get_fq_header:
    """
    Download the first 400000 lines of R1.fastq.gz file. Pipeline uses these lines for 10x version definition.
    """
    output: temp("{accession}_1.fastq.gz")
    conda: "{{ PathToCondaYml }}"
    log: "logs/prepare_loading_{accession}.log"
    benchmark: "benchmarks/prepare_loading_{accession}.txt"
    params: link=lambda wildcards, output: barcodes[f'{wildcards.accession}']
    {% if not test_mode %}
shell: "wget -q -O - ftp://{params.link}| zcat | head -400000 | gzip > {output} || true"
    {% elif test_mode %}
shell: "wget -q -O - {params.link}| zcat | head -400000 | gzip > {output} || true"
    {% endif %}

{% elif cell_ranger and bam and db == 'MTAB' and not fq_dump %}

rule prepare_loading:
    output: kallisto_script="download_the_beginning_of_bam.sh"
    conda: "{{ PathToCondaYml }}"
    log: "logs/prepare_loading.log"
    benchmark: "benchmarks/prepare_loading.txt"
    run:
        with open(output.kallisto_script, 'w') as out_file:
            bam_file = description[0]['bam']
            out_file.write(f"wget -q -O - {bam_file} | head -1000 > file.bam\n")

rule download_the_beginning_of_bam:
    input: rules.prepare_loading.output
    output: temp("file.bam")
    benchmark: "benchmarks/download_the_beginning_of_bam.txt"
    log: "logs/download_the_beginning_of_bam.log"
    conda: "{{ PathToCondaYml }}"
    run: shell("""
        set +o pipefail
        chmod +x  download_the_beginning_of_bam.sh
        ./download_the_beginning_of_bam.sh""")

{% elif cell_ranger and not bam and fq_dump %}

rule fastq_dump:
    output: sra_file="{accession}_sra.txt", tmp_fq=temp(directory('{accession}_tmp'))
    params: run_id=lambda wildcards: wildcards.accession
    conda: "{{ PathToCondaYml }}"
    threads: 1
    shell:
         """
         timeout 20s fastq-dump --outdir {output.tmp_fq} --split-files {params.run_id} || true
         esearch -db sra -query {params.run_id} | efetch -format metadata | grep -Po 'average="\K.*?(?=")' | head -$(ls {output.tmp_fq} | wc -l) > {params.run_id}_sra.txt
         """

{% elif not cell_ranger and not bam and not fq_dump %}

{% if Organism == "Mus musculus" %}
configfile: "{{ ConfMus }}"
{% elif Organism == "Homo sapiens" %}
configfile: "{{ ConfHomo }}"
{% elif Organism == "Rattus norvegicus" %}
configfile: "{{ ConfRat }}"
{% endif %}

index = config['index']
transcripts_to_genes = config['transcripts_to_genes']

rule prepare_kallisto:
    """
    Prepare kallisto using named pipes, defined technology and data from config file.
    """
    output: "kallisto.sh"
    conda: "{{ PathToCondaYml }}"
    log: "logs/kallisto.log"
    benchmark: "benchmarks/kallisto.txt"
    threads: {{ KallistoThreads }}
    shell: "python scripts/prepare_kallisto.py --threads {threads} --index {index} --transcripts_to_genes {transcripts_to_genes}"

{% endif %}