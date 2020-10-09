{% if cell_ranger and bam and db == 'GEO' and not fq_dump %}

rule download_bam_header:
    """
    Download the first ~1000 lines from bam file.
    Pipeline uses this file for definition of 10x version:
    """
    output: temp('file.bam') #todo script->bam and rm download_the_beginning_of_bam (rename prepare_lg to it); check all downstram rules are ok
    log: "logs/prepare_loading.log"
    benchmark: "benchmarks/prepare_loading.txt"
    shell: "wget -q -O - ftp://{bam_file} | head -1000 > file.bam"

{% elif cell_ranger and not bam and db == 'GEO' and not fq_dump %}

rule download_fq_header:
    """
    Download the first read for forward fastq file
    and gzipped it in to R1.gz file. Pipeline uses this file for definition of 10x version:
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
{% elif cell_ranger and not bam and db == 'MTAB' and not fq_dump %}

rule prepare_loading: #todo prepare_loading + down_beg_of_fq -> one rule
    """
    Download the first read for forward fastq file
    and gzipped it in to R1.gz file. Pipeline uses this file for definition of 10x version:
    """
    input: s_d="sample_description.csv"
    output: download_script="download_the_beginning_of_fqs.sh"
    run:
        import pandas as pd

        assert 'R1' in read_types.values() and 'R2' in read_types.values(), 'There are no R1 and / or R2 in given sample description'
        with open(output.download_script, 'w') as out_file:
            description = pd.read_csv(input.s_d).reset_index().to_dict("records")
            out_file.write(f"wget -q -O - {description[0][list(read_types.keys())[list(read_types.values()).index('R1')]]} | zcat | head -400000 | gzip > R1.gz\n")

rule download_the_beginning_of_fastqs:
    input: rules.prepare_loading.output
    output: "R1.gz"
    conda: "{{ PathToCondaYml }}"
    log: "logs/download_the_beginning_of_fqs.log"
    benchmark: "benchmarks/download_the_beginning_of_fqs.txt"
    run: shell("""
        set +o pipefail
        chmod +x  download_the_beginning_of_fqs.sh
        ./download_the_beginning_of_fqs.sh || true""")

{% elif cell_ranger and bam and db == 'MTAB' and not fq_dump %}

rule prepare_loading:
    """
    Download the first ~1000 lines from bam file.
    Pipeline uses this file for definition of 10x version:
    """
    input: "sample_description.csv"
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
    output: "{accession}_sra.txt"
    params: run_id=lambda wildcards: wildcards.accession
    conda: "{{ PathToCondaYml }}"
    threads: 1
    shell:
         """
         esearch -db sra -query {params.run_id} | efetch -format metadata | grep -Po 'average="(\d+)"' | awk '!x[$0]++' | sed 's/average=//g' | sed 's/"//g' > {params.run_id}_sra.txt
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