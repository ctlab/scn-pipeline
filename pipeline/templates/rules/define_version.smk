{% if Organism == "Mus musculus" and not test_mode %}
configfile: "{{ ConfMus }}"
{% elif Organism == "Homo sapiens" and not test_mode %}
configfile: "{{ ConfHomo }}"
{% elif Organism == "Rattus norvegicus" and not test_mode %}
configfile: "{{ ConfRat }}"
{% elif test_mode %}
configfile: "{{ ConfTest }}"
{% endif %}


index = config['index']
transcripts_to_genes = config['transcripts_to_genes']

{% if not bam and not fq_dump and cell_ranger %}

rule define_version:
    """
    Define version of 10x chemistry used for particular dataset. If dataset corresponds
    to 10xv1, then downstream analysis won't be started because kallisto needs
    index file for 10xv1 platform.
    If length of the reads doesn't correspond to any 10x versions, then analysis also
    won't be pushed, and dataset must be validated manually.
    If length of the reads corresponds to particular 10x version, then defined version (and
    its whitelist) will be used in pseudoalignment step.
    """
    input: reads=expand("{accession}_1.fastq.gz", accession=barcodes.keys())
    output: kallisto_script="kallisto.sh"
    conda: "{{ PathToCondaYml }}"
    log: "logs/define_version.log"
    benchmark: "benchmarks/define_version.txt"
    params: white_10xv2="{{ white_10xv2 }}", white_10xv3="{{ white_10xv3 }}"
    threads: {{ KallistoThreads }}
    shell: "python scripts/define_version.py --fq_r1 {input.reads[0]} --threads {threads} --index {index} --transcripts_to_genes {transcripts_to_genes} --white_10xv2 {params.white_10xv2} --white_10xv3 {params.white_10xv3}"

{% elif bam and not fq_dump %}

rule define_version:
    """
    Define version of 10x chemistry used for particular dataset using whitelists and CB:Z tag in the
    first 1e3 reads from bam file.
    In the output BAM file, the original uncorrected barcode is encoded in the CR tag, and the corrected
    barcode sequence is encoded in the CB tag.
    Reads that are not able to be assigned a corrected barcode will not have a CB tag.
    Source:
    https://kb.10xgenomics.com/hc/en-us/articles/115003822406-How-does-Cell-Ranger-correct-barcode-sequencing-errors-
    """
    input: bams=expand("{accession}_tmp.bam", accession=bam_files.keys()), s_d="sample_description.csv"
    output: kallisto_script="kallisto.sh"
    conda: "{{ PathToCondaYml }}"
    log: "logs/define_version.log"
    benchmark: "benchmarks/define_version.txt"
    params: white_10xv1="{{ white_10xv1 }}", white_10xv2="{{ white_10xv2 }}",
            white_10xv3="{{ white_10xv3 }}", threads = {{ KallistoThreads }}
    shell:
        """
        python scripts/define_version.py --s_d {input.s_d} --tmp_bam {input.bams[0]} \
        --threads {params.threads} --index {index} --transcripts_to_genes {transcripts_to_genes} \
        --white_10xv1 {params.white_10xv1} --white_10xv2 {params.white_10xv2} --white_10xv3 {params.white_10xv3}
        """

{% elif cell_ranger and not bam and fq_dump %}

rule define_version:
    input: "{accession}_sra.txt"
    output: temp("{accession}_1.fastq.gz"), temp("{accession}_2.fastq.gz"), sra=temp("sra/{accession}.sra")
    log: fq_dump="logs/fq_dump/{accession}_dump.log"
    conda: "{{ PathToCondaYml }}"
    params: run_id=lambda wildcards: wildcards.accession,
            white_10xv2="{{ white_10xv2 }}", white_10xv3="{{ white_10xv3 }}", max_size='50G'
    threads: {{ KallistoThreads }}
    shell:
         """
         python scripts/define_version.py --sra_file {input} --threads {threads} --index {index} --transcripts_to_genes {transcripts_to_genes} --white_10xv2 {params.white_10xv2} --white_10xv3 {params.white_10xv3}
         prefetch {params.run_id} -o {output.sra} --max-size {params.max_size}
         parallel-fastq-dump -s {output.sra} --split-files --threads {threads} -O . --tmpdir . --gzip
         """

{% endif %}