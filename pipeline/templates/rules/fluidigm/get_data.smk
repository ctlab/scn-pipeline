checkpoint get_sra:
    output: tmp_file=temp("out/{srr}.count")
    log: "logs/sra/{srr}.log"
    conda: "{{ PathToCondaYml }}"
    threads: {{ KallistoThreads }}
    params: sra_file="out/{srr}.sra"
    shell:
        """
        sleep 2m
        touch {output.tmp_file}
        prefetch {wildcards.srr} -o {params.sra_file} > {log} || true
        """

rule get_fastq:
    input: sra_file=prefetch_input
    output: fq_dir=temp(directory('out/fastq/{srr}'))
    log: "logs/sra/{srr}.log"
    conda: "{{ PathToCondaYml }}"
    threads: {{ KallistoThreads }}
    shell:
        """
        parallel-fastq-dump -s {input.sra_file} --split-files --threads {threads} -O {output.fq_dir} --tmpdir {output.fq_dir} --gzip
        """