rule get_fastq:
    output: fq_dir=temp(directory('out/fastq/{srr}')), sra_file=temp("out/{srr}.sra")
    log: "logs/sra/{srr}.log"
    conda: "{{ PathToCondaYml }}"
    threads: {{ KallistoThreads }}
    shell:
        """
        prefetch {wildcards.srr} -o {output.sra_file} > {log}
        parallel-fastq-dump -s {output.sra_file} --split-files --threads {threads} -O {output.fq_dir} --tmpdir {output.fq_dir} --gzip
        """