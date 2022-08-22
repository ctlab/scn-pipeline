rule get_bam_to_fastq:
    output:
        bam_to_fastq=config["resources"] + "/tools/bamtofastq"
    conda: "../../envs/git.yaml"
    log: config["logs_dir"] + "/resources/get_bam_to_fastq.log"
    benchmark: config["logs_dir"] + "/resources/get_bam_to_fastq.benchmark"
    shell: "wget -o {log} -O {output.bam_to_fastq} https://github.com/10XGenomics/bamtofastq/releases/download/v1.4.1/bamtofastq_linux"
