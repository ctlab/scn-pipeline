import os.path
from pathlib import Path
from workflow.scripts.DependencyDispatcher import DependencyDispatcher


dispatcher = DependencyDispatcher(config)

rule ncbi_prefetch:
    output:
        sra=temp(Path(config["ncbi_dir"], "sra/{run}.sra"))
    params:
        run="{run}",
        max_size="100g"
    log: "logs/{run}/ncbi_prefetch.log"
    benchmark: "logs/{run}/ncbi_prefetchbenchmark"
    conda: "entrez_direct_utils"
    shell: """
    prefetch --max-size {params.max_size} {params.run} 2>&1 > {log}
    """

rule get_file:
    """
    Download the file and check the MD5sum
    """
    output: "data/samples/{dataset}/{sample}/files/{run}/{filename}"
    conda: "define_technology"
    log: "logs/{dataset}/{sample}/{run}/{filename}_get_file.log"
    benchmark: "logs/{dataset}/{sample}/{run}/{filename}_get_file.benchmark"
    params:
        url=dispatcher.get_file_url,
        md5=dispatcher.get_file_md5,
        work_dir=lambda wildcards, output: os.path.split(output[0])[0],
        filename="{filename}"
    resources:
        mem_mb=4000
    shell: """
    mkdir -p {params.work_dir}
    wget -q -O {output:q} {params.url:q} 2> {log}
    cd {params.work_dir}
    echo "{params.md5}  {params.filename}" > {params.filename}.md5sum
    md5sum -c {params.filename}.md5sum
    """

use rule get_file as get_fastq_ftp with:
    params:
        url = dispatcher.get_file_url,
        md5 = dispatcher.get_file_md5,
        filename="{filename}",
        work_dir=lambda wildcards, output: os.path.split(output[0])[0]
    output:
        temp("data/samples/{dataset}/{sample}/files/{run}/fastq/{filename}")

use rule get_file as get_bam_ftp with:
    params:
        url = dispatcher.get_file_url,
        md5 = dispatcher.get_file_md5,
        filename = "{filename}",
        work_dir=lambda wildcards, output: os.path.split(output[0])[0]
    output:
        temp("data/samples/{dataset}/{sample}/files/{run}/bam/{filename}")


rule get_fastq_dump_files:
    input:
        sra=rules.ncbi_prefetch.output.sra
    output:
        outs=[
            temp("data/samples/{dataset}/{sample}/files/{run}/fastq_dump/{run}_1.fastq.gz"),
            temp("data/samples/{dataset}/{sample}/files/{run}/fastq_dump/{run}_2.fastq.gz"),
            temp("data/samples/{dataset}/{sample}/files/{run}/fastq_dump/{run}_3.fastq.gz"),
            temp("data/samples/{dataset}/{sample}/files/{run}/fastq_dump/{run}_4.fastq.gz")
        ]
    conda: "entrez_direct_utils"
    params:
        sra='{run}',
        work_dir=lambda wildcards, output: os.path.split(output.outs[0])[0]
    threads: 4
    log: "logs/{dataset}/{sample}/{run}/get_fastq_dump_files.log"
    benchmark: "logs/{dataset}/{sample}/{run}/get_fastq_dump_files.benchmark"
    resources:
        mem_mb=8000
    shell:
        """
        mkdir -p {params.work_dir}
        parallel-fastq-dump -s {input.sra} --split-files --threads {threads} --outdir {params.work_dir:q} --tmpdir {params.work_dir} --gzip 2>&1 > {log}
        touch {output.outs}
        """
