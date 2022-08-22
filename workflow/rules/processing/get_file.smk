from workflow.scripts.DependencyDispatcher import DependencyDispatcher

def get_file_param(wildcards):
    files = DependencyDispatcher(config).get_run(wildcards).files
    file = [file for file in files if file.filename == wildcards.filename][0]
    return file

def get_url(wildcards):
    return get_file_param(wildcards).url

def get_md5sum(wildcards):
    return get_file_param(wildcards).md5

rule ncbi_prefetch:
    output:
        sra=temp(config["ncbi_dir"] + "/sra/{run}.sra")
    params:
        run="{run}"
    conda: "../../envs/entrez_direct_utils.yaml"
    shell: """
    prefetch {params.run}
    """

rule get_file:
    """
    Download the file and check the MD5sum
    """
    ## output: temp(config["out_dir"] + "/data/samples/{dataset}/{sample}/files/{run}/{filename}")
    output: config["out_dir"] + "/data/samples/{dataset}/{sample}/files/{run}/{filename}"
    conda: "../../envs/define_technology.yaml"
    log: config['logs_dir'] + "/{dataset}/{sample}/{run}/{filename}_get_file.log"
    benchmark: config['logs_dir'] + "/{dataset}/{sample}/{run}/{filename}_get_file.benchmark"
    params: url=get_url, md5=get_md5sum, work_dir=config["out_dir"] + "/data/samples/{dataset}/{sample}/files/{run}", filename="{filename}"
    resources:
        mem_mb=4000
    shell: """
    mkdir -p {params.work_dir}
    cd {params.work_dir}
    echo "{params.md5}  {params.filename}" > {params.filename}.md5sum
    wget -q -O {output:q} {params.url:q} 2> {log}
    md5sum -c {params.filename}.md5sum
    """

use rule get_file as get_fastq_ftp with:
    params:
        url = get_url,
        md5 = get_md5sum,
        filename="{filename}",
        work_dir=config["out_dir"] + "/data/samples/{dataset}/{sample}/files/{run}/fastq"
    output:
        temp(config["out_dir"] + "/data/samples/{dataset}/{sample}/files/{run}/fastq/{filename}")

use rule get_file as get_bam_ftp with:
    params:
        url = get_url,
        md5 = get_md5sum,
        filename = "{filename}",
        work_dir=config["out_dir"] + "/data/samples/{dataset}/{sample}/files/{run}/bam"
    output:
        temp(config["out_dir"] + "/data/samples/{dataset}/{sample}/files/{run}/bam/{filename}")


rule get_fastq_dump_files:
    input:
        sra=rules.ncbi_prefetch.output.sra
    output:
        outs=[
            temp(config["out_dir"] + "/data/samples/{dataset}/{sample}/files/{run}/fastq_dump/{run}_1.fastq.gz"),
            temp(config["out_dir"] + "/data/samples/{dataset}/{sample}/files/{run}/fastq_dump/{run}_2.fastq.gz"),
            temp(config["out_dir"] + "/data/samples/{dataset}/{sample}/files/{run}/fastq_dump/{run}_3.fastq.gz"),
            temp(config["out_dir"] + "/data/samples/{dataset}/{sample}/files/{run}/fastq_dump/{run}_4.fastq.gz")
        ]
    conda: "../../envs/entrez_direct_utils.yaml"
    params:
        sra='{run}',
        work_dir=config["out_dir"] + "/data/samples/{dataset}/{sample}/files/{run}/fastq_dump"
    threads: 4
    resources:
        mem_mb=8000
    shell:
        """
        mkdir -p {params.work_dir}
        cd {params.work_dir}
        parallel-fastq-dump -s {input.sra} --split-files --threads {threads} --outdir {params.work_dir:q} --tmpdir . --gzip
        touch {output.outs}
        """