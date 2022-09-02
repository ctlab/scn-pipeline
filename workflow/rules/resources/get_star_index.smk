import os.path

include: "get_cellranger_index.smk"

rule get_star_index:
    input:
        fasta=rules.get_cellranger_fasta.output.fasta,
        gtf=rules.get_cellranger_fasta.output.gtf
    output:
        star_chr_length = "./resources/star/{species}/chrLength.txt",
        star_chr_name_length = "./resources/star/{species}/chrNameLength.txt",
        star_chr_name = "./resources/star/{species}/chrName.txt",
        star_chr_start = "./resources/star/{species}/chrStart.txt",
        star_exon_ge_tr_info = "./resources/star/{species}/exonGeTrInfo.tab",
        star_exon_info = "./resources/star/{species}/exonInfo.tab",
        star_gene_info = "./resources/star/{species}/geneInfo.tab",
        star_genome = "./resources/star/{species}/Genome",
        star_genome_parameters = "./resources/star/{species}/genomeParameters.txt",
        star_sa = "./resources/star/{species}/SA",
        star_sa_index = "./resources/star/{species}/SAindex",
        star_sjdb_info = "./resources/star/{species}/sjdbInfo.txt",
        star_sjdb_list_from_gtf_out = "./resources/star/{species}/sjdbList.fromGTF.out.tab",
        star_sjdb_list_out = "./resources/star/{species}/sjdbList.out.tab",
        star_transcript_info = "./resources/star/{species}/transcriptInfo.tab",
    params:
        out_dir=lambda wildcards, output: os.path.split(output.star_genome)[0]
    threads: 4
    resources:
        mem_mb=64000
    log: "./logs/resources/get_alevin_index_{species}.log"
    benchmark: "./logs/resources/get_alevin_index_{species}.benchmark"
    conda: "../../envs/star.yaml"
    shell: """
    
    STAR  --runMode genomeGenerate --runThreadN {threads} --genomeDir {params.out_dir} --genomeFastaFiles {input.fasta} --sjdbGTFfile {input.gtf} 2> {log}

    """