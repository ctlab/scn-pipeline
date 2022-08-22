include: "get_cellranger_index.smk"

rule get_star_index:
    input:
        fasta=rules.get_cellranger_fasta.output.fasta,
        gtf=rules.get_cellranger_fasta.output.gtf
    output:
        star_chr_length = config["resources"] + "/star/{species}/chrLength.txt",
        star_chr_name_length = config["resources"] + "/star/{species}/chrNameLength.txt",
        star_chr_name = config["resources"] + "/star/{species}/chrName.txt",
        star_chr_start = config["resources"] + "/star/{species}/chrStart.txt",
        star_exon_ge_tr_info = config["resources"] + "/star/{species}/exonGeTrInfo.tab",
        star_exon_info = config["resources"] + "/star/{species}/exonInfo.tab",
        star_gene_info = config["resources"] + "/star/{species}/geneInfo.tab",
        star_genome = config["resources"] + "/star/{species}/Genome",
        star_genome_parameters = config["resources"] + "/star/{species}/genomeParameters.txt",
        star_sa = config["resources"] + "/star/{species}/SA",
        star_sa_index = config["resources"] + "/star/{species}/SAindex",
        star_sjdb_info = config["resources"] + "/star/{species}/sjdbInfo.txt",
        star_sjdb_list_from_gtf_out = config["resources"] + "/star/{species}/sjdbList.fromGTF.out.tab",
        star_sjdb_list_out = config["resources"] + "/star/{species}/sjdbList.out.tab",
        star_transcript_info = config["resources"] + "/star/{species}/transcriptInfo.tab",
    params:
        out_dir=config["resources"] + "/star/{species}/"
    threads: 4
    resources:
        mem_mb=64000
    log: config["logs_dir"] + "/resources/get_alevin_index_{species}.log"
    benchmark: config["logs_dir"] + "/resources/get_alevin_index_{species}.benchmark"
    conda: "../../envs/star.yaml"
    shell: """
    
    STAR  --runMode genomeGenerate --runThreadN {threads} --genomeDir {params.out_dir} --genomeFastaFiles {input.fasta} --sjdbGTFfile {input.gtf} 2> {log}

    """