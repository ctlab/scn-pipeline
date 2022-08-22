def link_param(wildcards):
    links = {
        'hs': 'https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2020-A.tar.gz',
        'mm': 'https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-mm10-2020-A.tar.gz'
    }
    return links[wildcards.species]

def dir_name_param(wildcards):
    dir_names = {
        'hs': 'refdata-gex-GRCh38-2020-A',
        'mm': 'refdata-gex-mm10-2020-A'
    }
    return dir_names[wildcards.species]


rule get_cellranger_fasta:
    output:
        reference_json = config["resources"] + "/cellranger/{species}/reference.json",
        fasta=config["resources"] + "/cellranger/{species}/fasta/genome.fa",
        fasta_index=config["resources"] + "/cellranger/{species}/fasta/genome.fa.fai",
        gtf=config["resources"] + "/cellranger/{species}/genes/genes.gtf",
        pickle=config["resources"] + "/cellranger/{species}/pickle/genes.pickle",
        star_chr_length=config["resources"] + "/cellranger/{species}/star/chrLength.txt",
        star_chr_name_length=config["resources"] + "/cellranger/{species}/star/chrNameLength.txt",
        star_chr_name=config["resources"] + "/cellranger/{species}/star/chrName.txt",
        star_chr_start=config["resources"] + "/cellranger/{species}/star/chrStart.txt",
        star_exon_ge_tr_info=config["resources"] + "/cellranger/{species}/star/exonGeTrInfo.tab",
        star_exon_info=config["resources"] + "/cellranger/{species}/star/exonInfo.tab",
        star_gene_info=config["resources"] + "/cellranger/{species}/star/geneInfo.tab",
        star_genome=config["resources"] + "/cellranger/{species}/star/Genome",
        star_genome_parameters=config["resources"] + "/cellranger/{species}/star/genomeParameters.txt",
        star_sa=config["resources"] + "/cellranger/{species}/star/SA",
        star_sa_index=config["resources"] + "/cellranger/{species}/star/SAindex",
        star_sjdb_info=config["resources"] + "/cellranger/{species}/star/sjdbInfo.txt",
        star_sjdb_list_from_gtf_out=config["resources"] + "/cellranger/{species}/star/sjdbList.fromGTF.out.tab",
        star_sjdb_list_out=config["resources"] + "/cellranger/{species}/star/sjdbList.out.tab",
        star_transcript_info=config["resources"] + "/cellranger/{species}/star/transcriptInfo.tab",
        tmp_tar_gz=temp(config["resources"] + "/{species}.tar.gz")
    params:
        link=link_param,
        dir_name=dir_name_param,
        resource_dir=config["resources"],
        out_dir=config["resources"] + "/cellranger/",
        species="{species}"
    conda: "../../envs/git.yaml"
    log: config["logs_dir"] + "/resources/get_cellranger_fasta_{species}.log"
    benchmark: config["logs_dir"] + "/resources/get_cellranger_fasta_{species}.benchmark"
    shell:
        """
        cd {params.resource_dir}
        wget -o {log} -O {output.tmp_tar_gz} {params.link}
        tar -xzf {output.tmp_tar_gz} 2> {log}
        rm -rf {params.out_dir}/{params.species}
        mv {params.dir_name} {params.out_dir}/{params.species}
        """