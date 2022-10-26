import os.path

FASTA_TAR_NAMES = {
    'hs': 'https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2020-A.tar.gz',
    'mm': 'https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-mm10-2020-A.tar.gz'
}

FASTA_DIR_NAMES = {
    'hs': 'refdata-gex-GRCh38-2020-A',
    'mm': 'refdata-gex-mm10-2020-A'
}

rule get_cellranger_fasta:
    output:
        reference_json = "resources/cellranger/{species}/reference.json",
        fasta="resources/cellranger/{species}/fasta/genome.fa",
        fasta_index="resources/cellranger/{species}/fasta/genome.fa.fai",
        gtf="resources/cellranger/{species}/genes/genes.gtf",
        pickle="resources/cellranger/{species}/pickle/genes.pickle",
        star_chr_length="resources/cellranger/{species}/star/chrLength.txt",
        star_chr_name_length="resources/cellranger/{species}/star/chrNameLength.txt",
        star_chr_name="resources/cellranger/{species}/star/chrName.txt",
        star_chr_start="resources/cellranger/{species}/star/chrStart.txt",
        star_exon_ge_tr_info="resources/cellranger/{species}/star/exonGeTrInfo.tab",
        star_exon_info="resources/cellranger/{species}/star/exonInfo.tab",
        star_gene_info="resources/cellranger/{species}/star/geneInfo.tab",
        star_genome="resources/cellranger/{species}/star/Genome",
        star_genome_parameters="resources/cellranger/{species}/star/genomeParameters.txt",
        star_sa="resources/cellranger/{species}/star/SA",
        star_sa_index="resources/cellranger/{species}/star/SAindex",
        star_sjdb_info="resources/cellranger/{species}/star/sjdbInfo.txt",
        star_sjdb_list_from_gtf_out="resources/cellranger/{species}/star/sjdbList.fromGTF.out.tab",
        star_sjdb_list_out="resources/cellranger/{species}/star/sjdbList.out.tab",
        star_transcript_info="resources/cellranger/{species}/star/transcriptInfo.tab",
        tmp_tar_gz=temp("resources/{species}.tar.gz")
    params:
        link=lambda wildcards: FASTA_TAR_NAMES[wildcards.species],
        dir_name=lambda wildcards: FASTA_DIR_NAMES[wildcards.species],
        resource_dir=lambda wildcards, output: os.path.split(output.tmp_tar_gz)[0],
        out_dir=lambda wildcards, output: os.path.split(output.reference_json)[0]
    conda: "../../envs/git.yaml"
    priority: 1
    log: "logs/resources/get_cellranger_fasta_{species}.log"
    benchmark: "logs/resources/get_cellranger_fasta_{species}.benchmark"
    shell:
        """
        wget -o {log} -O {output.tmp_tar_gz} {params.link}
        tar -xzf {output.tmp_tar_gz} -C {params.resource_dir}
        rm -rf {params.out_dir}
        mv {params.resource_dir}/{params.dir_name} {params.out_dir}
        """