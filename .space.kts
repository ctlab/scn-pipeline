job("Test snakemake pipeline using test data: 10x fastq") {
    container("continuumio/miniconda3") {
        shellScript {
        workDir = "/mnt/space/work/pipeline"
            content = """
            conda env update -n base --file /mnt/space/work/scn.yml &&
            wget https://www.dropbox.com/s/3aa6x4jk0odfg6u/945396_mm.GRCm38.cdna.all.idx?dl=1 -O /mnt/space/work/945396_mm.GRCm38.cdna.all.idx &&
            wget https://www.dropbox.com/s/zsbrkioi7nocyrl/10xv1_whitelist.txt?dl=1 -O /mnt/space/work/10xv1_whitelist.txt &&
            wget https://www.dropbox.com/s/ptuhp5qlfs8pfro/10xv2_whitelist.txt?dl=1 -O /mnt/space/work/10xv2_whitelist.txt &&
            wget https://www.dropbox.com/s/xcxbnvhwfr6b02g/10xv3_whitelist.txt?dl=1 -O /mnt/space/work/10xv3_whitelist.txt &&
            wget https://www.dropbox.com/s/nacb098h1ubbekw/transcripts_to_genes_v2.txt?dl=1 -O /mnt/space/work/transcripts_to_genes_v2.txt &&
            python generate.py --yaml_path example &&
            mv sample_description.csv /mnt/space/work/GSE01/SRS01 &&
            bash /mnt/space/work/GSE01/SRS01/task.bash
            """
        }
    }
}