job("Test single, merge and bulk-like modes using prepared public available data") {
    parallel {
        container("continuumio/miniconda3") {
            workDir = "/mnt/space/work/pipeline"
            shellScript {
                content = """
                conda env update -n base --file /mnt/space/work/scn.yml
                wget https://www.dropbox.com/s/3aa6x4jk0odfg6u/945396_mm.GRCm38.cdna.all.idx?dl=1 -O /mnt/space/work/945396_mm.GRCm38.cdna.all.idx
                wget https://www.dropbox.com/s/zsbrkioi7nocyrl/10xv1_whitelist.txt?dl=1 -O /mnt/space/work/10xv1_whitelist.txt
                wget https://www.dropbox.com/s/ptuhp5qlfs8pfro/10xv2_whitelist.txt?dl=1 -O /mnt/space/work/10xv2_whitelist.txt
                wget https://www.dropbox.com/s/xcxbnvhwfr6b02g/10xv3_whitelist.txt?dl=1 -O /mnt/space/work/10xv3_whitelist.txt
                wget https://www.dropbox.com/s/nacb098h1ubbekw/transcripts_to_genes_v2.txt?dl=1 -O /mnt/space/work/transcripts_to_genes_v2.txt
                python generate.py --yaml_path example/single/fq
                mv example/single/fq/sample_description.csv /mnt/space/work/GSE01/SRS01
                bash /mnt/space/work/GSE01/SRS01/task.bash
                """
            }
        }
        container("continuumio/miniconda3") {
            workDir = "/mnt/space/work/pipeline"
            shellScript {
                content = """
                conda env update -n base --file /mnt/space/work/scn.yml
                python generate.py --yaml_path example/merge
                mv example/merge/SRS01/* /mnt/space/work/GSE02/SRS01
                mv example/merge/SRS02/* /mnt/space/work/GSE02/SRS02
                bash /mnt/space/work/GSE02/merged/task.bash
                """
            }
        }
        container("continuumio/miniconda3") {
            workDir = "/mnt/space/work/pipeline"
            shellScript {
                content = """
                conda env update -n base --file /mnt/space/work/scn.yml
                wget https://www.dropbox.com/s/3aa6x4jk0odfg6u/945396_mm.GRCm38.cdna.all.idx?dl=1 -O /mnt/space/work/945396_mm.GRCm38.cdna.all.idx
                wget https://www.dropbox.com/s/nacb098h1ubbekw/transcripts_to_genes_v2.txt?dl=1 -O /mnt/space/work/transcripts_to_genes_v2.txt
                python generate.py --yaml_path example/fluidigm
                mv example/fluidigm/sample_description.csv /mnt/space/work/GSE153313
                bash /mnt/space/work/GSE153313/task.bash
                """
            }
        }
        container("continuumio/miniconda3") {
            workDir = "/mnt/space/work/pipeline"
            shellScript {
                content = """
                conda env update -n base --file /mnt/space/work/scn.yml
                wget https://www.dropbox.com/s/3aa6x4jk0odfg6u/945396_mm.GRCm38.cdna.all.idx?dl=1 -O /mnt/space/work/945396_mm.GRCm38.cdna.all.idx
                wget https://www.dropbox.com/s/zsbrkioi7nocyrl/10xv1_whitelist.txt?dl=1 -O /mnt/space/work/10xv1_whitelist.txt
                wget https://www.dropbox.com/s/ptuhp5qlfs8pfro/10xv2_whitelist.txt?dl=1 -O /mnt/space/work/10xv2_whitelist.txt
                wget https://www.dropbox.com/s/xcxbnvhwfr6b02g/10xv3_whitelist.txt?dl=1 -O /mnt/space/work/10xv3_whitelist.txt
                wget https://www.dropbox.com/s/nacb098h1ubbekw/transcripts_to_genes_v2.txt?dl=1 -O /mnt/space/work/transcripts_to_genes_v2.txt
                python generate.py --yaml_path example/single/fq_unknown_10x
                mv example/single/fq_unknown_10x/sample_description.csv /mnt/space/work/GSE142345/SRS5863537
                bash /mnt/space/work/GSE142345/SRS5863537/task.bash
                """
            }
        }
    }
}
job("Test panglao mode using SRS1876792 as example dataset") {
    container("continuumio/miniconda3") {
        workDir = "/mnt/space/work/pipeline"
        shellScript {
            content = """
            conda env update -n base --file /mnt/space/work/scn.yml
            python generate.py --yaml_path example/single/panglao
            mv example/single/panglao/* /mnt/space/work/GSE092865/SRS1876792
            bash /mnt/space/work/GSE092865/SRS1876792/task.bash
            """
        }
    }
}