rule get_whitelists:
    output:
        whitelist_10x_v1 = "resources/10xv1_whitelist.txt",
        whitelist_10x_v2 = "resources/10xv2_whitelist.txt",
        whitelist_10x_v3 = "resources/10xv3_whitelist.txt"
    conda: "git"
    log: "logs/resources/get_whitelists.log"
    benchmark: "logs/resources/get_whitelists.benchmark"
    shell: """
    wget -o {log} -O {output.whitelist_10x_v1} https://github.com/10XGenomics/cellranger/raw/master/lib/python/cellranger/barcodes/737K-april-2014_rc.txt
    wget -o {log} -O {output.whitelist_10x_v2} https://github.com/10XGenomics/cellranger/raw/master/lib/python/cellranger/barcodes/737K-august-2016.txt
    wget -o {log}  https://github.com/10XGenomics/cellranger/raw/master/lib/python/cellranger/barcodes/3M-february-2018.txt.gz -O - | zcat > {output.whitelist_10x_v3}
    """
