from workflow.scripts.DependencyDispatcher import DependencyDispatcher

dispatcher = DependencyDispatcher(config)

def run_star_input(wildcards):
    sample = DependencyDispatcher(config).get_sample(wildcards)
    technology = sample.get_technology()
    processing_mode = sample.get_processing_mode()
    version = sample.get_version()

    result = {}

    if technology == "10X":
        if processing_mode == "fastq" or processing_mode == "fastq_dump":
            result["barcode"] = dispatcher.get_barcode_reads(wildcards)
            result["cdna"] = dispatcher.get_cdna_reads(wildcards)
            if version == 1:
                result["index"] = dispatcher.get_index_reads(wildcards)
    return result


include: 'star_jinja.smk'

rule run_star:
    input:
        unpack(run_star_input),
        star_script=rules.render_star_script.output,
        star_index=dispatcher.star_index,
    output:
        sam=temp(config["out_dir"] + "/data/samples/{dataset}/{sample}/star/solo/Aligned.out.bam"),
        sj_out_tab=config["out_dir"] + "/data/samples/{dataset}/{sample}/star/solo/SJ.out.tab",
        solo_barcode_stats=config["out_dir"] + "/data/samples/{dataset}/{sample}/star/solo/Solo.out/Barcodes.stats",
        solo_summary=config["out_dir"] + "/data/samples/{dataset}/{sample}/star/solo/Solo.out/Gene/Summary.csv",
        solo_features_stats=config["out_dir"] + "/data/samples/{dataset}/{sample}/star/solo/Solo.out/Gene/Features.stats",
        solo_umi_per_cell_sorted=config["out_dir"] + "/data/samples/{dataset}/{sample}/star/solo/Solo.out/Gene/UMIperCellSorted.txt",
        solo_raw_matrix=config["out_dir"] + "/data/samples/{dataset}/{sample}/star/solo/Solo.out/Gene/raw/matrix.mtx.gz",
        solo_raw_barcodes=config["out_dir"] + "/data/samples/{dataset}/{sample}/star/solo/Solo.out/Gene/raw/barcodes.tsv.gz",
        solo_raw_features=config["out_dir"] + "/data/samples/{dataset}/{sample}/star/solo/Solo.out/Gene/raw/features.tsv.gz",
        solo_filtered_matrix=config["out_dir"] + "/data/samples/{dataset}/{sample}/star/solo/Solo.out/Gene/filtered/matrix.mtx.gz",
        solo_filtered_barcodes=config["out_dir"] + "/data/samples/{dataset}/{sample}/star/solo/Solo.out/Gene/filtered/barcodes.tsv.gz",
        solo_filtered_features=config["out_dir"] + "/data/samples/{dataset}/{sample}/star/solo/Solo.out/Gene/filtered/features.tsv.gz",
        log_out=config["out_dir"] + "/data/samples/{dataset}/{sample}/star/solo/Log.out",
        log_progress_out=config["out_dir"] + "/data/samples/{dataset}/{sample}/star/solo/Log.progress.out",
        log_final_out=config["out_dir"] + "/data/samples/{dataset}/{sample}/star/solo/Log.final.out"
    resources:
        mem_mb=32000
    threads: 4
    log: config['logs_dir'] + "/{dataset}/{sample}/star/run_star.log"
    benchmark: config['logs_dir'] + "/{dataset}/{sample}/star/run_star.benchmark"
    conda: "../../../envs/star.yaml"
    shell: """
    # exec > {log} 2>&1
    bash {input.star_script:q}
    
    cd `dirname {output.solo_raw_matrix:q}`
    gzip *
    
    cd `dirname {output.solo_filtered_matrix:q}`
    gzip *
    """

