from workflow.scripts.DependencyDispatcher import DependencyDispatcher
from star_input import run_star_input

dispatcher = DependencyDispatcher(config)

include: 'star_jinja.smk'

rule run_star:
    input:
        unpack(run_star_input(dispatcher)),
        star_script=rules.render_star_script.output,
        star_index=dispatcher.star_index,
    output:
        sam=temp("data/samples/{dataset}/{sample}/star/solo/Aligned.out.bam"),
        sj_out_tab="data/samples/{dataset}/{sample}/star/solo/SJ.out.tab",
        solo_barcode_stats="data/samples/{dataset}/{sample}/star/solo/Solo.out/Barcodes.stats",
        solo_summary="data/samples/{dataset}/{sample}/star/solo/Solo.out/Gene/Summary.csv",
        solo_features_stats="data/samples/{dataset}/{sample}/star/solo/Solo.out/Gene/Features.stats",
        solo_umi_per_cell_sorted="data/samples/{dataset}/{sample}/star/solo/Solo.out/Gene/UMIperCellSorted.txt",
        solo_raw_matrix="data/samples/{dataset}/{sample}/star/solo/Solo.out/Gene/raw/matrix.mtx.gz",
        solo_raw_barcodes="data/samples/{dataset}/{sample}/star/solo/Solo.out/Gene/raw/barcodes.tsv.gz",
        solo_raw_features="data/samples/{dataset}/{sample}/star/solo/Solo.out/Gene/raw/features.tsv.gz",
        solo_filtered_matrix="data/samples/{dataset}/{sample}/star/solo/Solo.out/Gene/filtered/matrix.mtx.gz",
        solo_filtered_barcodes="data/samples/{dataset}/{sample}/star/solo/Solo.out/Gene/filtered/barcodes.tsv.gz",
        solo_filtered_features="data/samples/{dataset}/{sample}/star/solo/Solo.out/Gene/filtered/features.tsv.gz",
        log_out="data/samples/{dataset}/{sample}/star/solo/Log.out",
        log_progress_out="data/samples/{dataset}/{sample}/star/solo/Log.progress.out",
        log_final_out="data/samples/{dataset}/{sample}/star/solo/Log.final.out"
    resources:
        mem_mb=32000
    threads: 4
    log: "logs/{dataset}/{sample}/star/run_star.log"
    benchmark: "logs/{dataset}/{sample}/star/run_star.benchmark"
    conda: "../../../envs/star.yaml"
    shell: """
    # exec > {log} 2>&1
    bash {input.star_script:q}
    
    cd `dirname {output.solo_raw_matrix:q}`
    gzip *
    
    cd `dirname {output.solo_filtered_matrix:q}`
    gzip *
    """

