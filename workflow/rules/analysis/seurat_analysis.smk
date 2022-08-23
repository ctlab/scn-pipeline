RESOLUTIONS = [0.2, 0.4, 0.6, 0.8, 1.0]
DEFAULT_RESOLUTION = 0.6

include: "../processing/star/filter_counts.smk"

rule seurat_analysis:
    input:
        filtered_counts = rules.filter_counts_star.output.filtered_counts
    output:
        seurat = config["out_dir"] + "/data/samples/{dataset}/{sample}/seurat.rds",
        vln_feature_plots = report(config["out_dir"] + "/data/samples/{dataset}/{sample}/plots/vln_feature_plots.pdf"),
        umi_mt_plot= report(config["out_dir"] + "/data/samples/{dataset}/{sample}/plots/umi_mt_plot.pdf"),
        umi_features_log10_plot= report(config["out_dir"] + "/data/samples/{dataset}/{sample}/plots/umi_features_log10_plot.pdf"),
        umi_features_plot= report(config["out_dir"] + "/data/samples/{dataset}/{sample}/plots/umi_features_plot.pdf"),
        vln_feature_plots_after= report(config["out_dir"] + "/data/samples/{dataset}/{sample}/plots/vln_feature_plots_after.pdf"),
        vln_feature_plots_after_log= report(config["out_dir"] + "/data/samples/{dataset}/{sample}/plots/vln_feature_plots_after_log.pdf"),
        umi_mt_plot_after=report(config["out_dir"] + "/data/samples/{dataset}/{sample}/plots/umi_mt_plot_after.pdf"),
        umi_features_log10_plot_after=report(config["out_dir"] + "/data/samples/{dataset}/{sample}/plots/umi_features_log10_plot_after.pdf"),
        umi_features_plot_after=report(config["out_dir"] + "/data/samples/{dataset}/{sample}/plots/umi_features_plot_after.pdf"),
        miqc_plot_prob= report(config["out_dir"] + "/data/samples/{dataset}/{sample}/plots/miqc_plot_prob.pdf"),
        miqc_plot_keep= report(config["out_dir"] + "/data/samples/{dataset}/{sample}/plots/miqc_plot_keep.pdf"),
        elbow_plot= report(config["out_dir"] + "/data/samples/{dataset}/{sample}/plots/elbow_plot.pdf"),
        tsne_plot= report(config["out_dir"] + "/data/samples/{dataset}/{sample}/plots/tsne_plot.pdf"),
        umap_plot= report(config["out_dir"] + "/data/samples/{dataset}/{sample}/plots/umap_plot.pdf"),
        seurat_stats = config["out_dir"] + "/data/samples/{dataset}/{sample}/seurat_stats.json"
    params:
        resolutions = RESOLUTIONS,
        default_resolution = DEFAULT_RESOLUTION,
        sample="{sample}"
    threads: 4
    resources:
        mem_mb=16000
    log: config['logs_dir'] + "/{dataset}/{sample}/seurat_analysis.log"
    benchmark: config['logs_dir'] + "/{dataset}/{sample}/seurat_analysis.benchmark"
    conda: "../../envs/seurat_analysis.yaml"
    script: "../../scripts/seurat_analysis.R"


rule markers:
    input:
        seurat = rules.seurat_analysis.output.seurat
    output:
        marker_dir = directory(config["out_dir"] + "/data/samples/{dataset}/{sample}/markers")
    params:
        resolutions = RESOLUTIONS,
        sample="{sample}"
    threads: 4
    resources:
        mem_mb=16000
    log: config['logs_dir'] + "/{dataset}/{sample}/markers.log"
    benchmark: config['logs_dir'] + "/{dataset}/{sample}/markers.benchmark"
    conda: "../../envs/seurat_analysis.yaml"
    script: "../../scripts/markers.R"

rule markers_default:
    input:
        marker_dir = rules.markers.output.marker_dir
    output:
        markers = config["out_dir"] + "/data/samples/{dataset}/{sample}/markers.tsv",
        clusters = config["out_dir"] + "/data/samples/{dataset}/{sample}/clusters.tsv",
        pct = config["out_dir"] + "/data/samples/{dataset}/{sample}/pct.tsv",
    params:
        default_resolution = DEFAULT_RESOLUTION
    log: config['logs_dir'] + "/{dataset}/{sample}/markers_default.log"
    benchmark: config['logs_dir'] + "/{dataset}/{sample}/markers_default.benchmark"
    shell: """
    cp "{input.marker_dir}/clusters_{params.default_resolution}_pct.tsv" {output.pct}
    cp "{input.marker_dir}/clusters_{params.default_resolution}_average.tsv" {output.clusters}
    cp "{input.marker_dir}/markers_{params.default_resolution}.tsv" {output.markers}
    echo `date`> {log}
    echo "cp successfull" >> {log}
    """