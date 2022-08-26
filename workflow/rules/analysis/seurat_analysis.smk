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
