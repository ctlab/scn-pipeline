rule sample_report:
    output:
        report = config["out_dir"] + "/reports/{dataset}/{sample}.html"
    input:
        sample_description = config["out_dir"] + "/data/samples/{dataset}/{sample}/sample_description.csv",
        kallisto_run_info = config["out_dir"] + "/data/samples/{dataset}/{sample}/bus_out/run_info.json",
        bus_full_info = config["out_dir"] + "/data/samples/{dataset}/{sample}/bus_out/bus_full_info.json",
        bus_correct_info = config["out_dir"] + "/data/samples/{dataset}/{sample}/bus_out/bus_correct_info.json",
        filtering_stats = config["out_dir"] + "/data/samples/{dataset}/{sample}/filtering_stats.json",
        mtx = config["out_dir"] + "/data/samples/{dataset}/{sample}/bus_out/genes.mtx",
        genes = config["out_dir"] + "/data/samples/{dataset}/{sample}/bus_out/genes.genes.txt",
        barcodes = config["out_dir"] + "/data/samples/{dataset}/{sample}/bus_out/genes.barcodes.txt",
        filtered_counts = config["out_dir"] + "/data/samples/{dataset}/{sample}/filtered_counts.rds",
        seurat_stats = config["out_dir"] + "/data/samples/{dataset}/{sample}/seurat_stats.json",
        seurat_object = config["out_dir"] + "/data/samples/{dataset}/{sample}/seurat.rds"
    params:
        sample = "{sample}",
        dataset = "{dataset}"
    conda: "../../envs/reports.yaml"
    script: "../../scripts/sample_report.Rmd"

