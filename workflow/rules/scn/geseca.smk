from workflow.scripts.DependencyDispatcher import DependencyDispatcher

dispatcher = DependencyDispatcher(config)

def get_species_name(wildcards):
    species = dispatcher.get_species(wildcards)
    species_mapping =  {
        "mm": "mouse",
        "hs": "human"
    }
    return species_mapping[species]


include: "../analysis/seurat_analysis.smk"

rule geseca:
    input:
        seurat = rules.seurat_analysis.output.seurat
    output:
        geseca_full="data/samples/{dataset}/{sample}/geseca/" + f"geseca_full.tsv",
        geseca_var="data/samples/{dataset}/{sample}/geseca/" + f"geseca_var.tsv"
    params:
        sample="{sample}",
        species=get_species_name
    threads: 4
    resources: mem_mb=32000
    log: "logs/{dataset}/{sample}/geseca.log"
    benchmark: "logs/{dataset}/{sample}/geseca.benchmark"
    conda: "../../envs/geseca.yaml"
    script: "../../scripts/geseca.R"