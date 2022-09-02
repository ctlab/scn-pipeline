include: "get_star_index.smk"
include: "get_whitelists.smk"


species = ["hs", "mm"]

rule get_all_resources:
    input:
        rules.get_whitelists.output,
        expand(rules.get_cellranger_fasta.output, species=species),
        expand(rules.get_star_index.output, species=species)



