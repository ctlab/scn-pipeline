SPECIES = ['mm', 'hs']

rule find_single_cell_species:
    output:
        file_study=config['out_dir'] + "/meta/dump/{sp}_study.tsv",
        file_sample=config['out_dir'] + "/meta/dump/{sp}_sample.tsv",
        file_meta=config['out_dir'] + "/meta/dump/{sp}_sample_annotated.tsv"
    benchmark: config['logs_dir'] + "/dump/find_single_cell_species_{sp}.benchmark"
    # log: config['logs_dir'] + "/dump/find_single_cell_species_{sp}.log"
    params:
        species='{sp}',
        start_date="2022/01/01",
        end_date="2022/09/30"
    conda: "../../envs/define_technology.yaml"
    script: "../../scripts/FindSingleCell.py"

rule find_all_single_cell:
    input: expand(config['out_dir'] + "/meta/dump/{sp}_study.tsv", sp=SPECIES),
           expand(config['out_dir'] + "/meta/dump/{sp}_sample.tsv", sp=SPECIES),
           expand(config['out_dir'] + "/meta/dump/{sp}_sample_annotated.tsv", sp=SPECIES)