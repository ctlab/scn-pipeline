SPECIES = ['mm', 'rn', 'hs']

rule all:
    input: expand("out/{sp}_study.tsv", sp=SPECIES),
           expand("out/{sp}_sample.tsv", sp=SPECIES),
           expand("out/{sp}_sample_annotated.tsv", sp=SPECIES)

rule get_files:
    output: file_study="out/{sp}_study.tsv",
            file_sample="out/{sp}_sample.tsv",
            sample_meta="out/{sp}_sample_annotated.tsv"
    conda: "../find_sc.yml"
    benchmark: "benchmarks/get_files/{sp}.txt"
    params: sp_name=lambda wildcards: wildcards.sp, start_date="2019/12/02", end_date="2019/12/02"
    log: "logs/get_files/{sp}.txt"
    shell:
        """
        python scripts/find_sc.py --start_date {params.start_date} --end_date {params.end_date} \
        --sp_name {params.sp_name} --file_std {output.file_study} --file_smpl {output.file_sample} \
        --file_meta {output.sample_meta} --log_file {log}
        """