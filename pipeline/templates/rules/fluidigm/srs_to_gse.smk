rule srs_to_gse:
    input: srs_files=expand(rules.srr_to_srs.output.srs_file, srs=srs_ids)
    output: gse=protected("out/gse/{gse}.tsv")
    log: "logs/srs_to_{gse}.log"
    conda: "{{ PathToCondaYml }}"
    params: transcript_gene=config['transcripts_to_genes']
    shell: "Rscript scripts/srs_to_gse.R {params.transcript_gene} {output.gse} {input.srs_files} 2> {log}"