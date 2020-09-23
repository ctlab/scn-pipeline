rule srr_to_srs:
    input: srr_files=ancient(get_srr_files_for_srs)
    output: srs_file="out/srss/{srs}.tsv",
    log: "logs/srs/{srs}.log"
    conda: "{{ PathToCondaYml }}"
    params: transcript_gene=config['transcripts_to_genes']
    shell: "Rscript scripts/srr_to_srs.R {output.srs_file} {params.transcript_gene} {input.srr_files} 2> {log}"