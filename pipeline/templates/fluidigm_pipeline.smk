container: "docker://mfiruleva/scn:latest"

{% if Organism == "Mus musculus" and not test_mode %}
configfile: "{{ ConfMus }}"
{% elif Organism == "Homo sapiens" and not test_mode %}
configfile: "{{ ConfHomo }}"
{% elif Organism == "Rattus norvegicus" and not test_mode %}
configfile: "{{ ConfRat }}"
{% elif test_mode %}
configfile: "{{ ConfTest }}"
{% endif %}

import pandas as pd

description = pd.read_csv('sample_description.csv')
srs_ids = set(description['SRS'])


def get_srr_files_for_srs(wildcards) -> list:
    """
    :param wildcards: wildcards with srs ID
    :return: SRR IDs corresponding to the given srs ID
    """
    print(f"Getting SRR for {wildcards.srs}")
    srr_list = description[description['SRS'].values == wildcards.srs]['SRR'].to_list()
    print(f"SRR ids for {wildcards.srs}: {','.join(srr_list)}")
    srr_files = expand(rules.kallisto.output.tsv, srr=srr_list)
    return srr_files


include: "rules/get_data.smk"
include: "rules/kallisto.smk"
include: "rules/srr_to_srs.smk"
include: "rules/srs_to_gse.smk"

rule all:
    input: expand("out/gse/{gse}.tsv", gse="{{ RunName }}")