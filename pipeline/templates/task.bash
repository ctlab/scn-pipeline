#!/bin/bash

cd {{ AnalysisFolder }}


{% if Organism == "Mus musculus" and cell_ranger and not test_mode %}
snakemake -j {{ Threads }} --use-singularity --use-conda --conda-prefix {{ PathToCondaDir }} --singularity-prefix {{ PathToSingularityDir }} --singularity-args '--bind {{ IndMus }}:/home --bind {{ WhitelistDir }}:/files --bind {{ CondaFull }}:{{ CondaFull }}' --verbose
{% elif Organism == "Mus musculus" and not cell_ranger and not test_mode %}
snakemake -j {{ Threads }} --use-singularity --use-conda --conda-prefix {{ PathToCondaDir }} --singularity-prefix {{ PathToSingularityDir }} --singularity-args '--bind {{ IndMus }}:/home --bind {{ CondaFull }}:{{ CondaFull }}' --verbose
{% elif Organism == "Homo sapiens" and cell_ranger and not test_mode %}
snakemake -j {{ Threads }} --use-singularity --use-conda --conda-prefix {{ PathToCondaDir }} --singularity-prefix {{ PathToSingularityDir }} --singularity-args '--bind {{ IndHomo }}:/home --bind {{ WhitelistDir }}:/files --bind {{ CondaFull }}:{{ CondaFull }}' --verbose
{% elif Organism == "Homo sapiens" and not cell_ranger and not test_mode %}
snakemake -j {{ Threads }} --use-singularity --use-conda --conda-prefix {{ PathToCondaDir }} --singularity-prefix {{ PathToSingularityDir }} --singularity-args '--bind {{ IndHomo }}:/home --bind {{ CondaFull }}:{{ CondaFull }}' --verbose
{% elif Organism == "Rattus norvegicus" and cell_ranger and not test_mode %}
snakemake -j {{ Threads }} --use-singularity --use-conda --conda-prefix {{ PathToCondaDir }} --singularity-prefix {{ PathToSingularityDir }} --singularity-args '--bind {{ IndRat }}:/home --bind {{ WhitelistDir }}:/files --bind {{ CondaFull }}:{{ CondaFull }}' --verbose
{% elif Organism == "Rattus norvegicus" and not cell_ranger and not test_mode %}
snakemake -j {{ Threads }} --use-singularity --use-conda --conda-prefix {{ PathToCondaDir }} --singularity-prefix {{ PathToSingularityDir }} --singularity-args '--bind {{ IndRat }}:/home --bind {{ CondaFull }}:{{ CondaFull }}' --verbose
{% elif test_mode %}
snakemake -j 4 --verbose
{% endif %}