## Automated analysis of scRNA-seq data

This pipeline implements strategy to automated analysis scRNA-seq data.

Analysis builds on [Seurat v3.1](https://satijalab.org/seurat/) package using [SCTransform](https://www.biorxiv.org/content/10.1101/576827v2) normalization method and include a way to integrate different dataset with common biological features (e.g., datasets which have obtained using different sequencing technology or datasets from different donors) and batch effect removing using Seurat [integration](https://satijalab.org/seurat/v3.1/integration.html) mode.

## Before you start

Please, before you start, install all the requirements using pip `pip install -r requirements.txt`. Python3 is required.

After you did all that, please refer to _example.yaml_ as YAML-configuration example (check */path* to the real path).

If your data is droplet-based (Drop-Seq, inDrops, SureCell, CELSeq, CELSeq2), your command should be:
`python generate --yaml_path example.yaml`

Elif your data is 10x (10xv1, 10xv2, 10xv3) with fastq files as input, your command should be:
`python generate --yaml_path example.yaml --cell_ranger`

Else your data is 10x (10xv1, 10xv2, 10xv3) with bam file as input, your command should be:
`python generate --yaml_path example.yaml --cell_ranger --bam`