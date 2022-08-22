## How to get resources required to run the pipeline

To properly process the public dataset we will need following files:

* Whitelists for different versions of 10x chemistry (v1, v2 and v3)
* Index files (idx) and transcript to gene (t2g) mapping for kallisto bustools for all species that you want to analyze


### The instructions

To download the whitelists one can run:

```bash
wget -O 10xv1_whitelist.txt https://github.com/10XGenomics/cellranger/raw/master/lib/python/cellranger/barcodes/737K-april-2014_rc.txt
wget -O 10xv2_whitelist.txt https://github.com/10XGenomics/cellranger/raw/master/lib/python/cellranger/barcodes/737K-august-2016.txt
wget https://github.com/10XGenomics/cellranger/raw/master/lib/python/cellranger/barcodes/3M-february-2018.txt.gz -O - | zcat > 10xv3_whitelist.txt
```

In this pipeline we are mostly focused on mouse and human, to obtain corresponding kallisto `idx` and `t2g` files one can run

```bash
kb ref -d mouse -i mm.idx -g mm.t2g.txt
kb ref -d human -i hs.idx -g hs.t2g.txt
```

### We also have a snakemake rule to generate these files

Rule `get_resources` puts these files into file locations provided in the config file.

```bash
snakemake --use-conda -j 1 get_resources
```
