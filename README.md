# scNavigator processing pipeline

This pipeline is implemented using snakemake, 
so please make sure to install it using mamba.

```bash
conda install -n base -c conda-forge mamba
mamba create -c conda-forge -c bioconda -n snakemake snakemake
```

For all of the steps in this pipeline we have specified 
the minimum environment required to run the step, so please consider running this pipeline using

```bash
snakemake --use-conda ...
```

## Configuration

You first have to configure the project and provide paths to relevant 
files and folder. The configuration is storen in `configs/config.yaml` file and
consists of several fields. The only two necessary fields that are required to fill are `out_dir` and `ncbi_dir`.

`out_dir` is a directory that will be used to store results and preliminary results (and by default resources and logs).

`ncbi_dir` is a directory that configure via `vdb-config` (see [README.vdb-config](https://github.com/ncbi/sra-tools/blob/master/README-vdb-config) in https://github.com/ncbi/sra-tools)
and path to this directory is usually stored in `~/.ncbi/user-settings.mkfg`

```yaml
out_dir: '/path/to/out/dir' # output directory with all the results
ncbi_dir: '/path/to/configured/ncbi/folder' # directory from vdb-config

## Fields below are optional

datasets: '' # default: /path/to/this/repo/configs/datasets.json
samples: '' # default: /out/dir/meta/sample_descriptions.json
resources: '' # default: /out/dir/resources
logs_dir: '' # default: /out/dir/logs
```

## Running pipeline

Once pipeline is configured, fill the datasets in `/config/dataset.yaml`. 
Contents of the file should be just a list of dataset IDs (see example below)

```yaml
["GSE145241", "GSE116240"]
```

Once datasets are specified pipeline consists of two main steps:

* Acquiring meta information 
(we use FFQ + custom scripts to detect single-cell technology and version of chemistry)
* Processing (we use STAR + Seurat for processing and further analysis of the dataset)

To acquire all the meta information run

```bash
snakemake --use-conda get_all_meta
```

To process the datasets for which you already have meta information simply run

```angular2html
snakemake --use-conda process_all
```