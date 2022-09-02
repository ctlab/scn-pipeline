# scNavigator processing pipeline

This repository is a Snakemake pipeline that we use to process public scRNA-seq data in scNavigator.
The scope of the main branch is to process the raw data and run basic seurat analysis.

Read the full documentation at https://scn-pipeline.readthedocs.io/.

To install the pipeline, please, clone the main branch and install snakemake
```bash
$ git clone https://github.com/ctlab/scn-pipeline.git 
$ cd scn-pipeline
$ conda install -n base -c conda-forge mamba
$ mamba create -c conda-forge -c bioconda -n snakemake snakemake
$ conda activate snakemake
```

For all the steps in this pipeline we have specified 
the minimum environment required to run the step, so please consider running this pipeline using:

```bash
$ snakemake --use-conda ...
```

## Configuration

You first have to configure the project and provide paths to relevant 
files and folder. The configuration is storen in `configs/config.yaml` file and
consists of several fields. The only two necessary fields that are required to fill are `out_dir` and `ncbi_dir`.

`out_dir` is a directory that will be used to store results and preliminary results (and by default resources and logs).

`ncbi_dir` is a directory that is configured via `vdb-config` 
(see [README.vdb-config](https://github.com/ncbi/sra-tools/blob/master/README-vdb-config) in https://github.com/ncbi/sra-tools)
and path to this directory is usually stored in `~/.ncbi/user-settings.mkfg`

`datasets` is a path to json file that should contain list of IDs

`samples` is the file describing all the meta information collected for the datasets

`resources` is the directory where all resources (whitelists, genome indexes) will be stored

`logs_dir` is the directory to keep the executions logs

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

Once pipeline is configured, fill the datasets in `./config/dataset.yaml`
(or other file listed in `dataset` if you decided to override this value). 
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
$ snakemake --use-conda get_all_meta
```

To process the datasets for which you already have meta information simply run

```bash
$ snakemake --use-conda process_all
```