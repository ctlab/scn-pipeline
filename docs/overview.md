# Overview of the pipeline

When designing this pipeline we wanted to have easy-to-use pipeline that

* takes ID(s) of the dataset(s)
* returns analyzed data with dimensionality reductions,
clusters, markers, and other standard results that
are routinely generated in scientific publications describing novel scRNA-seq data 

The data that goes into such analysis exists on several levels: SRA Runs containing raw sequences,
SRS Samples (and corresponding GEO samples) that are present in the study, and datasets that contain
cells from different samples.

We decided to use two-step approach when developing this pipeline:

1. We first obtain all required meta information about the dataset:
how many samples are present in the dataset, how many runs are in each sample,
what technology was used to generate the sample 
(datasets can and do contain scRNA-seq samples performed with different scRNA-seq technologies),
what version of chemistry is used (in case of 10x) and so on.
2. Once we know all this information, and we know that our pipeline can run this data, we run the processing part.

This is done using two separate rules `get_all_meta` and `process_all`:

```bash
snakemake --use-conda get_all_meta ## Use this to obtain meta information
snakemake --use-conda process_all ## Use this to process and analyze scRNA-seq datasets
```