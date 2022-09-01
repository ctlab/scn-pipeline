# Wildcards and dependency dispatcher

The wildcards that are used in the pipeline are:

* `dataset` - dataset ID
* `sample` - sample ID
* `run` - SRR ID
* `species` - currently only `mm` and `hs` are supported

There are usually several runs in a sample, and several samples in a dataset.
There can be multi-species datasets, but sample typically comes from a single species.

This relationships are kept within the resulting file of
`get_all_meta` (see [how we obtain meta information](meta.md)).
This many-to-one relationships means that in many cases we have to use 
[input functions](https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#input-functions)
for snakemake rules.


## Dependency dispatcher

To facilitate writing input functions we isolated most of the often used functionality into files
`DependencyDispatcher.py` and `Classes.py`. 

Within `Classes.py` we define python classes (`Run`, `Sample`, `Dataset`) and implement serialization 
of FFQ results (which are in JSON) into these classes.

And within `DependencyDispatcher.py` we implement many functions that we use as input functions,
or within input functions. Methods of an instance of `DependencyDispatcher` take wildcards
as an input and return requested properties of the dataset or the sample.

Example

``` py
DependencyDispatcher(config).get_all_datasets() # returns list of dataset IDs
DependencyDispatcher(config).get_sample_names(wildcards) # will search for `dataset` wildcard, and return sample names
DependencyDispatcher(config).get_species(wildcards) # will search for `sample` wildcard, and return species
```