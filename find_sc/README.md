## scn-pipeline/find_sc

## Repository structure

* `get_files.smk` -- snakemake pipeline
* `scripts/functions.py` -- functions for geo parsing
* `scripts/find_sc.py` -- script for geo parsing

## Pipeline structure

![](images/code_structure.png)

## Set some arguments

In params section of `get_files` rule of `get_files.smk` file,
specify both the `start_date` and `end_date`.

If you want to exclude some organism from analysis, remove its notation from `SPECIES` list.

## Code execution

```commandline
snakemake -j 3 -s get_files.smk
```

## Examine output

`sp` -- wildcard for organism, e.g.:
```python
SPECIES =   {
            "mm": "Mus musculus",
            "rn": "Rattus norvegicus",
            "hs": "Homo sapiens"
            }
```

Output files:

```commandline
out/{sp}_study.tsv -- INFO AT STUDY LEVEL
out/{sp}_sample.tsv -- INFO AT SAMPLE LEVEL
out/{sp}_sample_annotated.tsv -- INFO AT METADATA LEVEL
```