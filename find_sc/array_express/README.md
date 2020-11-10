# Find scRNA-seq experiments: E-MTAB

## Description

The main aim of this pipeline is to find scRNA-seq experiments in arrayexpress database,
annotate and prepare them for downstream analysis.

## Dependencies

All dependencies are listed in the *../find_sc.yml* file.

## Run pipeline

```
snakemake -j 1 --use-conda --conda-prefix $(pwd) --config start_date=2020.10.01 end_date=2020.10.15
```

## Results

All result will be saved to *results/{start_date}_{end_date}* directory:
for example, *results/2020-10-01_2020-10-15/E-MTAB.tsv*
Additionally, you can read the log from *results/{start_date}_{end_date}.txt* file and
benchmark information from the *benchmarks/{start_date}_{end_date}.txt* file.