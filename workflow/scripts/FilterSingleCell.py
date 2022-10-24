import json
import os
import random

from ConcatTables import *
from pprint import pprint
from Constants import *
from random import shuffle

pd.options.display.width = 100000
random.seed(0)

os.chdir("/home/askmebefore/dump")

studies = read_studies("~/dump/20220101_20220930.study.tsv")
print(studies)

samples = read_samples("~/dump/20220101_20220930.sample.tsv")
print(samples)

any_single_cell = studies[list(map(lambda x: len(x) > 0, studies["technology"]))]
print(any_single_cell)

pprint(any_single_cell.values)

simple_cases = []
total_datasets = 0
total_samples = 0
for index, dataset in any_single_cell.iterrows():
    dataset_id = dataset["dataset"]
    if len(dataset["taxon_id"]) == 1 and len(dataset["technology"]) == 1:
        tech = dataset["technology"][0]
        species = dataset["taxon_id"][0]
        if tech in SUPPORTED_TECHS:
            subset = samples[samples["dataset"] == dataset_id]
            if all(map(lambda x: tech in x, subset["technology"])):
                simple_cases.append((dataset_id, dataset["n"], tech, species))
                total_datasets += 1
                total_samples += dataset["n"]

shuffle(simple_cases)
pprint(simple_cases)
print(total_datasets)
print(total_samples)


# chunking
chunk_size = 50
chunk_no = 0
i = 0
while i < len(simple_cases):
    chunk = simple_cases[i:i + chunk_size]
    chunk_datasets = list(map(lambda x: x[0], chunk))
    with open(f"datasets_chunk_{chunk_no}.json", "w") as f:
        json.dump(chunk_datasets, f)
    i += chunk_size
    chunk_no += 1



# any_single_cell = studies[list(map(lambda x: len(x) > 0, studies["tech"]))]
# print(any_single_cell)
#
# simple_cases = any_single_cell[list(map(lambda x: len(x) == 1, any_single_cell["tech"]))]
# simple_cases = simple_cases[list(map(lambda x: len(x) == 1, simple_cases["species"]))]
#
# print(simple_cases)