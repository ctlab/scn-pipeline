from ConcatTables import *
from pprint import pprint
from Constants import *

pd.options.display.width = 100000


studies = read_studies("/mnt/MyPassport/work/new_scn_pipeline/results/meta/dump/20220910_20220910.study.tsv")
print(studies)

samples = read_samples("/mnt/MyPassport/work/new_scn_pipeline/results/meta/dump/20220910_20220910.sample.tsv")
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
                simple_cases.append((dataset_id, dataset["N"], tech, species))
                total_datasets += 1
                total_samples += dataset["N"]

pprint(simple_cases)
print(total_datasets)
print(total_samples)

# any_single_cell = studies[list(map(lambda x: len(x) > 0, studies["tech"]))]
# print(any_single_cell)
#
# simple_cases = any_single_cell[list(map(lambda x: len(x) == 1, any_single_cell["tech"]))]
# simple_cases = simple_cases[list(map(lambda x: len(x) == 1, simple_cases["species"]))]
#
# print(simple_cases)