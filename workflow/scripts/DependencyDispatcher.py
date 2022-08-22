import os
from .Classes import *
from .Constants import *


class Object(object):
    pass


def if_empty_return(ret_val):
    def decorator(func):
        def wrapper(self, *args, **kwargs):
            if self.get_datasets() is None:
                return ret_val
            else:
                return func(self, *args, **kwargs)
        return wrapper
    return decorator


class DependencyDispatcher(object):

    @staticmethod
    def get_wildcards_placeholder():
        return Object()

    TAX_IDS = {
        9606: 'hs',
        10090: 'mm'
    }

    ORGANISM_MAPPING = {
        'Mus musculus': 'mm',
        'mus musculus': 'mm',
        'Homo sapiens': 'hs',
        'homo sapiens': 'hs',
    }

    def __init__(self, config):
        self.samples_file = config['samples']
        self.resources = config["resources"]
        self.out_dir = config["out_dir"]

    def get_datasets(self) -> Dict[str, Dataset]:
        dataset_list = {}
        if os.path.exists(self.samples_file):
            dataset_list = parse_sample_descriptions(json.load(open(self.samples_file, "r")))
        return dataset_list

    def get_dataset(self, wildcard) -> Dataset:
        return self.get_datasets()[wildcard.dataset]

    @if_empty_return([])
    def get_all_datasets(self) -> List[str]:
        return list(self.get_datasets().keys())

    @if_empty_return([])
    def get_all_samples(self) -> List[str]:
        samples = set()
        for dataset in self.get_datasets().values():
            s = dataset.samples.keys()
            samples = samples.union(s)
        return list(samples)

    @if_empty_return([])
    def get_samples(self, wildcards) -> Dict[str, Sample]:
        dataset = wildcards.dataset
        return self.get_datasets()[dataset].samples

    @if_empty_return([])
    def get_sample_names(self, wildcards) -> List[str]:
        return list(self.get_samples(wildcards).keys())

    @if_empty_return(None)
    def get_sample(self, wildcards) -> Sample:
        samples = self.get_samples(wildcards)
        sample = samples[wildcards.sample]
        return sample

    @if_empty_return([])
    def get_runs(self, wildcards) -> List[Run]:
        sample = self.get_sample(wildcards)
        return sample.get_all_runs()

    @if_empty_return(None)
    def get_run(self, wildcards) -> Run:
        runs = self.get_runs(wildcards)
        run = [run for run in runs if run.accession == wildcards.run][0]
        return run

    @if_empty_return([])
    def get_run_names(self, wildcards) -> List[str]:
        return [run.accession for run in self.get_runs(wildcards)]

    @if_empty_return(None)
    def get_species(self, wildcards) -> str:
        sample = self.get_sample(wildcards)
        return DependencyDispatcher.ORGANISM_MAPPING[sample.organism]

    @if_empty_return(None)
    def get_db(self, wildcards):
        return self.get_dataset(wildcards).db

    @if_empty_return(None)
    def kallisto_idx(self, wildcards):
        species = self.get_species(wildcards)
        return self.resources + f"/kallisto/{species}/index.idx"

    @if_empty_return(None)
    def kallisto_t2g(self, wildcards):
        species = self.get_species(wildcards)
        return self.resources + f"/kallisto/{species}/t2g.txt"

    @if_empty_return(None)
    def alevin_index(self, wildcards):
        species = self.get_species(wildcards)
        return self.resources + f"/salmon/{species}/index"

    @if_empty_return(None)
    def alevin_t2g(self, wildcards):
        species = self.get_species(wildcards)
        return self.resources + f"/salmon/{species}/ref/splici_fl86_t2g_3col.tsv"

    @if_empty_return(None)
    def star_index(self, wildcards):
        species = self.get_species(wildcards)
        return self.resources + f"/star/{species}/SA"

    @if_empty_return(None)
    def get_technology(self, wildcards):
        sample = self.get_sample(wildcards)
        return sample.get_technology()

    @if_empty_return(None)
    def get_processing_mode(self, wildcards):
        sample = self.get_sample(wildcards)
        return sample.get_processing_mode()

    @if_empty_return(None)
    def get_version(self, wildcards):
        sample = self.get_sample(wildcards)
        return sample.get_version()

    @if_empty_return(None)
    def get_whitelist(self, wildcards):
        sample = self.get_sample(wildcards)
        return sample.get_whitelist()

    @if_empty_return(None)
    def get_kallisto_tech_name(self, wildcards):
        sample = self.get_sample(wildcards)
        if sample.get_technology() == TECH_10X:
            version = sample.get_version()
            return f"10xv{version}"

    @if_empty_return(None)
    def get_alevin_tech_name(self, wildcards):
        sample = self.get_sample(wildcards)
        if sample.get_technology() == TECH_10X:
            version = sample.get_version()
            return f"v{version}"

    @if_empty_return([])
    def get_barcode_reads(self, wildcards):
        sample = self.get_sample(wildcards)
        runs = sample.get_all_runs()
        barcodes = [
            run.get_barcode_read().get_path(self.out_dir + f"/data/samples/{wildcards.dataset}/{wildcards.sample}")
            for run in runs if run.get_barcode_read() is not None
        ]
        return barcodes

    @if_empty_return([])
    def get_cdna_reads(self, wildcards):
        sample = self.get_sample(wildcards)
        runs = sample.get_all_runs()
        cdna = [
            run.get_cdna_read().get_path(self.out_dir + f"/data/samples/{wildcards.dataset}/{wildcards.sample}")
            for run in runs if run.get_cdna_read() is not None
        ]
        return cdna

    @if_empty_return([])
    def get_index_reads(self, wildcards):
        sample = self.get_sample(wildcards)
        runs = sample.get_all_runs()
        index = [
            run.get_index_read().get_path(self.out_dir + f"/data/samples/{wildcards.dataset}/{wildcards.sample}")
            for run in runs if run.get_index_read() is not None
        ]
        return index

    @if_empty_return([])
    def get_bam(self, wildcards):
        sample = self.get_sample(wildcards)
        runs = sample.get_all_runs()
        bam = [
            run.get_bam().get_path(self.out_dir + f"/data/samples/{wildcards.dataset}/{wildcards.sample}")
            for run in runs if run.get_bam() is not None
        ]
        return bam
