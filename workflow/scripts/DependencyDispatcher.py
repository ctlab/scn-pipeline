from workflow.scripts.Classes import *
from snakemake.io import Wildcards
from typing import Callable


def if_empty_return(ret_val) -> Callable:
    def decorator(func) -> Callable:
        def wrapper(self, *args, **kwargs):
            if self.get_datasets() is None:
                return ret_val
            else:
                return func(self, *args, **kwargs)
        return wrapper
    return decorator


class DependencyDispatcher(object):

    ORGANISM_MAPPING = {
        'Mus musculus': 'mm',
        'mus musculus': 'mm',
        'Homo sapiens': 'hs',
        'homo sapiens': 'hs',
    }

    def __init__(self, config: dict):
        self.samples_file = config['samples']
        self.resources = config["resources"]
        self.out_dir = config["out_dir"]

    def get_datasets(self) -> Dict[str, Dataset]:
        dataset_map = {}
        if os.path.exists(self.samples_file):
            dataset_map = parse_sample_descriptions(json.load(open(self.samples_file, "r")))
        return dataset_map

    def get_dataset(self, wildcards: Wildcards) -> Dataset:
        return self.get_datasets()[wildcards.get('dataset')]

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
    def get_samples(self, wildcards: Wildcards) -> Dict[str, Sample]:
        return self.get_datasets()[wildcards.get("dataset")].samples

    @if_empty_return([])
    def get_sample_names(self, wildcards: Wildcards) -> List[str]:
        return list(self.get_samples(wildcards).keys())

    @if_empty_return(None)
    def get_sample(self, wildcards: Wildcards) -> Sample:
        samples = self.get_samples(wildcards)
        sample = samples[wildcards.get("sample")]
        return sample

    @if_empty_return([])
    def get_runs(self, wildcards: Wildcards) -> List[Run]:
        sample = self.get_sample(wildcards)
        return sample.get_all_runs()

    @if_empty_return(None)
    def get_run(self, wildcards: Wildcards) -> Run:
        runs = self.get_runs(wildcards)
        run = [run for run in runs if run.accession == wildcards.get("run")][0]
        return run

    @if_empty_return([])
    def get_run_names(self, wildcards: Wildcards) -> List[str]:
        return [run.accession for run in self.get_runs(wildcards)]

    @if_empty_return(None)
    def get_species(self, wildcards: Wildcards) -> str:
        sample = self.get_sample(wildcards)
        return DependencyDispatcher.ORGANISM_MAPPING[sample.organism]

    @if_empty_return(None)
    def get_common_species(self, wildcards: Wildcards) -> str:
        dataset = self.get_dataset(wildcards)
        species = set()
        for sample in dataset.samples.values():
            species.add(sample.organism)
        return DependencyDispatcher.ORGANISM_MAPPING[list(species)[0]]

    @if_empty_return(None)
    def get_db(self, wildcards: Wildcards) -> str:
        return self.get_dataset(wildcards).db

    @if_empty_return(None)
    def star_index(self, wildcards: Wildcards) -> str:
        species = self.get_species(wildcards)
        return f"resources/star/{species}/SA"

    @if_empty_return(None)
    def get_technology(self, wildcards: Wildcards) -> str:
        sample = self.get_sample(wildcards)
        return sample.get_technology()

    @if_empty_return(None)
    def get_processing_mode(self, wildcards: Wildcards) -> str:
        sample = self.get_sample(wildcards)
        return sample.get_processing_mode()

    @if_empty_return(None)
    def get_version(self, wildcards: Wildcards) -> int:
        sample = self.get_sample(wildcards)
        return sample.get_version()

    @if_empty_return(None)
    def get_whitelist(self, wildcards: Wildcards) -> str:
        sample = self.get_sample(wildcards)
        return sample.get_whitelist()

    @if_empty_return([])
    def get_barcode_reads(self, wildcards: Wildcards) -> List[str]:
        sample = self.get_sample(wildcards)
        runs = sample.get_all_runs()
        barcodes = [
            run.get_barcode_read().get_path(f"data/samples/{wildcards.get('dataset')}/{wildcards.get('sample')}")
            for run in runs if run.get_barcode_read() is not None
        ]
        return barcodes

    @if_empty_return([])
    def get_cdna_reads(self, wildcards: Wildcards) -> List[str]:
        sample = self.get_sample(wildcards)
        runs = sample.get_all_runs()
        cdna = [
            run.get_cdna_read().get_path(f"data/samples/{wildcards.get('dataset')}/{wildcards.get('sample')}")
            for run in runs if run.get_cdna_read() is not None
        ]
        return cdna

    @if_empty_return([])
    def get_index_reads(self, wildcards: Wildcards) -> List[str]:
        sample = self.get_sample(wildcards)
        runs = sample.get_all_runs()
        index = [
            run.get_index_read().get_path(f"data/samples/{wildcards.get('dataset')}/{wildcards.get('sample')}")
            for run in runs if run.get_index_read() is not None
        ]
        return index

    @if_empty_return([])
    def get_bam(self, wildcards: Wildcards) -> List[str]:
        sample = self.get_sample(wildcards)
        runs = sample.get_all_runs()
        bam = [
            run.get_bam().get_path(f"data/samples/{wildcards.get('dataset')}/{wildcards.get('sample')}")
            for run in runs if run.get_bam() is not None
        ]
        return bam

    @if_empty_return([])
    def get_seurat_objects(self, wildcards: Wildcards) -> List[str]:
        samples = self.get_sample_names(wildcards)

        return [
            f"data/samples/{wildcards.get('dataset')}/{sample}/seurat.rds" for sample in samples
        ]

    @if_empty_return([])
    def get_seurat_objects_forced(self, wildcards: Wildcards) -> List[str]:
        return [file for file in self.get_seurat_objects(wildcards) if os.path.exists(file)]

    @if_empty_return(None)
    def get_file(self, wildcards: Wildcards) -> FFQFile:
        files = self.get_run(wildcards).files
        filename = wildcards.get("filename")
        file = [file for file in files if file.filename == filename][0]
        return file

    @if_empty_return(None)
    def get_file_url(self, wildcards: Wildcards) -> str:
        file = self.get_file(wildcards)
        return file.url

    @if_empty_return(None)
    def get_file_md5(self, wildcards: Wildcards) -> str:
        file = self.get_file(wildcards)
        return file.md5
