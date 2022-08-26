import yaml
import json
import os
from functools import reduce
from typing import List, Dict, Union, Any
from workflow.scripts.Constants import *


class FFQJSonEncoder(json.JSONEncoder):
    def default(self, o: Any) -> Any:
        if isinstance(o, JsonYamlSerializable):
            return o.__dict__
        elif isinstance(o, set):
            return list(o)
        else:
            return super().default(o)


class JsonYamlSerializable(object):
    def __repr__(self):
        return yaml.dump(self.__dict__)

    def yaml(self) -> str:
        return yaml.dump(self.__dict__)

    def json(self) -> str:
        return FFQJSonEncoder().encode(self.__dict__)

    @classmethod
    def load(cls, values):
        return cls.__new__(**values)

    @classmethod
    def load_yaml(cls, data: str):
        values = yaml.safe_load(data)
        return cls.load(values)

    @classmethod
    def load_json(cls, data: str):
        values = json.loads(data)
        return cls.load(values)


class FFQFile(JsonYamlSerializable):
    def __init__(self, **kwargs):
        self.accession: str = kwargs.get("accession")
        self.filename: str = kwargs.get("filename")
        self.filetype: str = kwargs.get("filetype")
        self.filesize: int = kwargs.get("filesize")
        self.filenumber: int = kwargs.get("filenumber")
        self.md5: str = kwargs.get("md5", "")
        self.urltype: str = kwargs.get("urltype", "")
        self.url: str = kwargs.get("url", "")
        self.barcode_length: int = kwargs.get("barcode_length", 0)
        self.type: str = kwargs.get("type", "")

    def get_path(self, prefix):
        return os.path.join(prefix, f"files/{self.accession}/{self.filetype}/{self.filename}")


class Run(JsonYamlSerializable):
    def __init__(self, **kwargs):
        self.accession: str = kwargs.get("accession")
        self.experiment: str = kwargs.get("experiment")
        self.study: str = kwargs.get("study", "")
        self.sample: str = kwargs.get("sample", "")
        self.title: str = kwargs.get("title", "")
        self.attributes: Dict[str, str] = kwargs.get("attributes", dict())

        if isinstance(kwargs["files"], dict):
            self.files: List[FFQFile] = [FFQFile(**file) for file in kwargs["files"]["ftp"]]
        elif isinstance(kwargs["files"], list):
            self.files: List[FFQFile] = [FFQFile(**file) for file in kwargs["files"]]
        else:
            self.files = None

        self.processing_type = kwargs.get("processing_type", "")
        self.technology: List[str] = kwargs.get("technology", "")
        self.version: int = kwargs.get("version", -1)
        self.whitelist: str = kwargs.get("whitelist", "")

    def get_processing_mode(self):
        return self.processing_type

    def get_bam(self):
        for file in self.files:
            if file.type == TYPE_BAM:
                return file

    def get_barcode_read(self):
        for file in self.files:
            if file.type == TYPE_BARCODE:
                return file

    def get_cdna_read(self):
        for file in self.files:
            if file.type == TYPE_CDNA:
                return file

    def get_index_read(self):
        index_files = [file for file in self.files if file.type == TYPE_INDEX]
        if len(index_files) == 0:
            return None
        elif len(index_files) == 1:
            return index_files[0]
        else:
            # select max read length
            index_files_ordered = sorted(index_files, key=lambda file: -file.read_length)
            return index_files_ordered[0]


class Experiment(JsonYamlSerializable):
    def __init__(self, **kwargs):
        self.accession: str = kwargs.get("accession")
        self.title: str = kwargs.get("title", "")
        self.platform: str = kwargs.get("platform", "")
        self.instrument: str = kwargs.get("instrument", "")
        self.runs: Dict[str, Run] = {key: Run(**run) for key, run in kwargs["runs"].items()}


class Sample(JsonYamlSerializable):
    def __init__(self, **kwargs):
        self.accession: str = kwargs.get("accession")
        self.title: str = kwargs.get("title", "")
        self.organism: str = kwargs.get("organism", "")
        self.attributes: Dict[str, str] = kwargs.get("attributes", dict())
        self.experiments: Dict[str, Experiment] = \
            {key: Experiment(**experiment) for key, experiment in kwargs['experiments'].items()}

    def get_all_runs(self):
        runs = []
        for experiment in self.experiments.values():
            for run in experiment.runs.values():
                runs.append(run)
        return runs

    def get_processing_mode(self):
        modes = list(set([run.get_processing_mode() for run in self.get_all_runs()]))
        if len(modes) > 1:
            raise Exception(f"Sample {self.accession} is inconsistent - multiple modes found: {modes}")
        return modes[0]

    def get_technology(self):
        technologies = list(set([run.technology for run in self.get_all_runs()]))
        if len(technologies) > 1:
            raise Exception(f"Sample {self.accession} is inconsistent - multiple technologies found: {technologies}")
        return technologies[0]

    def get_version(self):
        versions = list(set([run.version for run in self.get_all_runs()]))
        if len(versions) > 1:
            raise Exception(f"Sample {self.accession} is inconsistent - multiple versions found: {versions}")
        return versions[0]

    def get_whitelist(self):
        whitelists = list(set([run.whitelist for run in self.get_all_runs()]))
        if len(whitelists) > 1:
            raise Exception(f"Sample {self.accession} is inconsistent - multiple whitelists found: {whitelists}")
        return whitelists[0]

    @staticmethod
    def parse_geo_sample(values):
        geo_meta = {
            'accession': values["accession"],
            'supplementary_files': [
                FFQFile(**file) for file in values.get("supplementary_files", [])
            ],
            'platform': values.get("platform", {})
        }
        sample_id = list(values["samples"].keys())[0]
        sample = Sample(**values["samples"][sample_id])
        sample.geo_meta = geo_meta
        return sample


class Dataset(JsonYamlSerializable):
    GEO_DB_REGEX = re.compile("GSE\d+")
    ENA_DB_REGEX = re.compile("ERP\d+")
    GEO = "GEO"
    ENA = "ENA"

    def __init__(self, **kwargs):
        self.accession: str = kwargs.get("accession")
        self.title: str = kwargs.get("title", "")
        self.db: Union[str, None] = None
        self.samples: Dict[str, Sample] = {}
        self.organisms: set[str] = set()
        self.attributes: Dict[str, str] = {}

        if Dataset.GEO_DB_REGEX.match(self.accession):
            self.db = Dataset.GEO
        elif Dataset.ENA_DB_REGEX.match(self.accession):
            self.db = Dataset.ENA
        else:
            raise Exception(f"Dataset {self.accession} is from unknown db")

        if self.db == Dataset.GEO and "geo_samples" in kwargs:
            for key, sample in kwargs["geo_samples"].items():
                self.samples[key] = Sample.parse_geo_sample(sample)
        else:
            for key, sample in kwargs["samples"].items():
                self.samples[key] = Sample(**sample)

        self.organisms = set(
            [sample.organism for sample in self.samples.values()]
        )

        self.attributes = reduce(
            lambda attr1, attr2: {key: value for key, value in attr1.items() if attr2.get(key, None) == value},
            [sample.attributes for sample in self.samples.values()]
        )


def parse_sample_descriptions(values: dict) -> dict[str, Dataset]:
    return {key: Dataset(**value) for key, value in values.items()}

