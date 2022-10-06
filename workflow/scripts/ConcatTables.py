import pandas as pd
from ast import literal_eval


def read_studies(file_path: str) -> pd.DataFrame:

    studies = pd.read_table(file_path)
    print(studies)
    studies["technology"] = list(map(literal_eval, studies["technology"]))
    studies["taxon_id"] = list(map(literal_eval, studies["taxon_id"]))
    return studies


def read_samples(file_path: str) -> pd.DataFrame:
    samples = pd.read_table(file_path)
    print(samples)
    samples["technology"] = list(map(lambda x: list(literal_eval(x)), samples["technology"]))
    return samples


if __name__ == "__main__":
    studies = pd.concat(map(read_studies, snakemake.input["study"]))
    samples = pd.concat(map(read_samples, snakemake.input["sample"]))
    studies.to_csv(snakemake.output["study"], sep='\t', index=False, header=True)
    samples.to_csv(snakemake.output["sample"], sep='\t', index=False, header=True)
