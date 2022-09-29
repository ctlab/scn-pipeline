from DefineTechnologyUtils import *
from pathlib import Path

SPECIES = {
    "mm": "Mus musculus",
    "rn": "Rattus norvegicus",
    "hs": "Homo sapiens"
}


def find_geo_meta(dataset: str, file_study: str, file_sample: str):
    results = parse_gse_meta(dataset)
    if len(results["experiments"]) > 0:
        results = pd.DataFrame.from_dict(results["experiments"], orient="index")
        results["dataset"] = dataset

        print(results)

        aggregated = results.groupby(['dataset']).agg(
            n=pd.NamedAgg(column='alias', aggfunc="count"),
            sample=pd.NamedAgg(column='alias', aggfunc=lambda samples: samples),
            technology=pd.NamedAgg(column='technology', aggfunc=lambda techs: list(set().union(*techs))),
            taxon_id=pd.NamedAgg(column='taxon_id', aggfunc=lambda taxons_ids: list(set(taxons_ids)))
        )

        aggregated.to_csv(file_study, sep='\t', index=False, header=False)
        results.to_csv(file_sample, sep="\t", index=False, header=False)
    else:
        logging.warning(f"No SRX are present in {dataset}: most likely an array")
        Path(file_study).touch()
        Path(file_sample).touch()


if __name__ == '__main__':
    dataset = snakemake.params["dataset"]
    file_study = snakemake.output["file_study"]
    file_sample = snakemake.output["file_sample"]

    find_geo_meta(dataset, file_study, file_sample)

