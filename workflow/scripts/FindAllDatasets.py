import json
from DefineTechnologyUtils import *

SPECIES = {
    "mm": "Mus musculus",
    "rn": "Rattus norvegicus",
    "hs": "Homo sapiens"
}


def find_all_geo(start_date: str, end_date: str, species: list[str], output: str):
    # GET GSE IDS LOADED IN TO GEO DB DURING [START_DATE ; END_DATE]
    geo_ids = []
    for sp in species:
        geo_ids.extend(get_geo_ids(start_date, end_date, SPECIES[sp]))

    with open(output, "w") as out_file:
        json.dump(geo_ids, out_file)


if __name__ == '__main__':
    start_date = snakemake.params["start_date"]
    end_date = snakemake.params["end_date"]
    species = snakemake.params["species"]
    output = snakemake.output["datasets"]

    find_all_geo(
        start_date,
        end_date,
        species,
        output
    )

