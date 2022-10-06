import datetime
import json
from DefineTechnologyUtils import *

SPECIES = {
    "mm": "Mus musculus",
    "rn": "Rattus norvegicus",
    "hs": "Homo sapiens"
}


def find_all_geo(start_date: datetime.date, end_date: datetime.date, species: list[str], output: str):
    # GET GSE IDS LOADED IN TO GEO DB DURING [START_DATE ; END_DATE]
    geo_ids = []
    for sp in species:
        geo_ids.extend(get_geo_ids(start_date, end_date, SPECIES[sp]))

    with open(output, "w") as out_file:
        json.dump(geo_ids, out_file)


def parse_date(date: str) -> datetime.date:
    # assuming date in YYYYMMDD
    date_format = re.compile(r'^(\d\d\d\d)(\d\d)(\d\d)$')
    groups = date_format.match(date).groups()
    return datetime.date(year=int(groups[0]), month=int(groups[1]), day=int(groups[2]))


if __name__ == '__main__':
    start_date = parse_date(snakemake.params["start_date"])
    end_date = parse_date(snakemake.params["end_date"])
    species = snakemake.params["species"]
    output = snakemake.output["datasets"]

    find_all_geo(
        start_date,
        end_date,
        species,
        output
    )

