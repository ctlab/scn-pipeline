import json
from DefineVersionTenX import prepare_files_10x
from DefineTechnologyUtils import *

SPECIES = {
    "mm": "Mus musculus",
    "rn": "Rattus norvegicus",
    "hs": "Homo sapiens"
}

with open(snakemake.input["ffq_json"], "r") as in_file:
    datasets = json.load(in_file)

    for dataset_accession, dataset in datasets.items():
        meta = get_sra_ids(dataset_accession)
        if meta['super_series']:
            raise Exception(f"Dataset {dataset_accession} is Super Series")

        for gsm_accession, gsm_sample in dataset["geo_samples"].items():
            for srs_accession, srs_sample in gsm_sample["samples"].items():
                for srx_accession, srx_experiment in srs_sample["experiments"].items():
                    for srr_accession, srr_run in srx_experiment["runs"].items():
                        technology = parse_srs_for_tech(srs_accession)
                        srr_run["technology"] = technology

                        if technology == Constants.TECH_10X:

                            whitelists = {
                                1: snakemake.input["whitelist_10x_v1"],
                                2: snakemake.input["whitelist_10x_v2"],
                                3: snakemake.input["whitelist_10x_v3"]
                            }

                            files, processing_type, version = prepare_files_10x(
                                srr_accession=srr_accession,
                                files=srr_run["files"],
                                ncbi_dir=snakemake.params["ncbi_dir"],
                                header_dir=snakemake.params["header_dir"],
                                threads=snakemake.threads,
                                white_list_paths=whitelists
                            )

                            if version == -1:
                                raise Exception(f"Could not identify version for {srr_accession}")

                            srr_run["processing_type"] = processing_type
                            srr_run["version"] = version
                            srr_run["whitelist"] = whitelists[version]
                            srr_run["files"] = files

                        else:
                            raise Exception(f"Technology {technology} is yet unsupported for {srr_accession}")

                        # elif technology == Constants.TECH_DROPSEQ:
                        #     pass

    with open(snakemake.output["ffq_full"], "w") as out_file:
        json.dump(datasets, out_file)
