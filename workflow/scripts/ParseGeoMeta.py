from DefineTechnologyUtils import *

def find_geo_meta(dataset: str, file_study: str, file_sample: str):
    results = parse_gse_meta(dataset)
    df = pd.DataFrame.from_dict(results["experiments"], orient="index")
    df["dataset"] = dataset

    aggregated = df.groupby(['dataset']).agg(
        dataset=pd.NamedAgg(column='dataset', aggfunc=lambda x: x[0]),
        n=pd.NamedAgg(column='alias', aggfunc="count"),
        sample=pd.NamedAgg(column='alias', aggfunc=lambda samples: samples),
        technology=pd.NamedAgg(column='technology', aggfunc=lambda techs: list(set().union(*techs))),
        taxon_id=pd.NamedAgg(column='taxon_id', aggfunc=lambda taxons_ids: list(set(taxons_ids)))
    )

    aggregated["title"] = results["title"]
    aggregated["description"] = results["description"]
    aggregated["link"] = f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={dataset}"
    aggregated.to_csv(file_study, sep='\t', index=False, header=True)
    df.to_csv(file_sample, sep="\t", index=False, header=True)


if __name__ == '__main__':
    dataset = snakemake.params["dataset"]
    file_study = snakemake.output["file_study"]
    file_sample = snakemake.output["file_sample"]

    find_geo_meta(dataset, file_study, file_sample)

