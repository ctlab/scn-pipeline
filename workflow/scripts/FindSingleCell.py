import logging
from DefineTechnologyUtils import *

SPECIES = {
    "mm": "Mus musculus",
    "rn": "Rattus norvegicus",
    "hs": "Homo sapiens"
}

def find_all_single_cells_geo(
        start_date: str,
        end_date: str,
        species: str,
        file_study: str,
        file_sample: str,
        file_meta: str):

    # GET GSE IDS LOADED IN TO GEO DB DURING [START_DATE ; END_DATE]
    geo_ids = get_geo_ids(start_date, end_date, SPECIES[species])

    # PARSE EACH GSE ID
    for geo in geo_ids:
        print(geo)
        try:
            res = get_sra_ids(geo)
            print(res)
            # SKIP SUPER_SERIES
            if res['super_series']:
                continue
            # GET ENA LINK USING SRP OR PRJNA ID
            if res['srp'] != 'NA':
                url = f'http://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?save=efetch&db=sra&rettype=runtable&term={res["srp"]}'
                ena = requests.get(
                    f"https://www.ebi.ac.uk/ena/portal/api/filereport?accession={res['srp']}&result=read_run&fields"
                    f"=study_accession,secondary_sample_accession,sample_accession,experiment_accession,run_accession,tax_id,scientific_name,"
                    f"fastq_ftp,submitted_ftp,sra_ftp&format=tsv&download=false")
            elif res['prjna'] != 'NA':
                url = f'http://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?save=efetch&db=sra&rettype=runtable&term={res["prjna"]}'
                ena = requests.get(
                    f"https://www.ebi.ac.uk/ena/portal/api/filereport?accession={res['prjna']}&result=read_run&fields"
                    f"=study_accession,secondary_sample_accession,sample_accession,experiment_accession,run_accession,tax_id,scientific_name,"
                    f"fastq_ftp,submitted_ftp,sra_ftp&format=tsv&download=false")
            else:
                continue

            ## LOAD XML TREE WITH BOTH STUDY/SAMPLE INFORMATION LEVELS
            response = urlopen(url).read()
            print(response)
            tree = ET.fromstring(response)
            try:
                if is_rna_seq(tree):
                    for idx, el in enumerate(tree.findall('EXPERIMENT_PACKAGE/STUDY/DESCRIPTOR')):
                        if el.find('STUDY_TYPE').attrib['existing_study_type'] == 'Transcriptome Analysis':
                            tax = {el.attrib['tax_id'] for el in tree.findall('EXPERIMENT_PACKAGE/Pool/Member')}
                            print("Transcriptome ")
                            # FIND PATTERNS
                            seen_patterns = find_patterns(el)
                            if seen_patterns:

                                print(f"Seen patterns: {seen_patterns}")
                                # ADD SOME INFO
                                res['tax'] = ','.join(map(str, tax))
                                tech_pattern = get_tech(','.join(seen_patterns))

                                # IF YOU CAN'T FIND ANY TECH PATTERN, TRY TO PARSE ARTICLE AT PUBMED
                                if tech_pattern:
                                    res['pattern'] = tech_pattern
                                else:
                                    res['pattern'] = parse_article(tree)

                                # WRITE INFO AT STYDY, SAMPLE AND META-DATA (IDS, LINKS, ETC) LEVELS
                                study = pd.DataFrame.from_dict({k: v for k, v in res.items()}, orient='index').T
                                sample = pd.DataFrame.from_dict(get_sample_info(res, tree), orient='index')
                                study.to_csv(file_study, sep='\t', mode='a', index=False, header=False)
                                sample.to_csv(file_sample, sep="\t", mode='a', index=False, header=False)
                                meta = get_full_meta(ena)
                                meta.insert(0, 'GSE', f'{res["geo_id"]}')
                                meta.to_csv(file_meta, sep="\t", mode='a', index=False, header=False)
            except Exception as e:
                logging.info(f'{res["geo_id"]}: {e}')
        except Exception as e:
                logging.info(f'{geo}: {e}')


if __name__ == '__main__':
    start_date = snakemake.params["start_date"]
    end_date = snakemake.params["end_date"]
    species = snakemake.params["species"]
    file_study = snakemake.output["file_study"]
    file_sample = snakemake.output["file_sample"]
    file_meta = snakemake.output["file_meta"]

    find_all_single_cells_geo(
        start_date,
        end_date,
        species,
        file_study,
        file_sample,
        file_meta
    )

