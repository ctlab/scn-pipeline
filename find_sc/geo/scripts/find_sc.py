import argparse
import re
import xml.etree.ElementTree as ET
from urllib.request import urlopen

import pandas as pd
import requests
from functions import *

SPECIES = {
    "mm": "Mus musculus",
    "rn": "Rattus norvegicus",
    "hs": "Homo sapiens"
          }

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Find scRNA-seq data in GEO database")
    parser.add_argument('--start_date', type=str, help='Specify the date from which you want to search for data, e.g., 2020/08')
    parser.add_argument('--end_date', type=str, help='Specify the date after which you do not want to search for data, e.g., 2020/10')
    parser.add_argument('--sp_name', type=str, help='Organism name, e.g., Mus musculus')
    parser.add_argument('--file_std', type=str, help='Path to output at study level')
    parser.add_argument('--file_smpl', type=str, help='Path to output at sample level')
    parser.add_argument('--file_meta', type=str, help='Path to output at metadata level')
    parser.add_argument('--log_file', type=str, help='Path to log file')
    args = parser.parse_args()

    ## GET GSE IDS LOADED IN TO GEO DB DURING [START_DATE ; END_DATE]
    geo_ids = get_geo_ids(args.start_date, args.end_date, SPECIES[args.sp_name])

    ## PARSE EACH GSE ID

    for geo in geo_ids:
        try:
            res = get_sra_ids(geo)
            ## SKIP SUPER_SERIES
            if res['super_series']:
                continue
            ## GET ENA LINK USING SRP OR PRJNA ID
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
            response = url_open(url).read()
            tree = ET.fromstring(response)
            try:
                if is_rna_seq(tree):
                    for idx, el in enumerate(tree.findall('EXPERIMENT_PACKAGE/STUDY/DESCRIPTOR')):
                        if not idx:
                            if el.find('STUDY_TYPE').attrib['existing_study_type'] == 'Transcriptome Analysis':
                                tax = {el.attrib['tax_id'] for el in tree.findall('EXPERIMENT_PACKAGE/Pool/Member')}

                                ## FIND PATTERNS
                                seen_patterns = find_patterns(el)
                                if seen_patterns:

                                    ## ADD SOME INFO
                                    res['tax'] = ','.join(map(str, tax))
                                    tech_pattern = get_tech(','.join(seen_patterns))

                                    ## IF YOU CAN'T FIND ANY TECH PATTERN, TRY TO PARSE ARTICLE AT PUBMED
                                    if tech_pattern:
                                        res['pattern'] = tech_pattern
                                    else:
                                        res['pattern'] = parse_article(tree)

                                    ## WRITE INFO AT STYDY, SAMPLE AND META-DATA (IDS, LINKS, ETC) LEVELS
                                    study = pd.DataFrame.from_dict({k: v for k, v in res.items()}, orient='index').T
                                    sample = pd.DataFrame.from_dict(get_sample_info(res, tree), orient='index')
                                    study.to_csv(args.file_std, sep='\t', mode='a', index=False, header=False)
                                    sample.to_csv(args.file_smpl, sep="\t", mode='a', index=False, header=False)
                                    meta = get_full_meta(ena)
                                    meta.insert(0, 'GSE', f'{res["geo_id"]}')
                                    meta.to_csv(args.file_meta, sep="\t", mode='a', index=False, header=False)
            except:
                with open(args.log_file, 'a+') as f:
                    f.write(f'{geo}\n')
        except:
            with open(args.log_file, 'a+') as f:
                f.write(f'{geo}\n')


