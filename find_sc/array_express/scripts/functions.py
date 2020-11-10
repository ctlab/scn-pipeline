import re
from io import StringIO
from urllib.request import urlopen, Request

import pandas as pd
import requests

pattern = re.compile(
    r"drop[\s-]*seq|cell[\s-]*ranger|chromium|cel[\s-]*seq|scrb[\s-]*seq|sure[\s-]*cell|in[\s-]*drops|fluidigm|smart["
    r"\s-]*seq|10x[\s-]*genomics|mars[\s-]*seq",
    re.IGNORECASE)

drop_seq = re.compile(r"drop[\s-]*seq", re.IGNORECASE)
cel_seq = re.compile(r"cel[\s-]*seq", re.IGNORECASE)
scrb_seq = re.compile(r"scrb[\s-]*seq", re.IGNORECASE)
sure_cell = re.compile(r"sure[\s-]*cell", re.IGNORECASE)
in_drops = re.compile(r"in[\s-]*drops", re.IGNORECASE)
cell_ranger = re.compile(r"cell[\s-]*ranger|chromium|10x|[\s-]*genomics", re.IGNORECASE)
fluidigm = re.compile(r"fluidigm", re.IGNORECASE)
smart_seq = re.compile(r"smart", re.IGNORECASE)
mars_seq = re.compile(r"mars", re.IGNORECASE)


def extract_table(start_date: str, end_date: str) -> pd.DataFrame:
    """
    Get data loaded at arrayexpress db till specified release date
    for Mm, Hs and Rn species
    :param start_date: year-month-day, e.g., 2020-10-01
    :param end_date: year-month-day, e.g., 2020-10-20
    :return: dataframe with samplse uploaded in arrayexpress db during specified time period
    """
    tab_url = 'https://www.ebi.ac.uk/arrayexpress/ArrayExpress-Experiments.txt?keywords=&organism=&exptype%5B%5D' \
              '=%22rna+assay%22&exptype%5B%5D=%22sequencing+assay%22&array=&directsub=on '
    link = Request(tab_url, headers={'User-Agent': 'Mozilla/5.0'})
    tab = urlopen(link)
    tab_content = tab.read()
    data = StringIO(tab_content.decode('latin-1'))
    df = pd.read_csv(data, sep='\t')
    df = df.loc[df.Organism.str.contains('Mus musculus|Homo sapiens|Rattus Norvegicus'), :]
    df = df.loc[df.Type.str.contains('RNA-seq of coding RNA from single cells'), :]
    after_start_date = df['Release Date'] >= start_date
    before_end_date = df['Release Date'] <= end_date
    df = df[after_start_date & before_end_date]
    df = df[['Accession', 'Organism', 'Processed Data', 'Raw Data', 'ArrayExpress URL', 'Release Date']]
    return df


def find_tech(description: str) -> list:
    """
    Find technology using experiment description
    :param description: text description of mtab expriment
    :return: list with matched technologies
    """
    result = []
    if drop_seq.findall(description):
        result.append("DropSeq")
    if cel_seq.findall(description):
        result.append("CELSeq")
        result.append("CELSeq2")
    if scrb_seq.findall(description):
        result.append("SCRBSeq")
    if sure_cell.findall(description):
        result.append("SureCell")
    if smart_seq.findall(description):
        result.append("SmartSeq")
    if in_drops.findall(description):
        result.append("inDrops")
    if cell_ranger.findall(description):
        result.append("10x")
    if fluidigm.findall(description):
        result.append("Fluidigm")
    if mars_seq.findall(description):
        result.append("Mars-Seq")
    return result


def parse_table(df: pd.DataFrame) -> pd.DataFrame:
    """
    Extract meta-information from arrayexpress db table
    :param df: dataframe with information per sample for mtab studies
    :return: meta-info per sample for mtab studies
    """
    all_sp = pd.DataFrame(
        columns=['E-MTAB', 'Comment[ENA_SAMPLE]', 'Comment[BioSD_SAMPLE]', 'Characteristics[organism]',
                 'Comment[ENA_RUN]', 'Comment[BAM_URI]', 'Comment[FASTQ_URI]', 'technology'])
    for idx, row in df.iterrows():
        files = requests.get(
            f"https://www.ebi.ac.uk/arrayexpress/files/{row['Accession']}/{row['Accession']}.sdrf.txt")
        tab_text = StringIO(files.text)
        tab_df = pd.read_csv(tab_text, sep='\t').filter(regex='SAMPLE|organism|ENA_RUN|URI')
        tab_df = tab_df[tab_df.columns.drop(list(tab_df.filter(regex='part')))]
        description = requests.get(
            f"https://www.ebi.ac.uk/arrayexpress/files/{row['Accession']}/{row['Accession']}.idf.txt")
        tab_df['E-MTAB'] = row['Accession']
        tab_df['technology'] = ','.join(set(pattern.findall(description.text)))
        if tab_df['technology'][0]:
            tab_df['technology'] = ','.join(find_tech(tab_df['technology'][0]))
        all_sp = pd.concat([all_sp, tab_df], sort=True)
    return all_sp


def save_results(df: pd.DataFrame, start_date: str, end_date: str) -> None:
    """
    Save full result table and table for each technology
    :param df: dataframe which contains meta-info about each MTAB sample
    :param start_date: year-month-day, e.g., 2020-10-01
    :param end_date: year-month-day, e.g., 2020-10-20
    :return: None
    """
    df.to_csv(f'result/{start_date}_{end_date}/E-MTAB.tsv', sep='\t', index=False)
    df.loc[df['technology'] == '10x'].to_csv(f'result/{start_date}_{end_date}/per_tech/cell_ranger.csv', sep='\t',
                                             index=False)
    df.loc[df['technology'] == 'DropSeq'].to_csv(f'result/{start_date}_{end_date}/per_tech/drop_seq.csv', sep='\t',
                                                 index=False)
    df.loc[df['technology'] == 'inDrops'].to_csv(f'result/{start_date}_{end_date}/per_tech/in_drops.csv', sep='\t',
                                                 index=False)
    df.loc[df['technology'] == 'SCRBSeq'].to_csv(f'result/{start_date}_{end_date}/per_tech/scrb.csv', sep='\t',
                                                 index=False)
    df.loc[df['technology'] == 'Fluidigm'].to_csv(f'result/{start_date}_{end_date}/per_tech/fluidigm.csv', sep='\t',
                                                  index=False)
    df.loc[df['technology'] == 'SureCell'].to_csv(f'result/{start_date}_{end_date}/per_tech/sure_cell.csv', sep='\t',
                                                  index=False)
    df.loc[df['technology'] == 'Mars-Seq'].to_csv(f'result/{start_date}_{end_date}/per_tech/mars_seq.csv', sep='\t',
                                                  index=False)
    df.loc[df['technology'] == 'SmartSeq'].to_csv(f'result/{start_date}_{end_date}/per_tech/smart_seq.csv', sep='\t',
                                                  index=False)
    df.loc[df['technology'].isnull()].to_csv(f'result/{start_date}_{end_date}/per_tech/unknown.csv', sep='\t',
                                             index=False)
