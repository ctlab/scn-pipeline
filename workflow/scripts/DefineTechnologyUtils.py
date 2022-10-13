import logging
import os
import re
import subprocess
import time
import xml.etree.ElementTree as ET
import Constants
import datetime
from io import StringIO
from urllib.request import urlopen, Request, URLError

import pandas as pd
import requests


def get_tech(found_patterns: str) -> set[str]:
    results = []
    for key, value in Constants.ALL_TECH.items():
        if value.findall(found_patterns):
            results.append(key)
    return set(results)


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
        tab_df['technology'] = ','.join(set(Constants.SAMPLE_PATTERN.findall(description.text)))
        if tab_df['technology'][0]:
            tab_df['technology'] = ','.join(get_tech(tab_df['technology'][0]))
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


def get_geo_ids(start_date: datetime.date, end_date: datetime.date, org: str) -> list:
    """
    Find GSE which were published during given period for given organism
    :param start_date: date from which you want to search gse ids (e.g., "2019/12/02")
    :param end_date: date until you want to search gse ids (e.g., "2019/12/05")
    :param org: organism name (e.g., Mus musculus)
    :return: GSE ids corresponding to request
    """
    start_date_str = start_date.strftime("%Y/%m/%d")
    end_date_str = end_date.strftime("%Y/%m/%d")
    xml_url = f'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gds&term={start_date_str}:{end_date_str}[PDAT' \
              f']+AND+{re.sub(" ", "%20", org)}[orgn]AND+(gse[ETYP]+OR+gds[ETYP])&retmax=50000&usehistory=y '
    gds_ids = url_open(xml_url).read()
    gds_tree = ET.fromstring(gds_ids)
    gds_pattern = re.compile(r'^20+')
    gse_ids = list()
    for elem in gds_tree.findall('IdList/Id'):
        gse_ids.append(gds_pattern.sub('GSE', elem.text))
    return gse_ids


def get_sample_info(df_info: dict, xml_tree) -> dict:
    """
    Parse xml from ENA a sample level
    :param df_info: dict with study information (contains GSE id, tax name, sample id, sample name)
    :param xml_tree: xml tree from ENA
    :return: dict with meta-data about each sample in given xml file
    """
    sample_data = dict()
    for index, elem in enumerate(xml_tree.findall('EXPERIMENT_PACKAGE/Pool/Member')):
        if index not in sample_data:
            sample_data[index] = list()
        sample_data[index].append(df_info['geo_id'])
        sample_data[index].append(elem.attrib['organism'])
        sample_data[index].append(elem.attrib['accession'])
        sample_data[index].append(elem.attrib['sample_name'])
    for index, elem in enumerate(xml_tree.findall('LIBRARY_STRATEGY')):
        sample_data[index].append(elem.text)
    for index, elem in enumerate(xml_tree.iter('LIBRARY_SOURCE')):
        sample_data[index].append(elem.text)
    for index, elem in enumerate(xml_tree.findall('EXPERIMENT_PACKAGE/EXPERIMENT/TITLE')):
        tit_pat = get_tech(','.join(Constants.SAMPLE_PATTERN.findall(elem.text)))
        if tit_pat:
            sample_data[index].append(tit_pat)
        else:
            sample_data[index].append('NA')
    for index, elem in enumerate(xml_tree.iter('LIBRARY_CONSTRUCTION_PROTOCOL')):
        lib_pat = get_tech(','.join(Constants.SAMPLE_PATTERN.findall(elem.text)))
        if lib_pat:
            sample_data[index].append(lib_pat)
        else:
            sample_data[index].append('NA')
    return sample_data


def get_sra_ids(geo_id: str) -> dict:
    meta_data = dict()
    url = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=" + geo_id
    article = Request(url, headers={'User-Agent': 'Mozilla/5.0'})
    geo_page = url_open(article).read()
    geo_id = re.findall('GSE\d*', geo_page.decode('utf-8'))
    prjna = re.findall('PRJNA\d*', geo_page.decode('utf-8'))
    srp = re.findall('SRP\d*', geo_page.decode('utf-8'))
    super_ser = re.findall('This SuperSeries is composed of the following SubSeries', geo_page.decode('utf-8'))
    meta_data['geo_id'] = geo_id[0] if geo_id else 'NA'
    meta_data['prjna'] = prjna[0] if prjna else 'NA'
    meta_data['srp'] = srp[0] if srp else 'NA'
    meta_data['super_series'] = True if super_ser else False
    return meta_data


def get_full_meta(ena_descr) -> pd.DataFrame:
    """
    Transform text from ENA to pandas dataframe
    Example:
    https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJNA353098&result=\
    read_run&fields=study_accession,secondary_sample_accession,sample_accession,experiment_accession,\
    run_accession,tax_id,scientific_name,fastq_ftp,submitted_ftp,sra_ftp&format=tsv&download=false
    :param ena_descr: text with all sample ids / run ids / data links for all runs for certain GSE
    :return: input text converted to dataframe
    """
    df = pd.DataFrame([x.split('\t') for x in ena_descr.text.split('\n')])
    df = df.rename(columns=df.iloc[0]).drop(df.shape[0] - 1).drop(df.index[0])
    return df


def find_patterns(experiment_package) -> set:
    """
    Find patterns correspond to scRNAseq
    :param experiment_package: description of GSE study
    :return: set with found patterns
    """
    seen_pattern = set()
    where_to_look = ['EXPERIMENT/DESIGN/LIBRARY_DESCRIPTOR/LIBRARY_CONSTRUCTION_PROTOCOL',
                     'STUDY/DESCRIPTOR/STUDY_TITLE',
                     'STUDY/DESCRIPTOR/STUDY_ABSTRACT']

    for section in where_to_look:
        elements = experiment_package.findall(section)
        for elem in elements:
            if elem.text:
                found_patterns = Constants.DF_PATTERN.findall(elem.text)
                seen_pattern.update(found_patterns)
    return seen_pattern


def is_rna_seq(study_tree) -> bool:
    """
    Check whether GSE is RNA-seq or not
    :param study_tree: XML-tree from ENA
    :return: True if GSE contains RNA-seq, else None (=False)
    """
    for elem in study_tree.findall('EXPERIMENT_PACKAGE/STUDY/DESCRIPTOR'):
        if elem.find('STUDY_TYPE').attrib['existing_study_type'] == 'Transcriptome Analysis':
            return True
    return False


def parse_article(study_tree) -> set[str]:
    """
    Parse pubmed article in order to find tech patterns
    :param study_tree: XML-tree from ENA
    :return: tech pattern from the pubmed article text
    """
    try:
        url = Request(
            f"https://www.ncbi.nlm.nih.gov/pmc/articles/pmid/{study_tree.find('EXPERIMENT_PACKAGE/STUDY/STUDY_LINKS/STUDY_LINK/XREF_LINK/ID').text}",
            headers={'User-Agent': 'Mozilla/5.0'})
        link = url_open(url)
        article = Request(link.url, headers={'User-Agent': 'Mozilla/5.0'})
        article_text = url_open(article).read()
        return get_tech(article_text.decode('utf-8'))
    except:
        return set()


def url_open(url, timeout=2, n_trials=5, sleep_time=2):
    """
    Try open given url during specified timeout
    :param url: url string
    :param timeout: in seconds (e.g., 2 = 2 seconds)
    :param n_trials: number of trials
    :param sleep_time: in seconds (e.g., 2 = 2 seconds to sleep)
    :return: http.client.HTTPResponse object in case of success
    """
    latest_exception = None
    for trial in range(0, n_trials):
        try:
            response = urlopen(url, timeout=timeout)
            return response
        except URLError as e:
            latest_exception = e
            time.sleep(sleep_time)
    print(url)
    raise latest_exception


def get_from_entrez(accession: str, db: str, rettype: str) -> str:
    command = ["efetch", "-db", db, "-id", accession, "-format", rettype]
    logging.warning(command)
    env = os.environ.copy()
    output = subprocess.check_output(" ".join(command), shell=True, env=env)
    return output.decode("utf-8")


def get_xml_from_entrez(accession: str) -> str:
    return get_from_entrez(accession, "sra", "xml")


EFETCH_SUMMARY_RECORD = re.compile(
    r"\n"
    r"(?P<number>\d+)\. (?P<title>.+)\n"
    r"(\(Submitter supplied\) (?P<description>.+)\n)?"
    r"(.+\n)*?"
    r"((SRA Run Selector:\s+https://www\.ncbi\.nlm\.nih\.gov/Traces/study/\?acc=(?P<experiment>SRX\d+)\n)(.+\n)*)?"
    r"(?P<type>Platform|Series|Sample)\s+Accession:\s(?P<accession>\w+)\s+ID:\s(?P<id>\d+)\n"
)


def get_gse_meta(accession) -> str:
    return get_from_entrez(accession, "gds", "summary")


def parse_gse_meta(accession: str):
    meta = get_gse_meta(accession)

    result = {
        'accession': accession,
        'error': '',
        'title': '',
        'description': '',
        'experiments': {},
        'super_series': False
    }

    records = EFETCH_SUMMARY_RECORD.findall(meta)

    platforms = []
    series = None
    samples = []

    for m in EFETCH_SUMMARY_RECORD.finditer(meta):
        result_dict = m.groupdict()
        if result_dict['type'] == 'Platform':
            platforms.append(result_dict)
        elif result_dict['type'] == 'Sample':
            samples.append(result_dict)
        elif result_dict['type'] == "Series":
            series = result_dict

    if series is None:
        result['error'] = "Could not find Series info"
        return result

    result['title'] = series['title']
    result['description'] = series['description']
    super_ser = re.findall('This SuperSeries is composed', series['description'])
    result['super_series'] = True if super_ser else False
    if super_ser:
        for sample in samples:
            result['experiments'][sample['experiment']] = get_tech_from_srx(None, sample['accession'])
    else:
        for sample in samples:
            result['experiments'][sample['experiment']] = get_tech_from_srx(sample['experiment'], sample['accession'])

    return result


def get_tech_from_srx(srr_accession, alias):
    result = {
        'accession': srr_accession,
        'alias': alias,
        'link': f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={alias}",
        'title': '',
        'description': '',
        'error': "",
        'transcriptomic': False,
        'technology': [],
        'taxon_id': 0
    }

    if srr_accession is None:
        result["accession"] = "NO SRX"
        result["error"] = "NO SRX"
        return result

    xml_content = get_xml_from_entrez(srr_accession)
    try:
        tree = ET.fromstring(xml_content)
        error = tree.find("ERROR")

        if error is not None:
            result['error'] = error.text
            return result

        result['taxon_id'] = int(
            tree.find("EXPERIMENT_PACKAGE")
                .find("SAMPLE")
                .find("SAMPLE_NAME")
                .find("TAXON_ID").text
        )

        titles = [
            elem.text
            for index, elem in enumerate(tree.findall('EXPERIMENT_PACKAGE/SAMPLE/TITLE'))
            if elem.text
        ]
        titles.extend([
            elem.text
            for index, elem in enumerate(tree.findall('EXPERIMENT_PACKAGE/EXPERIMENT/TITLE'))
            if elem.text
        ])
        if len(titles) > 0:
            result['title'] = titles[0]

        descriptions = [
            elem.text for index, elem in enumerate(tree.findall('EXPERIMENT_PACKAGE/STUDY/DESCRIPTOR/STUDY_ABSTRACT'))
        ]

        descriptions.extend([
            elem.text for index, elem in enumerate(
                tree.findall('EXPERIMENT_PACKAGE/EXPERIMENT/DESIGN/DESIGN_DESCRIPTION')
            )
        ])

        if len(descriptions) > 0:
            result['description'] = descriptions[0]

        if not is_rna_seq(tree):
            return result
        else:
            result['transcriptomic'] = True
            package = tree.find('EXPERIMENT_PACKAGE')
            seen_patterns = find_patterns(package)
            tech_pattern = get_tech(','.join(seen_patterns))
            result['technology'] = tech_pattern
            return result
    except ET.ParseError as e:
        result["error"] = "Empty response"
        return result


def parse_srs_for_tech(srr_accession):
    xml_content = get_xml_from_entrez(srr_accession)
    tree = ET.fromstring(xml_content)
    if not is_rna_seq(tree):
        raise Exception(f"Run {srr_accession} does not seem to be transcriptomic")
    package = tree.find('EXPERIMENT_PACKAGE')
    seen_patterns = find_patterns(package)
    tech_pattern = get_tech(','.join(seen_patterns))
    if len(tech_pattern) == 0:
        tech_pattern = parse_article(tree)
    if len(tech_pattern) > 1:
        raise Exception(f"Run {srr_accession} is inconsistent and has several techs identified")
    return list(tech_pattern)[0]
