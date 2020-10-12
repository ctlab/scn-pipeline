import re
import xml.etree.ElementTree as ET
from typing import List
from urllib.request import urlopen, Request

import pandas as pd

## DEFINE PATTERNS: DATASET LEVEL, SAMPLE LEVEL

df_pattern = re.compile(
    r"single[\s-]*cell[\s-]*transcriptomic|single[\s-]*cell[\s-]*rna[\s-]*seq|single[\s-]*cell[\s-]*rna[\s-]*sequencing|sc[\s-]*rna[\s-]*seq|drop[\s-]*seq|cell[\s-]*ranger|chromium|cel[\s-]*seq|scrb[\s-]*seq|sure[\s-]*cell|in[\s-]*drops|seurat|fluidigm|smart[\s-]*seq|10x[\s-]*genomics|mars[\s-]*seq",
    re.IGNORECASE)

sample_pattern = re.compile(
    r"drop[\s-]*seq|cell[\s-]*ranger|chromium|cel[\s-]*seq|scrb[\s-]*seq|sure[\s-]*cell|in[\s-]*drops|fluidigm|10x[\s-]*genomics|smart[\s-]*seq",
    re.IGNORECASE)

## COMPILE PATTERNS PER TECHNOLOGY

drop_seq = re.compile(r"drop[\s-]*seq", re.IGNORECASE)
cel_seq = re.compile(r"cel[\s-]*seq", re.IGNORECASE)
scrb_seq = re.compile(r"scrb[\s-]*seq", re.IGNORECASE)
sure_cell = re.compile(r"sure[\s-]*cell", re.IGNORECASE)
in_drops = re.compile(r"in[\s-]*drops", re.IGNORECASE)
cell_ranger = re.compile(r"cell[\s-]*ranger|chromium|10x|[\s-]*genomics", re.IGNORECASE)
fluidigm = re.compile(r"fluidigm", re.IGNORECASE)
smart_seq = re.compile(r"smart", re.IGNORECASE)
mars_seq = re.compile(r"mars", re.IGNORECASE)


def get_tech(found_patterns: str) -> str:
    """
    Transform found patterns to technique-specific format
    :param found_patterns: string of found patterns
    :return: string with found techs in pattern
    """
    techs: List[str] = list()
    if cell_ranger.findall(found_patterns):
        techs.append('10x')
    if drop_seq.findall(found_patterns):
        techs.append('DropSeq')
    if cel_seq.findall(found_patterns):
        techs.append('CELSeq/CELSeq2')
    if in_drops.findall(found_patterns):
        techs.append('inDrops')
    if sure_cell.findall(found_patterns):
        techs.append('SureCell')
    if scrb_seq.findall(found_patterns):
        techs.append('SCRBSeq')
    if smart_seq.findall(found_patterns):
        techs.append('SmartSeq')
    if mars_seq.findall(found_patterns):
        techs.append('MARS-seq')
    if fluidigm.findall(found_patterns):
        techs.append('Fluidigm')
    return ','.join(techs)


def get_geo_ids(start_date: str, end_date: str, org: str) -> list:
    """
    Find GSE which were published during given period for given organism
    :param int_period: number of years / months / days
    :param name_period: name of period (years / moths / days)
    :param org: organism name (e.g., Mus musculus)
    :return: GSE ids corresponding to request
    """
    xml_url = f'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gds&term={start_date}:{end_date}[PDAT]+AND+{re.sub(" ", "%20", org)}[orgn]AND+(gse[ETYP]+OR+gds[ETYP])&retmax=50000&usehistory=y'
    gds_ids = urlopen(xml_url).read()
    gds_tree = ET.fromstring(gds_ids)
    gds_pattern = re.compile(r'^200')
    gse_ids = list()
    for elem in gds_tree.findall('IdList/Id'):
        gse_ids.append(gds_pattern.sub('https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE', elem.text))
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
        tit_pat = get_tech(','.join(sample_pattern.findall(elem.text)))
        if tit_pat:
            sample_data[index].append(tit_pat)
        else:
            sample_data[index].append('NA')
    for index, elem in enumerate(xml_tree.iter('LIBRARY_CONSTRUCTION_PROTOCOL')):
        lib_pat = get_tech(','.join(sample_pattern.findall(elem.text)))
        if lib_pat:
            sample_data[index].append(lib_pat)
        else:
            sample_data[index].append('NA')
    return sample_data


def get_sra_ids(geo_id: str) -> dict:
    meta_data = dict()
    article = Request(geo_id, headers={'User-Agent': 'Mozilla/5.0'})
    geo_page = urlopen(article).read()
    geo_id = re.findall('GSE\d*', geo_page.decode('utf-8'))
    prjna = re.findall('PRJNA\d*', geo_page.decode('utf-8'))
    srp = re.findall('SRP\d*', geo_page.decode('utf-8'))
    super_ser = re.findall('This SuperSeries is composed of the following SubSeries', geo_page.decode('utf-8'))
    meta_data['geo_id'] = geo_id[0] if geo_id else 'NA'
    meta_data['prjna'] = prjna[0] if prjna else 'NA'
    meta_data['srp'] = srp[0] if srp else 'NA'
    meta_data['super_series'] = True if super_ser else False
    return meta_data


def is_rna_seq(study_tree) -> bool:
    """
    Check whether GSE is RNA-seq or not
    :param study_tree: XML-tree from ENA
    :return: True if GSE contains RNA-seq, else None (=False)
    """
    for elem in study_tree.findall('EXPERIMENT_PACKAGE/STUDY/DESCRIPTOR'):
        if elem.find('STUDY_TYPE').attrib['existing_study_type'] == 'Transcriptome Analysis':
            return True


def find_patterns(study_desc) -> set:
    """
    Find patterns correspond to scRNAseq
    :param study_desc: description of GSE study
    :return: set with found patterns
    """
    seen_pattern = set()
    title_pat = df_pattern.findall(study_desc.find('STUDY_TITLE').text)
    abs_pat = df_pattern.findall(study_desc.find('STUDY_ABSTRACT').text)
    if title_pat:
        seen_pattern.add(','.join(title_pat))
    if abs_pat:
        seen_pattern.add(','.join(abs_pat))
    return seen_pattern


def parse_article(study_tree) -> str:
    """
    Parse pubmed article in order to find tech patterns
    :param study_tree: XML-tree from ENA
    :return: tech pattern from the pubmed article text
    """
    try:
        url = Request(
            f"https://www.ncbi.nlm.nih.gov/pmc/articles/pmid/{study_tree.find('EXPERIMENT_PACKAGE/STUDY/STUDY_LINKS/STUDY_LINK/XREF_LINK/ID').text}",
            headers={'User-Agent': 'Mozilla/5.0'})
        link = urlopen(url)
        article = Request(link.url, headers={'User-Agent': 'Mozilla/5.0'})
        article_text = urlopen(article).read()
        return ','.join(find_tech(article_text.decode('utf-8')))
    except:
        return "NA"


def get_full_meta(ena_descr: str) -> pd.DataFrame:
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
