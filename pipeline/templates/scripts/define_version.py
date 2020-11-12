import argparse
import os.path
import pandas as pd
import re

{% if not bam and not fq_dump and cell_ranger %}
import gzip


def tech_version(fq_r1: str, whielist_10xv2: str, whielist_10xv3: str) -> dict:
    res = {'barcode_len': '', 'technology': 'NA', 'whitelist': 'NA'}
    with gzip.open(fq_r1, 'r') as in_file:
        for idx, line in enumerate(in_file, start=1):
            if idx == 2:
                barcode_len = len(line.decode('ascii').strip())
                res['barcode_len'] = barcode_len
                if barcode_len == 26:
                    res['technology'] = "10xv2"
                    res['whitelist'] = whielist_10xv2
                elif barcode_len == 28:
                    res['technology'] = "10xv3"
                    res['whitelist'] = whielist_10xv3
                else:
                    res['technology'] = "NA"
    if res['technology'] == "10xv3":
        tt_count = 0
        nn_count = 0
        with gzip.open(fq_r1, 'r') as in_file:
            for idx, line in enumerate(in_file, start=1):
                if idx > 1e5:
                    break
                if idx and not idx % 2 and idx % 4:
                    if line.decode('ascii')[26:28] == 'TT':
                        tt_count += 1
                    if line.decode('ascii')[26:28] == 'NN':
                        nn_count += 1
        if ((tt_count + nn_count) / 1e5) > 0.5:  # 4e5 / 4 (lines per read in fq) = 1e5
            res['technology'] = "NA"
    return res


def get_info(fq_r1: str, threads: int, index: str, transcripts_to_genes: str,
             whielist_10xv2: str, whielist_10xv3: str, kallisto_script='kallisto.sh'):
    description = pd.read_csv('sample_description.csv').reset_index().to_dict("records")
    {% if db == 'GEO' %}

    r1 = [re.sub('ftp://', '', x['fastq_ftp'].split(';')[0]) for x in description]
    r2 = [re.sub('ftp://', '', x['fastq_ftp'].split(';')[1]) for x in description]

    {% elif db == 'MTAB' %}

    r1 = [x['fastq_1'] for x in description]
    r2 = [x['fastq_2'] for x in description]

    {% endif %}

    read_info = tech_version(fq_r1, whielist_10xv2, whielist_10xv3)
    if read_info['technology'] == 'NA':
        raise Exception("Technology wasn't defined")
    if not os.path.isfile(kallisto_script):
        with open(kallisto_script, 'w') as out_file:
            out_file.write("mkfifo R1.gz R2.gz\n")
        {% if not test_mode %}

            out_file.write("curl -Ls ftp://{" + ','.join(r1) + "} > R1.gz &\n")
            out_file.write("curl -Ls ftp://{" + ','.join(r2) + "} > R2.gz &\n")
        {% elif test_mode %}

            out_file.write("curl -Ls {" + ','.join(r1) + "} > R1.gz &\n")
            out_file.write("curl -Ls {" + ','.join(r2) + "} > R2.gz &\n")
        {% endif %}

            out_file.write(
                f'kallisto bus -i {index} -x {read_info["technology"]} -t {threads} -o bus_out/ R1.gz R2.gz \n')
            out_file.write('mkdir bus_out/tmp\n')
            out_file.write(
                f'bustools correct -w {read_info["whitelist"]} -o bus_out/correct_output.bus bus_out/output.bus\n')
            out_file.write(
                f'bustools sort -t {threads} -T bus_out/tmp/ -p bus_out/correct_output.bus | bustools count -o bus_out/genes -g {transcripts_to_genes} -e bus_out/matrix.ec -t bus_out/transcripts.txt --genecounts -')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Define 10x version")
    parser.add_argument('--fq_r1', type=str, required=True,
                        help='fastq file header')
    parser.add_argument('--threads', type=int, required=True,
                        help='Number of threads for parallel-fastq-dump program')
    parser.add_argument('--index', type=str, required=True,
                        help='Path to kallisto index')
    parser.add_argument('--transcripts_to_genes', type=str, required=True,
                        help='Path to transcripts_to_genes file')
    parser.add_argument('--white_10xv2', type=str, required=True,
                        help='Path to whitelist for 10xv2')
    parser.add_argument('--white_10xv3', type=str, required=True,
                        help='Path to whitelist for 10xv3')
    args = parser.parse_args()
    get_info(args.fq_r1, args.threads, args.index, args.transcripts_to_genes, args.white_10xv2, args.white_10xv3)

{% elif not bam and fq_dump and cell_ranger %}

import glob
import subprocess

from itertools import product

read_info = dict()
read_info['technology'] = 'NA'


def def_type(read_len: int, postfix: int, whitelist_10xv2: str, whitelist_10xv3: str) -> None:
    if read_len < 12:
        read_info['index'] = postfix
    if read_len == 26:
        read_info['barcode_len'] = read_len
        read_info['barcode'] = postfix
        read_info['technology'] = '10xv2'
        read_info['whitelist'] = whitelist_10xv2
    if read_len == 28:
        read_info['barcode_len'] = read_len
        read_info['barcode'] = postfix
        read_info['whitelist'] = whitelist_10xv3
        read_info['technology'] = '10xv3'
    if read_len > 50:
        read_info['bio'] = postfix

def match_barcodes(white_list: str, read_file: str) -> float:
    """
    Find match percentage of sliced bp seq from reads and barcodes from whitelist
    :param white_list: path to whitelist file
    :param read_file: path to unzipped fq file (only header)
    :return:
    """
    barcodes = set(line.strip() for line in open(white_list))
    with open(read_file, 'r') as file:
        count = 0
        no_of_seqs = 0
        for idx, line in enumerate(file, start=1):
            if idx and not idx % 2 and idx % 4:
                no_of_seqs += 1
                if line[0:16] in barcodes:
                    count += 1
    return 100 * (count / no_of_seqs)

def get_read_info(sra_file: str, whitelist_10xv2: str, whitelist_10xv3: str, fq_dir: str) -> None:
    """
    Define 10x version: find matches between 10x barcodes and reads[0:16] header for each fq read for certain SRR
    from fq header (timeout 30s fastq-dump, see rule fastq_dump)
    If more than 95% of sequences were matched to the certain whitelist, pipeline will use corresponding
    technology in downstream analysis
    :param sra_file: path to sra file
    :param whitelist_10xv2: path to 10xv2 whitelist
    :param whitelist_10xv3: path to 10xv3 whitelist
    :param fq_dir: path to dir with fq headers
    :return: None
    """
    lengths = dict()
    with open(sra_file) as file:
        for idx, line in enumerate(file):
            lengths[idx + 1] = int(float(line.strip()))
    whitelists = [whitelist_10xv2, whitelist_10xv3]
    reads = sorted(glob.glob(f'{fq_dir}/*'))
    res = {pair: match_barcodes(pair[0], pair[1]) for pair in product(whitelists, reads)}
    if max(res.values()) > 95:
        read_info['technology'] = re.findall('10xv\d', max(res, key=res.get)[0])[0]
        read_info['barcode'] = re.sub('_', '', re.findall('_\d', max(res, key=res.get)[1])[0])
        read_info['whitelist'] = max(res, key=res.get)[0]
        other_files = {k: v for k, v in lengths.items() if str(k) != read_info['barcode']}
        read_info['bio'] = max(other_files, key=other_files.get)

def get_info(sra_file, threads: int, index: str, transcripts_to_genes: str, whitelist_10xv2: str,
             whitelist_10xv3: str, tmp_fq: str, fq_dir: str, kallisto_script='kallisto.sh') -> None:
    number_of_files = int(subprocess.getoutput(f"wc -l {sra_file}")[0])
    sra = pd.read_csv(f'{sra_file}', names=['length']).transpose()
    sra.columns = range(1, number_of_files + 1)
    for postfix in range(1, number_of_files + 1):
        def_type(sra[postfix][0], postfix=postfix, whitelist_10xv2=whitelist_10xv2,
                 whitelist_10xv3=whitelist_10xv3)
    if read_info['technology'] == 'NA':
        get_read_info(sra_file, whitelist_10xv2, whitelist_10xv3, tmp_fq)
    if read_info['technology'] == 'NA':
        raise Exception("Technology wasn't defined")
    with open(kallisto_script, 'w') as out_file:
        out_file.write("mkfifo R1.gz R2.gz\n")
        out_file.write(f"cat {fq_dir}/*/*{read_info['barcode']}.fastq.gz > R1.gz &\n")
        out_file.write(f"cat {fq_dir}/*/*{read_info['bio']}.fastq.gz > R2.gz &\n")
        out_file.write(
            f'kallisto bus -i {index} -x {read_info["technology"]} -t {threads} -o bus_out/ R1.gz R2.gz \n')
        out_file.write('mkdir bus_out/tmp\n')
        out_file.write(
            f'bustools correct -w {read_info["whitelist"]} -o bus_out/correct_output.bus bus_out/output.bus\n')
        out_file.write(
            f'bustools sort -t {threads} -T bus_out/tmp/ -p bus_out/correct_output.bus | bustools count -o bus_out/genes -g {transcripts_to_genes} -e bus_out/matrix.ec -t bus_out/transcripts.txt --genecounts -')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Define 10x version")
    parser.add_argument('--sra_file', type=str, required=True,
                        help='File with average read length per each fastq file for corresponding run')
    parser.add_argument('--threads', type=int, required=True,
                        help='Number of threads for parallel-fastq-dump program')
    parser.add_argument('--index', type=str, required=True,
                        help='Path to kallisto index')
    parser.add_argument('--transcripts_to_genes', type=str, required=True,
                        help='Path to transcripts_to_genes file')
    parser.add_argument('--white_10xv2', type=str, required=True,
                        help='Path to whitelist for 10xv2')
    parser.add_argument('--white_10xv3', type=str, required=True,
                        help='Path to whitelist for 10xv3')
    parser.add_argument('--tmp_fq', type=str, required=True,
                        help='Path to temporary fq directory with fq headers')
    parser.add_argument('--fq_dir', type=str, required=True,
                        help='Path to directory with fq files')
    args = parser.parse_args()
    get_info(args.sra_file, args.threads, args.index, args.transcripts_to_genes, args.white_10xv2,
             args.white_10xv3, args.tmp_fq, args.fq_dir)

{% elif bam and not fq_dump %}

import argparse
import re
import warnings

warnings.simplefilter("ignore", UserWarning)
import numpy as np
import operator
import pysam
import pandas as pd


def get_info(s_d: str, tmp_bam: str, index: str, thread: str, transcripts_to_genes: str,
             white_10xv1: str, white_10xv2: str, white_10xv3: str, kallisto_script='kallisto.sh') -> None:
    whitelist_pathes = {"10xv1": white_10xv1,
                        "10xv2": white_10xv2,
                        "10xv3": white_10xv3}
    with open(white_10xv1) as f:
        whitelist_v1 = np.array(f.read().splitlines())
    with open(white_10xv2) as f:
        whitelist_v2 = np.array(f.read().splitlines())
    with open(white_10xv3) as f:
        whitelist_v3 = np.array(f.read().splitlines())
    with pysam.AlignmentFile(tmp_bam, "rb", ignore_truncation=True) as file:
        result = {"10xv1": 0, "10xv2": 0, "10xv3": 0}
        for idx, line in enumerate(file):
            try:
                if idx > 1e3:
                    break
                barcode = re.sub('-[0-9]*', '', line.get_tag('CB'))
                if barcode.strip() in whitelist_v1:
                    result["10xv1"] += 1
                if barcode.strip() in whitelist_v2:
                    result["10xv2"] += 1
                if barcode.strip() in whitelist_v3:
                    result["10xv3"] += 1
            except:
                continue
    max_n = max(result.items(), key=operator.itemgetter(1))[1]
    res = [x[0] for x in result.items() if x[1] == max_n]
    if not max_n or len(res) > 1:
        raise Exception("Technology wasn't defined")
    technology = res[0]
    whitelist = whitelist_pathes[technology]
    description = pd.read_csv(s_d).reset_index().to_dict("records")
    with open(kallisto_script, 'w') as out_file:
        bam = [f"ftp://{bam_file['submitted_ftp'].split(';')[0]}" for bam_file in description]
        file_names = [f'file{idx}.bam' for idx in range(0, len(bam))]
        if len(file_names) > 30:
            raise Exception(f"Too many files (there are {len(file_names)} files), break the process")
        out_file.write(f"mkfifo {' '.join(file_names)}\n")
        if len(file_names) > 1:
            out_file.write("mkfifo merged.bam\n")
        for idx, file in enumerate(file_names):
            out_file.write(f"curl -Ls {bam[idx]} > {file} &\n")
        if technology == "10xv1":
            if len(file_names) > 1:
                out_file.write(f"samtools cat -o merged.bam file* &\n")
                out_file.write("{{ PathToBamToFastq }} --reads-per-fastq=100000000000000 merged.bam out_bam\n")
            else:
                out_file.write(
                    f"{{ PathToBamToFastq }} --reads-per-fastq=100000000000000 {' '.join(file_names)} out_bam\n")
            out_file.write(f"cat out_bam/*/*R1*.gz > R1.gz\n")
            out_file.write(f"cat out_bam/*/*R2*.gz > R2.gz\n")
            out_file.write(f"cat out_bam/*/*R3*.gz > R3.gz\n")
            out_file.write(
                f'kallisto bus -i {index} -x 2,0,14:1,0,10:0,0,0 -t {thread} -o bus_out/ R1.gz R3.gz R2.gz\n')
        else:
            if len(file_names) > 1:
                out_file.write(f"samtools cat -o merged.bam file* &\n")
                out_file.write(
                    f'kallisto bus --bam -i {index} -x {technology} -t {thread} -o bus_out/ merged.bam\n')
            else:
                out_file.write(
                    f'kallisto bus --bam -i {index} -x {technology} -t {thread} -o bus_out/ {" ".join(file_names)}\n')
        out_file.write('mkdir bus_out/tmp\n')
        out_file.write(f'bustools correct -w {whitelist} -o bus_out/correct_output.bus bus_out/output.bus\n')
        out_file.write(
            f'bustools sort -t {thread} -T bus_out/tmp/ -p bus_out/correct_output.bus | bustools count -o bus_out/genes -g {transcripts_to_genes} -e bus_out/matrix.ec -t bus_out/transcripts.txt --genecounts -\n')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Define 10x version")
    parser.add_argument('--s_d', type=str, required=True,
                        help='Path to sample description')
    parser.add_argument('--tmp_bam', type=str, required=True,
                        help='Downloaded header of bam file')
    parser.add_argument('--threads', type=int, required=True,
                        help='Number of threads for kallisto')
    parser.add_argument('--index', type=str, required=True,
                        help='Path to kallisto index')
    parser.add_argument('--transcripts_to_genes', type=str, required=True,
                        help='Path to transcripts_to_genes file')
    parser.add_argument('--white_10xv1', type=str, required=True,
                        help='Path to whitelist for 10xv1')
    parser.add_argument('--white_10xv2', type=str, required=True,
                        help='Path to whitelist for 10xv2')
    parser.add_argument('--white_10xv3', type=str, required=True,
                        help='Path to whitelist for 10xv3')
    args = parser.parse_args()
    get_info(args.s_d, args.tmp_bam , args.index, args.threads, args.transcripts_to_genes,
             args.white_10xv1, args.white_10xv2, args.white_10xv3)
{% endif %}
