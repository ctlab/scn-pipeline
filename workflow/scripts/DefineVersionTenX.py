import re
import warnings
import gzip
import os
import logging
import json
import subprocess
import pysam
import operator
import numpy as np
from Constants import *

warnings.simplefilter("ignore", UserWarning)


def match_barcodes_10x_fq_r1(fastq_path: str, white_list_paths: str) -> float:
    """
    Find match percentage of sliced bp seq from reads and barcodes from whitelist
    :param fastq_path: path to gzipped fq file (only header)
    :param white_list_paths: path to whitelist file
    :return: ratio of whitelisted reads
    """
    barcode_list = [line.strip() for line in open(white_list_paths)]
    barcode_length = len(barcode_list[0])
    barcode_set = set(barcode_list)
    with gzip.open(fastq_path, 'r') as in_file:
        count = 0
        no_of_seqs = 0
        for idx, line in enumerate(in_file, start=1):
            if idx % 2 == 0 and idx % 4 != 0:
                no_of_seqs += 1
                if line.decode('utf-8')[0:barcode_length] in barcode_set:
                    count += 1
    return count / no_of_seqs


def get_read_info_10x_fq_r1(fastq_path: str, white_list_paths: dict[int, str]) -> int:
    """
    Define 10x version:
    find matches between 10x barcodes and fastq file
    If more than `RATIO_THRESHOLD` of sequences were matched to the certain whitelist, pipeline will use corresponding
    technology in downstream analysis
    :param fastq_path: path to gzipped fq file (only header)
    :param white_list_paths: dict pointing version number to whitelists
    :return: dict
    """

    current_max = RATIO_THRESHOLD
    current_max_version = -1

    results = [(version, match_barcodes_10x_fq_r1(fastq_path, whitelist))
               for version, whitelist in white_list_paths.items()]
    stats_string = ", ".join([f"v{version}: {percent:.2f}" for version, percent in results])
    logging.warning(f"Stats for {fastq_path}: {stats_string}")

    for version, match_ratio in results:
        if match_ratio > current_max:
            current_max = match_ratio
            current_max_version = version
    return current_max_version


def get_read_length_fastq(fastq_path: str) -> int:
    """
    Get read length from the fastq file. Reads second line in the file
    :param fastq_path: path to (gzipped) fastq file to check read length
    :type fastq_path: str
    """
    with gzip.open(fastq_path, 'r') as fq_file:
        fq_file.readline()  # skipping fastq header
        line = fq_file.readline()  # read line
        read_length = len(line.decode('ascii').strip())
        return read_length


def get_read_length_sam(sam_path: str) -> int:
    barcode = ""
    with pysam.AlignmentFile(sam_path) as file:
        for idx, line in enumerate(file):
            try:
                barcode = re.sub('-[0-9]*', '', line.get_tag('CB')).strip()
                break
            except KeyError:
                continue
            except (OSError, StopIteration) as e:
                break
    if barcode:
        return len(barcode)
    else:
        raise KeyError(f"Tag CB is not found in {sam_path}")


def tech_version_10x_fq_r1(fastq_path: str, white_list_paths: dict[int, str]) -> dict:
    """
    Define 10x version
    :param fastq_path: path to gzipped fq file (only header)
    :param white_list_paths dict from version number to path to whitelist
    :return: dict result
    """
    res = {
        'barcode_length': 0,
        'version': -1,
        'whitelist': None
    }

    read_length = get_read_length_fastq(fastq_path)
    version = get_read_info_10x_fq_r1(fastq_path, white_list_paths)
    res['barcode_length'] = read_length
    res['version'] = version
    if version != -1:
        res['whitelist'] = white_list_paths[version]
    return res


def get_read_type_fastq(fastq_paths: list[str], white_list_paths: dict[int, str]):
    """
    For 10x technology identify which read file is which (cdna, barcode and index) and identify chemistry version
    :type fastq_paths: list[str]
    :param fastq_paths: list of fastq files (or headers) to check
    :type white_list_paths: dict[int, str]
    :param white_list_paths: dict from version (int) to path to whitelist (str)
    """
    sorted_by_length = sorted(fastq_paths, key=get_read_length_fastq, reverse=True)
    results = [
        tech_version_10x_fq_r1(file, white_list_paths)
        for file in sorted_by_length
    ]

    barcode_idx = -1
    version = -1
    results_barcode = list(filter(lambda res: res["version"] != -1, results))

    if len(results_barcode) > 1:
        Exception(f"Too many files match whitelists for single run. Files {results_barcode}")

    files_n = len(fastq_paths)
    for index in range(files_n):
        if results[index]["version"] != -1:
            barcode_idx = index
            version = results[index]["version"]
            break

    results = {
        "files": {
            fastq_file: {
                "barcode_length": get_read_length_fastq(fastq_file)
            } for fastq_file in fastq_paths
        },
        "version": version
    }

    if barcode_idx == 0:
        results["files"][sorted_by_length[0]]["type"] = TYPE_BARCODE
        results["files"][sorted_by_length[1]]["type"] = TYPE_CDNA
    elif barcode_idx == 1:
        results["files"][sorted_by_length[0]]["type"] = TYPE_CDNA
        results["files"][sorted_by_length[1]]["type"] = TYPE_BARCODE
    else:
        raise Exception(f"Barcode file is not first or second. This is weird. Files {results_barcode}")
    if files_n > 2:
        results["files"][sorted_by_length[2]]["type"] = TYPE_INDEX
    if files_n > 3:
        results["files"][sorted_by_length[3]]["type"] = TYPE_INDEX

    return results


def get_read_type_sam(sam_path: str,
                      white_list_paths: dict[int, str]):

    whitelists = dict()
    results = dict()
    for version, white_list_path in white_list_paths.items():
        whitelists[version] = set(line.strip() for line in open(white_list_path))
        results[version] = 0
    total_count = 0

    with pysam.AlignmentFile(sam_path) as file:
        iterator = enumerate(file)
        while total_count < READS_TO_CHECK:
            try:
                idx, line = next(iterator)
                if idx > READS_TO_CHECK:
                    break
                barcode = re.sub('-[0-9]*', '', line.get_tag('CB')).strip()
                total_count += 1
                for version, whitelist in whitelists.items():
                    if barcode in whitelist:
                        results[version] += 1
            except KeyError:
                continue
            except (OSError, StopIteration) as e:
                break

    results = {
        key: count / total_count
        for key, count in results.items()
    }
    stats_string = f"reads tested: {total_count}, " + \
                   ", ".join([f"v{version}: {percent:.2f}" for version, percent in results.items()])
    logging.warning(f"Stats for {sam_path}: {stats_string}")

    current_max = RATIO_THRESHOLD
    version = -1
    for ver, res in results.items():
        if res > current_max:
            current_max = res
            version = ver

    return version


def get_files_from_fastq_dump(srr_accession,
                              ncbi_dir,
                              header_dir,
                              threads,
                              white_list_paths):
    if not os.path.exists(header_dir):
        os.makedirs(header_dir)

    command = f"parallel-fastq-dump " \
              f"--sra-id {srr_accession} " \
              f"--threads {threads} " \
              f"--outdir {header_dir}/ " \
              f"--tmpdir . " \
              f"--split-files " \
              f"--gzip " \
              f"-X {READS_TO_CHECK}"

    output = subprocess.check_output(command, shell=True)
    files_n = int(subprocess.check_output(f"find '{header_dir}' | grep {srr_accession} | wc -l", shell=True))

    if files_n < 2 or files_n > 4:
        return [], -1

    files = [{
        "accession": srr_accession,
        "filename": f"{srr_accession}_{index}.fastq.gz",
        "filetype": FILETYPE_FQDUMP,
        "filesize": 0,
        "filenumber": index,
        "md5": None,
        "urltype": "",
        "url": command,
    } for index in range(1, files_n + 1)]

    headers = list()
    mapping = dict()

    for file in files:
        out_file = os.path.join(header_dir, file["filename"])
        mapping[out_file] = file
        headers.append(out_file)

    results = get_read_type_fastq(headers, white_list_paths)
    for key, res in results["files"].items():
        mapping[key]["barcode_length"] = res["barcode_length"]
        mapping[key]["type"] = res["type"]

    return files, results["version"]


def get_files_from_ftp_fastq(files,
                             header_dir,
                             white_list_paths):
    if not os.path.exists(header_dir):
        os.makedirs(header_dir)

    headers = []
    mapping = {}
    for file in files:
        url = file["url"]
        filename = file["filename"]
        out_file = os.path.join(header_dir, filename)
        output = subprocess.check_output(
            f"wget -q -O - {url} | "
            f"zcat | "
            f"head -{READS_TO_CHECK} | "
            f"gzip > {out_file}",
            shell=True)
        mapping[out_file] = file
        headers.append(out_file)

    results = get_read_type_fastq(headers, white_list_paths)
    for key, res in results["files"].items():
        mapping[key]["barcode_length"] = res["barcode_length"]
        mapping[key]["type"] = res["type"]

    return files, results["version"]


def get_files_from_ftp_bam(file: dict,
                           header_dir: str,
                           white_list_paths: dict[int, str]):
    if not os.path.exists(header_dir):
        os.makedirs(header_dir)

    url = file["url"]
    filename = file["filename"]
    tmp_file = re.sub(r'\.bam', '.sam', filename)
    out_file = os.path.join(header_dir, tmp_file)

    # header first
    output = subprocess.check_output(
        f"wget -q -O - {url} | "
        f"samtools view -H > {out_file} ",
        shell=True)

    # now alignments
    output = subprocess.check_output(
        f"wget -q -O - {url} | "
        f"samtools view --no-header |"
        f"head -{READS_TO_CHECK} >> {out_file} ",
        shell=True)

    version = get_read_type_sam(out_file, white_list_paths)

    file["barcode_length"] = get_read_length_sam(out_file)
    file["type"] = TYPE_BAM
    return [file], version


def prepare_files_10x(srr_accession,
                      files,
                      ncbi_dir,
                      header_dir,
                      threads,
                      white_list_paths):
    logging.warning(f"Checking fastq-dump files first")
    result_files, version = get_files_from_fastq_dump(srr_accession, ncbi_dir, header_dir, threads, white_list_paths)
    if version != -1:
        return result_files, FILETYPE_FQDUMP, version

    logging.error(f"For {srr_accession} - fastq-dump didn't work, will try ftp files")

    # we check ftp and if we cant use ftp - we go for sra

    ftp_files = files["ftp"]
    fastq_files = list(filter(lambda x: x["filetype"] == FILETYPE_FASTQ, ftp_files))
    bam_files = list(filter(lambda x: x["filetype"] == FILETYPE_BAM, ftp_files))

    # number of files that are ok to use is from 2 to 4
    logging.warning(f"Found {len(fastq_files)} fastq files, "
                    f"files are {', '.join([file['filename'] for file in fastq_files])}")
    if 2 <= len(fastq_files) <= 4:
        logging.warning(f"Checking {len(fastq_files)} FTP fastq files")
        result_files, version = get_files_from_ftp_fastq(fastq_files, header_dir, white_list_paths)
        if version != -1:
            return result_files, FILETYPE_FASTQ, version

    total_bam_files = len(bam_files)
    logging.warning(f"Found {total_bam_files} bam files, "
                    f"files are {', '.join([file['filename'] for file in bam_files])}")

    bam_files = list(filter(lambda x: not x["filename"].endswith(".bai"), bam_files))
    if total_bam_files > len(bam_files):
        logging.warning(f"Ignoring {total_bam_files - len(bam_files)} bam index files")

    # number of files that are ok to use is 1
    if len(bam_files) == 1:
        logging.warning(f"Checking {len(bam_files)} FTP bam file")
        result_files, version = get_files_from_ftp_bam(bam_files[0], header_dir, white_list_paths)
        if version != -1:
            return result_files, FILETYPE_BAM, version

    logging.error(f"For {srr_accession} - could not identify version")
    return [], None, -1



