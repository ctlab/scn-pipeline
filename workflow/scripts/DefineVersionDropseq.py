import re
import warnings
import gzip
import os
import logging
import subprocess

from Constants import *

warnings.simplefilter("ignore", UserWarning)


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


def get_read_type_fastq(fastq_paths: list[str]):
    """
    For 10x technology identify which read file is which (cdna, barcode and index) and identify chemistry version
    :type fastq_paths: list[str]
    :param fastq_paths: list of fastq files (or headers) to check
    """
    sorted_by_length = sorted(fastq_paths, key=get_read_length_fastq, reverse=True)
    version = 0

    files_n = len(fastq_paths)
    results = {
        "files": {
            fastq_file: {
                "barcode_length": get_read_length_fastq(fastq_file)
            } for fastq_file in fastq_paths
        },
        "version": version
    }

    if files_n >= 2:
        results["files"][sorted_by_length[0]]["type"] = TYPE_CDNA
        results["files"][sorted_by_length[1]]["type"] = TYPE_BARCODE
    if files_n >= 3:
        results["files"][sorted_by_length[2]]["type"] = TYPE_INDEX

    return results


def get_files_from_fastq_dump(srr_accession,
                              ncbi_dir,
                              header_dir,
                              threads):
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

    if files_n < 2 or files_n > 3:
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

    results = get_read_type_fastq(headers)
    for key, res in results["files"].items():
        mapping[key]["barcode_length"] = res["barcode_length"]
        mapping[key]["type"] = res["type"]

    return files, results["version"]


def get_files_from_ftp_fastq(files,
                             header_dir):
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

    results = get_read_type_fastq(headers)
    for key, res in results["files"].items():
        mapping[key]["barcode_length"] = res["barcode_length"]
        mapping[key]["type"] = res["type"]

    return files, results["version"]


def prepare_files_dropseq(srr_accession,
                          files,
                          ncbi_dir,
                          header_dir,
                          threads):
    logging.warning(f"Checking fastq-dump files first")
    result_files, version = get_files_from_fastq_dump(srr_accession, ncbi_dir, header_dir, threads)
    if version != -1:
        return result_files, FILETYPE_FQDUMP, version

    logging.error(f"For {srr_accession} - fastq-dump didn't work, will try ftp files")

    # we check ftp and if we cant use ftp - we go for sra

    ftp_files = files["ftp"]
    fastq_files = list(filter(lambda x: x["filetype"] == FILETYPE_FASTQ, ftp_files))

    # number of files that are ok to use is from 2 to 3
    if 2 <= len(fastq_files) <= 3:
        logging.warning(f"Checking {len(fastq_files)} FTP fastq files")
        result_files, version = get_files_from_ftp_fastq(fastq_files, header_dir)
        if version != -1:
            return result_files, FILETYPE_FASTQ, version

    logging.error(f"For {srr_accession} - could not identify version")
    return [], None, -1



