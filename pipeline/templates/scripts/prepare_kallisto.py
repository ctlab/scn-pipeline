import argparse
import re

import pandas as pd


def prepare_kallisto(threads: int, index: str, transcripts_to_genes: str):
    description = pd.read_csv('sample_description.csv').reset_index().to_dict("records")
    R1 = [re.sub('ftp.sra.ebi.ac.uk/vol1/fastq/', '', R1['fastq_ftp'].split(';')[0]) for R1 in description]
    R2 = [re.sub('ftp.sra.ebi.ac.uk/vol1/fastq/', '', R2['fastq_ftp'].split(';')[1]) for R2 in description]
    technology = description[0]['technology']
    with open('kallisto.sh', 'w') as out_file:
        out_file.write("mkfifo R1.gz R2.gz\n")
        out_file.write("curl -Ls ftp://ftp.sra.ebi.ac.uk/vol1/fastq/{" + ','.join(R1) + "} > R1.gz &\n")
        out_file.write("curl -Ls ftp://ftp.sra.ebi.ac.uk/vol1/fastq/{" + ','.join(R2) + "} > R2.gz &\n")
        out_file.write(
            f'kallisto bus -i {index} -x {technology} -t {threads} -o bus_out/ R1.gz R2.gz\n')
        out_file.write('mkdir bus_out/tmp\n')
        out_file.write(
            f'bustools sort -t {threads} -T bus_out/tmp/ -p bus_out/output.bus | bustools count -o bus_out/genes -g {transcripts_to_genes} -e bus_out/matrix.ec -t bus_out/transcripts.txt --genecounts -')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Prepare kallisto script")
    parser.add_argument('--threads', type=int, required=True,
                        help='Number of threads for parallel-fastq-dump program')
    parser.add_argument('--index', type=str, required=True,
                        help='Path to kallisto index')
    parser.add_argument('--transcripts_to_genes', type=str, required=True,
                        help='Path to transcripts_to_genes file')
    args = parser.parse_args()
    prepare_kallisto(args.threads, args.index, args.transcripts_to_genes)