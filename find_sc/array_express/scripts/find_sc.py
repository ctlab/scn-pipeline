import argparse
from functions import *

import pandas as pd

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Find scRNA-seq data in GEO database")
    parser.add_argument('--start_date', type=str,
                        help='Specify the date from which you want to search for data, e.g., 2020-10-01')
    parser.add_argument('--end_date', type=str,
                        help='Specify the date after which you do not want to search for data, e.g., 2020/10')
    args = parser.parse_args()

    all_sp = pd.DataFrame(
        columns=['E-MTAB', 'Comment[ENA_SAMPLE]', 'Comment[BioSD_SAMPLE]', 'Characteristics[organism]',
                 'Comment[ENA_RUN]', 'Comment[BAM_URI]', 'Comment[FASTQ_URI]', 'technology'])
    df = extract_table(args.start_date, args.end_date)
    meta_info = parse_table(df)
    save_results(meta_info, args.start_date, args.end_date)