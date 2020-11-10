import re

start_date = re.sub('\.', '-', config['start_date'])
end_date = re.sub('\.', '-', config['end_date'])

rule all:
    input: expand("result/{start_date}_{end_date}/per_tech/cell_ranger.csv", start_date=start_date, end_date=end_date)

rule get_files:
    """
    Extract and annotate scRNA-seq experiments from  https://www.ebi.ac.uk/arrayexpress/experiments/
    """
    output: "result/{start_date}_{end_date}/per_tech/cell_ranger.csv"
    conda: "../find_sc.yml"
    log: "logs/{start_date}_{end_date}.txt"
    benchmark: "benchmarks/{start_date}_{end_date}.txt"
    params: start_date=start_date, end_date=end_date
    shell: "python scripts/find_sc.py --start_date {params.start_date} --end_date {params.end_date}"
