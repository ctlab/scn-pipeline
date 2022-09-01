# Dataset processing

![Rule graph for dataset processing](images/process_rulegraph.svg)

*Above: Rule graph for complete data processing*

## Obtaining raw data

Raw data can come from two supported sources:

* From `parallel-fastq-dump`
* From FTP (fastq files and bams)

The rule is determined based on the determined meta-information obtained from 
`get_all_meta` (see [how we obtain meta information](meta.md)).

Rules to obtain raw data are: `get_fastq_dump_files` for `parallel-fastq-dump`, 
and `get_fastq_ftp`, `get_bam_ftp` for FTP files.

## STARSolo

We then generate required STAR index (`get_star_index`) using FASTA and GTF files from cellranger index
(if index is not generated yet) for desired species and generate bash script (`render_star_script`) 
which will combine  all raw files from runs that correspond to single sample.

To see which parameters are used when we generate star_script, please,
refer to `worfklow/templates/star.bash` which is Jinja2 template we use to render the script.

We then run the generated script (`run_star`) to obtain gene/cell expression matrix for downstream analysis.

We are currently using unmodified filtered output matrix from STARSolo 
as an input for seurat, but we might change this in future in `filter_counts_star` if
we decide that there are better approaches to filter empty droplets from the matrix.