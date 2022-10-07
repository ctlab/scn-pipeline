import re

TECH_10X = '10X'
TECH_DROPSEQ = 'DROPSEQ'
TECH_CELSEQ = 'CELSEQ'
TECH_INDROPS = 'INDROPS'
TECH_SURECELL = 'SURECELL'
TECH_SCRBSEQ = 'SCRBSEQ'
TECH_SMARTSEQ = 'SMARTSEQ'
TECH_MARSSEQ = 'MARSSEQ'
TECH_FLUIDIGM = 'FLUIDIGM'
TECH_CITESEQ = 'CITESEQ'
TECH_SLIDESEQ = 'SLIDESEQ'
TECH_SPLITSEQ = 'SPLITSEQ'
TECH_SCIRNASEQ = 'SCIRNASEQ'
TECH_VISIUM = "VISIUM"
TECH_BDRHAPSODY = "BDRHAPSODY"


ALL_PATTERNS = {
    TECH_10X: r"cell[\s-]*ranger|chromium|10x[\s-]*genomics",
    TECH_VISIUM: r"space[\s-]*ranger|visium",
    TECH_DROPSEQ: r"drop[\s-]*seq",
    TECH_CELSEQ: r"cel[\s-]*seq",
    TECH_INDROPS: r"in[\s-]*drops",
    TECH_SURECELL: r"sure[\s-]*cell",
    TECH_SCRBSEQ: r"scrb[\s-]*seq",
    TECH_SMARTSEQ: r"smart",
    TECH_MARSSEQ: r"mars",
    TECH_FLUIDIGM: r"fluidigm",
    TECH_CITESEQ: r"cite[\s-]*seq",
    TECH_SLIDESEQ: r"slide[\s]*seq",
    TECH_SPLITSEQ: r"split[\s-]*seq",
    TECH_SCIRNASEQ: r"sci[\s-]*rna[\s-]*seq",
    TECH_BDRHAPSODY: r"bd[\s-]*rhapsody"
}

ALL_TECH = {
    tech: re.compile(pattern, re.IGNORECASE) for tech, pattern in ALL_PATTERNS.items()
}

ALL_TECH_COMBINED = r"|".join(ALL_PATTERNS.values())

SAMPLE_PATTERN = re.compile(
    ALL_TECH_COMBINED,
    re.IGNORECASE)

DF_PATTERN = re.compile(
    r"single[\s-]*cell[\s-]*transcriptomic|"
    r"single[\s-]*cell[\s-]*rna[\s-]*seq|"
    r"single[\s-]*cell[\s-]*rna[\s-]*sequencing|" +
    ALL_TECH_COMBINED,
    re.IGNORECASE)

FILETYPE_FASTQ = 'fastq'
FILETYPE_BAM = 'bam'
FILETYPE_FQDUMP = 'fastq_dump'

TYPE_BAM = "bam"
TYPE_CDNA = "cdna"
TYPE_BARCODE = "barcode"
TYPE_INDEX = "index"

ALLOWED_FILETYPES = [FILETYPE_FASTQ, FILETYPE_BAM, FILETYPE_FQDUMP]
SUPPORTED_TECHS = [TECH_10X, TECH_DROPSEQ]
AWS_TYPE = 'aws'
FTP_TYPE = 'ftp'
GCP_TYPE = 'gcp'
READS_TO_CHECK = 100000
RATIO_THRESHOLD = 0.70