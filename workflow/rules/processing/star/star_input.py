from snakemake.io import Wildcards
from workflow.scripts.DependencyDispatcher import DependencyDispatcher
from workflow.scripts.Constants import TECH_10X, TECH_DROPSEQ, \
    FILETYPE_FASTQ, FILETYPE_FQDUMP, FILETYPE_BAM, \
    TYPE_BARCODE, TYPE_BAM, TYPE_CDNA, TYPE_INDEX


def run_star_input(dispatcher: DependencyDispatcher):
    def bound_function(wildcards: Wildcards):
        sample = dispatcher.get_sample(wildcards)
        technology = sample.get_technology()
        processing_mode = sample.get_processing_mode()
        version = sample.get_version()

        result = {}

        if technology == TECH_10X:
            if processing_mode == FILETYPE_FASTQ or processing_mode == FILETYPE_FQDUMP:
                result[TYPE_BARCODE] = dispatcher.get_barcode_reads(wildcards)
                result[TYPE_CDNA] = dispatcher.get_cdna_reads(wildcards)
                if version == 1:
                    result[TYPE_INDEX] = dispatcher.get_index_reads(wildcards)
            if processing_mode == FILETYPE_BAM:
                result[TYPE_BAM] = dispatcher.get_bam(wildcards)

        if technology == TECH_DROPSEQ:
            if processing_mode == FILETYPE_FASTQ or processing_mode == FILETYPE_FQDUMP:
                result[TYPE_BARCODE] = dispatcher.get_barcode_reads(wildcards)
                result[TYPE_CDNA] = dispatcher.get_cdna_reads(wildcards)
            if processing_mode == FILETYPE_BAM:
                result[TYPE_BAM] = dispatcher.get_bam(wildcards)
        return result
    return bound_function
