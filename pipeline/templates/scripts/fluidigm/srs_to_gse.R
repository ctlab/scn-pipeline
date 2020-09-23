suppressMessages(library(data.table))
suppressMessages(library(dplyr))
suppressMessages(library(plyr))
suppressMessages(library(tidyverse))

args <- commandArgs(TRUE)
gene_anotation_file <- args[1]
gse_file <- args[2]
gse_id <- "{{ RunName }}"
gsm_files <- args[3:length(args)]

cat(sprintf("Working on: %s \n", gse_id))
cat(sprintf("Output file: %s \n", gse_file))
cat(sprintf("Gene annotation: %s \n", gene_anotation_file))
cat(sprintf("Aggregating %i GSM \n", length(gsm_files)))

#############################FUNCTIONS############################

# load gsm file

get_gsm <- function(gsm_file) {
  gsm <-
    gsm_file %>%
    read.csv(sep = "\t",
             stringsAsFactors = F) %>%
    tibble::column_to_rownames("GENE") %>%
    select(COUNTS) %>%
    magrittr::set_colnames(stringr::str_extract(gsm_file, "SRS\\d+"))
  gsm
}


aggregate_gse <- function(gsm_files, gene_mapping) {
  gse <- do.call(cbind, unname(lapply(gsm_files, get_gsm)))
  gse
}

##################EXECUTION##########################

print("Merge GSM tables to GSE table...")

gse <- aggregate_gse(gsm_files, gene_mapping)

print("Writing expression table...")
write.table(x=gse,
            gse_file,
            sep = "\t",
            row.names = F,
            quote = F
)

cat("Done.\n")


