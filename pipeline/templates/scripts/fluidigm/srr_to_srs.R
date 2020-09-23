suppressMessages(library(data.table))
suppressMessages(library(dplyr))
suppressMessages(library(plyr))
suppressMessages(library(tidyverse))

args <- commandArgs(TRUE)

cat(sprintf("Output GSM file: %s \n", args[1]))
cat(sprintf("Gene mapping: %s \n", args[2]))
cat(sprintf("Number of SRR files to aggregate: %i \n", length(args)-2))

# get gene_mapping

gene_mapping <-
  as.character(args[2]) %>%
  fread(stringsAsFactors = F,
           sep = "\t",
           header=F) %>%
  magrittr::set_colnames(c('TRANSCRIPT', 'GENE')) %>%
  mutate(GENE = gsub('.*\\|', '', GENE))


srr_list <- args[3:length(args)]

##############FUNCS#################

get_abundance <- function(file) {
  abundance <-
    file %>%
    read.table(sep = "\t",
               header = T,
               stringsAsFactors = F) %>%
    mutate(counts = round(est_counts, 0)) %>%
    select("target_id",	"counts")
  abundance
}


aggregate_gsm <- function(srr_list) {
  gsm <- lapply(srr_list, get_abundance) %>%
    bind_rows() %>%
    group_by(target_id) %>%
    summarise_all(funs(sum)) %>%
    as.data.table()
  gsm
}

collapse_transcripts <- function(gsm, gene_mapping) {
  gsm <-
    merge(gsm,
          gene_mapping,
          by.x = "target_id",
          by.y = "TRANSCRIPT",
          all = T) %>%
    select(-"target_id") %>%
    aggregate(. ~ GENE, ., FUN = sum)
  gsm
}

##################EXECUTION##################


gsm <- aggregate_gsm(srr_list)

gsm <- collapse_transcripts(gsm, gene_mapping) %>%
  magrittr::set_colnames(c("GENE", "COUNTS"))


write.table(
  gsm,
  args[1],
  row.names = F,
  quote = F,
  sep = "\t"
)
