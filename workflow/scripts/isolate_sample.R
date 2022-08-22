library(yaml)

sample_info <- read.csv2(snakemake@input[[1]], header=1, sep=",")
sample_info <- subset(sample_info, GSM==snakemake@params$sample)
write.table(sample_info, snakemake@output[[1]], sep=",", row.names = F, quote=F)

# SRS <- sample_info[1, "secondary_sample_accession"]
# yaml_anno <- read_yaml(snakemake@input[[2]])
# i <- which(yaml_anno$SampleIds == SRS)
# yaml_anno$Organism <- yaml_anno$Organism[i]
# yaml_anno$SampleId <- yaml_anno$SampleIds[i]
# yaml_anno$SampleName <- snakemake@params$sample
# yaml_anno$SampleIds <- NULL
# write_yaml(yaml_anno, snakemake@output[[2]])
