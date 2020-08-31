setwd('/scratch/mfiruleva/winter/sc_pipeline_v2/sc_pipeline/')
library(dplyr)
library(data.table)
library(yaml)

output_dir <- '/scratch/mfiruleva/winter/10x_new/'

preparedTable <- fread('/scratch/mfiruleva/winter/only_10x_data.tsv')
multiple <- fread('/scratch/mfiruleva/winter/10x/multiple_runs.txt', header = F) %>% set_colnames('SRS')
preparedTable <- filter(preparedTable, secondary_sample_accession %in% multiple$SRS)
gse_projects <- unique(preparedTable$GSE)

### Generate yaml for script generation ###

for (gse_proj in 1:length(gse_projects)) {
  task <- list()
  gsms <- filter(preparedTable, GSE == gse_projects[gse_proj])
  task$RunName <- gsms$GSE[1]
  task$Organism <- gsms$scientific_name[1]
  # task$Objects <- paste0(unique(gsms$secondary_sample_accession), '/counts.RData')
  task$Objects <- unique(paste0(paste0(output_dir, paste0(gsms$GSE, '/')), paste0(unique(gsms$secondary_sample_accession), '/counts.RData')))
  task$SampleIds <- unique(gsms$secondary_sample_accession)
  task$PathsToCounts <- FALSE
  task$PathToAnalysis <- paste0(output_dir, gsms$GSE[1])
  write_yaml(task, paste0(paste0(output_dir, "yamls/"), paste0(gsms$GSE[1], '.yaml')))
}


### Generate sample descriptions ###
for (gse_proj in 1:length(gse_projects)) {
  srs_common <- filter(preparedTable, GSE == gse_projects[gse_proj])
  for (srs_single in srs_common$secondary_sample_accession) {
    srs <- filter(srs_common, secondary_sample_accession == srs_single)
    write.csv(srs, file=unique(paste0(paste0(paste0(paste0(output_dir, unique(srs$GSE)), '/'), paste0(srs$secondary_sample_accession, '/sample_description.csv')))), quote = F) 
  }
}
