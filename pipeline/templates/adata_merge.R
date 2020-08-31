

library(reticulate)
loompy <- reticulate::import('loompy')

sample <- fread(paste("{{ PathToAnalysis }}", paste0(list.files(path = "{{ PathToAnalysis }}")[grepl("SRS", list.files(path = "{{ PathToAnalysis }}"))][1], '/sample_description.csv'), sep='/'))
sceasy:::seurat2anndata(whole.integrated, outFile=paste0(unique(sample$secondary_sample_accession), ".h5ad"),
                        assay="SCT", main_layer="counts", compression = "gzip")