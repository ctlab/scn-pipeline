

library(reticulate)
use_condaenv('anaconda3', required = T)
target_path <- gsub('counts.RData', '', "{{ Object }}")
sample <- fread(paste0(target_path, 'sample_description.csv'))
sceasy:::seurat2anndata(whole, outFile=paste0(unique(sample$biosd_sample), ".h5ad"),
                        assay="SCT", main_layer="counts", compression = "gzip")