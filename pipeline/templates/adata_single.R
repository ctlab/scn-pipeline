

library(reticulate)
target_path <- gsub('counts.RData', '', "{{ Object }}")
sample <- fread(paste0(target_path, 'sample_description.csv'))
adata <- sceasy:::seurat2anndata(whole, outFile=paste0(unique(sample$secondary_sample_accession), ".h5ad"),
                        assay="SCT", main_layer="counts")
adata$write(paste0(unique(sample$secondary_sample_accession), ".h5ad"), compression = "gzip")
