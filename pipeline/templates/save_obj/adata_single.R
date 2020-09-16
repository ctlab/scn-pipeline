

library(reticulate)

sceasy:::seurat2anndata(whole, outFile="{{ SampleId }}.h5ad",
                        assay="SCT", main_layer="counts")