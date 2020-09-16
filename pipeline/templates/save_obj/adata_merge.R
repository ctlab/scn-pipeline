

library(reticulate)

sceasy:::seurat2anndata(whole.integrated, outFile="{{ RunName }}.h5ad",
                        assay="SCT", main_layer="counts")