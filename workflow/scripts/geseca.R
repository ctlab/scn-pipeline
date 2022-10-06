library(data.table)
library(devtools)
library(msigdbr)
library(dplyr)
library(Seurat)
library(fgsea)

install_github(repo="ctlab/fgsea", ref="geseca")


data <- readRDS(snakemake@input$seurat)
data <- ProjectDim(data)

pca_loadings_var <- data@reductions$pca@feature.loadings
pca_loadings_full <- data@reductions$pca@feature.loadings.projected

pws <- msigdbr(species = snakemake@params$species, category = "H") %>%
  dplyr::distinct(gs_name, gene_symbol) %>%
  dplyr::group_by(gs_name)

pathway_names <- pws %>% group_keys() %>% pull(gs_name)
pathways_to_test <- list()
for (pathway in pathway_names) {
  pathways_to_test[[pathway]] <- pws %>% filter(gs_name==pathway) %>% pull(gene_symbol)
}

pcs <- snakemake@params$PCS
geseca_res_var <- lapply(pcs, function(k) geseca(pca_loadings_var[, 1:k], pathways_to_test, minSize = 15, maxSize = 500, eps=0))
geseca_res_var <- lapply(geseca_res_var, function(res) res[order(padj), ])
geseca_res_full <- lapply(pcs, function(k) geseca(pca_loadings_full[, 1:k], pathways_to_test, minSize = 15, maxSize = 500, eps=0))
geseca_res_full <- lapply(geseca_res_full, function(res) res[order(padj), ])

output_geseca_dir <- dirname(snakemake@output$geseca_full[1])
dir.create(output_geseca_dir, recursive = T)

for (i in 1:length(pcs)) {
  pc <- pcs[i]

  write.table(geseca_res_var[[i]],
              file.path(output_geseca_dir, paste0("geseca_var_", pc, ".tsv")),
              sep="\t",
              quote=F,
              row.names=F)

  write.table(geseca_res_full[[i]],
            file.path(output_geseca_dir, paste0("geseca_full_", pc, ".tsv")),
            sep="\t",
            quote=F,
            row.names=F)


}