library(Seurat)
library(SCNPrep)
library(data.table)

study <- as.data.frame(fread(snakemake@input$study_meta))
samples <- as.data.frame(fread(snakemake@input$sample_meta))

title <- ""
description <- ""
link <- ""
token <- snakemake@params$token
species <- snakemake@params$species


if (snakemake@params$level == "dataset") {
  title <- study$title
  description <- study$description
  link <- study$link
} else if (snakemake@params$level == "sample") {
  sample <- samples[samples$alias == token, ]
  title <- sample$title
  description <- sample$description
  link <- sample$link
}

file_name <- basename(snakemake@input$markers)
markers <- as.data.frame(fread(snakemake@input$markers))
markers$cluster <- as.factor(markers$cluster)
markers <- list(markers)
names(markers) <- c(file_name)

out_dir <- dirname(snakemake@output$descriptor)
seurat <- readRDS(snakemake@input$seurat)

migrateSeuratObject(seurat,
                    species=species,
                    name=title,
                    description=description,
                    link=link,
                    outdir=out_dir,
                    token=token,
                    markers=markers,
                    public = T,
                    curated = F,
                    generateGMTS = T)