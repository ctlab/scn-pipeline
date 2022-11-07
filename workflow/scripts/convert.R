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
  title <- snakemake@params$dataset
  description <- sprintf("%s | %s | %s", snakemake@params$dataset, study$title, study$description)
  link <- study$link
} else if (snakemake@params$level == "sample") {
  sample <- samples[samples$alias == token, ]
  title <- snakemake@params$sample
  description <- sprintf("%s | %s | %s | %s", snakemake@params$dataset, snakemake@params$sample, sample$title, sample$description)
  link <- sample$link
}

file_name <- basename(snakemake@input$markers)
markers <- as.data.frame(fread(snakemake@input$markers))
markers$cluster <- as.factor(markers$cluster)
markers <- list(markers)
names(markers) <- c(file_name)


markers <- list()
for (markers_file in snakemake@input$markers) {
  file_name <- basename(markers_file)

  table <- tryCatch({
    table <- as.data.frame(fread(markers_file))
    table$cluster <- as.factor(table$cluster)
    return(table)
  },
  error=function(cond) {
    message(paste("Error while fread:", markers_file))
    message("Here's the original error message:")
    message(cond)
    return(NULL)
  })
  markers[file_name] <- table
}



out_dir <- dirname(snakemake@output$descriptor)
seurat <- readRDS(snakemake@input$seurat)

migrateSeuratObject(seurat,
                    species=species,
                    name=title,
                    description=description,
                    link=link,
                    outdir=out_dir,
                    token=token,
                    generateMasks=F,
                    markers=markers,
                    public = T,
                    curated = F,
                    generateGMTS = T)

validateSCDataset(out_dir)