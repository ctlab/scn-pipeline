suppressMessages(library(DropletUtils))
suppressMessages(library(Matrix))
suppressMessages(library(data.table))
suppressMessages(library(ggplot2))
suppressMessages(library(magrittr))
suppressMessages(library(tidyverse))
suppressMessages(library(jsonlite))
suppressMessages(library(Seurat))

# Gene Mapping

data_dir <- dirname(snakemake@input$mtx)
data <- Read10X(data_dir)
data <- data[rowSums(data) > 0, ]
saveRDS(data, file=snakemake@output$filtered_counts)


