
suppressMessages(library(Seurat))
suppressMessages(library(DropletUtils))
suppressMessages(library(Matrix))
suppressMessages(library(data.table))
suppressMessages(library(ggplot2))
suppressMessages(library(magrittr))
suppressMessages(library(tidyverse))


# DEFINE FUNCTIONS

read_count_output <- function(dir, name = 'genes', tcc = FALSE) {
  dir <- normalizePath(dir, mustWork = TRUE)
  m <- readMM(paste0(dir, "/", name, ".mtx"))
  m <- t(m)
  m <- as(m, "dgCMatrix")
  ge <- if (tcc) ".ec.txt" else ".genes.txt"
  genes <- fread(paste0(dir, "/", name, ge), header = FALSE)$V2
  barcodes <- fread(paste0(dir, "/", name, ".barcodes.txt"), header = FALSE)$V1
  colnames(m) <- barcodes
  rownames(m) <- genes
  return(m)
}

createSeurat <- function(counts) {
  data <- read_count_output(counts, "genes", tcc = FALSE)
  data <- data[rowSums(data) > 0, ] # "We assume that any gene with counts of zero for all barcodes has already been filtered out, as this provides no information for distiguishing between barcodes"
  e.out <- emptyDrops(data, lower = {{ EmptyDropsLower }})
  e.out <- e.out %>% as.data.frame() %>%
    rownames_to_column(var = "cell") %>%
    select(cell, Total, LogProb, FDR)
  e.out[is.na(e.out$FDR),]$FDR <- 1
  e.out$sig <- ifelse(e.out$FDR <= 0.01, 'significant', 'not_significant')
  ggplot(data = e.out, aes(x = Total, y = -LogProb))+
    geom_point(aes(col=sig))+
    scale_color_manual(values=c("black", "red"))+
    theme_bw(base_size = 8)+
    xlab("Total UMI count")+
    ylab("-Log Probability")+
    theme(aspect.ratio = 1)
  ggsave("plots/emptyDrops.pdf")
  e.out <-
    e.out %>% select(cell, FDR) %>% 
    mutate(
      is_cell = FDR <= 0.01,
      FDR = NULL
    )
  data <- data[, e.out$is_cell]
  obj <- CreateSeuratObject(data, min.cells = 3, project = "{{ SampleId }}")
  obj
}

# CREATE COUNT MATRIX, FILTER CELLS AND FEATURES

dir.create('plots')
counts <- createSeurat('bus_out')

# SAVE COUNT MATRIX

save(counts, file='counts.RData')
