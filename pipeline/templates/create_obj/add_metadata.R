

add_metadata <- function(data) {
  mito.genes <-
    grep(pattern = "^Mt\\.|^MT\\.|^mt\\.|^Mt-|^MT-|^mt-",
         x = rownames(x = GetAssayData(object = data)),
         value = TRUE)
  percent.mito <-
    Matrix::colSums(GetAssayData(object = data, slot = "counts")[mito.genes, ]) /
    Matrix::colSums(GetAssayData(object = data, slot = "counts"))
  data[['percent.mito']] <- percent.mito
  data[['percent.mito_log10']] <- log10(data[['percent.mito']] + 1)
  data[['nCount_RNA_log10']] <- log10(data[['nCount_RNA']] + 1)
  data[['nFeature_RNA_log10']] <- log10(data[['nFeature_RNA']] + 1)
  data[['nCount_RNA_log2']] <- log2(data[['nCount_RNA']] + 1)
  data[['nFeature_RNA_log2']] <- log2(data[['nFeature_RNA']] + 1)
  data[['scaled_mito']] <- scale(percent.mito)
  data[['scaled_nCount_RNA']] <- scale(data[['nCount_RNA_log10']])
  attr(data$scaled_nCount_RNA, "scaled:center") <- NULL
  attr(data$scaled_nCount_RNA, "scaled:scale") <- NULL
  attr(data$scaled_mito, "scaled:center") <- NULL
  attr(data$scaled_mito, "scaled:scale") <- NULL
  data
}
