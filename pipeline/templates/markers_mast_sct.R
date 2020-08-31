library(Seurat)
library(dplyr)

load(paste0(stringr::str_extract("{{ PathToAnalysis }}", 'GSE[0-9]*'), '.RData'))

markers <- FindAllMarkers(object = whole.integrated,
                                   assay='SCT',
                                   only.pos = TRUE,
                                   min.pct = 0.10,
                                   test.use = 'MAST')


write.table(whole.markers, "markers.tsv", sep="\t", quote=F, row.names=F)

whole.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)

top50_log_fc <- whole.markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_logFC)
top100_log_fc <- whole.markers %>% group_by(cluster) %>% top_n(n = 100, wt = avg_logFC)
top200_log_fc <- whole.markers %>% group_by(cluster) %>% top_n(n = 200, wt = avg_logFC)

top50_adj_pval <- whole.markers %>% group_by(cluster) %>% top_n(n = 50, wt = p_val_adj)
top100_adj_pval <- whole.markers %>% group_by(cluster) %>% top_n(n = 100, wt = p_val_adj)
top200_adj_pval <- whole.markers %>% group_by(cluster) %>% top_n(n = 200, wt = p_val_adj)


# 
# allClusters <- sort(unique(whole.integrated$seurat_clusters))
# allClusters <- as.numeric(levels(allClusters))[allClusters]
# for (cluster in allClusters) {
#   fileName <- sprintf("conserved_markers_%02d.tsv", cluster)
#   conserved_markers <- FindConservedMarkers(object = whole.integrated,
#                                               ident.1 = cluster,
#                                               assay = 'SCT',
#                                               grouping.var = 'sample')
#   write.table(
#     conserved_markers,
#     fileName,
#     sep = "\t",
#     quote = F,
#     col.names = NA
#   )
# }