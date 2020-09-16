
analyze_object <- function(object, ident) {
  Idents(object) <- object[[ident]]
  out_dir <- paste0('markers/', ident)
  dir.create(out_dir, recursive = T)
  whole.markers <- FindAllMarkers(object = object,
                                  assay='SCT',
                                  only.pos = TRUE,
                                  min.pct = 0.10,
                                  test.use = 'MAST')
  write.table(whole.markers, paste(out_dir, "markers.tsv", sep = '/'), sep="\t", quote=F, row.names=F)
  top50_log_fc <- whole.markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_logFC)
  top100_log_fc <- whole.markers %>% group_by(cluster) %>% top_n(n = 100, wt = avg_logFC)
  top200_log_fc <- whole.markers %>% group_by(cluster) %>% top_n(n = 200, wt = avg_logFC)

  top50_adj_pval <- whole.markers %>% group_by(cluster) %>% top_n(n = 50, wt = p_val_adj)
  top100_adj_pval <- whole.markers %>% group_by(cluster) %>% top_n(n = 100, wt = p_val_adj)
  top200_adj_pval <- whole.markers %>% group_by(cluster) %>% top_n(n = 200, wt = p_val_adj)

  write.table(top50_log_fc, paste0(out_dir, "/top50_log_fc.tsv"), sep="\t", quote=F, row.names=F)
  write.table(top100_log_fc, paste0(out_dir, "/top100_log_fc.tsv"), sep="\t", quote=F, row.names=F)
  write.table(top200_log_fc, paste0(out_dir, "/top200_log_fc.tsv"), sep="\t", quote=F, row.names=F)
  write.table(top50_adj_pval, paste0(out_dir, "/top50_adj_pval.tsv"), sep="\t", quote=F, row.names=F)
  write.table(top100_adj_pval, paste0(out_dir, "/top100_adj_pval.tsv"), sep="\t", quote=F, row.names=F)
  write.table(top200_adj_pval, paste0(out_dir, "/top200_adj_pval.tsv"), sep="\t", quote=F, row.names=F)
}

{% if AnalysisType == "single" %}

idents <- paste0('SCT_snn_res.', c({{ Clustering.GraphBased.Resolution }}))
sapply(idents, function(ident) analyze_object(object = whole, ident = ident))


{% elif AnalysisType == "many" %}
idents <- paste0('integrated_snn_res.', c({{ Clustering.GraphBased.Resolution }}))
sapply(idents, function(ident) analyze_object(object = whole.integrated, ident = ident))

{% endif %}