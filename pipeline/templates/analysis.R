suppressMessages(library(data.table))
suppressMessages(library(dplyr))
suppressMessages(library(functools))
suppressMessages(library(ggplot2))
suppressMessages(library(gridExtra))
suppressMessages(library(Matrix))
suppressMessages(library(magrittr))
suppressMessages(library(sctransform))
suppressMessages(library(Seurat))
suppressMessages(library(reticulate))

setwd("{{ AnalysisFolder }}")
set.seed(1)

## FUNCTIONS
{% include 'create_obj/add_metadata.R' %}
{% include 'visualization/draw_plots.R' %}
{% if FilterMito %}
{% include 'filtering/get_conf_interval.R' %}
{% endif %}
{%  if FilterMito and AnalysisType == "many"%}
{% include 'filtering/filter_mito_func.R' %}
{% endif %}

{% if FilterUMI and not WholeUMI %}
{% include 'filtering/peakfinder.R' %}
{% endif %}
{%  if FilterUMI and not WholeUMI and AnalysisType == "many"%}
{% include 'filtering/filter_umi_func.R' %}
{% endif %}


## GATHERING DATA TOGETHER
{% if AnalysisType == "single" %}
    {% include 'create_obj/intro_single.R' %}
{% elif AnalysisType == "many" %}
    {% include 'create_obj/intro_merge.R' %}
{% endif %}


## CREATE DIRECTORY FOR PLOTS
{% include 'visualization/plot_dir_creation.R' %}


## Number of cells before
{% include 'create_obj/cells_before.R' %}


## FILTER MT CONTENT
{% if FilterMito and not test_mode %}
{% include 'filtering/mt_filtering.R' %}
{% endif %}


{% if not FilterUMI or WholeUMI %}
## NORMALIZATION
{% include 'normalization/normalization.R' %}


{% if AnalysisType == "many" %}
  {% include 'integration/integration.R' %}
{% endif %}


## PCA
{% include 'dim_reduction/pca.R' %}


## TSNE
{% include 'dim_reduction/tsne.R' %}


## UMAP
{% include 'dim_reduction/umap.R' %}


## CLUSTERING
{% include 'clustering/clustering.R' %}


## VISUALIZATION
{% include 'visualization/visualization.R' %}


## AVERAGING
{% include 'clustering/average.R' %}


## FINDING ANS SAVING MARKERS
{% include 'clustering/markers.R' %}


## SAVING
{% include 'save_obj/save.R' %}
{% if AnalysisType == "single" %}
  {% include 'save_obj/adata_single.R' %}
{% elif AnalysisType == "many" %}
  {% include 'save_obj/adata_merge.R' %}
{% endif %}


## Number of cells after
{% include 'save_obj/cells_after.R' %}



{% else %}
## FILTER nUMI
{% include 'filtering/umi_filtering.R' %}


## NORMALIZATION: FILTERED DATASET
{% include 'normalization/normalization.R' %}

{% if AnalysisType == "many" %}
  {% include 'integration/integration.R' %}
{% endif %}

## PCA: FILTERED DATASET
{% include 'dim_reduction/pca.R' %}


## TSNE: FILTERED DATASET
{% include 'dim_reduction/tsne.R' %}


## UMAP: FILTERED DATASET
{% include 'dim_reduction/umap.R' %}


## CLUSTERING: FILTERED DATASET
{% include 'clustering/clustering.R' %}


## VISUALIZATION: FILTERED DATASET
{% include 'visualization/visualization.R' %}


## AVERAGING: FILTERED DATASET
{% include 'clustering/average.R' %}


## FINDING MARKERS: FILTERED DATASET
{% include 'clustering/markers.R' %}


## SAVING: FILTERED DATASET
{% include 'save_obj/save.R' %}

{% if AnalysisType == "single" %}
  {% include 'save_obj/adata_single.R' %}
{% elif AnalysisType == "many" %}
  {% include 'save_obj/adata_merge.R' %}
{% endif %}


## Number of cells after
{% include 'save_obj/cells_after.R' %}
{% endif %}

