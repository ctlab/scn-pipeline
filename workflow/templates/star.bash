
export index_dir=`dirname {{ idx }}`


{% if technology == "10X" %}
  {% if processing_mode == "fastq" or processing_mode == "fastq_dump" %}
    {% if version == 2 or version == 3 %}

STAR --genomeDir $index_dir \
  --readFilesIn {{ cdna|join(",") }} {{ barcode|join(",") }} \
  --clipAdapterType CellRanger4 --outFilterScoreMin 30 \
  --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts --soloUMIfiltering MultiGeneUMI_CR --soloUMIdedup 1MM_CR \
  --soloCBlen {% if version==2 %} 16 {% elif version==3 %} 16 {% endif %} \
  --soloUMIlen {% if version==2 %} 10 {% elif version==3 %} 12 {% endif %} \
  --soloBarcodeReadLength 0 \
  --soloType CB_UMI_Simple \
  --soloCBwhitelist {{ whitelist }} \
  --outFileNamePrefix {{ out_dir }} \
  --runThreadN {{ threads }} \
  --readFilesCommand zcat \
  --outSAMtype BAM Unsorted

    {% endif %}
  {% endif %}
{% endif %}