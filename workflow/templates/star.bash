
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

  {% if processing_mode == "bam" %}

STAR --genomeDir $index_dir \
  --readFilesIn {{ bam|join(",") }} \
  --readFilesType SAM SE \
  --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts --soloUMIfiltering MultiGeneUMI_CR --soloUMIdedup 1MM_CR \
  --soloInputSAMattrBarcodeSeq CR UR --soloInputSAMattrBarcodeQual CY UY \
  --soloType CB_UMI_Simple \
  --soloCBwhitelist {{ whitelist }} \
  --outFileNamePrefix {{ out_dir }} \
  --runThreadN {{ threads }} \
  --readFilesCommand samtools view -F 0x100 \
  --outSAMtype BAM Unsorted

  {% endif %}

{% endif %}

{% if technology == "DROPSEQ" %}
  {% if processing_mode == "fastq" or processing_mode == "fastq_dump" %}

STAR --genomeDir $index_dir \
  --readFilesIn {{ cdna|join(",") }} {{ barcode|join(",") }} \
  --clipAdapterType CellRanger4 --outFilterScoreMin 30 \
  --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts --soloUMIfiltering MultiGeneUMI_CR --soloUMIdedup 1MM_CR \
  --soloCBstart 1 --soloCBlen 12 --soloUMIstart 13 --soloUMIlen 8 \
  --soloBarcodeReadLength 0 \
  --soloType CB_UMI_Simple \
  --soloCBwhitelist None \
  --outFileNamePrefix {{ out_dir }} \
  --runThreadN {{ threads }} \
  --readFilesCommand zcat \
  --outSAMtype BAM Unsorted

  {% endif %}
{% endif %}