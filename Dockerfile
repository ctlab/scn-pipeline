FROM condaforge/mambaforge:latest
LABEL io.github.snakemake.containerized="true"
LABEL io.github.snakemake.conda_env_hash="050534f6a96ab0b54775a3dba616ae04bc0b4b38c3158a275ac7bc99376ff429"

# Step 1: Retrieve conda environments

# Conda environment:
#   source: workflow/envs/define_technology.yaml
#   prefix: /conda-envs/a9f55bd4c09e265788f01728ff2d35a6
#   name: define_technology
#   channels:
#    - bioconda
#    - conda-forge
#   dependencies:
#    - python=3.9
#    - requests=2.27.1
#    - pysam=0.19.0
#    - numpy=1.22.3
#    - pandas=1.4.2
#    - jinja2=3.1.1
#    - entrez-direct=16.2
#    - parallel-fastq-dump=0.6.7
#    - samtools=1.15.1
RUN mkdir -p /conda-envs/define_technology
COPY workflow/envs/define_technology.yaml /conda-envs/define_technology/environment.yaml

# Conda environment:
#   source: workflow/envs/entrez_direct_utils.yaml
#   prefix: /conda-envs/4d48c8f47ae56aaa3c806c5a5a298a80
#   name: entrez_direct_utils
#   channels:
#    - bioconda
#    - conda-forge
#   dependencies:
#    - entrez-direct=16.2
#    - sra-tools=2.11.0
#    - parallel-fastq-dump=0.6.7
RUN mkdir -p /conda-envs/entrez_direct_utils
COPY workflow/envs/entrez_direct_utils.yaml /conda-envs/entrez_direct_utils/environment.yaml

# Conda environment:
#   source: workflow/envs/ffq.yaml
#   prefix: /conda-envs/b026d874309e542a75261a1d31c318d9
#   name: ffq
#   channels:
#    - bioconda
#    - conda-forge
#   dependencies:
#    - python=3.9
#    - ffq=0.3.0
RUN mkdir -p /conda-envs/ffq
COPY workflow/envs/ffq.yaml /conda-envs/ffq/environment.yaml

# Conda environment:
#   source: workflow/envs/filtering_counts.yaml
#   prefix: /conda-envs/ff80f7a8bea9de785025e5ff8d0053a0
#   name: filtering_counts
#   channels:
#    - conda-forge
#    - bioconda
#   dependencies:
#    - r-base=4.0.3
#    - r-seurat=4.1.1
#    - r-cairo=1.5_15
#    - r-matrix=1.4_1
#    - r-data.table=1.13.2
#    - r-ggplot2=3.3.2
#    - r-magrittr=2.0.1
#    - r-tidyverse=1.3.0
#    - r-biocmanager=1.30.10
#    - r-jsonlite=1.7.1
#    - bioconductor-dropletutils=1.10.0
RUN mkdir -p /conda-envs/filtering_counts
COPY workflow/envs/filtering_counts.yaml /conda-envs/filtering_counts/environment.yaml

# Conda environment:
#   source: workflow/envs/git.yaml
#   prefix: /conda-envs/6076c324ebc0e071433e9952599f3cb8
#   name: git
#
#   channels:
#     - anaconda
#
#   dependencies:
#     - anaconda::git=2.34.1
RUN mkdir -p /conda-envs/git
COPY workflow/envs/git.yaml /conda-envs/git/environment.yaml

# Conda environment:
#   source: workflow/envs/jinja.yaml
#   prefix: /conda-envs/1b3b6e79602055a6eea6e85b52bf5657
#   name: python-render-templates
#   channels:
#    - bioconda
#    - conda-forge
#   dependencies:
#    - python=3.9
#    - pandas=1.4.2
#    - jinja2=3.1.1
RUN mkdir -p /conda-envs/jinja
COPY workflow/envs/jinja.yaml /conda-envs/jinja/environment.yaml

# Conda environment:
#   source: workflow/envs/jq.yaml
#   prefix: /conda-envs/39fed979d0f9e53a7cc1510ae8428325
#   name: jq
#   channels:
#    - conda-forge
#   dependencies:
#    - python=3.9
#    - conda-forge::jq=1.6
RUN mkdir -p /conda-envs/jq
COPY workflow/envs/jq.yaml /conda-envs/jq/environment.yaml

# Conda environment:
#   source: workflow/envs/seurat_analysis.yaml
#   prefix: /conda-envs/ad273055aa94c52656f49bcd0d0665a3
#   name: seurat_analysis
#   channels:
#    - conda-forge
#   dependencies:
#    - fit-sne=1.2.1
#    - r-base=4.0.3
#    - r-cairo=1.5_15
#    - r-remotes=2.2.0
#    - r-ggplot2=3.3.2
#    - r-flexmix=2.3_17
#    - r-biocmanager=1.30.10
#    - r-matrix=1.4_1
#    - r-seurat=4.1.1
#    - r-jsonlite=1.7.1
#    - bioconda::bioconductor-glmgampoi=1.2.0
#    - bioconda::bioconductor-limma=3.46.0
#    - bioturing::r-seuratwrappers=0.3.0
RUN mkdir -p /conda-envs/seurat_analysis
COPY workflow/envs/seurat_analysis.yaml /conda-envs/seurat_analysis/environment.yaml

# Conda environment:
#   source: workflow/envs/star.yaml
#   prefix: /conda-envs/d88546f102dfebb502ed7945dcac76f6
#   name: star
#   channels:
#     - bioconda
#     - conda-forge
#   dependencies:
#     - star=2.7.10a
#     - samtools=1.15.1
RUN mkdir -p /conda-envs/star
COPY workflow/envs/star.yaml /conda-envs/star/environment.yaml

# Step 2 Installing Snakemake into base

RUN mamba install -c conda-forge -c bioconda -n base snakemake

# Step 3: Generate conda environments

RUN mamba env create -n define_technology --file /conda-envs/define_technology/environment.yaml && \
    mamba env create -n entrez_direct_utils --file /conda-envs/entrez_direct_utils/environment.yaml && \
    mamba env create -n ffq --file /conda-envs/ffq/environment.yaml && \
    mamba env create -n filtering_counts --file /conda-envs/filtering_counts/environment.yaml && \
    mamba env create -n git --file /conda-envs/git/environment.yaml && \
    mamba env create -n jinja --file /conda-envs/jinja/environment.yaml && \
    mamba env create -n jq --file /conda-envs/jq/environment.yaml && \
    mamba env create -n seurat_analysis --file /conda-envs/seurat_analysis/environment.yaml && \
    mamba env create -n star --file /conda-envs/star/environment.yaml && \
    mamba clean --all -y