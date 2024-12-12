FROM condaforge/miniforge3:latest
LABEL io.github.snakemake.containerized="true"
LABEL io.github.snakemake.conda_env_hash="2cce31bf26b8d746628308fecefbf599fd6c0c39a26c2cab2d528eed25f7b5c7"

# Step 2: Retrieve conda environments

# Conda environment:
#   source: https://github.com/snakemake/snakemake-wrappers/raw/v5.2.1/bio/bowtie2/build/environment.yaml
#   prefix: /conda-envs/4f56e195d46da93fbd8fc814617891ce
#   channels:
#     - conda-forge
#     - bioconda
#     - nodefaults
#   dependencies:
#     - bowtie2 =2.5.4
RUN mkdir -p /conda-envs/4f56e195d46da93fbd8fc814617891ce
ADD https://github.com/snakemake/snakemake-wrappers/raw/v5.2.1/bio/bowtie2/build/environment.yaml /conda-envs/4f56e195d46da93fbd8fc814617891ce/environment.yaml

# Conda environment:
#   source: https://github.com/snakemake/snakemake-wrappers/raw/v5.2.1/bio/cutadapt/se/environment.yaml
#   prefix: /conda-envs/7cffe387c7e7ef03de7e29fc51afc8f7
#   channels:
#     - conda-forge
#     - bioconda
#     - nodefaults
#   dependencies:
#     - cutadapt =4.9
RUN mkdir -p /conda-envs/7cffe387c7e7ef03de7e29fc51afc8f7
ADD https://github.com/snakemake/snakemake-wrappers/raw/v5.2.1/bio/cutadapt/se/environment.yaml /conda-envs/7cffe387c7e7ef03de7e29fc51afc8f7/environment.yaml

# Conda environment:
#   source: https://github.com/snakemake/snakemake-wrappers/raw/v5.2.1/bio/fastqc/environment.yaml
#   prefix: /conda-envs/24b8923f8e4abe077ffe95b01bfc1652
#   channels:
#     - conda-forge
#     - bioconda
#     - nodefaults
#   dependencies:
#     - fastqc =0.12.1
#     - snakemake-wrapper-utils =0.6.2
RUN mkdir -p /conda-envs/24b8923f8e4abe077ffe95b01bfc1652
ADD https://github.com/snakemake/snakemake-wrappers/raw/v5.2.1/bio/fastqc/environment.yaml /conda-envs/24b8923f8e4abe077ffe95b01bfc1652/environment.yaml

# Conda environment:
#   source: https://github.com/snakemake/snakemake-wrappers/raw/v5.2.1/bio/multiqc/environment.yaml
#   prefix: /conda-envs/9ab3161e947cc968897b9bb4302e9e45
#   channels:
#     - conda-forge
#     - bioconda
#     - nodefaults
#   dependencies:
#     - multiqc =1.25.2
#     - snakemake-wrapper-utils =0.6.2
RUN mkdir -p /conda-envs/9ab3161e947cc968897b9bb4302e9e45
ADD https://github.com/snakemake/snakemake-wrappers/raw/v5.2.1/bio/multiqc/environment.yaml /conda-envs/9ab3161e947cc968897b9bb4302e9e45/environment.yaml

# Conda environment:
#   source: workflow/envs/count.yaml
#   prefix: /conda-envs/f590ec8012fafe471c06e1cd52c78fc8
#   name: count
#   channels: 
#     - conda-forge
#     - bioconda
#     - defaults
#   dependencies:
#     - bioconda::bowtie2=2.5.4
#     - pandas=2.1.1
RUN mkdir -p /conda-envs/f590ec8012fafe471c06e1cd52c78fc8
COPY workflow/envs/count.yaml /conda-envs/f590ec8012fafe471c06e1cd52c78fc8/environment.yaml

# Conda environment:
#   source: workflow/envs/stats.yaml
#   prefix: /conda-envs/c82c98986f965bb549885eed808d81e3
#   name: stats
#   channels: 
#     - conda-forge
#     - bioconda
#     - defaults
#   dependencies:
#     - r-tidyverse=2.0.0
#     - r-cowplot=1.1.1
#     - r-viridis=0.6.4
#     - r-ggrepel=0.9.4
#     - r-ineq=0.2_13
#     - r-reshape2=1.4.4
#     - conda-forge::r-foreach=1.5.2
#     - conda-forge::r-doparallel=1.0.17
#     - python=3.10.12
#     - bioconda::mageck=0.5.9
#     - pandas=2.2.2
#     - scipy=1.14.1
RUN mkdir -p /conda-envs/c82c98986f965bb549885eed808d81e3
COPY workflow/envs/stats.yaml /conda-envs/c82c98986f965bb549885eed808d81e3/environment.yaml

# Step 3: Generate conda environments

RUN conda env create --prefix /conda-envs/4f56e195d46da93fbd8fc814617891ce --file /conda-envs/4f56e195d46da93fbd8fc814617891ce/environment.yaml && \
    conda env create --prefix /conda-envs/7cffe387c7e7ef03de7e29fc51afc8f7 --file /conda-envs/7cffe387c7e7ef03de7e29fc51afc8f7/environment.yaml && \
    conda env create --prefix /conda-envs/24b8923f8e4abe077ffe95b01bfc1652 --file /conda-envs/24b8923f8e4abe077ffe95b01bfc1652/environment.yaml && \
    conda env create --prefix /conda-envs/9ab3161e947cc968897b9bb4302e9e45 --file /conda-envs/9ab3161e947cc968897b9bb4302e9e45/environment.yaml && \
    conda env create --prefix /conda-envs/f590ec8012fafe471c06e1cd52c78fc8 --file /conda-envs/f590ec8012fafe471c06e1cd52c78fc8/environment.yaml && \
    conda env create --prefix /conda-envs/c82c98986f965bb549885eed808d81e3 --file /conda-envs/c82c98986f965bb549885eed808d81e3/environment.yaml && \
    conda clean --all -y