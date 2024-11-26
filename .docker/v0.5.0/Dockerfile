FROM condaforge/mambaforge:latest
LABEL io.github.snakemake.containerized="true"
LABEL io.github.snakemake.conda_env_hash="ce13dcdf719a1fa3ef1629bbb74b911cb8ffab3220e656fe3b9d6090bd472689"

# Step 1: Retrieve conda environments

# Conda environment:
#   source: https://github.com/snakemake/snakemake-wrappers/raw/v5.0.2/bio/bowtie2/build/environment.yaml
#   prefix: /conda-envs/84012557432533881e14141abce4094f
#   channels:
#     - conda-forge
#     - bioconda
#     - nodefaults
#   dependencies:
#     - bowtie2 =2.5.4
RUN mkdir -p /conda-envs/84012557432533881e14141abce4094f
ADD https://github.com/snakemake/snakemake-wrappers/raw/v5.0.2/bio/bowtie2/build/environment.yaml /conda-envs/84012557432533881e14141abce4094f/environment.yaml

# Conda environment:
#   source: https://github.com/snakemake/snakemake-wrappers/raw/v5.0.2/bio/cutadapt/se/environment.yaml
#   prefix: /conda-envs/afc705fc7d66595f67fba2816309d591
#   channels:
#     - conda-forge
#     - bioconda
#     - nodefaults
#   dependencies:
#     - cutadapt =4.9
RUN mkdir -p /conda-envs/afc705fc7d66595f67fba2816309d591
ADD https://github.com/snakemake/snakemake-wrappers/raw/v5.0.2/bio/cutadapt/se/environment.yaml /conda-envs/afc705fc7d66595f67fba2816309d591/environment.yaml

# Conda environment:
#   source: https://github.com/snakemake/snakemake-wrappers/raw/v5.0.2/bio/fastqc/environment.yaml
#   prefix: /conda-envs/90ea9a0a5f4de478eb6aa43932774175
#   channels:
#     - conda-forge
#     - bioconda
#     - nodefaults
#   dependencies:
#     - fastqc =0.12.1
#     - snakemake-wrapper-utils =0.6.2
RUN mkdir -p /conda-envs/90ea9a0a5f4de478eb6aa43932774175
ADD https://github.com/snakemake/snakemake-wrappers/raw/v5.0.2/bio/fastqc/environment.yaml /conda-envs/90ea9a0a5f4de478eb6aa43932774175/environment.yaml

# Conda environment:
#   source: https://github.com/snakemake/snakemake-wrappers/raw/v5.0.2/bio/multiqc/environment.yaml
#   prefix: /conda-envs/863b76ea04987efcef2c54ec3473f2c3
#   channels:
#     - conda-forge
#     - bioconda
#     - nodefaults
#   dependencies:
#     - multiqc =1.25.1
#     - snakemake-wrapper-utils =0.6.2
RUN mkdir -p /conda-envs/863b76ea04987efcef2c54ec3473f2c3
ADD https://github.com/snakemake/snakemake-wrappers/raw/v5.0.2/bio/multiqc/environment.yaml /conda-envs/863b76ea04987efcef2c54ec3473f2c3/environment.yaml

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
#   prefix: /conda-envs/6225fb395fc421ce212b701dc6c2cce9
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
#     - python=3.10.12
#     - bioconda::mageck=0.5.9
RUN mkdir -p /conda-envs/6225fb395fc421ce212b701dc6c2cce9
COPY workflow/envs/stats.yaml /conda-envs/6225fb395fc421ce212b701dc6c2cce9/environment.yaml

# Step 2: Generate conda environments

RUN mamba env create --prefix /conda-envs/84012557432533881e14141abce4094f --file /conda-envs/84012557432533881e14141abce4094f/environment.yaml && \
    mamba env create --prefix /conda-envs/afc705fc7d66595f67fba2816309d591 --file /conda-envs/afc705fc7d66595f67fba2816309d591/environment.yaml && \
    mamba env create --prefix /conda-envs/90ea9a0a5f4de478eb6aa43932774175 --file /conda-envs/90ea9a0a5f4de478eb6aa43932774175/environment.yaml && \
    mamba env create --prefix /conda-envs/863b76ea04987efcef2c54ec3473f2c3 --file /conda-envs/863b76ea04987efcef2c54ec3473f2c3/environment.yaml && \
    mamba env create --prefix /conda-envs/f590ec8012fafe471c06e1cd52c78fc8 --file /conda-envs/f590ec8012fafe471c06e1cd52c78fc8/environment.yaml && \
    mamba env create --prefix /conda-envs/6225fb395fc421ce212b701dc6c2cce9 --file /conda-envs/6225fb395fc421ce212b701dc6c2cce9/environment.yaml && \
    mamba clean --all -y