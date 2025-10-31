<p align="center">
  <img src="docs/_static/gpsw_logo9.png" width="200" alt="GPSW Logo" />
</p>

[![Snakemake](https://img.shields.io/badge/snakemake-≥8.25.5-brightgreen.svg)](https://snakemake.github.io)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)
[![Tests](https://github.com/niekwit/gps-orfeome/actions/workflows/main.yml/badge.svg?branch=main)](https://github.com/niekwit/gps-orfeome/actions/workflows/main.yml)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.15473715.svg)](https://doi.org/10.5281/zenodo.15473715)

[![Bioconda version](https://img.shields.io/conda/vn/bioconda/gpsw.svg)](https://anaconda.org/bioconda/gpsw)
[![Bioconda downloads](https://anaconda.org/bioconda/gpsw/badges/downloads.svg)](https://anaconda.org/bioconda/gpsw) 


# Description

`GPSW` is a tool for analysing Global Protein Stability Profiling data.

It can deal with two types of experiments:
1. **Protein stability profiling** using Protein Stability Index (PSI) as a metric, which is calculated from the proportion of reads across multiple bins. Publications that include this type of data are [Koren et al. Cell 2018](https://pubmed.ncbi.nlm.nih.gov/29779948/) and [Timms et al. Science 2019](https://pubmed.ncbi.nlm.nih.gov/31273098/).

2. **Pairwise comparison** of ORF counts between two conditions using MAGeCK/DrugZ 

# Citation

If you find this software helpful and use it in your publication, please cite it as follows:

Wit, N. GPSW [Computer software]. https://doi.org/10.5281/zenodo.15473715

# Documentation

Documentation is available [here](https://gps-orfeome.readthedocs.io/).

# Post-analysis scripts

We are also developing post-analysis scripts to further analyse the output of `GPSW` or to prepare publication quality figures. These scripts are available in a separate repository: [gps-orfeome-tools](https://github.com/niekwit/gps-orfeome-tools).
