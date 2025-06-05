.. meta::
   :description: GPSW is a tool for analysing Global Protein Stability Profiling data.
   :keywords: Snakemake, bioinformatics, protein stability

.. |rtdBadge| image:: https://readthedocs.org/projects/gpsw/badge/?version=latest
    :target: https://gpsw.readthedocs.io/en/latest/?badge=latest
    :alt: Documentation Status
.. |testsBadge| image:: https://github.com/niekwit/gps-orfeome/actions/workflows/main.yml/badge.svg?branch=main
   :target: https://github.com/niekwit/gps-orfeome/actions/workflows/main.yml
   :alt: Tests
.. |snakemakeBadge| image:: https://img.shields.io/badge/snakemake-%E2%89%A58.25.5-brightgreen.svg
   :target: https://snakemake.github.io
   :alt: Snakemake
.. |biocondaVersionBadge| image:: https://anaconda.org/bioconda/gpsw/badges/version.svg
   :target: https://anaconda.org/bioconda/gpsw
   :alt: Bioconda version
.. |biocondaDownloadsBadge| image:: https://anaconda.org/bioconda/gpsw/badges/downloads.svg
   :target: https://anaconda.org/bioconda/gpsw
   :alt: Bioconda downloads
.. |blackStyleBadge| image:: https://img.shields.io/badge/code%20style-black-000000.svg
   :target: https://github.com/psf/black
   :alt: Code style: black
.. |githubStars| image:: https://img.shields.io/github/stars/niekwit/gpsw?style=social
    :alt: GitHub stars
.. |zenodoDOI| image:: https://zenodo.org/badge/DOI/10.5281/zenodo.15473715.svg
  :target: https://doi.org/10.5281/zenodo.15473715
  :alt: DOI

==========
GPSW
==========

|rtdBadge| |testsBadge| |snakemakeBadge| |biocondaVersionBadge| |biocondaDownloadsBadge| |blackStyleBadge| |githubStars| |zenodoDOI|

Description
===========

`GPSW` is a Python package for analysing Global Protein Stability Profiling data as described in `Koren et al. Cell 2018 <https://pubmed.ncbi.nlm.nih.gov/29779948/>`_ and `Timms et al. Science 2019 <https://pubmed.ncbi.nlm.nih.gov/31273098/>`_.

It can deal with two types of experiments:

1. **Protein stability profiling** using Protein Stability Index (PSI) as a metric, which is calculated from the proportion of reads across multiple bins.

   .. figure:: images/psi-flow-plot.png
      :alt: Sorting for PSI analysis

      Cell sort strategy for PSI analysis

2. **Pairwise comparison** of ORF counts between two conditions (or populations) using MAGeCK/DrugZ

.. toctree::
   :maxdepth: 2
   :caption: Contents:
   
   background
   user_guide
   citation
   report_issues

.. toctree::
   :hidden:

   Changelog <https://github.com/niekwit/gps-orfeome/releases>

Quick Start
===========

To install `GPSW`, you can use `conda`:


::

   conda install -c bioconda gpsw

To download the workflow code:


::

   gpsw fetch

Configure your analysis in `config/config.yaml`, place your sequencing data in `reads`, and provide ORF metadata in a csv file in `resources`. Then run the workflow with:


::

   gpsw run