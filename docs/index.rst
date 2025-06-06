.. meta::
   :description: GPSW is a tool for analysing Global Protein Stability Profiling data.
   :keywords: Snakemake, bioinformatics, protein stability


==================
GPSW documentation
==================

**Version**: |version|

Description
===========

:mod:`GPSW` is a Python package for analysing Global Protein Stability Profiling data as described in `Koren et al. Cell 2018 <https://pubmed.ncbi.nlm.nih.gov/29779948/>`_ and `Timms et al. Science 2019 <https://pubmed.ncbi.nlm.nih.gov/31273098/>`_. It employs the `Snakemake <https://snakemake.readthedocs.io/>`_ workflow management system to provide a reproducible and scalable analysis workflow.

It can deal with two types of experiments:

1. **Protein stability profiling** using Protein Stability Index (PSI) as a metric, which is calculated from the proportion of reads across multiple bins.

2. **Pairwise comparison** of ORF counts between two conditions (or populations) using MAGeCK/DrugZ.


.. grid:: 1 2 2 2
   :gutter: 4
   :padding: 2 2 0 0
   :class-container: sd-text-center

   .. grid-item-card:: Background
      :class-card: intro-card
      :shadow: md

      The background section provides an overview of the Global Protein Stability Profiling (GPS) method and the output of the *GPSW* workflow.

      +++

      .. button-ref:: background
         :ref-type: ref
         :click-parent:
         :color: secondary
         :expand:

         To the background guides

   .. grid-item-card::  User guide
      :class-card: intro-card
      :shadow: md

      The user guide provides in-depth information on how to install and use *GPSW*.

      +++

      .. button-ref:: user_guide
         :ref-type: ref
         :click-parent:
         :color: secondary
         :expand:

         To the user guide

   .. grid-item-card::  About GPSW
      :class-card: intro-card
      :shadow: md

      Information about the GPSW project, its authors, and how to cite it.

      +++

      .. button-ref:: about
         :ref-type: ref
         :click-parent:
         :color: secondary
         :expand:

         To the about section

   .. grid-item-card::  Report issues
      :class-card: intro-card
      :shadow: md

      Saw a typo in the documentation? Want to improve
      existing functionalities? The contributing guidelines will guide
      you through the process of improving GPSW.

      +++

      .. button-ref:: report_issues
         :ref-type: ref
         :click-parent:
         :color: secondary
         :expand:

         To the issues guide


.. toctree::
   :maxdepth: 2
   :hidden:
   :titlesonly:
   
   background
   user_guide
   about
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