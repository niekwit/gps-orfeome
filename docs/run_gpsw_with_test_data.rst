Run ``GPSW`` with test data
================================================================================

First activate the `GPSW` Conda environment:

.. code-block:: shell

   $ conda activate gpsw

Download workflow code and small test data set:

.. code-block:: shell

   $ cd /path/to/test/dir
   $ gpsw fetch --test-data

This will download the workflow code, configuration files and test data, and should look as follows:

.. code-block:: text

   .
   ├── config
   │   └── config.yml
   ├── reads
   │   ├── Control_1.fastq.gz
   │   ├── Control_2.fastq.gz
   │   ├── Control_3.fastq.gz
   │   ├── Control_4.fastq.gz
   │   ├── Control_5.fastq.gz
   │   ├── Control_6.fastq.gz
   │   ├── Test_1.fastq.gz
   │   ├── Test_2.fastq.gz
   │   ├── Test_3.fastq.gz
   │   ├── Test_4.fastq.gz
   │   ├── Test_5.fastq.gz
   │   └── Test_6.fastq.gz
   ├── resources
   │   └── orfs.csv
   └── workflow
       ├── envs
       │   └── stats.yaml
       ├── report
       │   ├── alignment-rates.rst
       │   ├── barcoderank.rst
       │   ├── drugz.rst
       │   ├── lfc_neg.rst
       │   ├── lfc_pos.rst
       │   ├── missed-barcodes.rst
       │   ├── multiqc.rst
       │   ├── pca.rst
       │   └── plot-coverage.rst
       ├── rules
       │   ├── count.smk
       │   ├── qc.smk
       │   └── stats.smk
       ├── schemas
       │   └── config.schema.yaml
       ├── scripts
       │   ├── calculate_proportion_of_reads_in_bins.py
       │   ├── calculate_psi.py
       │   ├── count_barcodes.sh
       │   ├── create_count_table.py
       │   ├── csv_to_fasta.py
       │   ├── general_functions.smk
       │   ├── plot_alignment_rate.R
       │   ├── plot_barcode_multi_conditions_profiles.R
       │   ├── plot_barcode_profiles.R
       │   ├── plot_barcoderank.R
       │   ├── plot_coverage.R
       │   ├── plot_dotplot.R
       │   ├── plot_lfc.R
       │   ├── plot_missed_barcodes.R
       │   ├── plot_pca.R
       │   └── rename_to_barcode.py
       └── Snakefile

   10 directories, 45 files

To start the workflow on the test data:

.. code-block:: shell

   $ cd /path/to/test/dir
   $ gpsw run --profile $HOME/.config/snakemake/standard/


.. note::
   The ``--profile`` argument will only have to be provided on the first run of `GPSW`, as it will create a config file (``~/.gpsw/config.ini``) that will store the profile path. To run `GPSW` without a profile use ``--profile None``.