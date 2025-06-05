Workflow settings
================================================================================

Workflow setting are in `config.yaml`:

.. code-block:: yaml

   orfeome_name: uORFbarcodes

   # Conditions in experiment
   # Sample files should be in reads/ directory
   # and have this format condition_1.fastq.gz
   # _1 represents the bin number of that condition
   conditions:
     test: [Test]
     control: [Control]

   # Number of bins set during sorting of cells
   # (should be the same for each sample)
   # With bin_number higher than 1, it is assumed that the user wants to perform
   # a protein stability analysis using PSI (Protein Stability Index) as a metric.
   # With bin_number = 1, the user wants to perform a pairwise comparison
   # of ORF counts between two conditions using MAGeCK/DrugZ.
   bin_number: 6

   cutadapt:
     # Sequence of an adapter ligated to the 5' end.
     # The adapter and any preceding bases are trimmed.
     five_prime_adapter: CCAGTAGGTCCACTATGAGT

     # Sequence of an adapter ligated to the 3' end.
     # The adapter and subsequent bases are trimmed.
     three_prime_adapter: AGCTGTGTAAGCGGAACTAG

     # Length of barcode sequence
     # This option will override 3'adapter sequence trimming if > 0
     barcode_length: 20

     # Extra cutadapt arguments
     extra: "--discard-untrimmed"

   csv:
     # CSV file with the gene/ORF/barcode information
     # 0-indexed column numbers (First column is 0)
     gene_column: 4 # Column number with gene names
     orf_column: 2 # Column number with unique ORF names
     barcode_id_column: 0 # Column with unique barcode IDs
     sequence_column: 1 # Column number with barcode sequences

   # Alignment settings
   bowtie2:
     mismatch: 0 # mismatches allowed in the alignment
     extra: "" # Extra arguments for bowtie2

   # MAGeCK/DrugZ can be used when bin_number is set to 1
   mageck:
     run: True # Run MAGeCK analysis
     extra_mageck_arguments: "--sort-criteria pos"
     mageck_control_barcodes: all # All or file with control barcodes
     fdr: 0.25 # FDR threshold for downstream mageck analysis

   # DrugZ can be used when bin_number is set to 1
   drugz:
     run: True # Run DrugZ analysis
     extra: "" # Extra DrugZ arguments

   # Settings for the protein stability analysis when bin number is higher than 1
   psi:
     # Minimum value of sum of barcode counts across all bins to keep
     sob_threshold: 100

     # deltaPSI thresholds for hits
     hit_threshold: [0.75, 1.0, 1.25]

     # Exclude barcode with twin peaks
     exclude_twin_peaks: True
     # Proportion threshold for second peak of first peak
     proportion_threshold: [0.5, 0.4, 0.35]

     # Penalty factor for having less than median number of good barcodes
     penalty_factor: [4, 4, 4]

     # Barcode threshold for hits
     # Keep ORFs with at least bc_threshold barcodes
     bc_threshold: 2

     # SD threshold for most stringent hits
     # mean deltaPSI > sd_threshold * SD
     sd_threshold: [2, 2, 2.5]

Sample names
--------------------------------------------------------------------------------
The `conditions` section defines the conditions in the experiment. The sample files should be placed in the `reads/` directory and should follow the naming convention ``<condition>_<bin_number>.fastq.gz``, where ``<condition>`` is one of the conditions defined in the `config.yaml` file (e.g. `Test_1.fastq.gz`, `Control_1.fastq.gz`, etc.).

.. code-block:: yaml

   conditions:
     test: [Test]
     control: [Control]

`GPSW` calculates :math:`\Delta PSI` values by pairing corresponding test and control conditions. For instance, if you define tests as `[Test1, Test2]` and controls as `[Control1, Control2]`, `GPSW` will compare Test1 with Control1 and Test2 with Control2. It's essential that all samples share the same `bin_number`.

Bin number
--------------------------------------------------------------------------------
If the `bin_number` is set to 1, the workflow will perform a pairwise comparison of ORF counts between two conditions using MAGeCK/DrugZ. If the `bin_number` is greater than 1, the workflow will perform a protein stability analysis using Protein Stability Index (PSI) as a metric.

.. code-block:: yaml

   bin_number: 6

Cutadapt settings
--------------------------------------------------------------------------------
The `cutadapt` section defines the settings for trimming the raw reads. The `five_prime_adapter` and `three_prime_adapter` are the sequences of the adapters ligated to the 5' and 3' ends of the reads, respectively. The `barcode_length` is the length of the barcode sequence, which will override the 3' adapter sequence trimming if set to a value greater than 0.

Extra arguments for `cutadapt` can be specified in the `extra` field. For example, ``--discard-untrimmed`` will discard reads that were not trimmed (recommended).

.. code-block:: yaml

   cutadapt:
     # Sequence of an adapter ligated to the 5' end.
     # The adapter and any preceding bases are trimmed.
     five_prime_adapter: CCAGTAGGTCCACTATGAGT

     # Sequence of an adapter ligated to the 3' end.
     # The adapter and subsequent bases are trimmed.
     three_prime_adapter: AGCTGTGTAAGCGGAACTAG

     # Length of barcode sequence
     # This option will override 3'adapter sequence trimming if > 0
     barcode_length: 20

     # Extra cutadapt arguments
     extra: "--discard-untrimmed"

ORF library information
--------------------------------------------------------------------------------
Provide a CSV file with the ORF library information in `resources/` directory. The CSV file should contain the following columns: `ID`, `sequence`, `IOH_ID`, and `Gene_ID`. See the example below:

.. list-table::
   :header-rows: 1
   :widths: 25 25 15 15

   * - ID
     - sequence
     - IOH_ID
     - Gene_ID
   * - 1_IOH10003_2802_PLD2
     - ATCCGAGTATAGAGACGTAAACTA
     - IOH10003
     - PLD2
   * - 2_IOH10003_2802_PLD2
     - AACTACGTCATGAGCCGGATACCG
     - IOH10003
     - PLD2
   * - 3_IOH10003_2802_PLD2
     - TTGCGCGCTGTGTTGTAACGTTAT
     - IOH10003
     - PLD2
   * - 4_IOH10003_2802_PLD2
     - GACTAGGATGACTACGGAGTTTGC
     - IOH10003
     - PLD2
   * - 5_IOH10003_2802_PLD2
     - GCGTCCTGTTATTCGTGATTGCGC
     - IOH10003
     - PLD2
   * - 6_IOH10004_585_RAB22A
     - ATACAGAGTAAGTTTCTCAAAATA
     - IOH10004
     - RAB22A
   * - 7_IOH10004_585_RAB22A
     - CGGAGCATCTATTACAGAAAGGTA
     - IOH10004
     - RAB22A

In `config/config.yaml` set the columns for this info as follows:

.. code-block:: yaml

   csv:
     # CSV file with the gene/ORF/barcode information
     # 0-indexed column numbers (First column is 0)
     gene_column: 3 # Column number with gene names
     orf_column: 2 # Column number with unique ORF names
     barcode_id_column: 0 # Column with unique barcode IDs
     sequence_column: 1 # Column number with barcode sequences

Alignment settings
--------------------------------------------------------------------------------
`GPSW` uses `Bowtie2` for aligning the reads to the ORF library. The `bowtie2` section defines the settings for the alignment. The `mismatch` is the number of mismatches allowed in the alignment, and `extra` can be used to specify additional arguments for `Bowtie2`.

.. code-block:: yaml

   bowtie2:
     mismatch: 0 # mismatches allowed in the alignment
     extra: "" # Extra arguments for bowtie2

MAGeCK/DrugZ settings
--------------------------------------------------------------------------------
When `bin_number` is set to 1, the workflow runs MAGeCK/DrugZ. The `mageck` section defines the settings for the MAGeCK analysis. The `run` field specifies to run MAGeCK/DrugZ analysis, and `extra_mageck_arguments` can be used to specify additional arguments for MAGeCK. The `mageck_control_barcodes` field specifies whether to use all control barcodes or a file with control barcodes. The `fdr` field specifies the FDR threshold for downstream MAGeCK analysis.

.. code-block:: yaml

   mageck:
     run: True # Run MAGeCK analysis
     extra_mageck_arguments: "--sort-criteria pos"
     mageck_control_barcodes: all # All or file with control barcodes
     fdr: 0.25 # FDR threshold for downstream mageck analysis

   drugz:
     run: True # Run DrugZ analysis
     extra: "" # Extra DrugZ arguments

PSI settings
--------------------------------------------------------------------------------
When `bin_number` is greater than 1, the workflow performs a protein stability analysis using PSI as a metric. The `psi` section defines the settings for the PSI analysis. The `sob_threshold` is the minimum value of the sum of barcode counts across all bins to keep an ORF (100 is recommended). The `hit_threshold`, `proportion_threshold`, and `penalty_factor` are lists of values that define the thresholds for hits, proportion of reads in bins, and penalty factor for having less than median number of good barcodes, respectively. The `bc_threshold` is the minimum number of 'good' barcodes required to keep an ORF (examples in Note below), and the `sd_threshold` is the SD threshold for marking high confidence hits (i.e is :math:`\Delta PSI` > n :math:`\times` SD, with n the `sd_threshold`). More on the PSI analysis can be found in the :ref:`background` section below.

.. code-block:: yaml

   psi:
     # Minimum value of sum of barcode counts across all bins to keep
     sob_threshold: 100

     # deltaPSI thresholds for hits
     hit_threshold: [0.75, 1.0, 1.25]

     # Exclude barcode with twin peaks
     exclude_twin_peaks: True
     # Proportion threshold for second peak of first peak
     proportion_threshold: [0.5, 0.4, 0.35]

     # Penalty factor for having less than median number of good barcodes
     penalty_factor: [4, 4, 4]

     # Barcode threshold for hits
     # Keep ORFs with at least bc_threshold barcodes
     bc_threshold: 2

     # SD threshold for most stringent hits
     # mean deltaPSI > sd_threshold * SD
     sd_threshold: [2, 2, 2.5]

.. note::
   Good barcodes are defined as those which do not have a twin peak in the distribution of their counts across bins. Barcodes with twin peaks are defined as having two peaks that are at least two bins apart (:math:`\Delta Bin > 1`) and the second peak has to be a minimum proportion of the highest peak. This proportion is defined by the user in the config.yaml file (`proportion_threshold`). See the example below for a visual representation of this. Not all twin peaks are marked in this example.

   .. figure:: images/twin_peak_example.png
      :alt: Twin peaks example

      Twin peaks example