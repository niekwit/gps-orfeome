# Snakemake workflow: `gps-orfeome`

[![Snakemake](https://img.shields.io/badge/snakemake-≥8.25.5-brightgreen.svg)](https://snakemake.github.io)
[![GitHub actions status](https://github.com/niekwit/gps-orfeome/workflows/Tests/badge.svg?branch=main)](https://github.com/niekwit/gps-orfeome/actions?query=branch%3Amain+workflow%3ATests)


A Snakemake workflow for `gps-orfeome screen analysis`

## Installation of required software 

Make sure you have [Conda](https://docs.conda.io/projects/conda/en/latest/index.html) installed.

First, install [Snakemake](https://snakemake.readthedocs.io/en/stable/):

```shell
$ conda create -n snakemake bioconda::snakemake=8.25.5
```

It is highly recommended to install [Apptainer](https://apptainer.org):

```shell
$ conda install conda-forge::apptainer=1.3.6
```

To setup a profile for custom `Snakemake` command line arguments, create a new profile (`config.yaml`) in `$HOME/.config/snakemake/standard/`:

```yaml
cores: 40
latency-wait: 10
use-conda: True
rerun-incomplete: True
printshellcmds: True
show-failed-logs: True
use-apptainer: True
```

## Analysis preparation

Prepare an analysis directory as follows:

```shell
.
├── config
│   ├── config.yml
│   └── stats.csv
├── reads
│   ├── BirA_1.fastq.gz
│   ├── BirA_2.fastq.gz
│   ├── BirA_3.fastq.gz
│   ├── BirA_4.fastq.gz
│   ├── BirA_5.fastq.gz
│   ├── BirA_6.fastq.gz
│   ├── LbNOX_1.fastq.gz
│   ├── LbNOX_2.fastq.gz
│   ├── LbNOX_3.fastq.gz
│   ├── LbNOX_4.fastq.gz
│   ├── LbNOX_5.fastq.gz
│   ├── LbNOX_6.fastq.gz
│   ├── TPNOX_1.fastq.gz
│   ├── TPNOX_2.fastq.gz
│   ├── TPNOX_3.fastq.gz
│   ├── TPNOX_4.fastq.gz
│   ├── TPNOX_5.fastq.gz
│   └── TPNOX_6.fastq.gz
├── resources
│   └── uORF_ORF81_plusHIF1a.csv
└── workflow
    ├── envs
    │   └── stats.yaml
    ├── report
    ├── rules
    │   ├── count.smk
    │   ├── qc.smk
    │   └── stats.smk
    ├── schemas
    │   └── config.schema.yaml
    ├── scripts
    │   ├── calculate_psi.py
    │   ├── create_count_table.py
    │   ├── csv_to_fasta.py
    │   ├── general_functions.smk
    │   ├── plot_alignment_rate.R
    │   ├── plot_barcode_profiles.R
    │   ├── plot_barcoderank.R
    │   ├── plot_coverage.R
    │   ├── plot_lfc.R
    │   ├── plot_missed_barcodes.R
    │   └── rename_to_barcode.py
    └── Snakefile

```


The analysis settings are in `config.yaml`:

```yaml
orfeome_name: uORFbarcodes

# Number of bins set during sorting of cells 
# (should be the same for each sample)
# With bin_number higher than 1, it is assumed that the user wants to perform
# a protein stability analysis using PSI (Protein Stability Index) as a metric.
# With bin_number = 1, the user wants to perform a pairwise comparison 
# of ORF counts between two conditions using MAGeCK.
bin_number: 1

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
  gene_column: 5 # Column number with gene names
  orf_column: 3 # Column number with unique ORF names
  barcode_id_column: 0 # Column with unique barcode IDs
  sequence_column: 1 # Column number with barcode sequences

# Mismatches allowed during alignment
mismatch: 0 

# MAGeCK/drugZ can be used when bin_number is set to 1
mageck:
  run: True # Run MAGeCK analysis
  extra_mageck_arguments: "--sort-criteria pos" 
  mageck_control_barcodes: all # All or file with control barcodes
  fdr: 0.25 # FDR threshold for downstream mageck analysis

# drugZ can be used when bin_number is set to 1
drugz:
  run: True # Run drugZ analysis
  extra: "" # Extra drugZ arguments

# Settings for the protein stability analysis when bin number is higher than 1
psi:
  # Minimum value of sum of barcode counts across all bins to keep
  sob_threshold: 100

  # deltaPSI thresholds for hits
  hit_threshold: 0.75

  # Barcode threshold for hits
  # Only keep ORFs with at least bc_th barcodes
  bc_th: 2

  # Correction factor for creating a random third barcode
  # for orfs with only two barcodes
  # (otherwise an SD cannot be calculated)
  correction_factor: 1.2
  
  # SD threshold for most stringent hits
  # mean Euclidian distance test & control > sd_th * SD
  sd_th: 4

  # Ranking penalty value to correct for barcode count differences
  # This is a fraction of the signal-to-noise ratio
  # SNR - (penalty * SNR) * (barcode_number_median - number_of_barcodes)
  penalty: 0.1
```

The pairwise compatisons for the calculation of deltaPSI values should be defined in `stats.csv`:

| test_condition | control_condition |
|----------------|-------------------|
| LbNOX	         |      BirA         |
| TPNOX	         |      BirA         |


The sample names in `stats.csv` should match the file names in the reads/ directory and have the extension `fastq.gz`. A sample for a specific bin should be marked with `_[0-9]`. 


## ORF library data

Information on the ORF library should be stored in a csv file located in `resources/`.

See below for an example:

| ID                    | sequence                 | IOH_ID    | Gene_ID    |
|-----------------------|--------------------------|-----------|------------|
|1_IOH10003_2802_PLD2	  | ATCCGAGTATAGAGACGTAAACTA | IOH10003	 | PLD2       |
|2_IOH10003_2802_PLD2	  | AACTACGTCATGAGCCGGATACCG | IOH10003	 | PLD2       |
|3_IOH10003_2802_PLD2	  | TTGCGCGCTGTGTTGTAACGTTAT | IOH10003	 | PLD2       |
|4_IOH10003_2802_PLD2	  | GACTAGGATGACTACGGAGTTTGC | IOH10003	 | PLD2       |
|5_IOH10003_2802_PLD2	  | GCGTCCTGTTATTCGTGATTGCGC | IOH10003	 | PLD2       |
|6_IOH10004_585_RAB22A	| ATACAGAGTAAGTTTCTCAAAATA | IOH10004	 | RAB22A     |
|7_IOH10004_585_RAB22A	| CGGAGCATCTATTACAGAAAGGTA | IOH10004	 | RAB22A     |

In `config/config.yaml` set the columns for this info as follows:

```yaml
csv: 
  # CSV file with the gene/ORF/barcode information
  # 0-indexed column numbers (First column is 0)
  gene_column: 5 # Column number with gene names
  orf_column: 3 # Column number with unique ORF names
  barcode_id_column: 0 # Column with unique barcode IDs
  sequence_column: 1 # Column number with barcode sequences
```


## Usage

To test the workflow, execute a dry-run first:

```shell
$ snakemake -np
```

To create a rule graph:

```shell
$ mkdir images
$ snakemake --forceall --rulegraph | grep -v '\-> 0\|0\[label = \"all\"' | dot -Tpng > images/rule_graph.png
```

To run the workflow:

```shell
$ snakemake --profile $HOME/.config/snakemake/standard/
```


## Citation

If you use this workflow in a paper, don't forget to give credits to the authors by citing the URL of this (original) repository and its DOI (see above).

