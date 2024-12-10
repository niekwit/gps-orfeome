import datetime
import os
import glob
import sys
import re
import pandas as pd
from snakemake.utils import min_version, validate
from snakemake.logging import logger
from snakemake.shell import shell


def targets():
    TARGETS = [
        "results/qc/multiqc.html",
        "results/qc/alignment-rates.pdf",
        "results/qc/sequence-coverage.pdf",
        "results/qc/missed-barcodes.pdf",
    ]
    if config["bin_number"] == 1:
        if config["mageck"]["run"]:
            TARGETS.extend([
                expand("results/mageck/temp/{comparison}/{comparison}.gene_summary.txt", comparison=COMPARISONS),
                "results/count/barcode-counts-aggregated.tsv",
                expand("results/mageck/{comparison}/{comparison}.gene_summary.txt", comparison=COMPARISONS),
                expand("results/mageck/{comparison}/{comparison}.barcode_summary.txt", comparison=COMPARISONS),
                expand("results/mageck/{comparison}/{comparison}.normalized.txt", comparison=COMPARISONS),
                expand("results/mageck_plots/{comparison}/barcode_rank.pdf", comparison=COMPARISONS),
                expand("results/mageck_plots/{comparison}/{comparison}.lfc_{val}.pdf", comparison=COMPARISONS, val=["pos", "neg"]),
            ])
        if config["drugz"]["run"]:
            TARGETS.extend([
                expand("results/drugz/{comparison}.txt", comparison=COMPARISONS),
            ])
    else:
        TARGETS.extend([
            expand("results/psi/{comparison}_{threshold}.csv", comparison=COMPARISONS, threshold=THRESHOLD),
            expand("results/psi/{comparison}_{threshold}_ranked.csv", comparison=COMPARISONS, threshold=THRESHOLD),
            #expand("results/psi_plots/{comparison}/", comparison=COMPARISONS),
        ])
    return TARGETS


def csv():
    """
    Locate CSV file
    """
    csv = glob.glob("resources/*.csv")

    # Check if csv file is present
    assert len(csv) == 1, "Only one CSV file should be present in resources directory"
    return csv[0], csv[0].replace(".csv", ".fasta")


def cut_adapt_arg():
    """
    Generates Cutadapt argument for removing vector sequences 
    and extra arguments specified in config.yml
    """
    five_prime = f"-g {config['cutadapt']['five_prime_adapter']}"
    three_prime = f"-a {config['cutadapt']['three_prime_adapter']}"
    length = config['cutadapt']['barcode_length']
    extra = config["cutadapt"]["extra"]
    
    if length == 0:
        return f"{five_prime} {three_prime} {extra}"
    else:
        length = f"--length {length}"
        return f"{five_prime} {length} {extra}"


def sample_names():
    """
    Get sample names from fastq files and check for invalid characters
    """
    fastq = glob.glob("reads/*.fastq.gz")
    
    # Check if fastq files are present
    assert len(fastq) != 0, "No fastq files (.fastq.gz) found in reads directory"
    
    sample_names = [os.path.basename(x).replace(".fastq.gz","") for x in fastq]

    # Check for invalid characters in sample names (non-alphanumeric characters or _)
    invalid_sample_names = []
    for name in sample_names:
        if not re.search("^[a-zA-Z0-9_]*$", name):
            invalid_sample_names.append(name)
    
    if len(invalid_sample_names) > 0:
        invalid_sample_names = "\n".join(invalid_sample_names)
        print("ERROR: Only alphanumeric characters and _ are allowed in sample names")
        print(f"Invalid character(s) found in sample name(s):\n{invalid_sample_names}")
        sys.exit(1)

    return sample_names


def comparisons():
    """
    Load comparisons for MAGeCK/PSI
    """
    COMPARISONS = pd.read_csv("config/stats.csv")
    return COMPARISONS[["test_condition","control_condition"]].agg('_vs_'.join, axis=1).tolist()
   
  
def mageck_control():
    """
    Load control genes for MAGeCK
    """
    if config["mageck"]["mageck_control_barcodes"] == "all": # Use all genes as controls
        control = ""
    else: # Use genes from file set in config
        file = config["mageck"]["mageck_control_barcodes"]

        # Check if file exists
        assert os.path.exists(file), f"Control gene file ({file}) does not exist"
        control = f"--control-gene {file}" 
        
    return control