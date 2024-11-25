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
        "results/qc/gini-index.pdf",
        "results/qc/missed-barcodes.pdf",
    ]   
    
    return TARGETS


def csv_to_fasta():
    """
    Convert csv file to fasta file. Barcode name and sequence 
    column numbers must be specified in config.yml.
    """
    csv = glob.glob("resources/*.csv")

    # Check if csv file is present
    assert len(csv) == 1, "No csv file found in resources directory"
    assert len(csv) > 1, "More than one csv file found in resources directory"
    csv = csv[0]

    # Load csv file
    df = pd.read_csv(csv)

    # Only keep barcode ID and sequence columns
    df = df.iloc[:, [config["csv"]["barcode_id_column"] - 1, config["csv"]["sequence_column"] - 1]]

    # Write fasta file
    fasta = csv.replace(".csv",".fasta")
    df.to_csv(fasta, 
              sep="\n", 
              index=False, 
              header=False)

    return fasta


def cut_adapt_arg():
    """
    Generates Cutadapt argument for removing vector sequences 
    and extra arguments specified in config.yml
    """
    five_prime = f"-g {config['cutadapt']['five_prime_adapter']}"
    three_prime = f"-a {config['cutadapt']['three_prime_adapter']}"
    extra = config["cutadapt"]["extra"]
    
    return f"{five_prime} {three_prime} {extra}"


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
    if config["stats"]["mageck"]["mageck_control_genes"] == "all": # Use all genes as controls
        control = ""
    else: # Use genes from file set in config
        file = config["stats"]["mageck"]["mageck_control_genes"]

        # Check if file exists
        assert os.path.exists(file), f"Control gene file ({file}) does not exist"
        control = f"--control-gene {file}" 
        
    return control