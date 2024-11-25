"""
1. Normalise reads to largest data set.
2. Compute sum of bins for each sample.
3. Compute PSI values.
"""
import logging
import pandas as pd

# Set up logging
log = snakemake.log_fmt_shell(stdout=True, stderr=True)
logging.basicConfig(filename=log, level=logging.INFO)
logger = logging.getLogger(__name__)

# Load Snakemake variables
counts = snakemake.input["counts"]
MIN_SOB_THRESHOLD = snakemake.config["psi"]["sob_threshold"]
comparisons = snakemake.wildcards["comparison"]
REFERENCE_CONDITIONS = [x.split("_vs_")[1] for x in comparisons]
TEST_CONDITIONS = [x.split("_vs_")[0] for x in comparisons]
stab_th = snakemake.config["psi"]["stringent_hit_threshold"]
destab_th = snakemake.config["psi"]["stringent_negative_hit_threshold"]
sd_th = snakemake.config["psi"]["sd_th"]
OUTPUT_FILES = snakemake.output


def compute_psi(row, condition, num_bins):
    """
    Compute PSI values for a row and one condition.
    """
    sob = row[f"SOB_{condition}"]
    psi_score = 0
    for i in range(1, num_bins + 1):
        bin_prop = row[f"{condition}_{i}"] / sob
        psi_score += bin_prop * i
    return psi_score


# Read count for all samples
df_all = pd.read_csv(counts, sep="\t")

# Process data for each comparison separately
for reference, test in zip(REFERENCE_CONDITIONS, TEST_CONDITIONS):
    logger.info(f"Processing {test} vs {reference}...")
    
    # Select columns that are part of the comparison
    df = df_all.filter(regex=f"{reference}|{test}|^barcode|^gene")
    nrows = df.shape[0]
    logger.info(f" Barcodes present pre-filtering: {nrows}")
        
    # Remove barcodes where there are no counts accross all samples
    df = df[df.iloc[:,2:].sum(axis=1) > 0].reset_index(drop=True)
    nrows_zero_counts = nrows - df.shape[0]
    logger.info(f" Barcodes with no counts in any sample: {nrows_zero_counts}")

    # Identify the sample with the most reads 
    largest_sum = df.iloc[:,2:].sum().max()
    largest_sample = df.iloc[:,2:].sum().idxmax()
    logger.info(f" Largest sample: {largest_sample} with {largest_sum} reads")

    # Normalise reads to largest data set
    for col in df.columns:
        # Only columns with count data
        if df[col].dtype == "int64":
            correction_factor = largest_sum / df[col].sum()
            logger.info(f" Normalising {col} by {correction_factor}")
            df[col] = df[col].multiply(correction_factor)
            df[col] = df[col].astype(int)

    # Compute sum of bins for each sample
    sample_sums = {}
    for sample in [reference, test]:
        sample_sums[f"SOB_{sample}"] = df.filter(regex=f"^{sample}_").sum(axis=1)
    df = pd.concat([df, pd.DataFrame(sample_sums)], axis=1)  
    
   # Remove barcodes where there are low counts in the reference sample
    df = df[df[f"SOB_{reference}"] > MIN_SOB_THRESHOLD].reset_index()
    nrows_low_counts = nrows_zero_counts - df.shape[0]
    logger.info(f" Barcodes with low counts in {reference}: {nrows_low_counts} ")

    # Remove barcodes with no test count reads in any bin (to avoid division by zero)
    df = df[df.filter(regex=f"{test}_").sum(axis=1) > 0].reset_index(drop=True)
    nrows_no_test_counts = nrows_low_counts - df.shape[0]
    
    # Compute PSI values:
    # For each sample bin, divide the bin count by the sum of all bins for that sample
    # Multiple by bin number and sum all values for each sample
    logger.info(" Computing PSI values")
    for sample in [reference, test]:
        sample_bins = df.filter(regex=f"^{sample}").columns
        sample_bins = [int(x.replace(f"{sample}_", "")) for x in sample_bins]
        # It is assumed that the bins are numbered from 1 to max_bin!
        max_bin = max(sample_bins)
        # Iterate over rows to compute PSI values
        df[f"PSI_{sample}"] = df.apply(lambda row: compute_psi(row, sample, max_bin), axis=1)
    
    # Calculate mean PSI values for each ORF
    df[f"PSI_{reference}_mean"] = df.groupby("orf")[f"PSI_{reference}"].transform("mean")
    df[f"PSI_{test}_mean"] = df.groupby("orf")[f"PSI_{test}"].transform("mean")
        
    # Calculate deltaPSI for each single ORF
    df["deltaPSI"] = df[f"PSI_{test}"] - df[f"PSI_{reference}"]
        
    # Calculate SD of PSI values for each condition of each ORF
    df["delta_PSI_SD"] = df.groupby("orf")["deltaPSI"].transform("std")
    
    # Caclulate mean deltaPSI values for each ORF
    df["delta_PSI_mean"] = df.groupby("orf")["deltaPSI"].transform("mean")
    
    # Add total number of barcodes for each ORF
    df["num_barcodes"] = df.groupby("orf")["barcode"].transform("count")
    
    # Remove ORFs with only one barcode
    df = df[df["num_barcodes"] > 1].reset_index(drop=True)
    nrows_single_barcode = nrows_no_test_counts - df.shape[0]
    logger.info(f" ORFs removed that have only one barcode after filtering: {nrows_single_barcode}")
    
    # Sum of read counts for each ORF and each bin
    for sample in [reference, test]:
        # Get columns with bin counts
        sample_bins = df.filter(regex=f"^{sample}_").columns
        # Get sum of read counts for each ORF and each bin
        for bin in sample_bins:
            df[bin] = df.groupby("orf")[bin].transform("sum")
    
    # Identify ORFs that are stabilised in test condition
    df[f"stabilised_in_{test}"] = df["delta_PSI_mean"] > stab_th
    sum_ = df[["orf", f"stabilised_in_{test}"]].drop_duplicates()
    sum_ = sum_[f"stabilised_in_{test}"].sum()
    logger.info(f" Number of stabilised ORFs in {test}: {sum_}")
    
    # Identify high confidence hits for stabilised ORFs
    df[f"stabilised_in_{test}_hc"] = (df["delta_PSI_mean"] > stab_th) & (df["delta_PSI_mean"] > sd_th * df["delta_PSI_SD"])
    sum_ = df[["orf", f"stabilised_in_{test}_hc"]].drop_duplicates()
    sum_ = sum_[f"stabilised_in_{test}_hc"].sum()
    logger.info(f" Number of high confidence stabilised ORFs in {test}: {sum_}")
    
    # Identify ORFs that are destabilised in test condition
    df[f"destabilised_in_{test}"] = df["delta_PSI_mean"] < destab_th
    sum_ = df[["orf", f"destabilised_in_{test}"]].drop_duplicates()
    sum_ = sum_[f"destabilised_in_{test}"].sum()
    logger.info(f" Number of destabilised ORFs in {test}: {sum_}")
    
    # Identify high confidence hits for destabilised ORFs
    df[f"destabilised_in_{test}_hc"] = (df["delta_PSI_mean"] < destab_th) & (df["delta_PSI_mean"] > sd_th * df["delta_PSI_SD"])
    sum_ = df[["orf", f"destabilised_in_{test}_hc"]].drop_duplicates()
    sum_ = sum_[f"destabilised_in_{test}_hc"].sum()
    logger.info(f" Number of high confidence destabilised ORFs in {test}: {sum_}")
    
    # Write to file
    file = [x for x in OUTPUT_FILES if f"{test}_vs_{reference}" in x][0]
    logger.info(f" Writing to {file}")
    df.to_csv(file, index=False)

logger.info("Done.")