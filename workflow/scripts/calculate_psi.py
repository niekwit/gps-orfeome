"""
1. Normalise reads to largest data set.
2. Compute sum of bins for each sample.
3. Compute PSI values.
"""

import logging
import pandas as pd
import numpy as np

# Set up logging
log = snakemake.log[0]
logging.basicConfig(
    format="%(levelname)s:%(asctime)s:%(message)s",
    level=logging.DEBUG,
    datefmt="%Y-%m-%d %H:%M:%S",
    handlers=[logging.FileHandler(log)],
)

# Load Snakemake variables
counts = snakemake.input["counts"]
MIN_SOB_THRESHOLD = snakemake.config["psi"]["sob_threshold"]
comparison = snakemake.wildcards["comparison"]
reference = comparison.split("_vs_")[1]
test = comparison.split("_vs_")[0]
exclude_dual_peaks = snakemake.config["psi"]["exclude_dual_peaks"]
hit_th = float(snakemake.wildcards["ht"])
sd_th = float(snakemake.wildcards["st"])
pr_th = float(snakemake.wildcards["pt"])
bc_threshold = snakemake.config["psi"]["bc_threshold"]
max_bin = snakemake.config["bin_number"]
penalty = float(snakemake.wildcards["p"])
output_file_csv = snakemake.output["csv"]
output_file_rank = snakemake.output["ranked"]


def identify_dual_peaks(row, condition, cutoff):
    # Get values and bin names
    values = list(row.filter(regex=f"^{condition}_").values)
    keys = list(row.filter(regex=f"^{condition}_").keys())
    keys = [int(x.replace(f"{condition}_", "")) for x in keys]

    # Convert to dictionary
    dict_ = dict(zip(keys, values))

    # Sort dictionary by key
    dict_ = {k: v for k, v in sorted(dict_.items(), key=lambda item: item[0])}

    # Get the highest and second highest values in dict_ and their keys
    max_val = max(dict_.values())
    max_key = [k for k, v in dict_.items() if v == max_val][0]
    dict_.pop(max_key)
    second_max_val = max(dict_.values())
    second_max_key = [k for k, v in dict_.items() if v == second_max_val][0]

    # Check if second_max_val is above the cutoff
    if second_max_val > max_val * cutoff:
        # Check if difference between max_key and second_max_key is at least 2
        if abs(max_key - second_max_key) >= 2:
            return True
        else:
            return False
    return False


def compute_psi(row, condition, num_bins):
    """
    Compute PSI values for a row and one condition.
    https://www.science.org/doi/10.1126/science.aaw4912#sec-11

    For each sample bin, divide the bin count by
    the sum of all bins for that sample.
    Multiple by bin number and sum all values for each sample.
    """
    # Check if the barcode has dual peaks
    if row["dual_peaks"]:
        return np.nan

    sob = row[f"SOB_{condition}"]
    psi_score = 0
    for i in range(1, num_bins + 1):
        bin_prop = row[f"{condition}_{i}"] / sob
        psi_score += bin_prop * i
    return psi_score


# Read counts for all samples
df = pd.read_csv(counts, sep="\t")

### Filtering of data
logging.info(f"Filtering data for {test} vs {reference}")

# Select columns that are part of the comparison
df = df.filter(regex=f"{reference}|{test}|^barcode_id|^orf_id|^gene")
nrows = df.shape[0]
logging.info(f"  Barcodes present pre-filtering: {nrows}")

# Remove barcodes where there are no counts accross all samples
df = df[df.iloc[:, 3:].sum(axis=1) > 0].reset_index(drop=True)
nrows_zero_counts = nrows - df.shape[0]
logging.info(f"  Barcodes with no counts in any sample: {nrows_zero_counts}")

# Identify the sample with the most reads
largest_sum = df.iloc[:, 3:].sum().max()
largest_sample = df.iloc[:, 3:].sum().idxmax()
logging.info(f"  Largest sample: {largest_sample} with {largest_sum} reads")

# Normalise reads to largest data set
for col in df.columns:
    # Only columns with count data
    if df[col].dtype == "int64":
        correction_factor = largest_sum / df[col].sum()
        logging.info(f" Normalising {col} by {correction_factor}")
        df[col] = df[col].multiply(correction_factor)
        df[col] = df[col].astype(int)

# Compute sum of bins for each sample
sample_sums = {}
for sample in [reference, test]:
    sample_sums[f"SOB_{sample}"] = df.filter(regex=f"^{sample}_").sum(axis=1)
df = pd.concat([df, pd.DataFrame(sample_sums)], axis=1)

# Remove barcodes where there are low counts in the reference sample
nrows = df.shape[0]
df = df[df[f"SOB_{reference}"] > MIN_SOB_THRESHOLD].reset_index(drop=True)
low_counts = nrows - df.shape[0]
logging.info(f"  Barcodes with low counts in {reference}: {low_counts}")

# Remove barcodes with no test count reads in any bin (to avoid division by zero)
nrows = df.shape[0]
df = df[df.filter(regex=f"{test}_").sum(axis=1) > 0].reset_index(drop=True)
nrows_no_test_counts = nrows - df.shape[0]
logging.info(f"  Barcodes with no counts for {test} in any bin: {nrows_no_test_counts}")

# Add total number of barcodes for each ORF
df["num_barcodes"] = df.groupby("orf_id")["barcode_id"].transform("count")

# Remove ORFs with only one barcode
df = df[df["num_barcodes"] > 1].reset_index(drop=True)
nrows = df.shape[0]
nrows_single_barcode = nrows - df.shape[0]  # CHECK THIS!!!
logging.info(
    f"  ORFs removed that have only one barcode after filtering: {nrows_single_barcode}"
)

# Remove ORFs with less than a specified number of barcodes
nrows = df.shape[0]
df = df[df["num_barcodes"] >= bc_threshold].reset_index(drop=True)
nrows_low_barcodes = nrows - df.shape[0]  # CHECK THIS!!!
logging.info(
    f"  ORFs removed with less than {bc_threshold} barcodes after filtering: {nrows_low_barcodes}"
)

logging.info(f"  Number of barcodes present post-filtering: {nrows}")

# Identify whether barcode distributions have dual peaks
# i.e. two peaks with at least one bin between them
# Check this for each barcode and condition and mark as True if so
# These barcodes are excluded when calculating PSI values
if exclude_dual_peaks:
    nrows = df.shape[0]
    logging.info("  Marking dual peaked barcodes in:")
    logging.info(f"    {test}")
    df[f"dual_peaks_{test}"] = df.apply(
        lambda row: identify_dual_peaks(row, test, pr_th), axis=1
    ).reset_index(drop=True)

    logging.info(f"    {reference}")
    df[f"dual_peaks_{reference}"] = df.apply(
        lambda row: identify_dual_peaks(row, reference, pr_th), axis=1
    ).reset_index(drop=True)

    df_no_dpeaks = df[~df[f"dual_peaks_{test}"] & ~df[f"dual_peaks_{reference}"]]
    nrows_dpeaks = nrows - df_no_dpeaks.shape[0]
    logging.info(f"  Barcodes marked as having dual peaks: {nrows_dpeaks}")

    # Remove ORFs that when dual peak barcodes are removed
    # have less than bc_threshold barcodes
    df_no_dpeaks = df_no_dpeaks.copy()
    df_no_dpeaks["num_barcodes"] = df_no_dpeaks.groupby("orf_id")[
        "barcode_id"
    ].transform("count")
    orfs_to_remove = df_no_dpeaks[df_no_dpeaks["num_barcodes"] < bc_threshold]["orf_id"]
    df = df[~df["orf_id"].isin(orfs_to_remove)].reset_index(drop=True)
    nrows_dpeaks_removed = nrows - df.shape[0]
    logging.info(
        f"  ORFs removed with less than {bc_threshold} barcodes after removing  barcodes with dual peaks: {nrows_dpeaks_removed}"
    )

    # Make one column for dual peak status
    # and remove the individual columns
    df["dual_peaks"] = df[f"dual_peaks_{test}"] | df[f"dual_peaks_{reference}"]
    df = df.drop(columns=[f"dual_peaks_{test}", f"dual_peaks_{reference}"])
else:
    df["dual_peaks"] = False


### Compute PSI values:
logging.info(f"Computing PSI values for {test} vs {reference}")
for sample in [reference, test]:
    sample_bins = df.filter(regex=f"^{sample}").columns
    sample_bins = [int(x.replace(f"{sample}_", "")) for x in sample_bins]
    # Iterate over rows to compute PSI values
    df[f"PSI_{sample}"] = df.apply(
        lambda row: compute_psi(row, sample, max_bin), axis=1
    ).reset_index(drop=True)

# Calculate mean PSI values for each ORF
df[f"PSI_{reference}_mean"] = df.groupby("orf_id")[f"PSI_{reference}"].transform("mean")
df[f"PSI_{test}_mean"] = df.groupby("orf_id")[f"PSI_{test}"].transform("mean")

# Calculate deltaPSI for each single ORF
df["deltaPSI"] = df[f"PSI_{test}"] - df[f"PSI_{reference}"]

# Calculate mean deltaPSI values for each ORF
df["delta_PSI_mean"] = df.groupby("orf_id")["deltaPSI"].transform("mean")

# Calculate SD of PSI values for each condition of each ORF
df["delta_PSI_SD"] = df.groupby("orf_id")["deltaPSI"].transform("std")

# IS THIS NEEDED?
# Sum of read counts for each ORF and each bin
# for sample in [reference, test]:
#    # Get columns with bin counts
#    sample_bins = df.filter(regex=f"^{sample}_").columns
#    # Get sum of read counts for each ORF and each bin
#    for bin in sample_bins:
#        df[bin] = df.groupby("orf_id")[bin].transform("sum")

### Hit identification
logging.info("Calling hits")

# Identify ORFs that are stabilised in test condition
df[f"stabilised_in_{test}"] = df["delta_PSI_mean"] > hit_th
sum_ = df[["orf_id", f"stabilised_in_{test}"]].drop_duplicates()
sum_ = sum_[f"stabilised_in_{test}"].sum()
logging.info(f"  Number of stabilised ORFs in {test}: {sum_}")

# Identify ORFs that are destabilised in test condition
df[f"destabilised_in_{test}"] = df["delta_PSI_mean"] < -hit_th
sum_ = df[["orf_id", f"destabilised_in_{test}"]].drop_duplicates()
sum_ = sum_[f"destabilised_in_{test}"].sum()
logging.info(f"  Number of destabilised ORFs in {test}: {sum_}")

# Identify high confidence hits
df["high_confidence"] = df["delta_PSI_mean"] > sd_th * df["delta_PSI_SD"]

hc_stabilised = len(df[(df[f"stabilised_in_{test}"]) & (df["high_confidence"])])
logging.info(f"  Number of high confidence stabilised ORFs in {test}: {hc_stabilised}")
hc_destabilised = len(df[(df[f"destabilised_in_{test}"]) & (df["high_confidence"])])
logging.info(
    f"  Number of high confidence destabilised ORFs in {test}: {hc_destabilised}"
)

# Make separate columns for high confidence hits (easier for plotting)
df[f"stabilised_in_{test}_hc"] = df[f"stabilised_in_{test}"] & df["high_confidence"]
df[f"destabilised_in_{test}_hc"] = df[f"destabilised_in_{test}"] & df["high_confidence"]
df = df.drop(columns=["high_confidence"])

# Save to file
logging.info(f"Writing pre-ranked results to {output_file_csv}")
df.to_csv(output_file_csv, index=False)

### Ranking of hits
logging.info("Ranking hits")

# Remove all non-hit ORFs
df_rank = df[
    (df[f"stabilised_in_{test}"]) | (df[f"destabilised_in_{test}"])
].reset_index(drop=True)

# Collapse data to ORF level
df_rank = (
    df_rank[
        [
            "orf_id",
            "gene",
            "delta_PSI_mean",
            "delta_PSI_SD",
            f"{test}_mean_distance",
            f"{test}_sd_distance",
            f"{reference}_mean_distance",
            f"{reference}_sd_distance",
            "num_barcodes",
            f"stabilised_in_{test}",
            f"stabilised_in_{test}_hc",
            f"destabilised_in_{test}",
            f"destabilised_in_{test}_hc",
        ]
    ]
    .drop_duplicates()
    .reset_index(drop=True)
)

# Calculate absolute rank based on signal-to-noise ratio
df_rank["SNR"] = abs(df_rank["delta_PSI_mean"]) / df_rank["delta_PSI_SD"]

# Move SNR column after num_barcodes
cols = list(df_rank.columns)
cols.insert(9, cols.pop(cols.index("SNR")))
df_rank = df_rank[cols]

# Correct SNR for the number of barcodes
# SNR - (penalty * SNR) * (expected_barcodes - num_barcodes)
# penalty is a value between 0-1
# Use median number of barcodes as expected number
df_rank["SNR"] = df_rank["SNR"] - (penalty * df_rank["SNR"]) * (
    df_rank["num_barcodes"].median() - df_rank["num_barcodes"]
)

df_rank = df_rank.sort_values(by="SNR", ascending=False).reset_index(drop=True)
df_rank["absolute_rank"] = df_rank.index + 1

# Create separate rankings for stabilised and destabilised hits
df_rank_stab = (
    df_rank[df_rank[f"stabilised_in_{test}"]]
    .sort_values(by="SNR", ascending=False)
    .reset_index(drop=True)
)
df_rank_stab["stabilised_rank"] = df_rank_stab.index + 1

df_rank_destab = (
    df_rank[df_rank[f"destabilised_in_{test}"]]
    .sort_values(by="SNR", ascending=False)
    .reset_index(drop=True)
)
df_rank_destab["destabilised_rank"] = df_rank_destab.index + 1

# Add these rankings to df_rank (NA for non-hits)
df_rank = pd.merge(
    df_rank, df_rank_stab[["orf_id", "stabilised_rank"]], on="orf_id", how="left"
)
df_rank = pd.merge(
    df_rank, df_rank_destab[["orf_id", "destabilised_rank"]], on="orf_id", how="left"
)

# Replace all missing values with NA
df_rank = df_rank.fillna("NA")

# Write to file
logging.info(f"Writing ranked results to {output_file_rank}")
df_rank.to_csv(output_file_rank, index=False)

logging.info("Done")
