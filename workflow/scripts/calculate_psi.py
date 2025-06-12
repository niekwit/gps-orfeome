import logging
import matplotlib.pyplot as plt
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
exclude_twin_peaks = snakemake.config["psi"]["exclude_twin_peaks"]
hit_th = float(snakemake.wildcards["ht"])
sd_th = float(snakemake.wildcards["st"])
pr_th = float(snakemake.wildcards["pt"])
penalty_factor = float(snakemake.wildcards["pnth"])
bc_threshold = snakemake.config["psi"]["bc_threshold"]
MAX_BIN = snakemake.config["bin_number"]
output_file_csv = snakemake.output["csv"]
output_file_rank = snakemake.output["ranked"]


def identify_twin_peaks(row, condition, cutoff):
    # Get values and bin names
    values = list(row.filter(regex=f"^{condition}_").values)
    keys = list(row.filter(regex=f"^{condition}_").keys())
    keys = [int(x.replace(f"{condition}_", "")) for x in keys]

    # Convert to dictionary
    dict_ = dict(zip(keys, values))

    # Get the highest and second highest values in dict_ and their keys
    max_val = max(dict_.values())
    max_key = [k for k, v in dict_.items() if v == max_val][0]
    dict_.pop(max_key)  # remove the max value
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


def compute_psi(row, condition):
    """
    Compute PSI values for a row and one condition.
    https://www.science.org/doi/10.1126/science.aaw4912#sec-11

    For each sample bin, divide the bin count by
    the sum of all bins for that sample.
    Multiple by bin number and sum all values for each sample.
    """
    # Check if the barcode has twin peaks
    if row["twin_peaks"]:
        return np.nan

    sob = row[f"SOB_{condition}"]
    psi_score = 0
    for i in range(1, MAX_BIN + 1):
        bin_prop = row[f"{condition}_{i}"] / sob
        psi_score += bin_prop * i
    return psi_score


# Read counts for all samples
df = pd.read_csv(counts, sep="\t")


# Order bin count columns so that they are in numerical order
# Bin count columns are assumed to be in the format f"{reference}_bin" and f"{test}_bin"
# where bin is an integer, sort first by reference and then by test
def sort_key(x):
    parts = x.split("_")
    if len(parts) > 1 and parts[1].isdigit():
        return (parts[0], int(parts[1]))
    return (x, 0)


df = df.reindex(sorted(df.columns, key=sort_key), axis=1)

# Move barcode_id, orf_id, gene to the front
df = df[
    ["barcode_id", "orf_id", "gene"]
    + [col for col in df.columns if col not in ["barcode_id", "orf_id", "gene"]]
]

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

# Raise error if no barcodes remain
if df.shape[0] == 0:
    raise ValueError(
        f"Error: No barcodes remaining after filtering for {reference} with sob_threshold {MIN_SOB_THRESHOLD}"
    )

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

# Identify whether barcode distributions have twin peaks
# i.e. two peaks with at least one bin between them
# Check this for each barcode and condition and mark as True if so
# These barcodes are excluded when calculating PSI values
if exclude_twin_peaks:
    nrows = df.shape[0]
    logging.info("  Marking twin peaked barcodes in:")
    logging.info(f"    {test}")
    df[f"twin_peaks_{test}"] = df.apply(
        lambda row: identify_twin_peaks(row, test, pr_th), axis=1
    ).reset_index(drop=True)

    logging.info(f"    {reference}")
    df[f"twin_peaks_{reference}"] = df.apply(
        lambda row: identify_twin_peaks(row, reference, pr_th), axis=1
    ).reset_index(drop=True)

    df_no_dpeaks = df[~df[f"twin_peaks_{test}"] & ~df[f"twin_peaks_{reference}"]]
    nrows_dpeaks = nrows - df_no_dpeaks.shape[0]
    logging.info(f"  Barcodes marked as having twin peaks: {nrows_dpeaks}")

    # Remove ORFs that when twin peak barcodes are removed
    # have less than bc_threshold barcodes
    df_no_dpeaks = df_no_dpeaks.copy()
    df_no_dpeaks["num_barcodes"] = df_no_dpeaks.groupby("orf_id")[
        "barcode_id"
    ].transform("count")
    orfs_to_remove = df_no_dpeaks[df_no_dpeaks["num_barcodes"] < bc_threshold]["orf_id"]
    df = df[~df["orf_id"].isin(orfs_to_remove)].reset_index(drop=True)
    nrows_dpeaks_removed = nrows - df.shape[0]
    logging.info(
        f"  ORFs removed with less than {bc_threshold} barcodes after removing  barcodes with twin peaks: {nrows_dpeaks_removed}"
    )

    # Make one column for twin peak status
    # and remove the indivitwin columns
    df["twin_peaks"] = df[f"twin_peaks_{test}"] | df[f"twin_peaks_{reference}"]
    df = df.drop(columns=[f"twin_peaks_{test}", f"twin_peaks_{reference}"])
else:
    df["twin_peaks"] = False

# Get number of "good barcodes" for each ORF, i.e. barcodes without twin peaks
df["good_barcodes"] = df.groupby("orf_id")["twin_peaks"].transform(
    lambda x: x.value_counts().get(False, 0)
)

# Remove ORFs with less than bc_threshold of good barcodes
df = df[df["good_barcodes"] >= bc_threshold].reset_index(drop=True)

### Compute PSI values:
logging.info(f"Computing PSI values for {test} vs {reference}")
for sample in [reference, test]:
    sample_bins = df.filter(regex=f"^{sample}").columns
    sample_bins = [int(x.replace(f"{sample}_", "")) for x in sample_bins]
    # Iterate over rows to compute PSI values
    df[f"PSI_{sample}"] = df.apply(
        lambda row: compute_psi(row, sample), axis=1
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

# Plot distribution of delta_PSI_mean
data = df["delta_PSI_mean"].dropna().unique()
data = data[~np.isinf(data)]

logging.info("Plotting distribution of delta_PSI_mean")
plt.hist(data, bins=150, color="#419179")
plt.axvline(data.mean(), color="red", linestyle="dashed", linewidth=1)
plt.xlabel("delta_PSI_mean")
plt.ylabel("Frequency")
plt.title("Distribution of delta_PSI_mean")
plt.savefig(snakemake.output["hist"])

logging.info("Calculating z-scores")
# Calculate robust z-score based of median and Median Absolute Deviation (MAD)
# https://en.wikipedia.org/wiki/Median_absolute_deviation
# https://en.wikipedia.org/wiki/Robust_measures_of_scale
mad = abs(df["delta_PSI_mean"] - df["delta_PSI_mean"].median())
df["z_score"] = (
    df["delta_PSI_mean"] - df["delta_PSI_mean"].median()
) / (
    mad.median() * 1.4826
)  # 1.4826 is a constant to make MAD comparable to standard deviation

logging.info("Correcting z-scores for number of barcodes")
# Correct for number of barcodes, but only if good barcode number is less than median
median = df["good_barcodes"].median()
logging.info(f"  Median number of good barcodes: {median}")
df["z_score_corr"] = df["z_score"] / np.where(
    df["good_barcodes"] < median,
    (np.sqrt(1 + ((median - df["good_barcodes"]) / penalty_factor))),
    1,  # No correction applied if good_barcodes >= median
)

logging.info("Correcting z-scores for intra ORF variability")
# Correct for delta_PSI_SD
# in rare cases this can be zero, so we need to avoid division by zero
# by adding a small value to the SD (epsilon): lowest delta_PSI_SD value of all ORFs
epsilon = min(df["delta_PSI_SD"][df["delta_PSI_SD"] > 0])
df["z_score_corr"] = np.where(
    df["delta_PSI_SD"] > 0,
    df["z_score_corr"] / df["delta_PSI_SD"],
    df["z_score_corr"] / epsilon,
)

logging.info("Correcting z-scores for deltaPSI")
# Multiply z-scores by absolute delta_PSI value
# First, scale delta_PSI_mean:
# delta_PSI that equal cutoff are scaled to 1 or above
# This will penalise ORFS with very low delta_PSI values
# ONLY DO THIS WITH DPSI < THRESHOLD?
df["delta_PSI_mean_scaled"] = abs(df["delta_PSI_mean"]) / hit_th
df["z_score_corr"] = df["z_score_corr"] * df["delta_PSI_mean_scaled"]
df = df.drop(columns=["delta_PSI_mean_scaled"])

logging.info("Scaling z-scores")
# Scale values between -128 and 128:
# Scale Negative and Positive Values Separately
col = "z_score_corr"
scaled_col = f"{col}_scaled"

df[scaled_col] = np.nan  # Initialize column with NaN

# Separate positive and negative values
pos_mask = df[col] > 0
neg_mask = df[col] < 0

# Scale positive values (from 2 to 128)
if df[pos_mask].shape[0] > 0:
    pos_max = df.loc[pos_mask, col].max()
    pos_min = df.loc[pos_mask, col].min()
    df.loc[pos_mask, scaled_col] = 2 + (
        (df.loc[pos_mask, col] - pos_min) / (pos_max - pos_min)
    ) * (128 - 2)

# Scale negative values (from -128 to -2)
if df[neg_mask].shape[0] > 0:
    neg_max = df.loc[neg_mask, col].max()
    neg_min = df.loc[neg_mask, col].min()
    df.loc[neg_mask, scaled_col] = -128 + (
        (df.loc[neg_mask, col] - neg_min) / (neg_max - neg_min)
    ) * (-2 + 128)

# Ensure zero remains zero if present
df.loc[df[col] == 0, scaled_col] = 0

# Replace original z_score_corr with scaled values
df[col] = df[scaled_col]
df = df.drop(columns=[scaled_col])

# Log2 transform z-scores, while preserving sign
# Highest/lowest value is -2/2 so log2 is safe when taking absolute value
# This transformation is done to make it easier to plot
df["z_score_corr"] = np.log2(abs(df["z_score_corr"])) * np.sign(df["z_score_corr"])

# Save to file
logging.info(f"Writing barcode-level results to {output_file_csv}")
df.to_csv(output_file_csv, index=False, na_rep="NA")

### Hit identification
logging.info("Calling hits")

# Identify ORFs that are stabilised in test condition
df[f"stabilised_in_{comparison}"] = df["delta_PSI_mean"] >= hit_th
stab = len(df[(df[f"stabilised_in_{comparison}"])]["orf_id"].unique())
logging.info(f"  Number of stabilised ORFs in {comparison}: {stab}")

# Identify ORFs that are destabilised in test condition
df[f"destabilised_in_{comparison}"] = df["delta_PSI_mean"] <= -hit_th
destab = len(df[(df[f"destabilised_in_{comparison}"])]["orf_id"].unique())
logging.info(f"  Number of destabilised ORFs in {comparison}: {destab}")

# Identify high confidence hits:
# ORFs with delta_PSI_mean >= sd_th * delta_PSI_SD &
# ORFs with delta_PSI_mean >= hit_th
df["high_confidence"] = (abs(df["delta_PSI_mean"]) >= sd_th * df["delta_PSI_SD"]) & (
    abs(df["delta_PSI_mean"]) >= hit_th
)

hc_stabilised = len(
    df[(df[f"stabilised_in_{comparison}"]) & (df["high_confidence"])]["orf_id"].unique()
)
logging.info(f"  Number of high confidence stabilised ORFs in {comparison}: {hc_stabilised}")
hc_destabilised = len(
    df[(df[f"destabilised_in_{comparison}"]) & (df["high_confidence"])]["orf_id"].unique()
)
logging.info(
    f"  Number of high confidence destabilised ORFs in {comparison}: {hc_destabilised}"
)

# Make separate columns for high confidence hits (easier for plotting)
df[f"stabilised_in_{comparison}_hc"] = df[f"stabilised_in_{comparison}"] & df["high_confidence"]
df[f"destabilised_in_{comparison}_hc"] = df[f"destabilised_in_{comparison}"] & df["high_confidence"]
df = df.drop(columns=["high_confidence"])


### Ranking of hits
logging.info("Ranking hits")

# Collapse data to ORF level
df_rank = (
    df[
        [
            "orf_id",
            "gene",
            "delta_PSI_mean",
            "delta_PSI_SD",
            "num_barcodes",
            "good_barcodes",
            f"stabilised_in_{comparison}",
            f"stabilised_in_{comparison}_hc",
            f"destabilised_in_{comparison}",
            f"destabilised_in_{comparison}_hc",
            "z_score",
            "z_score_corr",
        ]
    ]
    .drop_duplicates()
    .reset_index(drop=True)
)

# Create separate rankings for stabilised and destabilised hits
df_rank_stab = (
    df_rank[df_rank[f"stabilised_in_{comparison}"]]
    .sort_values(by="z_score_corr", ascending=False)
    .reset_index(drop=True)
)
df_rank_stab["stabilised_rank"] = df_rank_stab.index + 1

df_rank_destab = (
    df_rank[df_rank[f"destabilised_in_{comparison}"]]
    .sort_values(by="z_score_corr", ascending=True)
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
# Convert ranks to integer values
df_rank["stabilised_rank"] = df_rank["stabilised_rank"].astype("Int64")
df_rank["destabilised_rank"] = df_rank["destabilised_rank"].astype("Int64")

# Write to file
logging.info(f"Writing ranked results to {output_file_rank}")
df_rank.to_csv(output_file_rank, index=False, na_rep="NA")

logging.info("Done")
