import logging
import pandas as pd

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
pr_th = float(snakemake.wildcards["pt"])
bc_threshold = snakemake.config["psi"]["bc_threshold"]
MAX_BIN = snakemake.config["bin_number"]
output_file_csv = snakemake.output["csv"]
output_file_rank = snakemake.output["ranked"]


def identify_twin_peaks(row, condition, cutoff):
    """
    Identify whether a given row contains "twin peaks" for a specified condition.

    A "twin peak" is defined as the presence of two highest values (peaks) in the row's columns
    corresponding to the given condition, where:
    - The second highest value is above a specified fraction (`cutoff`) of the highest value.
    - The positions (keys) of these two peaks differ by at least 2.

    Parameters
    ----------
    row : pandas.Series
        A row from a DataFrame containing columns named with the pattern '{condition}_<bin>'.
    condition : str
        The prefix for the columns to consider (e.g., 'sample1').
    cutoff : float
        The minimum fraction of the highest value that the second highest value must exceed
        to be considered a "twin peak".

    Returns
    -------
    bool
        True if the row contains twin peaks for the specified condition, False otherwise.
    """
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
    Calculate the Protein Stability Index (PSI) value for a given row and condition.

    The PSI value is computed by:
    1. For each bin (from 1 to MAX_BIN), dividing the count in that bin by the sum of all bins for the sample (SOB).
    2. Multiplying the resulting proportion by the bin number.
    3. Summing these values across all bins to obtain the PSI score.

    Args:
        row (pd.Series): A pandas Series containing bin counts and the sum of bins (SOB) for a sample under a specific condition.
        condition (str): The condition name used to access the relevant columns in the row.

    Returns:
        float: The computed PSI score for the given row and condition.

    Reference:
        https://www.science.org/doi/10.1126/science.aaw4912#sec-11
    """
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

# Compute sum of bins for each condition
sample_sums = {}
for sample in [reference, test]:
    sample_sums[f"SOB_{sample}"] = df.filter(regex=f"^{sample}_").sum(axis=1)
df = pd.concat([df, pd.DataFrame(sample_sums)], axis=1)

# Write csv with just the sums of bins for each condition of all orfs
# This is for plotting histogram
df[["orf_id", "gene"] + [f"SOB_{s}" for s in [reference, test]]].to_csv(
    snakemake.output["sums"], index=False
)

# Remove barcodes where there are low counts in the reference sample
nrows = df.shape[0]
df = df[df[f"SOB_{reference}"] > MIN_SOB_THRESHOLD].reset_index(drop=True)

# Do the same for the test sample
df = df[df[f"SOB_{test}"] > MIN_SOB_THRESHOLD].reset_index(drop=True)

# Raise error if no barcodes remain
if df.shape[0] == 0:
    raise ValueError(
        f"Error: No barcodes remaining after filtering for {reference} with sob_threshold {MIN_SOB_THRESHOLD}"
    )

low_counts = nrows - df.shape[0]
logging.info(
    f"  Barcodes with low counts in both reference and test condition: {low_counts}"
)

# Remove barcodes with no test count reads in any bin (to avoid division by zero)
nrows = df.shape[0]
df = df[df.filter(regex=f"{test}_").sum(axis=1) > 0].reset_index(drop=True)
nrows_no_test_counts = nrows - df.shape[0]
logging.info(f"  Barcodes with no counts for {test} in any bin: {nrows_no_test_counts}")

# Add total number of barcodes for each ORF
df["num_barcodes"] = df.groupby("orf_id")["barcode_id"].transform("count")

# Remove ORFs with less than a specified number of barcodes
norfs = df["orf_id"].nunique()
df = df[df["num_barcodes"] >= bc_threshold].reset_index(drop=True)
norfs_low_barcodes = norfs - df["orf_id"].nunique()
logging.info(
    f"  ORFs removed with less than {bc_threshold} barcodes after filtering: {norfs_low_barcodes}"
)

logging.info(f"  Number of barcodes present post-filtering: {df.shape[0]}")

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
# exclude barcodes with twin peaks (PSI value are still calculated for
# these barcodes, but should be excluded from deltaPSI calculation)
# Compute mean deltaPSI for each ORF, excluding barcodes with twin peaks
delta_psi_mean = df[df["twin_peaks"] == False].groupby("orf_id")["deltaPSI"].mean()
df["delta_PSI_mean"] = df["orf_id"].map(delta_psi_mean)

# Calculate SD of deltaPSI values for each ORF,
# excluding barcodes with twin peaks (consistent with delta_PSI_mean)
delta_psi_sd = df[df["twin_peaks"] == False].groupby("orf_id")["deltaPSI"].std()
df["delta_PSI_SD"] = df["orf_id"].map(delta_psi_sd)

# Convert normalised counts to proportions of reads in bins
# Do this in new df that will contain barcode-level results
logging.info("Calculating proportions of reads in bins")
df_barcodes = df.copy()

test_cols = [f"{test}_{i}" for i in range(1, MAX_BIN + 1)]
test_sums = df_barcodes[test_cols].sum(axis=1)  # Sum across rows for Test columns
for col in test_cols:
    df_barcodes[col] = df_barcodes[col] / test_sums
    # Handle possible division by zero
    df_barcodes[col] = df_barcodes[col].fillna(0)

# Calculate proportions for Reference columns
ref_cols = [f"{reference}_{i}" for i in range(1, MAX_BIN + 1)]
ref_sums = df_barcodes[ref_cols].sum(axis=1)
for col in ref_cols:
    df_barcodes[col] = df_barcodes[col] / ref_sums
    # Handle possible division by zero
    df_barcodes[col] = df_barcodes[col].fillna(0)

# Add column with comparison name (move to first position)
# This saves time when plotting
df_barcodes.insert(0, "Comparison", comparison)

# Save barcode level results to file
logging.info(f"Writing barcode-level results to {output_file_csv}")
df_barcodes.to_csv(output_file_csv, index=False, na_rep="NA")

### Hit identification
logging.info("Calling hits")

# Identify ORFs that are stabilised in test condition
df["stabilised"] = df["delta_PSI_mean"] >= hit_th
stab = len(df[(df["stabilised"])]["orf_id"].unique())
logging.info(f"  Number of stabilised ORFs in {comparison}: {stab}")

# Identify ORFs that are destabilised in test condition
df["destabilised"] = df["delta_PSI_mean"] <= -hit_th
destab = len(df[(df["destabilised"])]["orf_id"].unique())
logging.info(f"  Number of destabilised ORFs in {comparison}: {destab}")


### Ranking of hits
logging.info("Ranking hits")

# Collapse data to ORF level
df_rank = (
    df[
        [
            "orf_id",
            "gene",
            "good_barcodes",
            "delta_PSI_mean",
            "stabilised",
            "destabilised",
        ]
    ]
    .drop_duplicates()
    .reset_index(drop=True)
)

# Round delta_PSI_mean to 3 decimal places
df_rank["delta_PSI_mean"] = df_rank["delta_PSI_mean"].round(3)

# Create separate rankings for stabilised and destabilised hits
df_rank_stab = (
    df_rank[df_rank["stabilised"]]
    .sort_values(by="delta_PSI_mean", ascending=False)
    .reset_index(drop=True)
)
df_rank_stab["stabilised_rank"] = df_rank_stab.index + 1

df_rank_destab = (
    df_rank[df_rank["destabilised"]]
    .sort_values(by="delta_PSI_mean", ascending=True)
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
