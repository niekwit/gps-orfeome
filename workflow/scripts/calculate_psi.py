"""
1. Normalise reads to largest data set.
2. Compute sum of bins for each sample.
3. Compute PSI values.
"""
import logging
import random
from itertools import combinations
import pandas as pd
import numpy as np

# Set up logging
log = snakemake.log[0]
logging.basicConfig(format='%(levelname)s:%(asctime)s:%(message)s', 
                    level=logging.DEBUG,
                    datefmt='%Y-%m-%d %H:%M:%S',
                    handlers=[logging.FileHandler(log)])

# Load Snakemake variables
counts = snakemake.input["counts"]
MIN_SOB_THRESHOLD = snakemake.config["psi"]["sob_threshold"]
comparison = snakemake.wildcards["comparison"]
reference = comparison.split("_vs_")[1] 
test = comparison.split("_vs_")[0]
hit_th = snakemake.config["psi"]["hit_threshold"]
sd_th = snakemake.config["psi"]["sd_th"]
bc_th = snakemake.config["psi"]["bc_th"]
bin_number = snakemake.config["bin_number"]
cf = snakemake.config["psi"]["correction_factor"]
penalty = snakemake.config["psi"]["penalty"]
output_file_csv = snakemake.output["csv"]
output_file_rank = snakemake.output["ranked"]


def compute_psi(row, condition, num_bins):
    """
    Compute PSI values for a row and one condition.
    https://www.science.org/doi/10.1126/science.aaw4912#sec-11
    """
    sob = row[f"SOB_{condition}"]
    psi_score = 0
    for i in range(1, num_bins + 1):
        bin_prop = row[f"{condition}_{i}"] / sob
        psi_score += bin_prop * i
    return psi_score


def compute_euclidean_distance(curves_):
    """
    Calculate the pairwise Euclidean distances between n line curves.
    """
    def euclidean_distance(curve1, curve2):
        # https://en.wikipedia.org/wiki/Euclidean_distance
        return np.sqrt(np.sum((np.array(curve1) - np.array(curve2))**2))
    
    # For lists with just two curves, add a third data set
    # This is picked randomly from the two existing curves
    # But this value is multiplied by a specified factor
    # This is to punish orfs with just two barcodes
    # And to avoid not being able to calculate the mean/SD with just one distance
    length = len(curves_)
    if length == 2:
        seeds = range(0, 25)[0:bin_number]
        new_curve = []
        for i in seeds:
            random.seed(i)
            data_set = curves_[random.randint(0, 1)]
            new_curve.append(data_set[i] * cf)
        curves_.append(new_curve)       
        
    # Calculate the pairwise distances
    distances = {}
    for i, j in combinations(range(len(curves_)), 2):
        distances[(i, j)] = euclidean_distance(curves_[i], curves_[j])

    # Calculate the mean distance
    mean = np.mean(list(distances.values()))
    
    # Calculate the standard deviation of the distances
    std = np.std(list(distances.values()))
    
    return length * [[mean, std]]


def stringent_hit(row, test, reference, th):
    """
    Determine if a hit is significant based on the mean and SD of two conditions.
    """
    mean1 = row[f"{test}_mean_distance"]
    mean2 = row[f"{reference}_mean_distance"]
    sd1 = row[f"{test}_sd_distance"]
    sd2 = row[f"{reference}_sd_distance"]
    
    if mean1 > th * sd1 and mean2 > th * sd2:
        return True
    else:
        return False


# Read count for all samples
df = pd.read_csv(counts, sep="\t")

### Filtering of data
logging.info(f"Filtering data for {test} vs {reference}")

# Select columns that are part of the comparison
df = df.filter(regex=f"{reference}|{test}|^barcode_id|^orf_id|^gene")
nrows = df.shape[0]
logging.info(f"  Barcodes present pre-filtering: {nrows}")

# Remove barcodes where there are no counts accross all samples
df = df[df.iloc[:,3:].sum(axis=1) > 0].reset_index(drop=True)
nrows_zero_counts = nrows - df.shape[0]
logging.info(f"  Barcodes with no counts in any sample: {nrows_zero_counts}")

# Identify the sample with the most reads 
largest_sum = df.iloc[:,3:].sum().max()
largest_sample = df.iloc[:,3:].sum().idxmax()
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
nrows_single_barcode = nrows - df.shape[0]
logging.info(f"  ORFs removed that have only one barcode after filtering: {nrows_single_barcode}")

# Remove ORFs with less than a specified number of barcodes
nrows = df.shape[0]
df = df[df["num_barcodes"] >= bc_th].reset_index(drop=True)
nrows_low_barcodes = nrows - df.shape[0]
logging.info(f"  ORFs removed with less than {bc_th} barcodes after filtering: {nrows_low_barcodes}")

logging.info(f"  Number of barcodes present post-filtering: {nrows}")

### Compute PSI values:
# For each sample bin, divide the bin count by the sum of all bins for that sample
# Multiple by bin number and sum all values for each sample
logging.info(f"Computing PSI values for {test} vs {reference}")
for sample in [reference, test]:
    sample_bins = df.filter(regex=f"^{sample}").columns
    sample_bins = [int(x.replace(f"{sample}_", "")) for x in sample_bins]
    # It is assumed that the bins are numbered from 1 to max_bin!
    max_bin = max(sample_bins)
    # Iterate over rows to compute PSI values
    df[f"PSI_{sample}"] = df.apply(lambda row: compute_psi(row, sample, max_bin), axis=1).reset_index(drop=True)

# Calculate mean PSI values for each ORF
df[f"PSI_{reference}_mean"] = df.groupby("orf_id")[f"PSI_{reference}"].transform("mean")
df[f"PSI_{test}_mean"] = df.groupby("orf_id")[f"PSI_{test}"].transform("mean")
    
# Calculate deltaPSI for each single ORF
df["deltaPSI"] = df[f"PSI_{test}"] - df[f"PSI_{reference}"]
    
# Calculate mean deltaPSI values for each ORF
df["delta_PSI_mean"] = df.groupby("orf_id")["deltaPSI"].transform("mean")

# Calculate SD of PSI values for each condition of each ORF
df["delta_PSI_SD"] = df.groupby("orf_id")["deltaPSI"].transform("std") 

#IS THIS NEEDED?
# Sum of read counts for each ORF and each bin
#for sample in [reference, test]:
#    # Get columns with bin counts
#    sample_bins = df.filter(regex=f"^{sample}_").columns
#    # Get sum of read counts for each ORF and each bin
#    for bin in sample_bins:
#        df[bin] = df.groupby("orf_id")[bin].transform("sum")

### Hit identification
# Calculate the Euclidean distances between the barcodes
# This will be used to identify high confidence hits later
for sample in [reference, test]:
    # Get columns with normalised bin counts
    sample_columns = list(df.filter(regex=f"^{sample}").columns)
    sample_columns.sort()
    
    # Create a list of lists with the normalised bin counts for each ORF
    curves = [group[sample_columns].values.tolist() for _, group in df.groupby("orf_id", sort=False)]
    
    # Calculate the Euclidean distance between the curves
    ed_results = list(map(compute_euclidean_distance, curves))
    # Flatten the list of lists
    ed_results = [item for sublist in ed_results for item in sublist]
    
    # Add the mean and SD to the dataframe
    df[f"{sample}_mean_distance"] = [x[0] for x in ed_results]
    df[f"{sample}_sd_distance"] = [x[1] for x in ed_results]

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

# Identify high confidence hits based on Ecclidean distances between barcodes
df["high_confidence"] = df.apply(lambda row: stringent_hit(row, test, reference, sd_th), axis=1).reset_index(drop=True)
hc_stabilised = len(df[(df[f"stabilised_in_{test}"]) & (df["high_confidence"])])
logging.info(f"  Number of high confidence stabilised ORFs in {test}: {hc_stabilised}")
hc_destabilised = len(df[(df[f"destabilised_in_{test}"]) & (df["high_confidence"])])
logging.info(f"  Number of high confidence destabilised ORFs in {test}: {hc_destabilised}")

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
df_rank = df[(df[f"stabilised_in_{test}"]) | (df[f"destabilised_in_{test}"])].reset_index(drop=True)

# Collapse data to ORF level
df_rank = df_rank[["orf_id", "gene", "delta_PSI_mean", "delta_PSI_SD",
                   f"{test}_mean_distance", f"{test}_sd_distance",
                   f"{reference}_mean_distance", f"{reference}_sd_distance",
                   "num_barcodes", f"stabilised_in_{test}", 
                   f"stabilised_in_{test}_hc", f"destabilised_in_{test}", 
                   f"destabilised_in_{test}_hc"]].drop_duplicates().reset_index(drop=True)

# Calculate absolute rank based on signal-to-noise ratio (combined deltaPSI and Euclidean distances)
df_rank["SNR"] = abs(df_rank["delta_PSI_mean"]) / df_rank["delta_PSI_SD"] + df_rank[f"{test}_mean_distance"] / df_rank[f"{test}_sd_distance"] + df_rank[f"{reference}_mean_distance"] / df_rank[f"{reference}_sd_distance"]

# Move SNR column after num_barcodes
cols = list(df_rank.columns)
cols.insert(9, cols.pop(cols.index("SNR")))
df_rank = df_rank[cols]

# Correct SNR for the number of barcodes 
# SNR - (penalty * SNR) * (expected_barcodes - num_barcodes)
# Use median number of barcodes as expected number
df_rank["SNR"] = df_rank["SNR"] - (penalty * df_rank["SNR"]) * (df_rank["num_barcodes"].median() - df_rank["num_barcodes"])

df_rank = df_rank.sort_values(by="SNR", ascending=False).reset_index(drop=True)
df_rank["absolute_rank"] = df_rank.index + 1

# Create separate rankings for stabilised and destabilised hits
df_rank_stab = df_rank[df_rank[f"stabilised_in_{test}"]].sort_values(by="SNR", ascending=False).reset_index(drop=True)
df_rank_stab["stabilised_rank"] = df_rank_stab.index + 1

df_rank_destab = df_rank[df_rank[f"destabilised_in_{test}"]].sort_values(by="SNR", ascending=False).reset_index(drop=True)
df_rank_destab["destabilised_rank"] = df_rank_destab.index + 1

# Add these rankings to df_rank (NA for non-hits)
df_rank = pd.merge(df_rank, df_rank_stab[["orf_id", "stabilised_rank"]], on="orf_id", how="left")
df_rank = pd.merge(df_rank, df_rank_destab[["orf_id", "destabilised_rank"]], on="orf_id", how="left")

# Replace all missing values with NA
df_rank = df_rank.fillna("NA")

# Write to file
logging.info(f"Writing ranked results to {output_file_rank}")
df_rank.to_csv(output_file_rank, index=False)

logging.info("Done")