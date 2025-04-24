import logging
import pandas

# Set up logging
log = snakemake.log[0]
logging.basicConfig(
    format="%(levelname)s:%(asctime)s:%(message)s",
    level=logging.DEBUG,
    datefmt="%Y-%m-%d %H:%M:%S",
    handlers=[logging.FileHandler(log)],
)

# Load Snakemake variables
csv = snakemake.input["csv"]
bins = snakemake.params["bin_number"] + 1
comparison = snakemake.wildcards["comparison"]

test = comparison.split("_vs_")[0]
ref = comparison.split("_vs_")[1]

# Load data
logging.info(f"Loading data from {csv}")
data = pandas.read_csv(csv)

# Calculate proportions
logging.info("Calculating proportions")

df_copy = data.copy()  # Create a copy to avoid modifying the original DataFrame

# Calculate proportions for Test columns
test_cols = [f"{test}_{i}" for i in range(1, bins)]
test_sums = df_copy[test_cols].sum(axis=1)  # Sum across rows for Test columns
for col in test_cols:
    df_copy[col] = df_copy[col] / test_sums
    # Handle possible division by zero
    df_copy[col] = df_copy[col].fillna(0)

# Calculate proportions for Reference columns
ref_cols = [f"{ref}_{i}" for i in range(1, bins)]
ref_sums = df_copy[ref_cols].sum(axis=1)
for col in ref_cols:
    df_copy[col] = df_copy[col] / ref_sums
    # Handle possible division by zero
    df_copy[col] = df_copy[col].fillna(0)

# Add column with comparison name (move to first position)
# This saves time when plotting
df_copy.insert(0, "Comparison", comparison)

# Save results
output = snakemake.output["csv"]
logging.info(f"Saving results to {output}")
df_copy.to_csv(output, index=False)
