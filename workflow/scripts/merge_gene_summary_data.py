"""
Merge gene summary data from multiple conditions.
"""

import os
import pandas as pd

rank_files = snakemake.input.ranks
output_file = snakemake.output[0]

# Read all rank files into a list of DataFrames
dfs = [pd.read_csv(file) for file in rank_files]

data = []
for i, df in enumerate(dfs):
    # Get the columns that start with stabilised_ or destabilised_
    cols = [
        col for col in df.columns if col.startswith(("stabilised_", "destabilised_"))
    ]

    # Extract the condition from the columns
    condition = os.path.basename(rank_files[i]).replace("_gene.summary.csv", "")

    # Append the condition to the columns named delta_PSI_mean and z_score_corr
    df.rename(
        columns={
            "delta_PSI_mean": f"delta_PSI_mean_{condition}",
            "z_score_corr": f"z_score_corr_{condition}",
        },
        inplace=True,
    )

    # Only keep the columns orf_id, gene and the ones that contain the condition in their name
    df = df[["orf_id", "gene"] + [col for col in df.columns if condition in col]]

    # Add the DataFrame to the data list
    data.append(df)

# Outer join all DataFrames on orf_id and gene
df = (
    pd.concat(data, axis=0, ignore_index=True)
    .groupby(["orf_id", "gene"], as_index=False)
    .first()
)

# Save to final output file
df.to_csv(output_file, index=False, na_rep="NA")
