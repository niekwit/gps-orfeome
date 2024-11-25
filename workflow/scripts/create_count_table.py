import os
import pandas as pd

''' 
Merge count files to create count table.
With bin_number = 1,
MAGeCK requires the first two columns to be sgRNA and gene.
Will be renamed to barcode and ORF, respectively post-MAGeCK.

Otherwise, the count table will be created with the first column as the barcode name
and the second the ORF name.
'''

if snakemake.config["bin_number"] == 1:
    # For use with MAGeCK
    name = "sgRNA"
    gene = "gene"
else:
    name = "barcode"
    gene = "gene"

# Read all barcode names from csv and create ORF column
df = pd.read_csv(snakemake.input["csv"], header=None)
df = df[snakemake.config["csv"]["barcode_column"] - 1]
df[gene] = df[0].str.rsplit(pat="_", n=1, expand=True)[0]

# Rename first column 
df = df.rename(columns={0: name})

# Merge all count files to df to create count table
for file in snakemake.input["files"]:
    # Read count file
    tmp = pd.read_csv(file, sep=" ", header=None)
    
    # Rename headers to wildcard value (for counts) and sgRNA 
    wildcard_value = os.path.basename(file).replace(".barcode.counts.txt", "")
    tmp.columns = [wildcard_value, name]
    
    # Merge to df
    df = pd.merge(df, tmp, on=name, how="left")

# Replace missing values with zero
df = df.fillna(0)

# Save data frame to file
df.to_csv(snakemake.output[0], sep='\t', index=False)
