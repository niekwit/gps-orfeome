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

# Load Snakemake variables
csv_file = snakemake.input["csv"]
count_files = snakemake.input["files"]
bin_number = snakemake.config["bin_number"]
barcode_id_column = snakemake.config["csv"]["barcode_id_column"]
orf_column = snakemake.config["csv"]["orf_column"]
gene_column = snakemake.config["csv"]["gene_column"]
out_file = snakemake.output[0]

if bin_number == 1:
    # For use with MAGeCK
    name = "sgRNA"
else:
    name = "barcode_id"

# Read all barcode names from csv and create ORF column
csv = pd.read_csv(csv_file)

# With MAGeCK, use ORF names as gene names
if bin_number == 1:
    df = csv.iloc[:, [barcode_id_column, orf_column]]
else:
    df = csv.iloc[:, [barcode_id_column, gene_column]]
df.columns = [name, "gene"]

# Merge all count files to df to create count table
for file in count_files:
    # Read count file
    tmp = pd.read_csv(file, sep=" ", header=None)
    
    # Rename headers to wildcard value (for counts) and sgRNA 
    wildcard_value = os.path.basename(file).replace(".barcode.counts.txt", "")
    tmp.columns = [wildcard_value, name]
    
    # Merge to df
    df = pd.merge(df, tmp, on=name, how="left")

# Replace missing values with zero
df = df.fillna(0)

# For non-MAGeCK data, add ORF column as 3rd column
if bin_number != 1:
    df.insert(1, "orf_id", csv[csv.columns[orf_column]])

# Save data frame to file
df.to_csv(out_file, sep='\t', index=False)
