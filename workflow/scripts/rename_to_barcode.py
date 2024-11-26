''' 
Rename all instances of sgRNA/gene to barcode/orf in all relevant files
'''
import pandas as pd

### MAGeCK gene summary
df = pd.read_csv(snakemake.input["gs"], sep='\t')

# Rename sgrna substring in columns to barcode
df.columns = df.columns.str.replace("sgrna", "barcode")
df.to_csv(snakemake.output["gs"], sep='\t', index=False)

### MAGeCK sgrna summary
df = pd.read_csv(snakemake.input["ss"], sep='\t')

# Rename sgrna column to barcode
df = df.rename(columns={"sgrna": "barcode"})
df.to_csv(snakemake.output["ss"], sep='\t', index=False)

### Normalised MAGeCK count file
df = pd.read_csv(snakemake.input["norm"], sep='\t')

# Rename sgRNA/Gene columns to barcode/orf
df = df.rename(columns={"sgRNA": "barcode", "Gene": "orf"})
df.to_csv(snakemake.output["norm"], sep='\t', index=False)

