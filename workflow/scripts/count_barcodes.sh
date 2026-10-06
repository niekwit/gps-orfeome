#!/usr/bin/env bash

# Counts reads per barcode, as in the original GPS analysis scripts:
# - only end-to-end alignments with at most `mismatch` mismatches are kept
#   (mismatch: 0 = perfect matches only, i.e. --score-min C,0,0); the
#   maximum bowtie2 mismatch penalty is 6, so the minimum score is -6 per
#   allowed mismatch
# - reads with more than one equally good hit are kept (bowtie2 reports one
#   of them at random) instead of discarded
# - unaligned reads (no reference, "*") are not counted

mm=${snakemake_params[mm]}
min_score=$((-6 * mm))
seed_mm=$((mm > 0 ? 1 : 0))

zcat ${snakemake_input[fq]} | \
bowtie2 ${snakemake_params[extra]} --no-hd -p ${snakemake[threads]} -t \
    --end-to-end --score-min C,${min_score},0 -N ${seed_mm} \
    -x ${snakemake_params[idx]} - 2> ${snakemake_log[0]} | \
awk -F '\t' '$3 != "*" { print $3 }' | \
sort | \
uniq -c | \
sed 's/^ *//' > ${snakemake_output[0]}
