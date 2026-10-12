import logging
import sys

import numpy as np
import pandas as pd

# Set up logging
log = snakemake.log[0]
logging.basicConfig(
    format="%(levelname)s:%(asctime)s:%(message)s",
    level=logging.DEBUG,
    datefmt="%Y-%m-%d %H:%M:%S",
    handlers=[logging.FileHandler(log)],
)

# Sensitivity of twin peak filtering (calculate_psi.identify_twin_peaks) to its
# proportion_threshold ("ratio"), for the comparison/hit_threshold combination
# given by the {comparison}/{ht} wildcards. The bin gap (2) stays fixed, as it
# is not configurable in the workflow.
#
# For each ratio threshold in a fixed sweep, barcodes are re-filtered for twin
# peaks from the aggregated counts (replicating the pre-twin-peak filtering of
# calculate_psi.py) and scored on:
# - cost: fraction of barcodes flagged, ORFs lost (< bc_threshold good
#   barcodes left after filtering)
# - precision: within-ORF dPSI SD after removing the flagged barcodes,
#   relative to removing the same number of random barcodes
# - hits: ORFs whose hit call (|dPSI| >= hit_threshold) changes between the
#   twin-peak-filtered and unfiltered dPSI

# Load Snakemake variables
counts = snakemake.input["counts"]
comparison = snakemake.wildcards["comparison"]
reference = comparison.split("_vs_")[1]
test = comparison.split("_vs_")[0]
MAX_BIN = snakemake.config["bin_number"]
sob_threshold = snakemake.config["psi"]["sob_threshold"]
bc_threshold = snakemake.config["psi"]["bc_threshold"]
hit_th = float(snakemake.wildcards["ht"])
configured_pt = float(snakemake.wildcards["pt"])
output_csv = snakemake.output["csv"]

GAP_MIN = 2  # fixed bin gap, matches calculate_psi.identify_twin_peaks
RATIOS = np.round(np.arange(0.05, 0.71, 0.05), 2)
N_PERM = 100  # random removals per ORF, for the precision null distribution
SEED = 42
# sd_gain_vs_random is a noisy median-of-medians, so ratios immediately after
# the best one are still treated as equivalent ("on the plateau") as long as
# each successive 0.05 step drops by no more than this many percentage
# points from the point before it; this is a per-step, not a cumulative,
# tolerance, so a real decline is not mistaken for a plateau just because it
# is still within some fixed distance of the original peak. The most
# permissive (highest) ratio still on the plateau is picked as the optimum,
# since it is strictly cheaper (fewer barcodes flagged, fewer ORFs lost) for
# essentially the same precision.
PLATEAU_STEP_TOLERANCE = 0.005

logging.info(f"Twin peak proportion_threshold sensitivity for {comparison}")
logging.info(
    f"  hit_threshold: {hit_th}, configured proportion_threshold: {configured_pt}"
)
logging.info(f"  Ratios swept: {RATIOS.tolist()}")

### Replicate calculate_psi.py filtering, up to (not including) twin peak
### identification, so the barcode set matches the production pipeline
df = pd.read_csv(counts, sep="\t")
df = df.filter(regex=f"{reference}|{test}|^barcode_id|^orf_id|^gene")

# Remove barcodes with no counts across all samples
df = df[df.iloc[:, 3:].sum(axis=1) > 0].reset_index(drop=True)

# Normalise reads to largest data set, as in calculate_psi.py
largest_sum = df.iloc[:, 3:].sum().max()
count_cols = [
    col for col in df.columns if col.startswith((f"{reference}_", f"{test}_"))
]
for col in count_cols:
    correction_factor = largest_sum / df[col].sum()
    df[col] = df[col].multiply(correction_factor)
    df[col] = df[col].round().astype(int)

# Sum of bins per condition
for sample in [reference, test]:
    bins = [f"{sample}_{i}" for i in range(1, MAX_BIN + 1)]
    df[f"SOB_{sample}"] = df[bins].sum(axis=1)

# Remove barcodes with low counts in either condition
df = df[df[f"SOB_{reference}"] > sob_threshold].reset_index(drop=True)
df = df[df[f"SOB_{test}"] > sob_threshold].reset_index(drop=True)

# Remove barcodes with no test counts in any bin (avoid division by zero)
df = df[df.filter(regex=f"{test}_").sum(axis=1) > 0].reset_index(drop=True)

# Remove ORFs with fewer than bc_threshold barcodes
df["num_barcodes"] = df.groupby("orf_id")["barcode_id"].transform("count")
df = df[df["num_barcodes"] >= bc_threshold].reset_index(drop=True)

logging.info(f"  Barcodes entering twin peak sweep: {df.shape[0]}")
logging.info(f"  ORFs entering twin peak sweep: {df['orf_id'].nunique()}")

if df.shape[0] == 0:
    raise ValueError(
        f"No barcodes remaining for {comparison} after filtering; "
        "cannot run the twin peak sensitivity sweep"
    )


### PSI / deltaPSI, as in calculate_psi.compute_psi
def psi(df, sample):
    bins = [f"{sample}_{i}" for i in range(1, MAX_BIN + 1)]
    m = df[bins].to_numpy(dtype=float)
    sob = df[f"SOB_{sample}"].to_numpy(dtype=float)
    weights = np.arange(1, MAX_BIN + 1)
    return (m * weights).sum(axis=1) / sob


df["deltaPSI"] = psi(df, test) - psi(df, reference)


### Peak ratio/gap per barcode and condition. These do not depend on the
### ratio threshold, so they are computed once and reused for every ratio in
### the sweep. Equivalent to calculate_psi.identify_twin_peaks, vectorised:
### the highest bin wins ties by being first in ascending bin order, same as
### the dict-based lookup there (verified against it on synthetic data).
def peak_stats(df, sample):
    bins = [f"{sample}_{i}" for i in range(1, MAX_BIN + 1)]
    m = df[bins].to_numpy(dtype=float)
    order = np.argsort(-m, axis=1, kind="stable")
    rows = np.arange(m.shape[0])
    top_idx = order[:, 0]
    second_idx = order[:, 1]
    top_val = m[rows, top_idx]
    second_val = m[rows, second_idx]
    with np.errstate(invalid="ignore", divide="ignore"):
        ratio = np.where(top_val > 0, second_val / top_val, 0.0)
    gap = np.abs(top_idx - second_idx)
    return ratio, gap


ratio_ref, gap_ref = peak_stats(df, reference)
ratio_test, gap_test = peak_stats(df, test)

delta_psi = df["deltaPSI"].to_numpy()
# Positional indices of each ORF's barcodes, computed once since orf_id does
# not change across the ratio sweep
groups = df.groupby("orf_id").indices
n_orfs_total = len(groups)

rng = np.random.default_rng(SEED)


def null_sd_median(vals, k, n_perm, rng):
    """Median SD after removing k random values from vals, over n_perm draws."""
    n = len(vals)
    keep = n - k
    rand_keys = rng.random((n_perm, n))
    perm_idx = np.argsort(rand_keys, axis=1)[:, :keep]
    samples = vals[perm_idx]
    return np.median(np.std(samples, axis=1, ddof=1))


records = []
for ratio_th in RATIOS:
    twin = ((ratio_ref > ratio_th) & (gap_ref >= GAP_MIN)) | (
        (ratio_test > ratio_th) & (gap_test >= GAP_MIN)
    )

    rows = []
    for orf, idx in groups.items():
        vals = delta_psi[idx]
        t = twin[idx]
        n = len(idx)
        k = int(t.sum())
        dpsi_unfiltered = vals.mean()
        dpsi_filtered = vals[~t].mean() if k < n else np.nan
        rows.append((orf, n, k, dpsi_unfiltered, dpsi_filtered))

    orfs = pd.DataFrame(
        rows, columns=["orf_id", "n", "k", "dpsi_unfiltered", "dpsi_filtered"]
    )
    # Keep ORFs with at least bc_threshold non-twin barcodes left (workflow rule)
    orfs = orfs[(orfs["n"] - orfs["k"]) >= bc_threshold]

    orfs_lost = n_orfs_total - len(orfs)
    hits_changed = int(
        (
            (orfs["dpsi_filtered"].abs() >= hit_th)
            != (orfs["dpsi_unfiltered"].abs() >= hit_th)
        ).sum()
    )

    # Precision: SD of dPSI after twin peak removal vs after random removal
    # of the same number of barcodes, for ORFs that actually had twin peaks
    # and still have at least 2 non-twin barcodes left (SD needs n >= 2)
    tp = orfs[(orfs["k"] > 0) & ((orfs["n"] - orfs["k"]) >= 2)]
    sd_twin = []
    sd_null_median = []
    for _, row in tp.iterrows():
        idx = groups[row["orf_id"]]
        vals = delta_psi[idx]
        t = twin[idx]
        sd_twin.append(np.std(vals[~t], ddof=1))
        sd_null_median.append(null_sd_median(vals, int(row["k"]), N_PERM, rng))

    if sd_twin:
        sd_gain_vs_random = 1 - np.median(sd_twin) / np.median(sd_null_median)
    else:
        sd_gain_vs_random = np.nan

    records.append(
        {
            "ratio": ratio_th,
            "barcodes_flagged": float(np.mean(twin)),
            "orfs_kept": len(orfs),
            "orfs_lost": int(orfs_lost),
            "orfs_with_twin": len(tp),
            "sd_gain_vs_random": sd_gain_vs_random,
            "hits_changed": hits_changed,
        }
    )
    logging.info(
        f"  ratio={ratio_th}: barcodes_flagged={np.mean(twin):.3f} "
        f"orfs_lost={orfs_lost} sd_gain_vs_random={sd_gain_vs_random} "
        f"hits_changed={hits_changed}"
    )

results = pd.DataFrame.from_records(records)
results.to_csv(output_csv, index=False)

### Optimum: starting from the ratio with the best SD reduction, extend
### towards higher (more permissive) ratios for as long as each next point
### drops by no more than PLATEAU_STEP_TOLERANCE from the point before it
### (not from the best value itself, so a real, sustained decline is not
### mistaken for a plateau just because it stays within some fixed distance
### of the peak), and take the last such ratio. sd_gain_vs_random is a noisy
### median-of-medians, so adjacent ratios are often statistically
### indistinguishable; among those, a higher ratio is cheaper (fewer
### barcodes flagged, fewer ORFs lost) for essentially the same precision.
### Warn if the configured proportion_threshold is still more permissive
### (higher) than this, i.e. it filters fewer real twin peaks than the
### sweep suggests is achievable without a real loss of precision.
valid = (
    results.dropna(subset=["sd_gain_vs_random"])
    .sort_values("ratio")
    .reset_index(drop=True)
)
if valid.empty:
    logging.warning(
        "No ORFs with twin peaks and enough remaining barcodes were found "
        "across the ratio sweep; cannot determine an optimal "
        "proportion_threshold for this comparison/hit_threshold."
    )
else:
    best_idx = int(valid["sd_gain_vs_random"].idxmax())
    best_gain = valid.loc[best_idx, "sd_gain_vs_random"]
    end_idx = best_idx
    for i in range(best_idx + 1, len(valid)):
        step = valid.loc[i - 1, "sd_gain_vs_random"] - valid.loc[i, "sd_gain_vs_random"]
        if step <= PLATEAU_STEP_TOLERANCE:
            end_idx = i
        else:
            break
    optimum = float(valid.loc[end_idx, "ratio"])
    plateau_ratios = valid.loc[best_idx:end_idx, "ratio"].tolist()
    logging.info(
        f"Best SD gain vs random: {best_gain:.4f} at ratio "
        f"{float(valid.loc[best_idx, 'ratio'])}; contiguous plateau (each "
        f"step dropping no more than {PLATEAU_STEP_TOLERANCE} from the "
        f"previous ratio): {plateau_ratios}"
    )
    logging.info(
        f"Optimal proportion_threshold (most permissive ratio on that "
        f"plateau): {optimum}"
    )
    if configured_pt > optimum:
        msg = (
            f"proportion_threshold ({configured_pt}) for {comparison} "
            f"(hit_threshold {hit_th}) is above the optimal threshold "
            f"({optimum}) found by the precision sweep: a lower "
            "proportion_threshold gives a larger reduction in within-ORF "
            f"dPSI SD from twin peak removal. See {output_csv} and the "
            "accompanying plot."
        )
        logging.warning(msg)
        print(f"WARNING: {msg}", file=sys.stderr)

logging.info("Done")
