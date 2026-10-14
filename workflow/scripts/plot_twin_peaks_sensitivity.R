# Redirect R output to log
log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

library(tidyverse)
library(cowplot)
source(file.path(snakemake@scriptdir, "theme_gpsw.R"))

# Plots the ratio sweep written by twin_peaks_sensitivity.py:
# A: precision, i.e. within-ORF dPSI SD reduction from twin peak removal,
#    relative to removing the same number of random barcodes
# B: cost, i.e. ORFs lost and ORFs whose hit call changes, for context
# The configured proportion_threshold (dashed grey) and the optimum (dashed
# red) are marked on both panels. The optimum is the most permissive (highest)
# ratio on the plateau starting at the best SD reduction found, walked
# towards higher ratios for as long as each next 0.05 step drops by no more
# than plateau_step_tolerance from the point before it (not from the best
# value itself, so a real, sustained decline is not mistaken for a plateau
# just because it stays within some fixed distance of the peak): since
# sd_gain_vs_random is a noisy median-of-medians, adjacent ratios are often
# statistically indistinguishable, and a higher ratio is cheaper (fewer
# barcodes flagged, fewer ORFs lost) for essentially the same precision.
# This must match the optimum picked by twin_peaks_sensitivity.py.

csv <- snakemake@input[["csv"]]
pdf <- snakemake@output[["pdf"]]
configured_pt <- as.numeric(snakemake@wildcards[["pt"]])
hit_th <- as.numeric(snakemake@wildcards[["ht"]])
comparison <- snakemake@wildcards[["comparison"]]
plateau_step_tolerance <- 0.005

results <- read_csv(csv, show_col_types = FALSE)

valid <- results %>%
  filter(!is.na(sd_gain_vs_random)) %>%
  arrange(ratio)
has_optimum <- nrow(valid) > 0
if (has_optimum) {
  best_idx <- which.max(valid$sd_gain_vs_random)
  end_idx <- best_idx
  if (best_idx < nrow(valid)) {
    for (i in seq(best_idx + 1, nrow(valid))) {
      step <- valid$sd_gain_vs_random[i - 1] - valid$sd_gain_vs_random[i]
      if (step <= plateau_step_tolerance) {
        end_idx <- i
      } else {
        break
      }
    }
  }
  optimum <- valid$ratio[end_idx]
} else {
  optimum <- NA_real_
}
exceeds_optimum <- has_optimum && configured_pt > optimum

mark_lines <- function(p) {
  p <- p +
    geom_vline(
      xintercept = configured_pt,
      linetype = "dashed",
      colour = "grey40"
    )
  if (has_optimum) {
    p <- p +
      geom_vline(
        xintercept = optimum,
        linetype = "dashed",
        colour = "firebrick"
      )
  }
  p
}

subtitle <- if (exceeds_optimum) {
  sprintf(
    "proportion_threshold (%.2f, grey) is above the optimum (%.2f, red)",
    configured_pt,
    optimum
  )
} else if (has_optimum) {
  sprintf(
    "proportion_threshold (%.2f, grey); optimum %.2f (red)",
    configured_pt,
    optimum
  )
} else {
  "No ORFs with twin peaks found across the sweep"
}

# A: precision
p_precision <- mark_lines(
  ggplot(results, aes(ratio, sd_gain_vs_random)) +
    geom_line(linewidth = 0.5) +
    geom_point(size = 1.5) +
    scale_y_continuous(labels = scales::label_percent(1)) +
    labs(
      title = sprintf("%s, hit_threshold %.2f", comparison, hit_th),
      subtitle = paste0(
        "Within-ORF dPSI SD reduction vs random barcode removal\n",
        subtitle
      ),
      x = "proportion_threshold (ratio)",
      y = "SD reduction vs random"
    ) +
    theme_gpsw() +
    theme(
      plot.title = element_text(size = 11),
      plot.subtitle = element_text(
        size = 9,
        colour = if (exceeds_optimum) "firebrick" else "grey30"
      )
    )
)

# B: cost (ORFs lost, hits changed), on a shared count axis
cost <- results %>%
  select(ratio, `ORFs lost` = orfs_lost, `Hits changed` = hits_changed) %>%
  pivot_longer(-ratio, names_to = "metric", values_to = "value")

p_cost <- mark_lines(
  ggplot(cost, aes(ratio, value, colour = metric)) +
    geom_line(linewidth = 0.5) +
    geom_point(size = 1.5) +
    labs(
      title = "Cost of twin peak filtering",
      x = "proportion_threshold (ratio)",
      y = "ORFs",
      colour = NULL
    ) +
    theme_gpsw() +
    theme(legend.position = "top")
)

p <- plot_grid(p_precision, p_cost, ncol = 1, align = "v")

ggsave(pdf, p, width = 8, height = 8)
