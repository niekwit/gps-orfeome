# Redirect R output to log
log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

library(tidyverse)
library(ggrepel)
source(file.path(snakemake@scriptdir, "theme_gpsw.R"))

csv_barcodes <- snakemake@input[["csv"]]
dpsi_cutoff <- as.numeric(snakemake@wildcards[["ht"]])
pdf <- snakemake@output[["pdf"]]
csv_out <- snakemake@output[["csv"]]
n_labels <- 10 # Number of flagged ORFs to label on either end

# ORF-level dPSI with twin peaks excluded (equivalent to delta_PSI_mean of
# calculate_psi.py) and included (mean over all barcodes). Both are computed
# from the (rounded) barcode-level deltaPSI so that ORFs without twin peak
# barcodes have a ddPSI of exactly 0
data <- read_csv(csv_barcodes, show_col_types = FALSE) %>%
  group_by(orf_id, gene) %>%
  summarise(
    twin_peaks_included_FALSE = mean(deltaPSI[!twin_peaks]),
    twin_peaks_included_TRUE = mean(deltaPSI),
    .groups = "drop"
  ) %>%
  filter(!is.na(twin_peaks_included_FALSE)) %>%
  mutate(
    deltaPSI_difference = twin_peaks_included_TRUE - twin_peaks_included_FALSE,
    # Order once here so every plot layer shares the same x positions
    gene = fct_reorder(gene, deltaPSI_difference),
    # |dPSI| passes the hit threshold in one analysis but not the other
    crosses_cutoff = (abs(twin_peaks_included_FALSE) >= dpsi_cutoff) !=
      (abs(twin_peaks_included_TRUE) >= dpsi_cutoff)
  ) %>%
  arrange(desc(deltaPSI_difference))

write_csv(mutate(data, across(where(is.double), ~ round(.x, 3))), csv_out)

# Label the flagged ORFs with the largest change at each end
flagged <- filter(data, crosses_cutoff)
labels <- bind_rows(
  slice_min(flagged, deltaPSI_difference, n = n_labels),
  slice_max(flagged, deltaPSI_difference, n = n_labels)
) %>%
  distinct()

# x positions of the contiguous block of ORFs with ddPSI = 0, i.e. without
# twin peaks (tolerance guards against floating point noise)
zero_pos <- which(
  levels(data$gene) %in% data$gene[abs(data$deltaPSI_difference) < 1e-8]
)

# Round limits out to the next 0.5 so the axis breaks span all points
y_limits <- c(
  floor(min(data$deltaPSI_difference) / 0.5) * 0.5,
  ceiling(max(data$deltaPSI_difference) / 0.5) * 0.5
)

# Create plot
p <- ggplot(data, aes(x = gene, y = deltaPSI_difference))

if (length(zero_pos) > 0) {
  zero_start <- min(zero_pos) - 0.5
  zero_end <- max(zero_pos) + 0.5
  p <- p +
    geom_vline(xintercept = c(zero_start, zero_end), linetype = "dotted") +
    annotate(
      "text",
      x = (zero_start + zero_end) / 2,
      y = max(data$deltaPSI_difference),
      label = paste0("ddPSI = 0\nn = ", length(zero_pos)),
      vjust = 1,
      size = 3
    )
}

p <- p +
  # One layer so all points share the x order; sorting draws red on top
  geom_point(
    data = arrange(data, crosses_cutoff),
    aes(colour = crosses_cutoff),
    size = 1.5
  ) +
  scale_colour_manual(
    values = c(`FALSE` = "grey70", `TRUE` = "firebrick"),
    guide = "none"
  ) +
  geom_text_repel(
    data = labels,
    aes(label = gene),
    size = 2.5,
    max.overlaps = Inf,
    min.segment.length = 0
  ) +
  labs(
    subtitle = paste0(
      "Red: |dPSI| crosses ",
      dpsi_cutoff,
      " in only one analysis (n = ",
      sum(data$crosses_cutoff),
      ")"
    ),
    x = "ORF",
    y = "ddPSI (twin peaks included - excluded)"
  ) +
  theme_gpsw() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),
    plot.subtitle = element_text(size = 12)
  ) +
  scale_x_discrete(expand = expansion(mult = c(0.02, 0.02))) +
  scale_y_continuous(
    limits = y_limits,
    breaks = seq(y_limits[1], y_limits[2], by = 0.5),
    expand = expansion(mult = c(0.05, 0))
  )

# Save plot
ggsave(pdf, p, width = 8, height = 6)
