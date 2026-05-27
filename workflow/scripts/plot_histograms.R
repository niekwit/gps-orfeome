# Redirect R output to log
log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

library(tidyverse)
source(file.path(snakemake@scriptdir, "theme_gpsw.R"))

csv <- snakemake@input[["csv"]]
comparison <- snakemake@wildcards[["comparison"]]

### Histogram for dPSI/dPSI SD values
plot_histogram <- function(data, column, x_label, pdf) {
  p <- ggplot(data, aes(x = .data[[column]])) +
    geom_histogram(
      bins = 60,
      colour = "black",
      fill = GPSW_COLOUR
    ) +
    theme_gpsw() +
    labs(x = x_label, y = "Count") +
    scale_y_continuous(expand = c(0, 0)) +
    geom_vline(
      xintercept = mean(data[[column]], na.rm = TRUE),
      linetype = "dashed",
      color = "red",
      linewidth = 0.5
    )

  # Save plot
  ggsave(pdf, p, width = 8, height = 6)
}

data <- read_csv(csv, show_col_types = FALSE) %>%
  dplyr::select(orf_id, gene, delta_PSI_mean, delta_PSI_SD) %>%
  distinct()

# Create histogram for dPSI values
print("Creating histogram for dPSI values...")
plot_histogram(
  data,
  "delta_PSI_mean",
  "deltaPSI",
  snakemake@output[["dpsi"]]
)

# Create histogram for dPSI SD values
print("Creating histogram for dPSI SD values...")
plot_histogram(
  data,
  "delta_PSI_SD",
  "deltaPSI SD",
  snakemake@output[["dpsi_sd"]]
)

### Histogram of PSI values for test vs control conditions
print("Creating histogram for PSI values...")
# Get column names that contain PSI values of test and control conditions
test <- str_split(comparison, "_vs_")[[1]][1]
test_column <- paste0("PSI_", test, "_mean")

control <- str_split(comparison, "_vs_")[[1]][2]
control_column <- paste0("PSI_", control, "_mean")

# Load data
data <- read_csv(csv, show_col_types = FALSE) %>%
  dplyr::select(
    orf_id,
    gene,
    all_of(test_column),
    all_of(control_column)
  ) %>%
  distinct() %>%
  pivot_longer(
    cols = c(all_of(test_column), all_of(control_column)),
    names_to = "condition",
    values_to = "psi"
  ) %>%
  mutate(
    condition = str_remove(condition, "PSI_"),
    condition = str_remove(condition, "_mean")
  )

# Create histogram
# Include histogram and smoothed density for each condition
p <- ggplot(data, aes(x = psi, fill = condition)) +
  geom_histogram(
    aes(y = ..density..),
    position = "identity",
    alpha = 0.35,
    bins = 60,
    colour = "black"
  ) +
  geom_density(
    aes(y = ..density.., colour = condition, fill = NULL),
    alpha = 0.35
  ) +
  theme_gpsw() +
  labs(
    x = "PSI",
    y = "Density",
  ) +
  scale_fill_manual(values = c("red", "grey")) +
  scale_color_manual(values = c("red", "grey")) +
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    legend.background = element_rect(color = NA),
  ) +
  scale_x_continuous(breaks = seq(0, 6, 1)) +
  scale_y_continuous(expand = c(0, 0))

# Save plot
ggsave(snakemake@output[["psi"]], p, width = 8, height = 6)


### Histogram of SOB values for test vs control conditions
print("Creating histogram for SOB values...")
csv <- snakemake@input[["sums"]]
data <- read_csv(csv, show_col_types = FALSE) %>%
  dplyr::select(orf_id, gene, starts_with("SOB")) %>%
  distinct() %>%
  pivot_longer(
    cols = starts_with("SOB"),
    names_to = "condition",
    values_to = "sob",
    names_prefix = "SOB_"
  )

# Create histogram
p <- ggplot(data, aes(x = sob)) +
  geom_histogram(
    aes(y = ..density.., fill = condition),
    position = "identity",
    alpha = 0.35,
    colour = "black"
  ) +
  geom_density(
    aes(y = ..density.., colour = condition),
    alpha = 0.35
  ) +
  theme_gpsw() +
  labs(
    x = "Sum of barcodes (SOB)",
    y = "Density",
  ) +
  scale_fill_manual(values = c("red", "grey")) +
  scale_color_manual(values = c("red", "grey")) +
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    legend.background = element_rect(color = NA),
  ) +
  scale_y_continuous(expand = c(0, 0)) +
  # set x-axis limits to include only 99.9th percentile of SOB values
  xlim(0, quantile(data$sob, 0.999, na.rm = TRUE)) +
  geom_vline(
    xintercept = snakemake@params[["sob_threshold"]],
    linetype = "dashed",
    color = "black",
    linewidth = 0.5
  )

# Save plot
ggsave(snakemake@output[["sob"]], p, width = 8, height = 6)
