# Redirect R output to log
log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

library(tidyverse)
library(cowplot)
library(ggrepel)

csv <- snakemake@input[["csv"]]
csv.rank <- snakemake@input[["ranked"]]
dpsi.cutoff <- as.numeric(snakemake@wildcards[["ht"]])
pdf <- snakemake@output[["pdf"]]

# Load data
data <- read_csv(csv)
data.rank <- read_csv(csv.rank) %>%
  select(orf_id, stabilised_rank, destabilised_rank)

# Add rank to data
data <- data %>%
  left_join(data.rank, by = "orf_id") %>%
  filter(!is.na(delta_PSI_mean))

# Determine x min/x masx for delta_PSI_mean
max <- ceiling(max(abs(data$delta_PSI_mean), na.rm = TRUE))

# Create values for color
data <- data %>%
  mutate(
    category = ifelse(delta_PSI_mean > 0, "Stabilised", "Destabilised"),
    colour_group = factor(
      case_when(
        category == "Stabilised" & delta_PSI_mean > dpsi.cutoff ~
          "Significant Stabilised",
        category == "Destabilised" & delta_PSI_mean < -dpsi.cutoff ~
          "Significant Destabilised",
        TRUE ~ "Other" # Catches all other cases
      ),
      levels = c(
        "Significant Stabilised",
        "Significant Destabilised",
        "Within Threshold",
        "Other"
      )
    ) # Define factor levels for consistent legend order
  )

my.colours <- c(
  "Significant Stabilised" = "green3",
  "Significant Destabilised" = "red",
  "Other" = "black"
)

# Create df for labels
stabilised <- data %>%
  filter(delta_PSI_mean > 0, z_score_corr > 0)

labels.stabilised <- stabilised %>%
  filter(stabilised_rank <= 7) %>%
  select(
    gene,
    z_score_corr,
    delta_PSI_mean,
    delta_PSI_SD,
    category
  ) %>%
  distinct()

destabilised <- data %>%
  filter(delta_PSI_mean < 0, z_score_corr < 0)

labels.destabilised <- destabilised %>%
  filter(destabilised_rank <= 7) %>%
  select(
    gene,
    z_score_corr,
    delta_PSI_mean,
    delta_PSI_SD,
    category
  ) %>%
  distinct()

labels <- rbind(labels.stabilised, labels.destabilised)

# Data for vertical lines
vline_data <- data.frame(
  category = factor(
    c("Destabilised", "Stabilised"),
    levels = c("Destabilised", "Stabilised")
  ),
  xintercept_value = c(-dpsi.cutoff, dpsi.cutoff)
)

# Create plot
p <- ggplot(data, aes(y = abs(z_score_corr), x = delta_PSI_mean)) +
  geom_point(aes(colour = colour_group)) +
  facet_wrap(~category, scales = "free_x") +
  theme_cowplot(15) +
  scale_color_manual(
    values = my.colours,
  ) +
  labs(y = "|z-score|", x = "dPSI") +
  guides(
    colour = "none"
  ) +
  geom_vline(
    data = vline_data,
    aes(xintercept = xintercept_value),
    color = "grey",
    linetype = "dashed"
  ) +
  geom_label_repel(
    data = labels,
    aes(label = gene, ),
    box.padding = 0.5,
    point.padding = 0.5,
    size = 2.5,
    max.overlaps = 20
  )

# Save plot
ggsave(pdf, p, width = 8, height = 5)
