# Redirect R output to log
log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

library(tidyverse)
source(file.path(snakemake@scriptdir, "theme_gpsw.R"))

csv.rank <- snakemake@input[["ranked"]]
dpsi.cutoff <- as.numeric(snakemake@wildcards[["ht"]])
pdf <- snakemake@output[["pdf"]]
N_LABELS <- 15 # Number of ORFs to label on either end

# Load ORF-level data and order by dPSI (lowest to highest)
data <- read_csv(csv.rank, show_col_types = FALSE) %>%
  filter(!is.na(delta_PSI_mean)) %>%
  arrange(delta_PSI_mean) %>%
  mutate(
    rank = row_number(),
    colour_group = factor(
      case_when(
        stabilised ~ "Stabilised",
        destabilised ~ "Destabilised",
        TRUE ~ "Other"
      ),
      levels = c("Stabilised", "Destabilised", "Other")
    )
  )

my.colours <- c(
  "Stabilised" = "green3",
  "Destabilised" = "red",
  "Other" = "black"
)

# Labels for the (up to) N_LABELS destabilised/stabilised ORFs with the lowest/highest dPSI values
# Labels are stacked in evenly spaced columns in the empty space above (lowest)
# and below (highest) the curve, in the same order as the points so that the
# connecting lines do not cross
# Labels use a fixed spacing and start at the end of the column closest to the
# points, so that a few labels stay together
n.orfs <- nrow(data)
median.dpsi <- median(data$delta_PSI_mean)
y.range <- diff(range(data$delta_PSI_mean))
low.start <- median.dpsi + 0.05 * y.range
low.step <- (max(data$delta_PSI_mean) - low.start) / (N_LABELS - 1)
high.end <- median.dpsi - 0.05 * y.range
high.step <- (high.end - min(data$delta_PSI_mean)) / (N_LABELS - 1)

labels.low <- data %>%
  filter(destabilised) %>%
  slice_head(n = N_LABELS) %>%
  mutate(
    label_x = n.orfs * 0.1,
    label_y = low.start + (row_number() - 1) * low.step,
    hjust = 0
  )
labels.high <- data %>%
  filter(stabilised) %>%
  slice_tail(n = N_LABELS) %>%
  mutate(
    label_x = n.orfs * 0.9,
    label_y = high.end - (n() - row_number()) * high.step,
    hjust = 1
  )
labels <- bind_rows(labels.low, labels.high)

# Create plot
p <- ggplot(data, aes(x = rank, y = delta_PSI_mean)) +
  geom_hline(
    yintercept = c(-dpsi.cutoff, dpsi.cutoff),
    color = "grey",
    linetype = "dashed"
  ) +
  geom_segment(
    data = labels,
    aes(xend = label_x, yend = label_y),
    colour = "grey50",
    linewidth = 0.2
  ) +
  geom_point(aes(colour = colour_group), size = 0.8) +
  geom_text(
    data = labels,
    aes(x = label_x, y = label_y, label = gene, hjust = hjust),
    size = 2.5
  ) +
  theme_gpsw() +
  scale_color_manual(values = my.colours) +
  labs(x = "ORF rank", y = "dPSI") +
  guides(colour = "none")

# Save plot
ggsave(pdf, p, width = 8, height = 6)
