# Redirect R output to log
log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

library(tidyverse)
library(ggrepel)
source(file.path(snakemake@scriptdir, "theme_gpsw.R"))

csv_rank <- snakemake@input[["ranked"]]
dpsi_cutoff <- as.numeric(snakemake@wildcards[["ht"]])
pdf <- snakemake@output[["pdf"]]
n_labels <- 10 # Number of ORFs to label on either end

# Load ORF-level data. ORFs with a single good barcode have no SD (NA) and
# are dropped, since they cannot be placed on this plot
data <- read_csv(csv_rank, show_col_types = FALSE) %>%
  filter(!is.na(delta_PSI_mean), !is.na(delta_PSI_SD)) %>%
  mutate(
    colour_group = factor(
      case_when(
        stabilised ~ "Stabilised",
        destabilised ~ "Destabilised",
        TRUE ~ "Other"
      ),
      levels = c("Stabilised", "Destabilised", "Other")
    )
  )

my_colours <- c(
  "Stabilised" = "#14bd14",
  "Destabilised" = "red",
  "Other" = "grey70"
)

# Label the top n_labels stabilised/destabilised hits (by existing rank)
labels <- data %>%
  filter(
    (stabilised & stabilised_rank <= n_labels) |
      (destabilised & destabilised_rank <= n_labels)
  )

# Create plot
p <- ggplot(data, aes(x = delta_PSI_mean, y = delta_PSI_SD)) +
  geom_vline(
    xintercept = c(-dpsi_cutoff, dpsi_cutoff),
    colour = "grey",
    linetype = "dotted"
  ) +
  geom_point(aes(colour = colour_group), size = 1) +
  geom_text_repel(
    data = labels,
    aes(label = gene),
    size = 3,
    max.overlaps = Inf,
    min.segment.length = 0
  ) +
  theme_gpsw() +
  scale_colour_manual(values = my_colours) +
  labs(x = "dPSI", y = "dPSI SD") +
  guides(colour = "none")

# Save plot
ggsave(pdf, p, width = 8, height = 6)
