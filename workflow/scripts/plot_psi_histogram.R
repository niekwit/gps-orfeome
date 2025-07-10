# Redirect R output to log
log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

library(tidyverse)
library(cowplot)

# Get column names that contain PSI values of test and control conditions
comparison <- snakemake@wildcards[["comparison"]]

test <- str_split(comparison, "_vs_")[[1]][1]
test_column <- paste0("PSI_", test, "_mean")

control <- str_split(comparison, "_vs_")[[1]][2]
control_column <- paste0("PSI_", control, "_mean")

# Get PSI values
data <- read_csv(snakemake@input[["csv"]], show_col_types = FALSE) %>%
    dplyr::select(orf_id, gene, {{ test_column }}, {{ control_column }}) %>%
    # remove duplciated rows (only keep 1 copy)
    distinct()

# Convert to long format
data_long <- data %>%
    pivot_longer(
        cols = c({{ test_column }}, {{ control_column }}),
        names_to = "condition",
        values_to = "psi"
    ) %>%
    mutate(
        condition = str_remove(condition, "PSI_"),
        condition = str_remove(condition, "_mean")
    )

# Create histogram plot
# Include histogram and smoothed density for each condition
p <- ggplot(data_long, aes(x = psi, fill = condition)) +
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
    theme_cowplot(16) +
    labs(
        x = "PSI",
        y = "Density",
    ) +
    scale_fill_manual(values = c("grey", "red")) +
    scale_color_manual(values = c("grey", "red")) +
    theme(
        legend.position = "top",
        legend.title = element_blank(),
        legend.background = element_rect(color = NA),
    ) +
    scale_x_continuous(breaks = seq(0, 6, 1)) +
    scale_y_continuous(expand = c(0, 0))


# Save plot
ggsave(snakemake@output[["pdf"]], p, width = 8, height = 6)
