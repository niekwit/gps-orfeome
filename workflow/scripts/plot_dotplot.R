# Redirect R output to log
log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

library(tidyverse)
library(cowplot)
library(ggrepel)

csv <- snakemake@input[["csv"]]
csv.rank <- snakemake@input[["ranked"]]
dpsi.cutoff <- snakemake@wildcards[["ht"]]
pdf <- snakemake@output[["pdf"]]

# Load data
data <- read_csv(csv)
data.rank <- read_csv(csv.rank) %>%
  select(orf_id, stabilised_rank, destabilised_rank)

# Add rank to data
data <- data %>%
  left_join(data.rank, by = "orf_id") %>%
  filter(!is.na(delta_PSI_mean))

# Split data in stabilised and destabilised genes/ORFs
# And log2 transform z_score_corr (preserve sign)
stabilised <- data %>%
  filter(delta_PSI_mean > 0,
         z_score_corr > 0)

labels.stabilised <- stabilised %>%
  filter(stabilised_rank <= 7) %>%
  select(gene, z_score_corr, delta_PSI_mean, delta_PSI_SD) %>%
  distinct()

destabilised <- data %>%
  filter(delta_PSI_mean < 0,
         z_score_corr < 0)

labels.destabilised <- destabilised %>%
  filter(destabilised_rank <= 10) %>%
  select(gene, z_score_corr, delta_PSI_mean, delta_PSI_SD) %>%
  distinct()

# Determine shared size limits for dot sizes
size_limits <- range(0 : ceiling(max(abs(data$delta_PSI_mean))))
range <- range(c(size_limits[1], size_limits[length(size_limits)]))

plotting <-function(data, labels, text, no.legend) {
  p <- ggplot(data, aes(x = z_score_corr, 
                             y = -log10(delta_PSI_SD))) +
  geom_point(aes(color = abs(delta_PSI_mean) > dpsi.cutoff,
                 size = abs(delta_PSI_mean))) +
  theme_cowplot(18) +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "black"),
                     name = paste("abs(dPSI) >", dpsi.cutoff)) +
  geom_label_repel(data = labels, 
                   aes(label = gene,),
                   box.padding = 0.5, 
                   point.padding = 0.5,
                   size = 5) +
  labs(title = text,
       x = "z-score",
       y = "-log10(dPSI SD)") +
  scale_size_continuous(limits = size_limits, 
                        range = range) +
  guides(
      color = guide_legend(title = paste("|dPSI| >", dpsi.cutoff)),
      size = guide_legend(title = "|dPSI|")
    )
  
  if (no.legend) {
    p <- p + theme(legend.position = "none")
  }
  return(p)
} 

# Create plots
p1 <- plotting(destabilised, labels.destabilised, "Destabilised proteins", TRUE)
p2 <- plotting(stabilised, labels.stabilised, "Stabilised proteins", FALSE)
# Extract legend and turn into plot
legend <- get_legend(p2)
p2 <- p2 + theme(legend.position = "none")
p3 <- plot_grid(p1, p2, legend, nrow = 1, rel_widths = c(1, 1, 0.25))

# Save plot
ggsave(pdf, p3, width = 20, height = 8)

# Close log file
sink(log, type = "output")
sink(log, type = "message")