# Redirect R output to log
log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

# ADAPTED (ANNOTATED AND SIMPLIFIED) FROM https://github.com/WubingZhang/MAGeCKFlute/blob/master/R/sgRankView.R

# Load required libraries
library(tidyverse)

# Load data
data <- read.delim(snakemake@input[[1]])

# Set parameters
select <- 5 # Number of genes to plot for enriched and depleted
binwidth <- 0.3 #
interval <- 0.1

# Create median LFC-based gene ranking
df <- data %>%
  group_by(Gene) %>%
  mutate(med.LFC = median(LFC)) %>%
  ungroup() %>%
  group_by(med.LFC) %>%
  mutate(rank = cur_group_id())

# Select top x genes with highest median lfc
df.top <- df[df$rank > max(df$rank) - select, ] %>%
  arrange(med.LFC)

# Select top x genes with lowest median lfc
df.bottom <- df %>%
  filter(rank <= select) %>%
  arrange(med.LFC)

# Add top and bottom genes data together and get genes to plot
df.sub <- rbind(df.top, df.bottom)
genes <- unique(df.sub$Gene)

# Add index, y values for plotting and colour param
df.sub <- df.sub %>%
  mutate(
    Gene = factor(Gene, levels = genes),
    index = rep(1:length(genes), as.numeric(table(Gene)[genes])),
    y = (binwidth + interval) * (index - 1),
    yend = (binwidth + interval) * index - interval,
    colour = ifelse(LFC > 0, "pos", "neg")
  ) %>%
  dplyr::select(c("sgrna", "Gene", "LFC", "y", "yend", "colour", "index")) %>%
  as.data.frame()

# Set the scale of x-axis
a <- -Inf
b <- Inf

# Set values for rectangle dimensions (x/y values)/colour fill to plot sgRNAs lfc inside
bgcol <- as.vector(sapply(seq(1, max(df.sub$index), 1), function(x) {
  rep(x, 4)
})) %>%
  as.data.frame() %>%
  rename("id" = 1) %>%
  mutate(
    x = rep(c(a, b, b, a), max(df.sub$index)),
    y = as.vector(sapply(seq(1, max(df.sub$index), 1), function(x) {
      c(
        (interval + binwidth) * (x - 1),
        (interval + binwidth) * (x - 1),
        (interval + binwidth) * x - interval,
        (interval + binwidth) * x - interval
      )
    }))
  )

# Create and save plot
p <- ggplot() +
  geom_polygon(
    aes_string("x", "y", group = "id"),
    fill = "#dedede",
    color = "gray20",
    data = bgcol
  ) +
  geom_segment(
    aes_string("LFC", "y", xend = "LFC", yend = "yend", color = "colour"),
    data = df.sub
  ) +
  scale_color_manual(values = c("pos" = "#e41a1c", "neg" = "#377eb8")) +
  scale_y_continuous(
    breaks = bgcol$y[seq(1, nrow(bgcol), 4)] + binwidth / 2,
    labels = genes,
    expand = c(0, 0)
  ) +
  labs(x = "Log2(Fold change)", y = NULL) +
  theme_bw(base_size = 14) +
  theme(legend.position = "none") +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank()
  )

ggsave(
  plot = p,
  filename = snakemake@output[[1]],
  units = "in",
  width = 10,
  height = 7
)

# Close redirection of output/messages
sink(log, type = "output")
sink(log, type = "message")
