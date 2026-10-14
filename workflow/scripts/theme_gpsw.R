library(cowplot)

gpsw_colour <- "#419179"
gpsw_base_size <- 16

theme_gpsw <- function(...) {
  theme_cowplot(gpsw_base_size, ...)
}

geom_bar_gpsw <- function(...) {
  geom_bar(stat = "identity", fill = gpsw_colour, colour = "black", ...)
}
