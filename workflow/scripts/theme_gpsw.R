library(cowplot)

GPSW_COLOUR <- "#419179"
GPSW_BASE_SIZE <- 16

theme_gpsw <- function() {
  theme_cowplot(GPSW_BASE_SIZE)
}

geom_bar_gpsw <- function(...) {
  geom_bar(stat = "identity", fill = GPSW_COLOUR, colour = "black", ...)
}
