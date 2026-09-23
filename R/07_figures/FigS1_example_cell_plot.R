# =============================================================================
# FIGURE S1, STEP 2 of 2 : TIME SERIES OF ONE EXAMPLE GRID CELL
# Original file name: Figure_pixel_drydown.R
# Methods 2.6
# -----------------------------------------------------------------------------
# Draws the example figure from the .rds file made by
# FigS1_example_cell_prepare.R: weekly CWD and EVI of one grid cell from 2000
# to 2023, with the largest drydown of each year shown as a shaded band and its
# CWD peak marked.
#
# Only needs data.table, ggplot2 and patchwork. Can be run on a laptop.
#
# Input  : FigureS1_data_<zone>.rds
# Output : FigureS1_example_timeseries.png
#
# Run    : set the two paths below and run the file.
# =============================================================================

library(data.table); library(ggplot2); library(patchwork)
setwd("/Net/Groups/BGI/scratch/mterristi/PhD/CWD_0_1/New_Climate_Zones_percentiles/Drydowns/Figures")
## ---- paths (change these lines) --------------------------------------------
data_file <- "FigureS1_data_temperate_humid.rds"          # the file from step 1
out_png   <- "FigureS1_example_timeseries.png"
## ---- read the prepared data ------------------------------------------------
d <- readRDS(data_file)
ts_cwd <- as.data.table(d$cwd)      # Date, CWD
evi_px <- as.data.table(d$evi)      # Date, EVI
bands  <- as.data.table(d$bands)    # year, start, end, peak_date, peak_cwd
peaks  <- as.data.table(d$peaks)    # Date, CWD, EVI
YEAR0  <- d$years[1]; YEAR1 <- d$years[2]
cat(sprintf("Pixel: lat %.4f, lon %.4f (%s) | %d drydown bands\n",
            d$lat, d$lon, d$zone, nrow(bands)))
## ---- figure ----------------------------------------------------------------
LINE <- "#E8663C"; DOT <- "#C0392B"; BAND <- "#E3E3E3"
xlim <- as.Date(c(paste0(YEAR0,"-01-01"), paste0(YEAR1,"-12-31")))
yr_breaks <- as.Date(paste0(seq(YEAR0, YEAR1, by = 2), "-01-01"))
base_ts <- theme_minimal(base_size = 15) +
  theme(panel.grid.minor   = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(colour = "grey90", linewidth = 0.5),
        axis.line          = element_line(colour = "grey40", linewidth = 0.4),
        axis.ticks         = element_line(colour = "grey40", linewidth = 0.4),
        axis.title.y       = element_text(size = 15),
        plot.margin        = margin(4, 12, 2, 6))
band_layer <- geom_rect(data = bands, aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf),
                        fill = BAND, inherit.aes = FALSE)
p_cwd <- ggplot() + band_layer +
  geom_line(data = ts_cwd, aes(Date, CWD), colour = LINE, linewidth = 0.4) +
  geom_point(data = peaks, aes(Date, CWD), colour = DOT, size = 2) +
  scale_x_date(limits = xlim, breaks = yr_breaks, date_labels = "%Y", expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
  labs(x = NULL, y = "CWD (mm)") +
  base_ts + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
                  axis.line.x = element_blank())
p_evi <- ggplot() + band_layer +
  geom_line(data = evi_px, aes(Date, EVI), colour = LINE, linewidth = 0.4) +
  geom_point(data = peaks, aes(Date, EVI), colour = DOT, size = 2) +
  scale_x_date(limits = xlim, breaks = yr_breaks, date_labels = "%Y", expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
  labs(x = NULL, y = "EVI") +
  base_ts + theme(axis.text.x = element_text(angle = 45, hjust = 1))
figS1 <- p_cwd / p_evi + plot_layout(heights = c(1, 1))
ggsave(out_png, figS1, width = 16, height = 6.5, dpi = 300, bg = "white", limitsize = FALSE)
cat("Saved:", out_png, "\n")