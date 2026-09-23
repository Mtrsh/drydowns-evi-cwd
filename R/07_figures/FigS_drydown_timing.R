# =============================================================================
# SUPPLEMENTARY FIGURES : TIMING, DURATION AND VARIABILITY OF THE DRYDOWNS
# Original file name: Timing_drydowns_2.R
# Methods 2.6 (timing of the drydowns, circular median, Fig. S2)
# -----------------------------------------------------------------------------
# Describes the largest drydown of each year for each grid cell, from 2000 to
# 2023, using the drydown tables of 04a_drydown_events.R. Makes three figures
# from one table with one row per grid cell:
#   duration     median duration of the largest drydown, in weeks: map and
#                boxplot by zone
#   timing       start and end of the drydown: season (DJF, MAM, JJA, SON) as
#                maps, and month as rose diagrams by zone
#   variability  variability of the start month between years (circular
#                standard deviation): map and boxplot by zone
#
# The largest drydown is the drydown with the highest maximum CWD in the grid
# cell and year, as in Methods 2.6. We use the maximum CWD of the drydown
# itself here, not CWD_EVI95 (the CWD at the peak state).
#
# For each grid cell, we use the circular median of the months and the median
# of the durations over the years, for grid cells with at least five
# drydowns. Months are calendar months, without any change between the two
# hemispheres. They show when the drydown starts and ends, not when the peak
# and low states happen. Arid Hot is removed.
#
# Circular statistics: months are placed on a circle. The circular median is
# the direction given by the median of the sines and cosines. The circular
# standard deviation is sqrt(-2 ln R), with R the length of the mean vector,
# converted from radians to months.
#
# Input  : <zone>_drydown_events.Rdata from 04a_drydown_events.R, six zones
# Output : Figure_drydown_duration.png
#          Figure_drydown_timing.png
#          Figure_drydown_variability.png
#          Figure_drydown_characteristics_pixel.csv, one row per grid cell
#          Figure_drydown_characteristics_summary.csv, medians by zone
#          The file names come from our working folder, see the README.
#
# Run    : Rscript FigS_drydown_timing.R
# =============================================================================

## ---- environment -----------------------------------------------------------

# Our R package folder on the cluster, used only if it exists.
custom_libs <- c(
  "/Net/Groups/BGI/scratch/mterristi/R_lodomeria/4.4",
  "/Net/Groups/BSI/work_scratch/quincy/model/software/r_packages/r_4.4.x"
)
.libPaths(c(custom_libs[dir.exists(custom_libs)], .libPaths()))

# Libraries needed on our cluster, skipped if they do not exist.
for (lib in c("/opt/ohpc/pub/libs/hwloc/lib/libhwloc.so.15",
              "/opt/ohpc/pub/libs/gnu9/openmpi4/hdf5/1.10.8/lib/libhdf5_hl.so.100",
              "/opt/ohpc/pub/apps/gdal/3.5.1/lib/libgdal.so.31")) {
  if (file.exists(lib)) dyn.load(lib)
}

library(data.table)
library(ggplot2)
library(patchwork)
library(sf)
library(rnaturalearth)
library(lubridate)
library(scales)
library(viridis)

## ---- paths (change these lines to run the script elsewhere) ----------------
drydowns_dir <- "/Net/Groups/BGI/scratch/mterristi/PhD/CWD_0_1/New_Climate_Zones/Drydowns"
out_dir      <- "/Net/Groups/BGI/scratch/mterristi/PhD/CWD_0_1/New_Climate_Zones_percentiles/Drydowns/Figures"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

## ---- drydown tables, six zones ---------------------------------------------
cat("Loading drydown events per zone (Arid Hot excluded)...\n")

zone_files <- list(
  "Arid Cold"       = "arid_cold_drydown_events.Rdata",
  "Tropical"        = "tropical_drydown_events.Rdata",
  "Temperate Dry"   = "temperate_dry_drydown_events.Rdata",
  "Temperate Humid" = "temperate_humid_drydown_events.Rdata",
  "Cold Dry"        = "cold_dry_drydown_events.Rdata",
  "Cold Humid"      = "cold_humid_drydown_events.Rdata"
)

all_events <- list()
for (zn in names(zone_files)) {
  f <- file.path(drydowns_dir, zone_files[[zn]])
  if (!file.exists(f)) { cat("  Missing:", f, "- skipping\n"); next }
  env <- new.env(); load(f, envir = env)
  ev <- as.data.table(get(ls(env)[1], envir = env))
  ev[, climate_zone := zn]
  all_events[[zn]] <- ev
  cat("  Loaded", zn, ":", nrow(ev), "events\n")
}

events <- rbindlist(all_events, fill = TRUE)
cat("Total events:", nrow(events), "\n")

stopifnot(all(c("lat", "lon", "start_date", "end_date",
                "max_deficit", "year") %in% names(events)))

events[, start_date  := as.Date(start_date)]
events[, end_date    := as.Date(end_date)]
events[, start_month := month(start_date)]
events[, end_month   := month(end_date)]

## ---- largest drydown per grid cell and year --------------------------------
cat("\nSelecting annual-maximum event per pixel x year (max max_deficit)...\n")

worst <- events[events[, .I[which.max(max_deficit)],
                       by = .(lat, lon, year)]$V1]
worst <- worst[year >= 2000 & year <= 2023]
worst[, duration_w := as.numeric(end_date - start_date) / 7 + 1]
cat("Annual-maximum events kept:", nrow(worst), "\n")

## ---- functions: circular statistics and season -----------------------------
# Circular median of months 1 to 12 (see header).
circular_median <- function(months) {
  m <- months[!is.na(months)]
  if (length(m) == 0) return(NA_real_)
  angles <- 2 * pi * (m - 1) / 12
  med_angle <- atan2(median(sin(angles)), median(cos(angles)))
  (((med_angle / (2 * pi)) * 12) %% 12) + 1
}

# Circular standard deviation of months (see header). When all months are
# the same (R = 1), the result is 0, not NA.
circular_sd_months <- function(months) {
  m <- months[!is.na(months)]
  if (length(m) < 2) return(NA_real_)
  angles <- 2 * pi * (m - 1) / 12
  R <- sqrt(mean(sin(angles))^2 + mean(cos(angles))^2)
  R <- min(R, 1)                 # guard floating-point R slightly above 1
  if (R <= 0) return(NA_real_)   # fully dispersed onset
  (sqrt(-2 * log(R)) / (2 * pi)) * 12
}

# Season of a month, the same for both hemispheres.
month_to_quarter <- function(month) {
  mm <- round(month) %% 12; mm[mm == 0] <- 12
  q <- fifelse(mm %in% c(12, 1, 2), "DJF",
               fifelse(mm %in% c(3, 4, 5),  "MAM",
                       fifelse(mm %in% c(6, 7, 8),  "JJA", "SON")))
  factor(q, levels = c("DJF", "MAM", "JJA", "SON"))
}

## ---- summary over the years for each grid cell -----------------------------
cat("Computing per-pixel statistics...\n")

pixel_stats <- worst[, .(
  n_years            = .N,
  median_duration_w  = median(duration_w, na.rm = TRUE),
  median_start_month = circular_median(start_month),
  median_end_month   = circular_median(end_month),
  sd_start_month     = circular_sd_months(start_month),
  climate_zone       = first(climate_zone)
), by = .(lat, lon)]

pixel_stats <- pixel_stats[n_years >= 5]
cat("Pixels with >= 5 events:", nrow(pixel_stats), "\n")

## ---- Robinson projection and coastlines ------------------------------------
robin_crs <- "+proj=robin +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"

proj_pts <- st_transform(st_as_sf(pixel_stats, coords = c("lon", "lat"),
                                  crs = 4326), crs = robin_crs)
coords <- st_coordinates(proj_pts)
pixel_stats[, X_robin := coords[, 1]]
pixel_stats[, Y_robin := coords[, 2]]

world_robin <- st_crop(
  st_transform(ne_countries(scale = "medium", returnclass = "sf"), robin_crs),
  xmin = -17000000, xmax = 17000000, ymin = -6000000, ymax = 8500000)

## ---- zones, months and seasons for the plots -------------------------------
zone_levels <- c("Tropical", "Arid Cold",
                 "Temperate Dry", "Temperate Humid",
                 "Cold Dry", "Cold Humid")
pixel_stats[, climate_zone := factor(climate_zone, levels = zone_levels)]

month_labels <- c("Jan","Feb","Mar","Apr","May","Jun",
                  "Jul","Aug","Sep","Oct","Nov","Dec")
month_colors <- c("#A6CEE3","#1F78B4","#B2DF8A","#33A02C",
                  "#FB9A99","#E31A1C","#FDBF6F","#FF7F00",
                  "#CAB2D6","#6A3D9A","#FFFF99","#B15928")

# Month numbers for the rose diagrams.
to_month_factor <- function(x) {
  i <- round(x) %% 12; i[i == 0] <- 12
  factor(i, levels = 1:12, labels = month_labels)
}
pixel_stats[, start_label     := to_month_factor(median_start_month)]
pixel_stats[, end_label       := to_month_factor(median_end_month)]
pixel_stats[, start_month_int := as.integer(start_label)]
pixel_stats[, end_month_int   := as.integer(end_label)]

# Season for the maps.
quarter_colors <- c(DJF = "#3288BD", MAM = "#66C2A5",
                    JJA = "#FDAE61", SON = "#8C510A")
pixel_stats[, start_quarter := month_to_quarter(median_start_month)]
pixel_stats[, end_quarter   := month_to_quarter(median_end_month)]

# Function for the season maps.
quarter_map <- function(fill_col, title_text) {
  ggplot() +
    geom_tile(data = pixel_stats,
              aes(x = X_robin, y = Y_robin, fill = .data[[fill_col]]),
              width = 11000, height = 11000) +
    geom_sf(data = world_robin, fill = NA, color = "gray30", linewidth = 0.2) +
    scale_fill_manual(values = quarter_colors, name = "Calendar quarter", drop = FALSE) +
    coord_sf(crs = robin_crs, expand = FALSE) +
    labs(title = title_text) +
    theme_void(base_size = 18) +
    theme(plot.title   = element_text(size = 18, face = "plain", hjust = 0),
          legend.position = "bottom",
          legend.title = element_text(size = 15),
          legend.text  = element_text(size = 13)) +
    guides(fill = guide_legend(nrow = 1, byrow = TRUE))
}

# Function for the rose diagrams: distribution of a month per zone.
make_rose <- function(month_col, title_text) {
  d <- pixel_stats[, .N, by = c("climate_zone", month_col)]
  setnames(d, month_col, "month_int")
  d[, month_int := factor(month_int, levels = 1:12, labels = month_labels)]
  ggplot(d, aes(x = month_int, y = N, fill = month_int)) +
    geom_bar(stat = "identity", width = 1, color = "white", linewidth = 0.2) +
    coord_polar(start = -pi / 12) +
    facet_wrap(~ climate_zone, ncol = 3) +
    scale_fill_manual(values = setNames(month_colors, month_labels),
                      guide = "none") +
    labs(title = title_text, x = NULL, y = NULL) +
    theme_bw(base_size = 14) +
    theme(strip.text   = element_text(size = 13, face = "plain"),
          axis.text.x  = element_text(size = 11),
          axis.text.y  = element_blank(),
          axis.ticks.y = element_blank(),
          panel.grid.minor = element_blank(),
          plot.title   = element_text(size = 16, hjust = 0, face = "plain"))
}

## ---- duration: map and boxplot by zone -------------------------------------
cat("\n=== FIGURE 1: DURATION ===\n")
pixel_stats[, duration_disp := pmin(median_duration_w, 52)]  # display cap, one year

p_map_dur <- ggplot() +
  geom_tile(data = pixel_stats,
            aes(x = X_robin, y = Y_robin, fill = duration_disp),
            width = 11000, height = 11000) +
  geom_sf(data = world_robin, fill = NA, color = "gray30", linewidth = 0.2) +
  scale_fill_viridis(option = "magma", direction = -1,
                     limits = c(0, 52), breaks = c(0, 13, 26, 39, 52),
                     labels = c("0","13 (3mo)","26 (6mo)","39 (9mo)",">=52 (12mo)"),
                     name = "Median duration\n(weeks)") +
  coord_sf(crs = robin_crs, expand = FALSE) +
  labs(title = "a) Median duration of the annual-maximum drydown (2000-2023)") +
  theme_void(base_size = 18) +
  theme(plot.title = element_text(size = 18, face = "plain", hjust = 0),
        legend.position = "bottom", legend.key.width = unit(2, "cm"),
        legend.title = element_text(size = 15),
        legend.text  = element_text(size = 13))

p_box_dur <- ggplot(pixel_stats,
                    aes(y = climate_zone, x = median_duration_w,
                        fill = climate_zone)) +
  geom_boxplot(outlier.size = 0.3, outlier.alpha = 0.4) +
  scale_fill_viridis_d(option = "viridis", guide = "none") +
  scale_y_discrete(limits = rev(zone_levels)) +
  scale_x_continuous(breaks = c(0, 13, 26, 39, 52, 65)) +
  labs(title = "b) Duration per climate zone",
       x = "Median duration (weeks)", y = NULL) +
  theme_bw(base_size = 16) +
  theme(plot.title  = element_text(size = 16, hjust = 0, face = "plain"),
        axis.text.y = element_text(size = 14),
        axis.text.x = element_text(size = 13),
        panel.grid.major.y = element_blank())

fig1 <- p_map_dur / p_box_dur + plot_layout(heights = c(1.7, 1))
ggsave(file.path(out_dir, "Figure_drydown_duration.png"),
       fig1, width = 16, height = 16, dpi = 300, bg = "white")
cat("Saved: Figure_drydown_duration.png\n")

## ---- timing: season maps and month rose diagrams ---------------------------
cat("\n=== FIGURE 2: TIMING (quarter maps + month roses) ===\n")

p_map_start  <- quarter_map("start_quarter", "a) Calendar quarter of drydown onset (2000-2023)")
p_map_end    <- quarter_map("end_quarter",   "b) Calendar quarter of drydown end (2000-2023)")
p_rose_start <- make_rose("start_month_int", "c) Distribution of median onset month per zone")
p_rose_end   <- make_rose("end_month_int",   "d) Distribution of median end month per zone")

fig2 <- (p_map_start | p_map_end) / (p_rose_start | p_rose_end) +
  plot_layout(heights = c(1.4, 1.2))
ggsave(file.path(out_dir, "Figure_drydown_timing.png"),
       fig2, width = 24, height = 18, dpi = 300, bg = "white")
cat("Saved: Figure_drydown_timing.png\n")

## ---- variability of the start month: map and boxplot by zone ---------------
cat("\n=== FIGURE 3: VARIABILITY OF ONSET MONTH ===\n")

p_map_var <- ggplot() +
  geom_tile(data = pixel_stats[!is.na(sd_start_month)],
            aes(x = X_robin, y = Y_robin, fill = sd_start_month),
            width = 11000, height = 11000) +
  geom_sf(data = world_robin, fill = NA, color = "gray30", linewidth = 0.2) +
  scale_fill_viridis(option = "viridis", direction = 1, limits = c(0, 3),
                     name = "Circular SD\nof onset month\n(months)") +
  coord_sf(crs = robin_crs, expand = FALSE) +
  labs(title = "a) Interannual variability of onset month (2000-2023)") +
  theme_void(base_size = 18) +
  theme(plot.title = element_text(size = 18, face = "plain", hjust = 0),
        legend.position = "bottom", legend.key.width = unit(2, "cm"),
        legend.title = element_text(size = 15),
        legend.text  = element_text(size = 13))

p_box_var <- ggplot(pixel_stats[!is.na(sd_start_month)],
                    aes(y = climate_zone, x = sd_start_month,
                        fill = climate_zone)) +
  geom_boxplot(outlier.size = 0.3, outlier.alpha = 0.4) +
  scale_fill_viridis_d(option = "viridis", guide = "none") +
  scale_y_discrete(limits = rev(zone_levels)) +
  labs(title = "b) Interannual variability per climate zone",
       x = "Circular SD of onset month (months)", y = NULL) +
  theme_bw(base_size = 16) +
  theme(plot.title  = element_text(size = 16, hjust = 0, face = "plain"),
        axis.text.y = element_text(size = 14),
        axis.text.x = element_text(size = 13),
        panel.grid.major.y = element_blank())

fig3 <- p_map_var / p_box_var + plot_layout(heights = c(1.7, 1))
ggsave(file.path(out_dir, "Figure_drydown_variability.png"),
       fig3, width = 16, height = 16, dpi = 300, bg = "white")
cat("Saved: Figure_drydown_variability.png\n")

## ---- table per grid cell and summary per zone ------------------------------
fwrite(pixel_stats, file.path(out_dir, "Figure_drydown_characteristics_pixel.csv"))

zone_summary <- pixel_stats[, .(
  n_pixels              = .N,
  median_duration_w     = round(median(median_duration_w, na.rm = TRUE), 1),
  median_sd_start_month = round(median(sd_start_month,    na.rm = TRUE), 2),
  median_start_month    = round(circular_median(median_start_month), 1),
  median_end_month      = round(circular_median(median_end_month),   1)
), by = climate_zone]
print(zone_summary)
fwrite(zone_summary, file.path(out_dir, "Figure_drydown_characteristics_summary.csv"))

cat("\n=== ALL DONE ===\n")