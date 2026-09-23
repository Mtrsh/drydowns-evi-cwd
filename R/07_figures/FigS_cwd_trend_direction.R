# =============================================================================
# SUPPLEMENTARY FIGURE : DIRECTION OF THE CWD TRENDS AT THE PEAK AND LOW STATES
# Original file name: Figures_CWD_trensholds_trends.R
# Methods 2.8
# -----------------------------------------------------------------------------
# Same as Fig. 4, but for CWD: each grid cell gets a class from the signs of
# the trends of CWD_EVI95 and CWD_EVI05, per regime, with the distribution of
# each class by latitude. Grid cells with a zero or missing trend for one of
# the two states are removed. Grid cells without a significant trend
# (p >= 0.05) are shown with transparency.
#
# The colours are the opposite of the EVI figure: an increasing deficit is
# brown, a decreasing deficit is teal.
#
# Arid Hot is removed, by keeping only the grid cells of the pooled analysis.
#
# Input  : sens_slope_CWD_atE95_<regime>_v2.Rdata, sens_slope_CWD_atE05_<regime>_v2.Rdata
#          combined_pooled_percentiles_summary_v2.Rdata (climate zones)
# Output : Figure_trend_direction_CWD_by_regime.png
#          The file name comes from our working folder, see the README.
#
# Run    : Rscript FigS_cwd_trend_direction.R
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
library(scales)

## ---- paths (change these lines to run the script elsewhere) ----------------
trends_dir <- "/Net/Groups/BGI/scratch/mterristi/PhD/CWD_0_1/New_Climate_Zones_percentiles/Trends/Trends_V2"
v2_dir     <- "/Net/Groups/BGI/scratch/mterristi/PhD/CWD_0_1/New_Climate_Zones_percentiles/Drydowns/Drydowns_percentiles_pooled_v2"
out_dir    <- "/Net/Groups/BGI/scratch/mterristi/PhD/CWD_0_1/New_Climate_Zones_percentiles/Trends/Trends_V2/Figures"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

p_threshold <- 0.05

## ---- climate zone of each grid cell, without Arid Hot ----------------------
cat("Loading V2 pooled data for climate-zone filter...\n")
load(file.path(v2_dir, "combined_pooled_percentiles_summary_v2.Rdata"))
setDT(all_summaries)

zone_lookup <- unique(all_summaries[
  climate_zone != "Arid Hot",
  .(lat = Latitude, lon = Longitude, climate_zone)
])
cat("Pixels in lookup (no Arid Hot):", nrow(zone_lookup), "\n")

## ---- trends ----------------------------------------------------------------
load_trend <- function(metric, regime) {
  f <- file.path(trends_dir,
                 paste0("sens_slope_", metric, "_", regime, "_v2.Rdata"))
  if (!file.exists(f)) { cat("  Missing:", f, "\n"); return(NULL) }
  env <- new.env(); load(f, envir = env)
  d <- as.data.table(get(ls(env)[1], envir = env))
  d[, metric := metric][, regime := regime]
  d
}

cat("\nLoading CWD_atE95 and CWD_atE05 trends...\n")
c95_b <- load_trend("CWD_atE95", "before")
c95_a <- load_trend("CWD_atE95", "after")
c05_b <- load_trend("CWD_atE05", "before")
c05_a <- load_trend("CWD_atE05", "after")

## ---- sign of the trends ----------------------------------------------------
classify <- function(d) {
  d[, sign := fcase(
    is.na(sens_slope),                       "flat",
    sens_slope > 0,                          "up",
    sens_slope < 0,                          "down",
    default = "flat"
  )]
  d[, sig := !is.na(sens_pvalue) & sens_pvalue < p_threshold]
  d
}

c95_b <- classify(c95_b)
c95_a <- classify(c95_a)
c05_b <- classify(c05_b)
c05_a <- classify(c05_a)

## ---- one table per regime --------------------------------------------------
build_regime <- function(d95, d05, regime_label) {
  m <- merge(d95[, .(lat, lon, sign_C95 = sign, sig_C95 = sig)],
             d05[, .(lat, lon, sign_C05 = sign, sig_C05 = sig)],
             by = c("lat", "lon"))
  m <- merge(m, zone_lookup, by = c("lat", "lon"))
  m[, regime := regime_label]
  
  # Class from the signs of the two trends.
  m[, category := fcase(
    sign_C95 == "flat" & sign_C05 == "flat",
    "CWD_atE95 \u2192 / CWD_atE05 \u2192",
    sign_C95 == "up"   & sign_C05 == "up",
    "CWD_atE95 \u2191 / CWD_atE05 \u2191",
    sign_C95 == "down" & sign_C05 == "down",
    "CWD_atE95 \u2193 / CWD_atE05 \u2193",
    sign_C95 == "down" & sign_C05 == "up",
    "CWD_atE95 \u2193 / CWD_atE05 \u2191",
    sign_C95 == "up"   & sign_C05 == "down",
    "CWD_atE95 \u2191 / CWD_atE05 \u2193",
    default = "CWD_atE95 \u2192 / CWD_atE05 \u2192"
  )]
  
  m[, both_sig := sig_C95 & sig_C05]
  m
}

cat("\nBuilding regime datasets...\n")
before <- build_regime(c95_b, c05_b, "before")
after  <- build_regime(c95_a, c05_a, "after")
cat("Before pixels:", nrow(before), "| After pixels:", nrow(after), "\n")

cat_levels <- c(
  "CWD_atE95 \u2191 / CWD_atE05 \u2191",
  "CWD_atE95 \u2193 / CWD_atE05 \u2191",
  "CWD_atE95 \u2193 / CWD_atE05 \u2193",
  "CWD_atE95 \u2191 / CWD_atE05 \u2193",
  "CWD_atE95 \u2192 / CWD_atE05 \u2192"
)
before[, category := factor(category, levels = cat_levels)]
after[,  category := factor(category, levels = cat_levels)]

# Share of each class per regime.
cat("\n=== Category distribution — BEFORE ===\n")
print(before[, .(N = .N, pct = round(100 * .N / nrow(before), 1)),
             by = category][order(category)])
cat("\n=== Category distribution — AFTER ===\n")
print(after[, .(N = .N, pct = round(100 * .N / nrow(after), 1)),
            by = category][order(category)])

cat("\nBoth significant (BEFORE): ",
    sum(before$both_sig, na.rm = TRUE), "/", nrow(before),
    sprintf(" (%.1f%%)\n", 100 * mean(before$both_sig, na.rm = TRUE)))
cat("Both significant (AFTER):  ",
    sum(after$both_sig,  na.rm = TRUE), "/", nrow(after),
    sprintf(" (%.1f%%)\n", 100 * mean(after$both_sig,  na.rm = TRUE)))

## ---- Robinson projection ---------------------------------------------------
robin_crs <- "+proj=robin +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"

project_robin <- function(d) {
  pts <- st_as_sf(d, coords = c("lon", "lat"), crs = 4326)
  pts <- st_transform(pts, crs = robin_crs)
  coords <- st_coordinates(pts)
  d[, X_robin := coords[, 1]]
  d[, Y_robin := coords[, 2]]
  d
}

before <- project_robin(before)
after  <- project_robin(after)

world <- ne_countries(scale = "medium", returnclass = "sf")
world_robin <- st_transform(world, robin_crs)
world_robin <- st_crop(world_robin,
                       xmin = -17000000, xmax = 17000000,
                       ymin = -6000000,  ymax = 8500000)

## ---- colours, opposite of the EVI figure -----------------------------------
# CWD increasing (more deficit): brown
# CWD decreasing (less deficit): teal
cat_colors <- c(
  "CWD_atE95 \u2191 / CWD_atE05 \u2191" = "#3D2E20",   # dark brown: both up = worsening
  "CWD_atE95 \u2193 / CWD_atE05 \u2191" = "#C77F2A",   # orange: mismatch
  "CWD_atE95 \u2193 / CWD_atE05 \u2193" = "#2C7A7B",   # dark teal: both down = improving
  "CWD_atE95 \u2191 / CWD_atE05 \u2193" = "#76C7C0",   # light teal: mismatch
  "CWD_atE95 \u2192 / CWD_atE05 \u2192" = "#D9D9D9"    # gray: flat
)

## ---- functions for the panels ----------------------------------------------
make_map_panel <- function(d, title_text) {
  ggplot() +
    geom_tile(data = d,
              aes(x = X_robin, y = Y_robin, fill = category),
              width = 11000, height = 11000) +
    geom_sf(data = world_robin, fill = NA, color = "gray30", linewidth = 0.2) +
    scale_fill_manual(values = cat_colors, name = "CWD trend direction",
                      drop = FALSE) +
    coord_sf(crs = robin_crs, expand = FALSE) +
    labs(title = title_text) +
    theme_void(base_size = 20) +
    theme(plot.title       = element_text(size = 22, face = "bold", hjust = 0),
          legend.position  = "bottom",
          legend.title     = element_text(size = 20),
          legend.text      = element_text(size = 18),
          legend.key.size  = unit(1.4, "cm"),
          legend.spacing.x = unit(0.6, "cm")) +
    guides(fill = guide_legend(nrow = 1, byrow = TRUE,
                               override.aes = list(size = 8)))
}

make_lat_panel <- function(d) {
  ggplot(d, aes(y = lat, color = category,
                weight = cos(lat * pi / 180))) +
    geom_density(linewidth = 0.9) +
    scale_color_manual(values = cat_colors, guide = "none") +
    scale_y_continuous(limits = c(-60, 80),
                       breaks = seq(-60, 80, 20),
                       labels = function(x) paste0(x, "\u00b0")) +
    scale_x_continuous(limits = c(0, 0.04),
                       breaks = c(0, 0.01, 0.02, 0.03, 0.04),
                       labels = c("0", "0.01", "0.02", "0.03", "0.04"),
                       oob    = scales::squish,
                       expand = c(0, 0)) +
    labs(x = NULL, y = NULL) +
    theme_bw(base_size = 18) +
    theme(panel.grid.minor = element_blank(),
          axis.text.x      = element_text(size = 16),
          axis.text.y      = element_blank(),
          axis.ticks.y     = element_blank())
}

## ---- remove the class without trend ----------------------------------------
flat_label <- "CWD_atE95 \u2192 / CWD_atE05 \u2192"
n_before_flat <- sum(before$category == flat_label, na.rm = TRUE)
n_after_flat  <- sum(after$category  == flat_label, na.rm = TRUE)
cat(sprintf("\nRemoving 'flat' (NA) pixels: %d before, %d after\n",
            n_before_flat, n_after_flat))
before <- before[category != flat_label]
after  <- after[category  != flat_label]
before[, category := droplevels(category)]
after[,  category := droplevels(category)]

## ---- panels ----------------------------------------------------------------
cat("\nBuilding map + lat panels for BEFORE...\n")
p_map_b <- make_map_panel(before, "a) Before regime")
p_lat_b <- make_lat_panel(before)

cat("Building map + lat panels for AFTER...\n")
p_map_a <- make_map_panel(after,  "b) After regime")
p_lat_a <- make_lat_panel(after)

row_before <- p_map_b + p_lat_b + plot_layout(widths = c(5, 1))
row_after  <- p_map_a + p_lat_a + plot_layout(widths = c(5, 1))

fig <- row_before / row_after +
  plot_layout(guides = "collect") +
  plot_annotation(theme = theme(legend.position = "bottom"))

ggsave(file.path(out_dir, "Figure_trend_direction_CWD_by_regime.png"),
       fig, width = 22, height = 18, dpi = 300, bg = "white")
cat("Saved PNG\n")

cat("\nDone.\n")