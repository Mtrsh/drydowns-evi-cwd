# =============================================================================
# SUPPLEMENTARY MAPS : TRENDS OF THE FOUR STATES, PER REGIME
# Original file name: Maps_trends.R
# Methods 2.8
# -----------------------------------------------------------------------------
# Maps of the significant trends (p < 0.05) of EVI95, CWD_EVI95, EVI05 and
# CWD_EVI05, as four maps per regime and as single maps. Slopes are converted
# to values per decade, with fixed colour scales: +/- 0.05 EVI per decade and
# +/- 50 mm per decade (higher values are shown at the limit). Increasing EVI
# is green and decreasing EVI is brown; increasing CWD is red and decreasing
# CWD is blue.
#
# Input  : sens_slope_<metric>_<regime>_v2.Rdata
# Output : Maps_4panel_trends_after.png, Maps_4panel_trends_before.png, and one
#          map per variable and regime
#          The file names come from our working folder, see the README.
#
# Run    : Rscript FigS_trend_slope_maps.R
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
library(scales)
library(viridis)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)

## ---- paths (change these lines to run the script elsewhere) ----------------
input_path  <- "/Net/Groups/BGI/scratch/mterristi/PhD/CWD_0_1/New_Climate_Zones_percentiles/Trends/Trends_V2"
output_path <- file.path(input_path, "Figures")
dir.create(output_path, showWarnings = FALSE, recursive = TRUE)

## ---- Robinson projection ---------------------------------------------------
robin_crs <- "+proj=robin +datum=WGS84 +no_defs"

cat("Loading country borders...\n")
world <- ne_countries(scale = "medium", returnclass = "sf")
world_robin <- st_transform(world, robin_crs)
world_robin <- st_crop(world_robin,
                       st_transform(st_as_sfc(st_bbox(c(xmin = -180, xmax = 180,
                                                        ymin = -60, ymax = 90),
                                                      crs = 4326)),
                                    robin_crs))

project_pixels <- function(dt, value_col) {
  setDT(dt)
  # Find the names of the lat and lon columns
  lon_col <- intersect(c("Longitude", "longitude", "lon", "x"), names(dt))[1]
  lat_col <- intersect(c("Latitude",  "latitude",  "lat", "y"), names(dt))[1]
  if (is.na(lon_col) || is.na(lat_col)) {
    stop("Could not find lat/lon columns. Available: ", paste(names(dt), collapse = ", "))
  }
  d <- dt[!is.na(get(value_col)) & !is.na(get(lon_col)) & !is.na(get(lat_col))]
  pts <- st_as_sf(d, coords = c(lon_col, lat_col), crs = 4326)
  pts <- st_transform(pts, robin_crs)
  coords <- st_coordinates(pts)
  data.table(x = coords[, 1], y = coords[, 2], val = d[[value_col]])
}

map_theme <- theme_void(base_size = 22) +
  theme(
    legend.position   = "right",
    legend.title      = element_text(size = 22),
    legend.text       = element_text(size = 20),
    legend.key.height = unit(1.3, "cm"),
    legend.key.width  = unit(0.55, "cm"),
    plot.title        = element_text(size = 20, hjust = 0,
                                     margin = margin(0, 0, 10, 0)),
    plot.margin       = margin(10, 10, 10, 10)
  )

## ---- function for the trend maps -------------------------------------------
make_trend_map <- function(dt_proj, fill_label, palette_type = "EVI",
                           limits = NULL, title = "") {
  
  if (palette_type == "EVI") {
    # brown (browning) -> white -> green (greening)
    fill_scale <- scale_fill_gradient2(
      name      = fill_label,
      low       = "#8B4513",
      mid       = "white",
      high      = "#2E7D32",
      midpoint  = 0,
      limits    = limits,
      oob       = squish
    )
  } else {
    # blue (wetting) -> white -> red (drying)
    fill_scale <- scale_fill_gradient2(
      name      = fill_label,
      low       = "#1976D2",
      mid       = "white",
      high      = "#C62828",
      midpoint  = 0,
      limits    = limits,
      oob       = squish
    )
  }
  
  ggplot() +
    geom_tile(data = dt_proj, aes(x = x, y = y, fill = val),
              width = 11000, height = 11000) +
    geom_sf(data = world_robin, fill = NA, color = "grey40", linewidth = 0.2) +
    fill_scale +
    coord_sf(crs = robin_crs, expand = FALSE) +
    labs(title = title) +
    map_theme
}

## ---- trends: significant only, per decade ----------------------------------
load_trend <- function(metric, regime) {
  fp <- file.path(input_path, paste0("sens_slope_", metric, "_", regime, "_v2.Rdata"))
  cat("  Loading", basename(fp), "...\n")
  e <- new.env()
  load(fp, envir = e)
  d <- e$metric_result
  setDT(d)
  if (metric == "EVI_E95" && regime == "after") {
    cat("    First file structure:\n")
    cat("    columns:", paste(names(d), collapse = ", "), "\n")
    cat("    nrow:", nrow(d), "\n")
  }
  
  # From per year to per decade.
  d[, slope_decade := sens_slope * 10]
  
  # Significant trends only.
  d_sig <- d[!is.na(sens_pvalue) & sens_pvalue < 0.05]
  cat("    Total:", nrow(d), "| Significant (p<0.05):", nrow(d_sig), "\n")
  d_sig
}

metrics <- c("EVI_E95", "CWD_atE95", "EVI_E05", "CWD_atE05")
regimes <- c("after", "before")

trend_data <- list()
for (m in metrics) {
  for (r in regimes) {
    key <- paste0(m, "_", r)
    cat("Processing", key, "\n")
    trend_data[[key]] <- load_trend(m, r)
  }
}

## ---- colour limits ---------------------------------------------------------
# One legend per variable, the same for both regimes.
get_lims <- function(metric) {
  combined <- rbind(trend_data[[paste0(metric, "_after")]],
                    trend_data[[paste0(metric, "_before")]])
  v <- combined$slope_decade
  q <- quantile(v, c(0.02, 0.98), na.rm = TRUE)
  m <- max(abs(q))
  c(-m, m)
}

# Fixed limits (see header).
lim_EVI_E95   <- c(-0.05, 0.05)  # EVI/decade
lim_EVI_E05   <- c(-0.05, 0.05)
lim_CWD_atE95 <- c(-50, 50)      # mm/decade
lim_CWD_atE05 <- c(-50, 50)

cat("\nLimits used (capped for visibility):\n")
cat("  EVI_E95, EVI_E05: ±0.05 / decade\n")
cat("  CWD_atE95, CWD_atE05: ±50 mm / decade\n")

## ---- maps per regime -------------------------------------------------------
build_regime_panel <- function(regime) {
  cat("\nBuilding maps for regime:", regime, "\n")
  
  cat("  Projecting points...\n")
  p_E95   <- project_pixels(trend_data[[paste0("EVI_E95_",   regime)]], "slope_decade")
  p_atE95 <- project_pixels(trend_data[[paste0("CWD_atE95_", regime)]], "slope_decade")
  p_E05   <- project_pixels(trend_data[[paste0("EVI_E05_",   regime)]], "slope_decade")
  p_atE05 <- project_pixels(trend_data[[paste0("CWD_atE05_", regime)]], "slope_decade")
  
  cat("  Building maps...\n")
  m1 <- make_trend_map(p_E95,   "ΔEVI_E95\n(decade⁻¹)",   "EVI", lim_EVI_E95,
                       paste0("a)  Sen's slope — EVI_E95 (", regime, ")"))
  m2 <- make_trend_map(p_atE95, "ΔCWD_atE95\n(mm/decade)", "CWD", lim_CWD_atE95,
                       paste0("b)  Sen's slope — CWD_atE95 (", regime, ")"))
  m3 <- make_trend_map(p_E05,   "ΔEVI_E05\n(decade⁻¹)",   "EVI", lim_EVI_E05,
                       paste0("c)  Sen's slope — EVI_E05 (", regime, ")"))
  m4 <- make_trend_map(p_atE05, "ΔCWD_atE05\n(mm/decade)", "CWD", lim_CWD_atE05,
                       paste0("d)  Sen's slope — CWD_atE05 (", regime, ")"))
  
  # Save the single maps
  ggsave(file.path(output_path, paste0("Map_trend_EVI_E95_",   regime, ".png")),
         m1, width = 18, height = 9, dpi = 200, bg = "white")
  ggsave(file.path(output_path, paste0("Map_trend_CWD_atE95_", regime, ".png")),
         m2, width = 18, height = 9, dpi = 200, bg = "white")
  ggsave(file.path(output_path, paste0("Map_trend_EVI_E05_",   regime, ".png")),
         m3, width = 18, height = 9, dpi = 200, bg = "white")
  ggsave(file.path(output_path, paste0("Map_trend_CWD_atE05_", regime, ".png")),
         m4, width = 18, height = 9, dpi = 200, bg = "white")
  
  # Four maps on one figure
  p_combined <- (m1 | m2) / (m3 | m4)
  ggsave(file.path(output_path, paste0("Maps_4panel_trends_", regime, ".png")),
         p_combined, width = 32, height = 16, dpi = 200, bg = "white")
  
  cat("  Saved 4-panel:", paste0("Maps_4panel_trends_", regime, ".png"), "\n")
}

## ---- run -------------------------------------------------------------------
build_regime_panel("after")
build_regime_panel("before")

cat("\nAll done! Maps saved to:", output_path, "\n")