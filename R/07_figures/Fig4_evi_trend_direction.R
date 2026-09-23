# =============================================================================
# FIGURE 4, FIGURES S3 AND S4 : EVI TRENDS AT THE PEAK AND LOW STATES
# Original file name: Figure_EVI_abs_trends.R
# Methods 2.8
# -----------------------------------------------------------------------------
# Uses the trends of EVI95 and EVI05 (06_trends_theil_sen.R), separately for
# energy-limited and water-limited grid cells.
#
# Fig. 4, direction of the trends. Each grid cell gets a class from the signs
# of its two trends: both increasing, peak decreasing and low increasing, both
# decreasing, or peak increasing and low decreasing. Only the sign is used,
# whether the trend is significant or not. Grid cells with a zero or missing
# trend for one of the two states are removed. One map per regime, with the
# distribution of each class by latitude next to it.
#
# Fig. S4, significance. Mann-Kendall p-value of each trend, in classes, as
# four maps (regime x state).
#
# Fig. S3, size of the trends. The slopes themselves, as four maps, with one
# colour scale centred on zero and limited to +/- COLLIM.
#
# Console: share of each class per regime, with and without weighting by
# cos(latitude), for all grid cells and for grid cells with a significant trend
# (p < 0.05) for at least one state, and the same by zone. These numbers are
# used in the supplementary tables on trend classes.
#
# Arid Hot is removed, by keeping only the grid cells of the pooled analysis.
#
# Input  : sens_slope_EVI_E95_<regime>_v2.Rdata, sens_slope_EVI_E05_<regime>_v2.Rdata
#          combined_pooled_percentiles_summary_v2.Rdata (climate zones)
# Output : Figure_evi_trend_direction_by_regime.png     Fig. 4
#          Figure_EVI_MK_pvalue_maps.png                 Fig. S4
#          FigureS_EVI_trend_magnitude_4maps.png         Fig. S3
#          The file names come from our working folder, see the README.
#
# Run    : Rscript Fig4_evi_trend_direction.R
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
# Labels with subscripts.
EVI95 <- "EVI₉₅"
EVI05 <- "EVI₀₅"
## ---- climate zone of each grid cell, without Arid Hot ----------------------
cat("Loading V2 pooled data for climate-zone filter...\n")
load(file.path(v2_dir, "combined_pooled_percentiles_summary_v2.Rdata"))
setDT(all_summaries)
zone_lookup <- unique(all_summaries[climate_zone != "Arid Hot",
                                    .(lat = Latitude, lon = Longitude, climate_zone)])
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
cat("\nLoading EVI_E95 and EVI_E05 trends...\n")
ev95_b <- load_trend("EVI_E95", "before")
ev95_a <- load_trend("EVI_E95", "after")
ev05_b <- load_trend("EVI_E05", "before")
ev05_a <- load_trend("EVI_E05", "after")
## ---- sign of the trends ----------------------------------------------------
classify <- function(d) {
  d[, sign := fcase(
    is.na(sens_slope), "flat",
    sens_slope > 0,    "up",
    sens_slope < 0,    "down",
    default = "flat"
  )]
  d[, sig := !is.na(sens_pvalue) & sens_pvalue < p_threshold]
  d
}
ev95_b <- classify(ev95_b); ev95_a <- classify(ev95_a)
ev05_b <- classify(ev05_b); ev05_a <- classify(ev05_a)
## ---- world map and projection ----------------------------------------------
robin_crs <- "+proj=robin +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
project_robin <- function(d) {
  pts <- st_as_sf(d, coords = c("lon", "lat"), crs = 4326)
  pts <- st_transform(pts, crs = robin_crs)
  coords <- st_coordinates(pts)
  d[, X_robin := coords[, 1]]
  d[, Y_robin := coords[, 2]]
  d
}
world <- ne_countries(scale = "medium", returnclass = "sf")
world_robin <- st_transform(world, robin_crs)
world_robin <- st_crop(world_robin,
                       xmin = -17000000, xmax = 17000000,
                       ymin = -6000000,  ymax = 8500000)
## ---- latitude and longitude lines and labels -------------------------------
grat_lat <- c(-60, -30, 0, 30, 60)                     # parallels (latitude)
grat_lon <- c(-120, -60, 0, 60, 120)                   # meridians (longitude)
LAT_LO   <- -60; LAT_HI <- 84                          # plotted vertical frame
# vertical limits of the map in Robinson metres (so that -60 and 60 are inside)
yb <- st_coordinates(st_transform(
  st_as_sf(data.frame(lon = 0, lat = c(LAT_LO, LAT_HI)),
           coords = c("lon", "lat"), crs = 4326), robin_crs))[, "Y"]
# lines at grat_lat and grat_lon, cut to the map limits
grat <- st_transform(st_graticule(lon = grat_lon, lat = grat_lat,
                                  crs = st_crs(4326)), robin_crs)
grat <- st_crop(grat, xmin = -17e6, xmax = 17e6, ymin = yb[1], ymax = yb[2])
# latitude labels on the left edge
lat_xy  <- st_coordinates(st_transform(
  st_as_sf(data.frame(lon = -180, lat = grat_lat),
           coords = c("lon", "lat"), crs = 4326), robin_crs))
lat_lab <- data.table(x = lat_xy[, "X"], y = lat_xy[, "Y"],
                      label = paste0(abs(grat_lat), "°",
                                     ifelse(grat_lat > 0, "N",
                                            ifelse(grat_lat < 0, "S", ""))))
# longitude labels along the bottom edge
lon_xy  <- st_coordinates(st_transform(
  st_as_sf(data.frame(lon = grat_lon, lat = LAT_LO),
           coords = c("lon", "lat"), crs = 4326), robin_crs))
lon_lab <- data.table(x = lon_xy[, "X"], y = lon_xy[, "Y"],
                      label = paste0(abs(grat_lon), "°",
                                     ifelse(grat_lon > 0, "E",
                                            ifelse(grat_lon < 0, "W", ""))))
## ---- Fig. 4: trend direction maps ------------------------------------------
# The four classes.
lab_up_up <- paste0(EVI95, " ↑ / ", EVI05, " ↑")
lab_dn_up <- paste0(EVI95, " ↓ / ", EVI05, " ↑")
lab_dn_dn <- paste0(EVI95, " ↓ / ", EVI05, " ↓")
lab_up_dn <- paste0(EVI95, " ↑ / ", EVI05, " ↓")
lab_flat  <- paste0(EVI95, " → / ", EVI05, " →")
build_regime <- function(d95, d05, regime_label) {
  m <- merge(d95[, .(lat, lon, sign_E95 = sign, sig_E95 = sig, p_E95 = sens_pvalue)],
             d05[, .(lat, lon, sign_E05 = sign, sig_E05 = sig, p_E05 = sens_pvalue)],
             by = c("lat", "lon"))
  m <- merge(m, zone_lookup, by = c("lat", "lon"))
  m[, regime := regime_label]
  m[, category := fcase(
    sign_E95 == "up"   & sign_E05 == "up",   lab_up_up,
    sign_E95 == "down" & sign_E05 == "up",   lab_dn_up,
    sign_E95 == "down" & sign_E05 == "down", lab_dn_dn,
    sign_E95 == "up"   & sign_E05 == "down", lab_up_dn,
    default = lab_flat
  )]
  m[, both_sig   := sig_E95 & sig_E05]
  m[, either_sig := sig_E95 | sig_E05]
  m
}
cat("\nBuilding regime datasets...\n")
before <- build_regime(ev95_b, ev05_b, "Energy-limited")
after  <- build_regime(ev95_a, ev05_a, "Water-limited")
cat("Energy-limited pixels:", nrow(before), "| Water-limited pixels:", nrow(after), "\n")
# Remove grid cells with a zero or missing trend for one of the states.
before <- before[category != lab_flat]
after  <- after[category  != lab_flat]
cat_levels <- c(lab_up_up, lab_dn_up, lab_dn_dn, lab_up_dn)
before[, category := factor(category, levels = cat_levels)]
after[,  category := factor(category, levels = cat_levels)]
# Share of each class, without weights.
cat("\n=== Category distribution — ENERGY-LIMITED ===\n")
print(before[, .(N = .N, pct = round(100 * .N / nrow(before), 1)), by = category][order(category)])
cat("\n=== Category distribution — WATER-LIMITED ===\n")
print(after[, .(N = .N, pct = round(100 * .N / nrow(after), 1)), by = category][order(category)])
cat(sprintf("\nEither significant — Energy-limited: %.1f%% | Water-limited: %.1f%%\n",
            100 * mean(before$either_sig, na.rm = TRUE),
            100 * mean(after$either_sig,  na.rm = TRUE)))
before <- project_robin(before)
after  <- project_robin(after)
cat_colors <- c("#2C7A7B", "#76C7C0", "#3D2E20", "#C77F2A")  # up/up, dn/up, dn/dn, up/dn
names(cat_colors) <- cat_levels
make_map_panel <- function(d, title_text) {
  ggplot() +
    geom_tile(data = d, aes(x = X_robin, y = Y_robin, fill = category),
              width = 11000, height = 11000) +
    geom_sf(data = grat, color = "grey85", linewidth = 0.2) +          # our own grid lines
    geom_sf(data = world_robin, fill = NA, color = "gray30", linewidth = 0.2) +
    geom_text(data = lat_lab, aes(x = x, y = y, label = label), inherit.aes = FALSE,
              hjust = 1, nudge_x = -4e5, size = 6, colour = "grey45") +   # latitude labels (left)
    geom_text(data = lon_lab, aes(x = x, y = y, label = label), inherit.aes = FALSE,
              vjust = 1, nudge_y = -4e5, size = 6, colour = "grey45") +   # longitude labels (bottom)
    scale_fill_manual(values = cat_colors, name = "EVI trend direction", drop = FALSE) +
    coord_sf(crs = robin_crs, ylim = yb, expand = FALSE, clip = "off") +
    labs(title = title_text) +
    theme_void(base_size = 24) +
    theme(plot.title       = element_text(size = 26, hjust = 0),
          plot.margin      = margin(6, 6, 22, 34),          # room for bottom/left labels
          legend.position  = "bottom",
          legend.title     = element_text(size = 26),
          legend.text      = element_text(size = 24),
          legend.key.size  = unit(1.4, "cm"),
          legend.spacing.x = unit(0.6, "cm")) +
    guides(fill = guide_legend(nrow = 1, byrow = TRUE, override.aes = list(size = 8)))
}
make_lat_panel <- function(d) {
  ggplot(d, aes(y = lat, color = category, weight = cos(lat * pi / 180))) +
    geom_density(linewidth = 0.9) +
    scale_color_manual(values = cat_colors, guide = "none") +
    scale_y_continuous(limits = c(-60, 80), breaks = seq(-60, 80, 20),
                       labels = function(x) paste0(x, "°")) +
    scale_x_continuous(limits = c(0, 0.05), breaks = seq(0, 0.05, 0.01),
                       labels = c("0","0.01","0.02","0.03","0.04","0.05"),
                       oob = scales::squish,
                       expand = expansion(mult = c(0, 0.04))) +   # a hair of room so "0.05" isn't clipped
    labs(x = NULL, y = NULL) +
    theme_bw(base_size = 22) +
    theme(panel.grid.minor = element_blank(),
          axis.text.x  = element_text(size = 22),
          axis.text.y  = element_blank(),
          axis.ticks.y = element_blank(),
          plot.margin  = margin(6, 20, 6, 2))                     # right margin for the last x label
}
cat("\nBuilding FIGURE 1...\n")
row_before <- make_map_panel(before, "a) Energy-limited") + make_lat_panel(before) +
  plot_layout(widths = c(4.5, 1))
row_after  <- make_map_panel(after,  "b) Water-limited")  + make_lat_panel(after) +
  plot_layout(widths = c(4.5, 1))
fig1 <- row_before / row_after +
  plot_layout(guides = "collect") +
  plot_annotation(theme = theme(legend.position = "bottom"))
ggsave(file.path(out_dir, "Figure_evi_trend_direction_by_regime.png"),
       fig1, width = 22, height = 18, dpi = 300, bg = "white")
cat("Saved FIGURE 1.\n")
## ---- Fig. S4: Mann-Kendall p-values, regime x state ------------------------
cat("\nBuilding FIGURE 2 (p-value maps)...\n")
mk <- function(d, state_label, regime_label) {
  d[, .(lat, lon, pval = sens_pvalue, state = state_label, regime = regime_label)]
}
dt_p <- rbindlist(list(
  mk(ev95_b, EVI95, "Energy-limited"),
  mk(ev05_b, EVI05, "Energy-limited"),
  mk(ev95_a, EVI95, "Water-limited"),
  mk(ev05_a, EVI05, "Water-limited")
), use.names = TRUE)
dt_p <- merge(dt_p, zone_lookup[, .(lat, lon)], by = c("lat", "lon"))
dt_p <- dt_p[!is.na(pval)]
dt_p[, pbin := cut(pval, breaks = c(-Inf, 0.01, 0.05, 0.10, Inf),
                   labels = c("p < 0.01", "0.01 ≤ p < 0.05",
                              "0.05 ≤ p < 0.10", "p ≥ 0.10"), right = FALSE)]
dt_p[, state  := factor(state,  levels = c(EVI95, EVI05))]
dt_p[, regime := factor(regime, levels = c("Energy-limited", "Water-limited"))]
dt_p <- project_robin(dt_p)
# Labels shown on some panels only:
#   latitude labels in the left column (state == EVI95), inside the map with
#     a white box (geom_label) so they can be read over the colours;
#   longitude labels in the bottom row (regime == "Water-limited").
lat_lab2 <- copy(lat_lab)[, state  := factor(EVI95, levels = c(EVI95, EVI05))]
lon_lab2 <- copy(lon_lab)[, regime := factor("Water-limited",
                                             levels = c("Energy-limited", "Water-limited"))]
# Viridis colours for the significant classes, light grey for p >= 0.10.
p_colors <- c(
  "p < 0.01"        = "#440154",
  "0.01 ≤ p < 0.05" = "#31688E",
  "0.05 ≤ p < 0.10" = "#35B779",
  "p ≥ 0.10"        = "#E8E8E8"
)
fig2 <- ggplot() +
  geom_tile(data = dt_p, aes(x = X_robin, y = Y_robin, fill = pbin),
            width = 11000, height = 11000) +
  geom_sf(data = grat, color = "grey80", linewidth = 0.15) +           # grid lines (every panel)
  geom_sf(data = world_robin, fill = NA, color = "gray30", linewidth = 0.15) +
  geom_label(data = lat_lab2, aes(x = x, y = y, label = label), inherit.aes = FALSE,
             hjust = 0, nudge_x = 1.5e5, size = 4, colour = "grey25",
             fill = "white", label.size = 0, label.padding = unit(0.08, "lines")) +  # lat (left col)
  geom_text(data = lon_lab2, aes(x = x, y = y, label = label), inherit.aes = FALSE,
            vjust = 1, nudge_y = -3e5, size = 4.5, colour = "grey45") +               # lon (bottom row)
  scale_fill_manual(values = p_colors, name = "Mann–Kendall p", drop = FALSE) +
  coord_sf(crs = robin_crs, ylim = yb, expand = FALSE, clip = "off") +
  facet_grid(regime ~ state, switch = "y") +
  theme_void(base_size = 22) +
  theme(strip.text      = element_text(size = 24),
        legend.position = "bottom",
        legend.title    = element_text(size = 24),
        legend.text     = element_text(size = 22),
        legend.key.size = unit(1.1, "cm"),
        panel.spacing   = unit(0.9, "lines"),
        plot.margin     = margin(6, 8, 20, 8)) +               # room for bottom lon labels
  guides(fill = guide_legend(nrow = 1, override.aes = list(size = 6)))
ggsave(file.path(out_dir, "Figure_EVI_MK_pvalue_maps.png"),
       fig2, width = 20, height = 12, dpi = 300, bg = "white")
cat("Saved FIGURE 2.\n")
## ---- Fig. S3: slopes, regime x state ---------------------------------------
# Rows: regime (energy- or water-limited). Columns: state (EVI95 or EVI05).
# Colour: Theil-Sen slope, one scale centred on 0 for all maps (brown for a
# decrease, teal for an increase).
cat("\nBuilding FIGURE S (4-map signed-trend layout)...\n")
# One table per state, with both regimes.
build_slope <- function(d_before, d_after) {
  s <- rbindlist(list(
    d_before[, .(lat, lon, slope = sens_slope, regime = "Energy-limited")],
    d_after[,  .(lat, lon, slope = sens_slope, regime = "Water-limited")]
  ))
  s <- merge(s, zone_lookup[, .(lat, lon)], by = c("lat", "lon"))   # drop Arid Hot
  s <- s[!is.na(slope)]                                             # keep computable trends
  project_robin(s)
}
sl_E95 <- build_slope(ev95_b, ev95_a)[, state := EVI95]
sl_E05 <- build_slope(ev05_b, ev05_a)[, state := EVI05]
sl_all <- rbindlist(list(sl_E95, sl_E05))
sl_all[, state  := factor(state,  levels = c(EVI95, EVI05))]
sl_all[, regime := factor(regime, levels = c("Energy-limited", "Water-limited"))]
# ---- colour limit, EVI per year ----------------------------------------------
# Most slopes are close to zero, with a few extreme values. A fixed limit
# (extreme values shown at the limit) shows the patterns better.
qs <- quantile(abs(sl_all$slope), c(0.90, 0.95, 0.98), na.rm = TRUE)
cat(sprintf("|slope| percentiles: p90=%.4f  p95=%.4f  p98=%.4f  (EVI/yr)\n",
            qs[1], qs[2], qs[3]))
COLLIM  <- 0.01     # symmetric colour limit, EVI per year
lim_all <- COLLIM
figS <- ggplot() +
  geom_tile(data = sl_all, aes(x = X_robin, y = Y_robin, fill = slope),
            width = 11000, height = 11000) +
  geom_sf(data = world_robin, fill = NA, color = "gray30", linewidth = 0.15) +
  scale_fill_gradient2(low = "#8C510A", mid = "grey95", high = "#01665E",
                       midpoint = 0, limits = c(-lim_all, lim_all),
                       breaks = c(-1, -0.5, 0, 0.5, 1) * lim_all,
                       oob = scales::squish, name = "EVI trend (EVI/yr)") +
  coord_sf(crs = robin_crs, ylim = yb, expand = FALSE) +   # no graticule / no lon-lat labels
  facet_grid(regime ~ state, switch = "y") +       # rows = regime (labels on LEFT), cols = state
  theme_void(base_size = 22) +
  theme(strip.text         = element_text(size = 24),   # horizontal, incl. the left regime labels
        legend.position    = "bottom",
        legend.title       = element_text(size = 22),
        legend.text        = element_text(size = 18),
        legend.key.width   = unit(3, "cm"),
        legend.key.height  = unit(0.7, "cm"),
        legend.box.spacing = unit(0.4, "cm"),
        panel.spacing.x    = unit(1.2, "lines"),
        panel.spacing.y    = unit(0.8, "lines"),
        plot.margin        = margin(8, 10, 10, 10)) +
  guides(fill = guide_colourbar(title.position = "top", title.hjust = 0.5))
ggsave(file.path(out_dir, "FigureS_EVI_trend_magnitude_4maps.png"),
       figS, width = 22, height = 13, dpi = 300, bg = "white")
cat("Saved FIGURE S: 4-map signed-trend layout.\n\nDone.\n")
aw_sig_class <- function(d) {
  d <- copy(d)[either_sig == TRUE]; d[, w := cos(lat*pi/180)]
  tot <- d[, sum(w)]
  d[, .(pct = round(100*sum(w)/tot,1)), by = category][order(category)]
}
cat("EL significant, area-weighted:\n"); print(aw_sig_class(before))
cat("WL significant, area-weighted:\n"); print(aw_sig_class(after))
## ---- share of each class, weighted by cos(latitude), printed only ----------
# All zones, per regime.
area_frac <- function(d) {
  d <- copy(d); d[, w := cos(lat * pi/180)]
  tot <- d[, sum(w)]
  d[, .(pct = round(100 * sum(w) / tot, 1)), by = category][order(category)]
}
cat("\n== Energy-limited — area-weighted class % ==\n"); print(area_frac(before))
cat("\n== Water-limited — area-weighted class % ==\n"); print(area_frac(after))
# Per climate zone.
area_frac_zone <- function(d) {
  d <- copy(d); d[, w := cos(lat * pi/180)]
  d[, .(w = sum(w)), by = .(climate_zone, category)][
    , pct := round(100 * w / sum(w), 1), by = climate_zone][order(climate_zone, category)]
}
cat("\n== Energy-limited — area-weighted % par zone ==\n"); print(area_frac_zone(before))
cat("\n== Water-limited — area-weighted % par zone ==\n"); print(area_frac_zone(after))
# Share of grid cells with a significant trend for at least one state,
# weighted by cos(latitude).
aw_sig <- function(d) { w <- cos(d$lat*pi/180); round(100*sum(w[d$either_sig], na.rm=TRUE)/sum(w),1) }
cat(sprintf("\nEither-sig (area-weighted): EL %.1f%% | WL %.1f%%\n", aw_sig(before), aw_sig(after)))