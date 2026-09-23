# =============================================================================
# FIGURE 2 : PEAK AND LOW STATES BY CLIMATE ZONE AND REGIME
# Original file name: Figure_boxplot_thresholds_distribution.R
# Methods 2.6
# -----------------------------------------------------------------------------
# Boxplots of the pooled states in the six climate zones, for energy-limited
# (EL) and water-limited (WL) grid cells:
#   a) EVI95        b) EVI05        c) dEVI (relative change, with sign)
#   d) CWD_EVI95    e) CWD_EVI05    f) dCWD
# and g) a map of the regimes. We also save, for the supplement, a panel with
# the share of each regime per zone and a summary table by zone and regime.
#
# How dEVI and dCWD are computed (the same in all figure scripts):
#   dEVI = 100 * (EVI05 - EVI95) / EVI95  in WL grid cells (negative: EVI
#          decreases as the deficit grows)
#   dEVI = 100 * (EVI95 - EVI05) / EVI95  in EL grid cells (positive: EVI
#          increases as the deficit grows)
#   dCWD = |CWD_EVI05 - CWD_EVI95|, in mm
# The regime comes from p05_position: "before" is EL, "after" is WL.
#
# We keep grid cells with both states and a regime ("equal" is removed).
# Arid Hot is removed from all figures and tables, because the yearly CWD
# reset at the wettest month does not work well in this zone.
#
# The output file names come from our working folder (Figure3_V2_*). This is
# Fig. 2 of the paper.
#
# Input  : combined_pooled_percentiles_summary_v2.Rdata
# Output : Figure3_V2_all_split_with_map_noAridHot.png       Fig. 2
#          Figure3_V2_panel_gh_only_noAridHot.png            supplement
#          Figure3_V2_all_split_summary_noAridHot.csv        supplement
#
# Run    : Rscript Fig2_states_by_zone_and_regime.R
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
library(sf)
library(rnaturalearth)
## ---- paths (change these lines to run the script elsewhere) ----------------
input_path  <- "/Net/Groups/BGI/scratch/mterristi/PhD/CWD_0_1/New_Climate_Zones_percentiles/Drydowns/Drydowns_percentiles_pooled_v2"
output_path <- "/Net/Groups/BGI/scratch/mterristi/PhD/CWD_0_1/New_Climate_Zones_percentiles/Drydowns/Figures"
dir.create(output_path, showWarnings = FALSE, recursive = TRUE)
## ---- pooled states, with the filters of the header -------------------------
cat("Loading V2 pooled data...\n")
load(file.path(input_path, "combined_pooled_percentiles_summary_v2.Rdata"))
setDT(all_summaries)
all_summaries <- all_summaries[!is.na(EVI_E05) & !is.na(CWD_atE05) &
                                 p05_position %in% c("before", "after")]
n_before_excl <- nrow(all_summaries)
all_summaries <- all_summaries[climate_zone != "Arid Hot"]
cat(sprintf("Pixels after Arid Hot exclusion: %d (removed %d, %.1f%%)\n",
            nrow(all_summaries), n_before_excl - nrow(all_summaries),
            100 * (n_before_excl - nrow(all_summaries)) / n_before_excl))
## ---- dEVI and dCWD (see header) --------------------------------------------
all_summaries[, Delta_EVI := 100 * fifelse(
  p05_position == "after",
  (EVI_E05 - EVI_E95) / EVI_E95,
  (EVI_E95 - EVI_E05) / EVI_E95)]
all_summaries[, Delta_CWD := abs(CWD_atE05 - CWD_atE95)]
## ---- order of the zones and names of the regimes ---------------------------
zone_levels <- c("Tropical", "Arid Cold", "Temperate Dry",
                 "Temperate Humid", "Cold Dry", "Cold Humid")
all_summaries[, climate_zone := factor(climate_zone, levels = zone_levels)]
all_summaries[, regime := factor(fifelse(p05_position == "before", "EL", "WL"),
                                 levels = c("EL", "WL"))]
print(all_summaries[, .N, by = regime])
col_EL <- "#2E7D32"; col_WL <- "#8B4513"
regime_colors <- c("EL" = col_EL, "WL" = col_WL)
# One colour scale for all panels, so that all legends are the same.
# key_glyph="rect" is set in the geoms below, not in the scale.
regime_scale <- scale_fill_manual(values = regime_colors, name = "Regime",
                                  breaks = c("EL", "WL"), limits = c("EL", "WL"))
## ---- theme -----------------------------------------------------------------
BASE <- 30
panel_theme <- theme_bw(base_size = BASE) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(color = "grey90"),
    panel.grid.minor   = element_blank(),
    axis.text.x        = element_text(angle = 32, hjust = 1, size = 28),
    axis.text.y        = element_text(size = 28),
    axis.title.y       = element_text(size = 30),
    plot.title         = element_text(size = 30, face = "plain", hjust = 0),
    legend.position    = "bottom",
    legend.direction   = "horizontal",
    legend.title       = element_text(size = 30),
    legend.text        = element_text(size = 30),
    legend.key.size    = unit(1.5, "cm")
  )
make_split_panel <- function(d, yvar, ytitle, title_text, ylim = NULL,
                             hline_at = NULL) {
  p <- ggplot(d, aes(x = climate_zone, y = .data[[yvar]], fill = regime)) +
    geom_boxplot(outlier.size = 0.4, outlier.alpha = 0.5,
                 position = position_dodge(width = 0.75), width = 0.65,
                 key_glyph = "rect") +                       # <- identical key
    regime_scale +
    labs(title = title_text, x = NULL, y = ytitle) +
    guides(fill = "none") +
    panel_theme
  if (!is.null(ylim))     p <- p + coord_cartesian(ylim = ylim)
  if (!is.null(hline_at)) p <- p + geom_hline(yintercept = hline_at,
                                              linetype = "dashed", linewidth = 0.4)
  p
}
## ---- panels a to f ---------------------------------------------------------
p_a <- make_split_panel(all_summaries, "EVI_E95", "EVI",
                        expression("a) " * EVI[95]),     ylim = c(0, 1))
p_b <- make_split_panel(all_summaries, "EVI_E05", "EVI",
                        expression("b) " * EVI["05"]),     ylim = c(0, 1))
p_c <- make_split_panel(all_summaries, "Delta_EVI", expression(Delta * EVI),
                        expression("c) " * Delta * EVI), ylim = c(-100, 100),
                        hline_at = 0) +
  scale_y_continuous(labels = function(x) paste0(x, "%"))
p_d <- make_split_panel(all_summaries, "CWD_atE95", "CWD (mm)",
                        expression("d) " * CWD[EVI95]),   ylim = c(0, 425))
p_e <- make_split_panel(all_summaries, "CWD_atE05", "CWD (mm)",
                        expression("e) " * CWD["EVI05"]),   ylim = c(0, 425))
p_f <- make_split_panel(all_summaries, "Delta_CWD", expression(Delta * CWD ~ "(mm)"),
                        expression("f) " * Delta * CWD),  ylim = c(0, 425))
## ---- panel g: regime map, Robinson projection ------------------------------
# The latitude and longitude lines and their labels are drawn by hand, from
# the same values, so they always match. Same method as in
# Fig6_regime_shifts_extreme_years.R.
cat("\nBuilding panel g (regime map)...\n")
robin_crs <- "+proj=robin +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
proj_pts <- st_transform(st_as_sf(all_summaries,
                                  coords = c("Longitude", "Latitude"), crs = 4326),
                         crs = robin_crs)
coords <- st_coordinates(proj_pts)
all_summaries[, X_robin := coords[, 1]]; all_summaries[, Y_robin := coords[, 2]]
world_robin <- st_crop(
  st_transform(ne_countries(scale = "medium", returnclass = "sf"), robin_crs),
  xmin = -17000000, xmax = 17000000, ymin = -6000000, ymax = 8500000)
# ---- graticule ---------------------------------------------------------------
grat_lat <- c(-60, -30, 0, 30, 60)                     # parallels (latitude)
grat_lon <- c(-180, -120, -60, 0, 60, 120, 180)        # meridians (longitude)
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
p_g <- ggplot() +
  geom_tile(data = all_summaries,
            aes(x = X_robin, y = Y_robin, fill = regime),
            width = 11000, height = 11000, key_glyph = "rect") +   # <- identical key
  geom_sf(data = grat, color = "grey85", linewidth = 0.2) +        # our own grid lines
  geom_sf(data = world_robin, fill = NA, color = "gray30", linewidth = 0.2) +
  geom_text(data = lat_lab, aes(x = x, y = y, label = label), inherit.aes = FALSE,
            hjust = 1, nudge_x = -4e5, size = 12, colour = "grey45") +  # latitude labels (left)
  geom_text(data = lon_lab, aes(x = x, y = y, label = label), inherit.aes = FALSE,
            vjust = 1, nudge_y = -4e5, size = 12, colour = "grey45") +  # longitude labels (bottom)
  regime_scale +
  coord_sf(crs = robin_crs, ylim = yb, expand = FALSE, clip = "off") +
  labs(title = "g) Regime classification") +
  theme_void(base_size = BASE) +
  theme(plot.title = element_text(size = 30, hjust = 0),
        plot.margin = margin(6, 12, 26, 34),           # room for bottom/left labels
        legend.position  = "bottom",
        legend.direction = "horizontal",
        legend.title = element_text(size = 30),
        legend.text  = element_text(size = 30),
        legend.key.size = unit(1.3, "cm"))
## ---- put the panels together, one legend at the bottom ---------------------
fig <- (p_a | p_b | p_c) /
  (p_d | p_e | p_f) /
  p_g /
  guide_area() +
  plot_layout(heights = c(1, 1, 1.5, 0.18), guides = "collect")
ggsave(file.path(output_path, "Figure3_V2_all_split_with_map_noAridHot.png"),
       fig, width = 32, height = 28, dpi = 300, bg = "white", limitsize = FALSE)
cat("Saved: Figure3_V2_all_split_with_map_noAridHot.png\n")
## ---- panel h: share of each regime per zone, for the supplement ------------
zone_summary <- all_summaries[, .N, by = .(climate_zone, regime)]
zone_summary[, pct := 100 * N / sum(N), by = climate_zone]
zone_summary[, label := sprintf("%.1f%%", pct)]
zone_summary[, regime := factor(regime, levels = c("EL", "WL"))]
zone_summary[, climate_zone := factor(climate_zone, levels = rev(zone_levels))]
p_h <- ggplot(zone_summary, aes(x = pct, y = climate_zone, fill = regime)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7) +
  geom_text(aes(label = label), position = position_stack(vjust = 0.5),
            color = "white", size = 7, fontface = "bold") +
  regime_scale +
  guides(fill = "none") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.02)),
                     breaks = seq(0, 100, 25), labels = function(x) paste0(x, "%")) +
  labs(title = " % of pixels per regime, per climate zone",
       x = "% of pixels", y = NULL) +
  theme_bw(base_size = BASE) +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_line(color = "grey90"),
        panel.grid.minor   = element_blank(),
        axis.text.y = element_text(size = 22), axis.text.x = element_text(size = 20),
        plot.title  = element_text(size = 30, face = "plain", hjust = 0))

fig_gh <-  p_h / guide_area() +
  plot_layout(heights = c(1.6, 0.8, 0.15), guides = "collect")
ggsave(file.path(output_path, "Figure3_V2_panel_gh_only_noAridHot.png"),
       fig_gh, width = 20, height = 16, dpi = 200, bg = "white", limitsize = FALSE)
cat("Saved: Figure3_V2_panel_gh_only_noAridHot.png\n")
## ---- summary table by zone and regime, for the supplement ------------------
stats <- all_summaries[, .(
  n = .N,
  EVI_95    = round(median(EVI_E95,   na.rm = TRUE), 3),
  EVI_05    = round(median(EVI_E05,   na.rm = TRUE), 3),
  CWD_EVI95 = round(median(CWD_atE95, na.rm = TRUE), 1),
  CWD_EVI05 = round(median(CWD_atE05, na.rm = TRUE), 1),
  Delta_EVI = round(median(Delta_EVI, na.rm = TRUE), 3),
  Delta_CWD = round(median(Delta_CWD, na.rm = TRUE), 1)
), by = .(climate_zone, regime)]
print(stats)
fwrite(stats, file.path(output_path, "Figure3_V2_all_split_summary_noAridHot.csv"))
cat("\nDone.\n")