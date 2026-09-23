# =============================================================================
# SUPPLEMENTARY GLOBAL MAPS : STATES, dEVI, dCWD AND REGIME
# Original file name: Figure_global_maps_thresholds.R
# Methods 2.6
# -----------------------------------------------------------------------------
# Maps (Robinson projection) of the pooled variables of Fig. 2, one map each:
# the four states (EVI95, EVI05, CWD_EVI95, CWD_EVI05), dEVI, dCWD and the
# regime. Then, printed in the console (not saved), the numbers cited in Sect. 3.2
# and a split of the variance of log S into the part due to dEVI and the part
# due to dCWD. The last part maps dEVI and dCWD from the yearly states for
# 2010, 2015 and 2018, one row per year.
#
# dEVI, dCWD and the filters are the same as in the other figure scripts
# (both states, EVI95 > 0, a regime, Arid Hot removed).
#
# Input  : combined_pooled_percentiles_summary_v2.Rdata
#          combined_per_year_percentiles_summary_v2.Rdata (last part)
# Output : Map_EVI_E95.png, Map_EVI_E05.png, Map_CWD_atE95.png, Map_CWD_atE05.png
#          Maps_4panel_thresholds_v2.png
#          Map_Delta_EVI_relative.png, Map_Delta_CWD.png, Maps_Delta_2panel_relative.png
#          Map_regime.png
#          Maps_Delta_allyears_panel.png (the three years)
#          The file names come from our working folder, see the README.
#
# Run    : Rscript FigS_global_maps.R
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

input_path  <- "/Net/Groups/BGI/scratch/mterristi/PhD/CWD_0_1/New_Climate_Zones_percentiles/Drydowns/Drydowns_percentiles_pooled_v2"
output_path <- "/Net/Groups/BGI/scratch/mterristi/PhD/CWD_0_1/New_Climate_Zones_percentiles/Drydowns/Figures"
dir.create(output_path, showWarnings = FALSE, recursive = TRUE)

cat("Loading combined data...\n")
load(file.path(input_path, "combined_pooled_percentiles_summary_v2.Rdata"))
setDT(all_summaries)

# Same filter for all maps. EVI_E95 must be > 0 because we divide by it in
# dEVI.
all_summaries <- all_summaries[!is.na(EVI_E05) & !is.na(CWD_atE05) &
                                 !is.na(EVI_E95) & EVI_E95 > 0 &
                                 p05_position %in% c("before", "after") &
                                 climate_zone != "Arid Hot"]
cat("Pixels with valid thresholds:", nrow(all_summaries), "\n")
print(all_summaries[, .N, by = p05_position])

## ---- dEVI and dCWD ---------------------------------------------------------
# dCWD: always positive.
all_summaries[, Delta_CWD := abs(CWD_atE05 - CWD_atE95)]

# dEVI: change relative to the peak, in %, with a sign that depends on the
# regime.
#   after  (water-limited) : 100 * (EVI_E05 - EVI_E95) / EVI_E95, negative
#   before (energy-limited): 100 * (EVI_E95 - EVI_E05) / EVI_E95, positive
all_summaries[, Delta_EVI := 100 * fifelse(
  p05_position == "after",
  (EVI_E05 - EVI_E95) / EVI_E95,       # negative
  (EVI_E95 - EVI_E05) / EVI_E95        # positive
)]

cat("\nDelta_EVI (%) summary by regime:\n")
print(all_summaries[, .(median = round(median(Delta_EVI), 3),
                        q05    = round(quantile(Delta_EVI, 0.05), 3),
                        q95    = round(quantile(Delta_EVI, 0.95), 3),
                        min    = round(min(Delta_EVI),    3),
                        max    = round(max(Delta_EVI),    3)),
                    by = p05_position])

cat("\nDelta_CWD summary by regime (both should be positive):\n")
print(all_summaries[, .(median = round(median(Delta_CWD), 1),
                        min    = round(min(Delta_CWD),    1),
                        max    = round(max(Delta_CWD),    1)),
                    by = p05_position])

# ============================================================================
# Robinson projection
# ============================================================================
robin_crs <- "+proj=robin +datum=WGS84 +no_defs"

cat("\nLoading country borders...\n")
world <- ne_countries(scale = "medium", returnclass = "sf")
world_robin <- st_transform(world, robin_crs)
world_robin <- st_crop(world_robin,
                       st_transform(st_as_sfc(st_bbox(c(xmin = -180, xmax = 180,
                                                        ymin = -60, ymax = 90),
                                                      crs = 4326)),
                                    robin_crs))

project_pixels <- function(dt, value_col) {
  d <- dt[!is.na(get(value_col))]
  pts <- st_as_sf(d, coords = c("Longitude", "Latitude"), crs = 4326)
  pts <- st_transform(pts, robin_crs)
  coords <- st_coordinates(pts)
  data.table(x = coords[, 1], y = coords[, 2], val = d[[value_col]])
}

map_theme <- theme_void(base_size = 26) +
  theme(
    legend.position   = "right",
    legend.title      = element_text(size = 26),
    legend.text       = element_text(size = 24),
    legend.key.height = unit(1.3, "cm"),
    legend.key.width  = unit(0.55, "cm"),
    plot.title        = element_text(size = 24, hjust = 0,
                                     margin = margin(0, 0, 10, 0)),
    plot.margin       = margin(10, 10, 10, 10)
  )

make_map <- function(dt_proj, fill_label, palette = "viridis",
                     limits = NULL, direction = 1, title = "") {
  ggplot() +
    geom_tile(data = dt_proj, aes(x = x, y = y, fill = val),
              width = 11000, height = 11000) +
    geom_sf(data = world_robin, fill = NA, color = "grey40", linewidth = 0.2) +
    scale_fill_viridis(option = palette, name = fill_label,
                       limits = limits, oob = squish, direction = direction) +
    coord_sf(crs = robin_crs, expand = FALSE) +
    labs(title = title) +
    map_theme
}

# ---- the four states on one figure ----
cat("\nProjecting threshold points...\n")
pc_evi_E95   <- project_pixels(all_summaries, "EVI_E95")
pc_evi_E05   <- project_pixels(all_summaries, "EVI_E05")
pc_cwd_atE95 <- project_pixels(all_summaries, "CWD_atE95")
pc_cwd_atE05 <- project_pixels(all_summaries, "CWD_atE05")

cat("Building threshold maps...\n")
# Labels with subscripts, as in Fig. 2: EVI[95], EVI[05], CWD[EVI95],
# CWD[EVI05]
p1 <- make_map(pc_evi_E95,   expression("EVI"["95"]),                  "viridis", c(0, 1),    1, expression("a)  EVI"["95"]))
p2 <- make_map(pc_evi_E05,   expression("EVI"["05"]),                  "viridis", c(0, 1),    1, expression("b)  EVI"["05"]))
p3 <- make_map(pc_cwd_atE95, expression(atop("CWD"["EVI95"], "(mm)")), "magma",   c(0, 400), -1, expression("c)  CWD"["EVI95"]))
p4 <- make_map(pc_cwd_atE05, expression(atop("CWD"["EVI05"], "(mm)")), "magma",   c(0, 400), -1, expression("d)  CWD"["EVI05"]))

ggsave(file.path(output_path, "Map_EVI_E95.png"),   p1, width = 18, height = 9, dpi = 200, bg = "white")
ggsave(file.path(output_path, "Map_EVI_E05.png"),   p2, width = 18, height = 9, dpi = 200, bg = "white")
ggsave(file.path(output_path, "Map_CWD_atE95.png"), p3, width = 18, height = 9, dpi = 200, bg = "white")
ggsave(file.path(output_path, "Map_CWD_atE05.png"), p4, width = 18, height = 9, dpi = 200, bg = "white")

p_combined <- (p1 | p2) / (p3 | p4)
ggsave(file.path(output_path, "Maps_4panel_thresholds_v2.png"),
       p_combined, width = 32, height = 16, dpi = 200, bg = "white")
cat("4-panel thresholds map saved\n")

# ============================================================================
# MAPS OF dEVI AND dCWD
# ============================================================================
cat("\nBuilding Delta maps (relative ΔEVI)...\n")

cat("Projecting Delta points...\n")
pc_dEVI <- project_pixels(all_summaries, "Delta_EVI")
pc_dCWD <- project_pixels(all_summaries, "Delta_CWD")

# dEVI: colours centred on 0, brown when negative (water-limited), green when
# positive (energy-limited). The scale goes from -evi_lim to +evi_lim.
evi_lim <- 100
cat("ΔEVI (%) symmetric limit:", evi_lim, "\n")

make_map_dEVI <- function(dt_proj, title = "") {
  ggplot() +
    geom_tile(data = dt_proj, aes(x = x, y = y, fill = val),
              width = 11000, height = 11000) +
    geom_sf(data = world_robin, fill = NA, color = "grey40", linewidth = 0.2) +
    scale_fill_gradient2(
      name      = "ΔEVI\n(%)",
      low       = "#8B4513",   # brown: greenness declines with CWD (water-limited)
      mid       = "#F5F5DC",   # beige/cream at zero
      high      = "#2E7D32",   # green: greenness rises with CWD (energy-limited)
      midpoint  = 0,
      limits    = c(-evi_lim, evi_lim),
      breaks    = seq(-evi_lim, evi_lim, by = 20),
      labels    = function(x) sprintf("%d", as.integer(x)),
      oob       = squish
    ) +
    coord_sf(crs = robin_crs, expand = FALSE) +
    labs(title = title) +
    map_theme +
    theme(legend.key.height = unit(2.2, "cm"))   # taller ΔEVI colourbar
}

# dCWD: one colour gradient with a maximum value; higher values are shown
# with the colour of the maximum.
make_map_dCWD <- function(dt_proj, title = "") {
  ggplot() +
    geom_tile(data = dt_proj, aes(x = x, y = y, fill = val),
              width = 11000, height = 11000) +
    geom_sf(data = world_robin, fill = NA, color = "grey40", linewidth = 0.2) +
    scale_fill_viridis(option = "magma", name = "ΔCWD\n(mm)",
                       limits = c(0, 100), oob = squish, direction = -1) +
    coord_sf(crs = robin_crs, expand = FALSE) +
    labs(title = title) +
    map_theme
}

p_dEVI <- make_map_dEVI(pc_dEVI, "a)  ΔEVI (relative change with CWD)")
p_dCWD <- make_map_dCWD(pc_dCWD, "b)  ΔCWD (deficit budget magnitude)")

ggsave(file.path(output_path, "Map_Delta_EVI_relative.png"), p_dEVI, width = 30, height = 15, dpi = 300, bg = "white")
ggsave(file.path(output_path, "Map_Delta_CWD.png"),          p_dCWD, width = 30, height = 15, dpi = 300, bg = "white")

# The two maps on one figure
p_delta_combined <- p_dEVI / p_dCWD
ggsave(file.path(output_path, "Maps_Delta_2panel_relative.png"),
       p_delta_combined, width = 32, height = 18, dpi = 300, bg = "white")
cat("Delta maps saved (relative ΔEVI).\n")

# ============================================================================
# MAP OF THE REGIMES
# ============================================================================
cat("\nBuilding regime map...\n")
regime_df <- all_summaries[p05_position %in% c("before", "after")]
regime_df[, regime := factor(p05_position, levels = c("before", "after"))]

pts_regime <- st_as_sf(regime_df, coords = c("Longitude", "Latitude"), crs = 4326)
pts_regime <- st_transform(pts_regime, robin_crs)
coords_regime <- st_coordinates(pts_regime)
regime_proj <- data.table(x = coords_regime[, 1], y = coords_regime[, 2],
                          regime = regime_df$regime)

p_regime <- ggplot() +
  geom_tile(data = regime_proj, aes(x = x, y = y, fill = regime),
            width = 11000, height = 11000) +
  geom_sf(data = world_robin, fill = NA, color = "grey40", linewidth = 0.2) +
  scale_fill_manual(values = c("before" = "#2E7D32", "after" = "#8B4513"),
                    name   = "Regime") +
  coord_sf(crs = robin_crs, expand = FALSE) +
  labs(title = "Regime classification: before vs after") +
  map_theme +
  theme(legend.key.height = unit(0.8, "cm"))

ggsave(file.path(output_path, "Map_regime.png"),
       p_regime, width = 18, height = 9, dpi = 200, bg = "white")
cat("Regime map saved\n")

cat("\nAll done! Maps saved to:", output_path, "\n")

## ---- numbers cited in Sect. 3.2, printed only ------------------------------
frac <- all_summaries[, .N, by = p05_position]; frac[, pct := round(100*N/sum(N),1)]; print(frac)
cat("n total =", nrow(all_summaries), "\n\n")
cat(sprintf("|dEVI| %% : median %.1f  IQR %.1f  range %.1f-%.1f\n",
            median(abs(all_summaries$Delta_EVI)), IQR(abs(all_summaries$Delta_EVI)),
            min(abs(all_summaries$Delta_EVI)), max(abs(all_summaries$Delta_EVI))))
cat(sprintf(" dCWD mm : median %.1f  IQR %.1f  range %.1f-%.1f\n\n",
            median(all_summaries$Delta_CWD), IQR(all_summaries$Delta_CWD),
            min(all_summaries$Delta_CWD), max(all_summaries$Delta_CWD)))

## ---- variance of log S = log|dEVI| - log dCWD, printed only ----------------
# Part of the variance of log|S| due to dEVI and due to dCWD.
dec <- all_summaries[abs(Delta_EVI) > 0 & Delta_CWD > 0]
dec[, A := log(abs(Delta_EVI))]; dec[, B := log(Delta_CWD)]; dec[, logS := A - B]
vS <- var(dec$logS)
cat(sprintf("Var(log|dEVI|)=%.3f  Var(log dCWD)=%.3f  cor=%.2f\n",
            var(dec$A), var(dec$B), cor(dec$A, dec$B)))
cat(sprintf("Share of Var(log|S|):  dEVI=%.0f%%   dCWD=%.0f%%   (n=%d)\n",
            100*cov(dec$logS,dec$A)/vS, 100*cov(dec$logS,-dec$B)/vS, nrow(dec)))



## ---- maps for single years -------------------------------------------------
input_path  <- "/Net/Groups/BGI/scratch/mterristi/PhD/CWD_0_1/New_Climate_Zones_percentiles/Drydowns/Drydowns_percentiles_per_year_v2"
output_path <- "/Net/Groups/BGI/scratch/mterristi/PhD/CWD_0_1/New_Climate_Zones_percentiles/Drydowns/Figures"
per_year_file <- file.path(input_path, "combined_per_year_percentiles_summary_v2.Rdata")
dir.create(output_path, showWarnings = FALSE, recursive = TRUE)

SELECT_YEARS <- c(2010L, 2015L, 2018L)   # one row of maps per year

## ---- yearly states, same filters as above ----------------------------------
cat("Loading per-year data...\n")
load(per_year_file)                       # -> all_summaries
setDT(all_summaries)
stopifnot("Year" %in% names(all_summaries))

all_summaries <- all_summaries[!is.na(EVI_E05) & !is.na(CWD_atE05) &
                                 !is.na(EVI_E95) & EVI_E95 > 0 &
                                 p05_position %in% c("before", "after") &
                                 climate_zone != "Arid Hot"]
cat("Pixel-years with valid thresholds:", nrow(all_summaries), "\n")

all_summaries[, Delta_CWD := abs(CWD_atE05 - CWD_atE95)]
all_summaries[, Delta_EVI := 100 * fifelse(
  p05_position == "after",
  (EVI_E05 - EVI_E95) / EVI_E95,       # negative
  (EVI_E95 - EVI_E05) / EVI_E95        # positive
)]

## ---- Robinson projection, as above -----------------------------------------
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
  d <- dt[!is.na(get(value_col))]
  pts <- st_as_sf(d, coords = c("Longitude", "Latitude"), crs = 4326)
  pts <- st_transform(pts, robin_crs)
  coords <- st_coordinates(pts)
  data.table(x = coords[, 1], y = coords[, 2], val = d[[value_col]])
}

map_theme <- theme_void(base_size = 26) +
  theme(
    legend.position   = "right",
    legend.title      = element_text(size = 24),
    legend.text       = element_text(size = 20),
    legend.key.height = unit(1.1, "cm"),
    legend.key.width  = unit(0.5, "cm"),
    plot.title        = element_text(size = 26, hjust = 0,
                                     margin = margin(0, 0, 8, 0)),
    plot.margin       = margin(8, 8, 8, 8)
  )

evi_lim <- 100

make_map_dEVI <- function(dt_proj, title = "") {
  ggplot() +
    geom_tile(data = dt_proj, aes(x = x, y = y, fill = val),
              width = 11000, height = 11000) +
    geom_sf(data = world_robin, fill = NA, color = "grey40", linewidth = 0.2) +
    scale_fill_gradient2(
      name      = "ΔEVI\n(%)",
      low       = "#8B4513",   # brown: greenness declines with CWD (water-limited)
      mid       = "#F5F5DC",   # beige/cream at zero
      high      = "#2E7D32",   # green: greenness rises with CWD (energy-limited)
      midpoint  = 0,
      limits    = c(-evi_lim, evi_lim),
      breaks    = seq(-evi_lim, evi_lim, by = 20),
      labels    = function(x) sprintf("%d", as.integer(x)),
      oob       = squish
    ) +
    coord_sf(crs = robin_crs, expand = FALSE) +
    labs(title = title) +
    map_theme +
    theme(legend.key.height = unit(1.8, "cm"))
}

make_map_dCWD <- function(dt_proj, title = "") {
  ggplot() +
    geom_tile(data = dt_proj, aes(x = x, y = y, fill = val),
              width = 11000, height = 11000) +
    geom_sf(data = world_robin, fill = NA, color = "grey40", linewidth = 0.2) +
    scale_fill_viridis(option = "magma", name = "ΔCWD\n(mm)",
                       limits = c(0, 100), oob = squish, direction = -1) +
    coord_sf(crs = robin_crs, expand = FALSE) +
    labs(title = title) +
    map_theme
}

## ---- one row (dEVI, dCWD) per year -----------------------------------------
letters_seq <- letters   # a, b, c, ...
k <- 0
rows <- list()
for (yr in SELECT_YEARS) {
  dy <- all_summaries[Year == yr]
  if (nrow(dy) == 0) { cat("No pixels for year", yr, "- skipped\n"); next }
  cat(sprintf("Year %d: %d pixels\n", yr, nrow(dy)))
  
  pc_dEVI <- project_pixels(dy, "Delta_EVI")
  pc_dCWD <- project_pixels(dy, "Delta_CWD")
  
  l1 <- letters_seq[k + 1]; l2 <- letters_seq[k + 2]; k <- k + 2
  p_evi <- make_map_dEVI(pc_dEVI, sprintf("%s)  ΔEVI: %d", l1, yr))
  p_cwd <- make_map_dCWD(pc_dCWD, sprintf("%s)  ΔCWD: %d", l2, yr))
  
  rows[[length(rows) + 1]] <- (p_evi | p_cwd)
}

panel <- wrap_plots(rows, ncol = 1)

ggsave(file.path(output_path, "Maps_Delta_allyears_panel.png"),
       panel, width = 30, height = 8 * length(rows), dpi = 200,
       bg = "white", limitsize = FALSE)

cat("\nDone. Saved Maps_Delta_allyears_panel.png (", length(rows), "years )\n")