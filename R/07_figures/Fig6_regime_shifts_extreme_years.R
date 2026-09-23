# =============================================================================
# FIGURE 6 : ENERGY-LIMITED GRID CELLS THAT BECAME WATER-LIMITED IN EXTREME YEARS
# Original file name: Figures_Switched_of_regimes.R
# Methods 2.6 (regimes from the pooled and yearly states)
# -----------------------------------------------------------------------------
# Which grid cells that are energy-limited over the whole period became
# water-limited in 2010, 2015 and 2018?
#
# Usual regime: from the pooled states (all years together). We keep only the
# grid cells that are energy-limited (EL).
# Regime of a given year: from the yearly states of that year, computed only
# from the drydown of that year. An EL grid cell is marked "EL to WL" if, that
# year, its low state is at a higher CWD than its peak state. Otherwise it is
# marked "EL to EL".
#
# A grid cell is used only if its regime is clear both over the whole period
# and in the given year: both states exist, EVI95 > 0, the regime is known, and
# the two states are more than 5 mm of CWD apart (CWD_FLOOR). This avoids
# counting a change when both states are at almost the same CWD. Arid Hot is
# removed.
#
# Part 1 (Fig. 6): one map per year of the EL grid cells, changed or not; a
# bar chart of the share of changed grid cells by zone and year; the same
# values in a CSV file.
# Part 2 (supplement): the share of changed grid cells for every year from
# 2000 to 2023, weighted by cos(latitude), per zone and for all zones, as a
# time series with the three extreme years marked. Part 2 reads the data again
# so it can be run alone.
#
# The label "EL → EL" is used in three places: in the fifelse() text, in the
# factor levels and in the colour names. It must be exactly the same in all
# three, otherwise these grid cells disappear from the map.
#
# Input  : combined_pooled_percentiles_summary_v2.Rdata
#          combined_per_year_percentiles_summary_v2.Rdata
# Output : FigureS_EL_to_WL_flips_map.png                    Fig. 6
#          FigureS_EL_to_WL_flips_by_zone.png                 supplement
#          TableS_EL_to_WL_flips_by_zone.csv
#          FigureS_Regime_flip_frequency_timeseries.png       supplement
#          TableS_EL_to_WL_flip_frequency_timeseries.csv
#          The file names come from our working folder, see the README.
#
# Run    : Rscript Fig6_regime_shifts_extreme_years.R
# =============================================================================

## ---- environment -----------------------------------------------------------

# Our R package folder on the cluster, used only if it exists.
custom_libs <- c("/Net/Groups/BGI/scratch/mterristi/R_lodomeria/4.4",
                 "/Net/Groups/BSI/work_scratch/quincy/model/software/r_packages/r_4.4.x")
.libPaths(c(custom_libs[dir.exists(custom_libs)], .libPaths()))
# Libraries needed on our cluster, skipped if they do not exist.
for (lib in c("/opt/ohpc/pub/libs/hwloc/lib/libhwloc.so.15",
              "/opt/ohpc/pub/libs/gnu9/openmpi4/hdf5/1.10.8/lib/libhdf5_hl.so.100",
              "/opt/ohpc/pub/apps/gdal/3.5.1/lib/libgdal.so.31")) {
  if (file.exists(lib)) dyn.load(lib)
}
library(data.table); library(ggplot2); library(sf); library(rnaturalearth); library(patchwork)
## ---- paths and settings (change these lines to run the script elsewhere) ---
per_year_file <- "/Net/Groups/BGI/scratch/mterristi/PhD/CWD_0_1/New_Climate_Zones_percentiles/Drydowns/Drydowns_percentiles_per_year_v2/combined_per_year_percentiles_summary_v2.Rdata"
pooled_file   <- "/Net/Groups/BGI/scratch/mterristi/PhD/CWD_0_1/New_Climate_Zones_percentiles/Drydowns/Drydowns_percentiles_pooled_v2/combined_pooled_percentiles_summary_v2.Rdata"
output_path   <- "/Net/Groups/BGI/scratch/mterristi/PhD/CWD_0_1/New_Climate_Zones_percentiles/Drydowns/Figures"
dir.create(output_path, showWarnings = FALSE, recursive = TRUE)
MAP_YEARS   <- c(2010L, 2015L, 2018L)                 # the extreme years shown
CWD_FLOOR   <- 5     # mm; regime well-defined only if |CWD_atE05 - CWD_atE95| > this
zone_levels <- c("Tropical","Arid Cold","Temperate Dry","Temperate Humid","Cold Dry","Cold Humid")
year_cols   <- c("2010"="#3182bd","2015"="#fd8d3c","2018"="#e31a1c",
                 "2020"="#6a3d9a","2023"="#1b9e77")   # blue / orange / red / purple / green
## ---- regime of each grid cell, with the filters of the header --------------
get_regime <- function(dt) {
  dt <- as.data.table(dt)
  dt <- dt[!is.na(EVI_E05) & !is.na(CWD_atE05) &
             !is.na(EVI_E95) & EVI_E95 > 0 &
             p05_position %in% c("before","after") &
             climate_zone != "Arid Hot"]
  dt[, Delta_CWD := abs(CWD_atE05 - CWD_atE95)]
  dt <- dt[Delta_CWD > CWD_FLOOR]
  dt[, regime := fifelse(p05_position=="before","EL","WL")]
  dt[, climate_zone := factor(climate_zone, levels=zone_levels)]
  dt[, cell := paste(round(Latitude,3), round(Longitude,3), sep="_")]
  dt
}
## ---- grid cells that are EL over the whole period --------------------------
load(pooled_file)                                    # -> all_summaries (pooled)
pool <- get_regime(all_summaries); rm(all_summaries)
pool <- unique(pool[, .(cell, climate_zone, base_regime = regime)], by="cell")
pool_EL <- pool[base_regime == "EL"]
cat("Climatologically-EL cells:", nrow(pool_EL), "\n")
## ---- yearly regimes --------------------------------------------------------
load(per_year_file)                                  # -> all_summaries (per year)
py <- get_regime(all_summaries); rm(all_summaries)
py <- py[, .(cell, Latitude, Longitude, Year, year_regime = regime)]
## ---- regime of these grid cells in each extreme year -----------------------
el <- merge(py[Year %in% MAP_YEARS], pool_EL, by="cell")   # inner: EL & defined that year
# The two cases; the labels must match the colour names below (see header).
el[, status := factor(fifelse(year_regime=="WL", "EL → WL", "EL → EL"),
                      levels=c("EL → EL","EL → WL"))]
el[, Year := factor(Year, levels=MAP_YEARS)]
el[, w := cos(Latitude * pi/180)]                          # area weight = cos(latitude)
## ---- share of changed grid cells by zone and year --------------------------
flip_zone <- el[, .(base_EL      = .N,
                    EL_to_WL_pct = round(100 * sum(w * (year_regime=="WL")) / sum(w), 1)),
                by=.(Year, climate_zone)][order(Year, climate_zone)]
cat("\n======  EL -> WL flip rate (% of climatologically-EL pixels)  ======\n")
print(flip_zone)
fwrite(flip_zone, file.path(output_path, "TableS_EL_to_WL_flips_by_zone.csv"))
cat("\n[TEXT] Global EL->WL flip rate by year:\n")
print(el[, .(EL_to_WL_pct = round(100*sum(w*(year_regime=="WL"))/sum(w),1), base_EL=.N), by=Year][order(Year)])
## ---- Figure 6: map per year ------------------------------------------------
robin <- "+proj=robin +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
# Colour names must match the labels above (see header).
pal   <- c("EL → EL"="#2166ac", "EL → WL"="#b2182b")
mp <- el[order(status)]                              # stayed first, flips drawn on top
pts <- st_transform(st_as_sf(mp, coords=c("Longitude","Latitude"), crs=4326), crs=robin)
xy  <- st_coordinates(pts)
dmap <- data.table(X=xy[,1], Y=xy[,2], status=mp$status, Year=mp$Year)
world <- st_crop(st_transform(ne_countries(scale="medium", returnclass="sf"), robin),
                 xmin=-17e6, xmax=17e6, ymin=-6e6, ymax=8.5e6)
# ---- graticule: lines and labels from the same lon/lat vectors --------------
grat_lat <- c(-60, -30, 0, 30, 60)                     # parallels (latitude)
grat_lon <- c(-180, -120, -60, 0, 60, 120, 180)        # meridians (longitude)
LAT_LO   <- -60; LAT_HI <- 84                          # plotted vertical frame
# vertical limits of the map in Robinson metres (so that -60 and 60 are inside)
yb <- st_coordinates(st_transform(
  st_as_sf(data.frame(lon = 0, lat = c(LAT_LO, LAT_HI)), coords = c("lon","lat"), crs = 4326),
  robin))[, "Y"]
# lines at grat_lat and grat_lon, cut to the map limits
grat <- st_transform(st_graticule(lon = grat_lon, lat = grat_lat, crs = st_crs(4326)), robin)
grat <- st_crop(grat, xmin = -17e6, xmax = 17e6, ymin = yb[1], ymax = yb[2])
# latitude labels on the left edge, all panels
lat_xy  <- st_coordinates(st_transform(
  st_as_sf(data.frame(lon = -180, lat = grat_lat), coords = c("lon","lat"), crs = 4326), robin))
lat_lab <- data.table(x = lat_xy[, "X"], y = lat_xy[, "Y"],
                      label = paste0(abs(grat_lat), "°", ifelse(grat_lat > 0, "N", ifelse(grat_lat < 0, "S", ""))))
# longitude labels along the bottom edge, last panel only
lon_xy  <- st_coordinates(st_transform(
  st_as_sf(data.frame(lon = grat_lon, lat = LAT_LO), coords = c("lon","lat"), crs = 4326), robin))
lon_lab <- data.table(x = lon_xy[, "X"], y = lon_xy[, "Y"],
                      label = paste0(abs(grat_lon), "°", ifelse(grat_lon > 0, "E", ifelse(grat_lon < 0, "W", ""))),
                      Year  = factor(MAP_YEARS[length(MAP_YEARS)], levels = MAP_YEARS))
p_map <- ggplot() +
  geom_tile(data=dmap, aes(X, Y, fill=status), width=11000, height=11000) +
  geom_sf(data=grat,  color="grey85", linewidth=0.2) +            # our own grid lines
  geom_sf(data=world, fill=NA, color="gray30", linewidth=.2) +
  geom_text(data=lat_lab, aes(x=x, y=y, label=label), inherit.aes=FALSE,
            hjust=1, nudge_x=-3e5, size=6, colour="grey45") +     # latitude labels (left)
  geom_text(data=lon_lab, aes(x=x, y=y, label=label), inherit.aes=FALSE,
            vjust=1, nudge_y=-3e5, size=6, colour="grey45") +     # longitude labels (bottom)
  scale_fill_manual(values=pal, name=NULL, drop=FALSE) +
  coord_sf(crs=robin, ylim=yb, expand=FALSE, clip="off") +
  facet_wrap(~Year, ncol=1) +                          # 2010 / 2015 / 2018 stacked
  guides(fill=guide_legend(nrow=1, override.aes=list(size=6))) +
  # white background, no frame; the lines drawn above replace the default grid
  theme_minimal(base_size=26) +
  theme(legend.position  = "bottom",
        legend.text      = element_text(size=26),
        strip.text       = element_text(size=28, face="bold", margin=margin(b=6)),
        panel.grid       = element_blank(),            # no auto graticule (we draw our own)
        axis.title       = element_blank(),
        axis.text        = element_blank(),            # no auto axis labels (we place our own)
        axis.ticks       = element_blank(),
        plot.margin      = margin(6, 10, 22, 30))
ggsave(file.path(output_path, "FigureS_EL_to_WL_flips_map.png"),
       p_map, width=13, height=6.5*length(MAP_YEARS), dpi=300, bg="white", limitsize=FALSE)
## ---- share of changed grid cells by zone and year, bar chart, supplement ---
p_bar <- ggplot(flip_zone, aes(climate_zone, EL_to_WL_pct, fill=Year)) +
  geom_col(position=position_dodge(width=0.8), width=0.7) +
  scale_fill_manual(values=year_cols, name="Year") +
  labs(x=NULL, y="% of climatologically-EL area that switched to WL (ΔEVI < 0)") +
  theme_bw(base_size=20) +
  theme(panel.grid.minor=element_blank(),
        axis.text.x=element_text(angle=25, hjust=1))
ggsave(file.path(output_path, "FigureS_EL_to_WL_flips_by_zone.png"),
       p_bar, width=13, height=7, dpi=300, bg="white", limitsize=FALSE)
cat("\nDone. Wrote:\n",
    " - FigureS_EL_to_WL_flips_map.png\n",
    " - FigureS_EL_to_WL_flips_by_zone.png\n",
    " - TableS_EL_to_WL_flips_by_zone.csv\n", sep="")

EXTREME_YEARS <- c(2010L, 2015L, 2018L)   # dashed markers
CWD_FLOOR     <- 5
zone_levels   <- c("Tropical","Arid Cold","Temperate Dry","Temperate Humid","Cold Dry","Cold Humid")
zone_cols     <- c("Tropical"="#1b9e77","Arid Cold"="#d95f02","Temperate Dry"="#7570b3",
                   "Temperate Humid"="#e7298a","Cold Dry"="#66a61e","Cold Humid"="#e6ab02",
                   "Global (all zones)"="black")

## ---- regime of each grid cell, same filters as part 1 ----------------------
get_regime <- function(dt) {
  dt <- as.data.table(dt)
  dt <- dt[!is.na(EVI_E05) & !is.na(CWD_atE05) &
             !is.na(EVI_E95) & EVI_E95 > 0 &
             p05_position %in% c("before","after") &
             climate_zone != "Arid Hot"]
  dt[, Delta_CWD := abs(CWD_atE05 - CWD_atE95)]
  dt <- dt[Delta_CWD > CWD_FLOOR]
  dt[, regime := fifelse(p05_position=="before","EL","WL")]
  dt[, climate_zone := factor(climate_zone, levels=zone_levels)]
  dt[, cell := paste(round(Latitude,3), round(Longitude,3), sep="_")]
  dt
}

## ---- grid cells that are EL over the whole period --------------------------
load(pooled_file)                                    # -> all_summaries (pooled)
pool <- get_regime(all_summaries); rm(all_summaries)
pool_EL <- unique(pool[regime=="EL", .(cell, climate_zone)], by="cell")
cat("Climatologically-EL cells:", nrow(pool_EL), "\n")

## ---- yearly regimes, all years ---------------------------------------------
load(per_year_file)                                  # -> all_summaries (per year)
py <- get_regime(all_summaries); rm(all_summaries)
py <- py[, .(cell, Latitude, Year, year_regime = regime)]

## ---- EL grid cells, all years ----------------------------------------------
el <- merge(py, pool_EL, by="cell")                  # inner: EL & defined that year
el[, w := cos(Latitude * pi/180)]                    # area weight = cos(latitude)

miss <- setdiff(EXTREME_YEARS, unique(el$Year))
if (length(miss)) cat("WARNING: extreme year(s) absent from the per-year file:",
                      paste(miss, collapse=", "), "\n")

## ---- share of changed grid cells per year, weighted by area ----------------
ts_zone <- el[!is.na(climate_zone),
              .(EL_to_WL_pct = 100*sum(w*(year_regime=="WL"))/sum(w), nEL = .N),
              by=.(Year, climate_zone)]
ts_zone[, series := as.character(climate_zone)]

ts_glob <- el[, .(EL_to_WL_pct = 100*sum(w*(year_regime=="WL"))/sum(w), nEL = .N), by=Year]
ts_glob[, series := "Global (all zones)"]

ts <- rbind(ts_zone[, .(Year, series, EL_to_WL_pct, nEL)],
            ts_glob[, .(Year, series, EL_to_WL_pct, nEL)])
ts[, series := factor(series, levels=c(zone_levels, "Global (all zones)"))]
setorder(ts, series, Year)

fwrite(ts, file.path(output_path, "TableS_EL_to_WL_flip_frequency_timeseries.csv"))
cat("\nGlobal EL->WL flip % in the extreme years:\n")
print(ts_glob[Year %in% EXTREME_YEARS][order(Year)])

## ---- time series figure, supplement ----------------------------------------
xbreaks <- sort(unique(c(seq(min(ts$Year), max(ts$Year), by = 4), EXTREME_YEARS)))

p_ts <- ggplot(ts, aes(Year, EL_to_WL_pct, colour=series)) +
  geom_vline(xintercept=EXTREME_YEARS, linetype="dashed", colour="grey55", linewidth=.4) +
  geom_line(aes(linewidth = series == "Global (all zones)")) +
  geom_point(size=1.5) +
  scale_linewidth_manual(values=c("FALSE"=0.8, "TRUE"=1.8), guide="none") +
  scale_colour_manual(values=zone_cols, name=NULL) +
  scale_x_continuous(breaks=xbreaks) +
  labs(x=NULL, y="% of climatologically-EL area that switched to WL (ΔEVI < 0)") +
  theme_bw(base_size=20) +
  theme(panel.grid.minor=element_blank(),
        legend.position="bottom",
        axis.text.x=element_text(angle=45, hjust=1))

ggsave(file.path(output_path, "FigureS_Regime_flip_frequency_timeseries.png"),
       p_ts, width=14, height=8, dpi=300, bg="white", limitsize=FALSE)

cat("\nDone. Wrote:\n",
    " - FigureS_Regime_flip_frequency_timeseries.png\n",
    " - TableS_EL_to_WL_flip_frequency_timeseries.csv\n", sep="")