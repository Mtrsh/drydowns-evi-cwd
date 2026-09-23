# =============================================================================
# CLIMATE ZONES : WEEKLY CWD FOR EACH ZONE
# Original file name: Merge_CWD_CZnew.R
# Methods 2.4
# -----------------------------------------------------------------------------
# Splits the weekly CWD maps by climate zone and saves each zone as a table,
# with one row per land grid cell and week:
#   lon, lat, time, cwd
# All the next steps (drydowns, EVI, states) use these tables.
#
# The zones come from grouped_koppen_1991_2020.nc, the Koeppen-Geiger map for
# 1991 to 2020 (Beck et al., 2023) grouped into eight classes:
#   1 Tropical, 2 Arid Hot, 3 Arid Cold, 4 Temperate Dry, 5 Temperate Humid,
#   6 Cold Dry, 7 Cold Humid, 8 Polar
# Polar and Arid Hot are kept here and removed later in the analysis.
#
# Grid: the zone map is cropped to the extent of a CWD file, because the CWD
# files have a southern edge at -60.00001 and not exactly -60. A difference of
# one row would shift all zones by one grid cell, so the script checks that
# both grids are the same (compareGeom) and stops if they are not.
#
# Table format: time is stored as a factor with the week dates as text, the
# same dates as in 01b_cwd_weekly.R. Coordinates are rounded to two decimals so
# that tables can be joined on lon and lat in the next scripts.
#
# Input  : cwd_weekly_<year>.nc, 1980 to 2023, from 01b_cwd_weekly.R
#          grouped_koppen_1991_2020.nc
# Output : cwd_<zone>.Rdata (object cwd_df), one file per zone, for example
#          cwd_arid_hot.Rdata
#
# Run    : one job, all zones and years one after the other.
# =============================================================================

## ---- paths (change these lines to run the script elsewhere) ----------------
cwd_path <- "/Net/Groups/BGI/scratch/mterristi/PhD/CWD_0_1/cwd_computation/"
out_path <- "/Net/Groups/BGI/scratch/mterristi/PhD/CWD_0_1/New_Climate_Zones/CWD_cz"
path_kg  <- "/Net/Groups/BGI/scratch/mterristi/PhD/Data/Koppen/grouped_koppen_1991_2020.nc"

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

library(terra)
library(data.table)
library(lubridate)

dir.create(out_path, recursive = TRUE, showWarnings = FALSE)

## ---- zone map on the CWD grid ----------------------------------------------
print("Loading climate zones...")
climate_map <- rast(path_kg)

# Extent taken from a CWD file (see header).
sample_cwd <- rast(file.path(cwd_path, "cwd_weekly_2023.nc"))[[1]]
cwd_extent <- ext(sample_cwd)

print("Cropping climate zones to exact CWD extent...")
climate_map <- crop(climate_map, cwd_extent)

# If a small difference remains after cropping, use the CWD extent.
if (!all(dim(climate_map) == dim(sample_cwd))) {
  set.ext(climate_map, cwd_extent)
}

# Stop if the grids still differ: otherwise grid cells would be put in the
# wrong zone without any error.
if (compareGeom(climate_map, sample_cwd, stopOnError = FALSE)) {
  cat("SUCCESS: Climate Map and CWD files are perfectly aligned.\n")
} else {
  stop("Grids still do not match. Check if the resolutions are identical.")
}

zones <- c("Tropical" = 1, "Arid Hot" = 2, "Arid Cold" = 3, "Temperate Dry" = 4,
           "Temperate Humid" = 5, "Cold Dry" = 6, "Cold Humid" = 7, "Polar" = 8)

years <- 1980:2023

## ---- one table per zone ----------------------------------------------------
for (z_name in names(zones)) {

  cat("\n========================================\n")
  cat("PROCESSING CLIMATE ZONE:", z_name, "\n")
  cat("========================================\n")

  z_val <- zones[z_name]

  zone_mask <- ifel(climate_map == z_val, 1, NA)

  dt_list <- list()

  for (yr in years) {
    infile <- file.path(cwd_path, paste0("cwd_weekly_", yr, ".nc"))
    if (!file.exists(infile)) { cat("  Year:", yr, "Skipped\n"); next }

    cat("  Year:", yr, "... ")
    cwd_yr <- rast(infile)

    cwd_masked <- mask(cwd_yr, zone_mask)

    # One row per grid cell, one column per week; ocean and other zones removed.
    dt_yr <- as.data.table(cwd_masked, xy = TRUE, na.rm = TRUE)

    if (nrow(dt_yr) > 0) {
      # One row per grid cell and week.
      dt_long <- melt(dt_yr, id.vars = c("x", "y"), variable.name = "layer", value.name = "cwd")

      n_weeks <- nlyr(cwd_yr)
      dates <- seq(as.Date(paste0(yr, "-01-01")), by = "week", length.out = n_weeks)

      # Layer names become week dates (time stays a factor, see header).
      levels(dt_long$layer) <- as.character(dates)
      setnames(dt_long, c("x", "y", "layer"), c("lon", "lat", "time"))

      # Round to two decimals so that lon and lat can be used to join tables.
      dt_long[, `:=`(lon = round(lon, 2), lat = round(lat, 2))]

      dt_list[[as.character(yr)]] <- dt_long
      cat("Done\n")
    } else {
      cat("No data\n")
    }

    rm(cwd_yr, cwd_masked, dt_yr, dt_long); gc()
  }

  ## ---- write the zone ------------------------------------------------------
  if (length(dt_list) > 0) {
    cwd_df <- as.data.frame(rbindlist(dt_list))
    safe_name <- tolower(gsub(" ", "_", z_name))
    out_file <- file.path(out_path, paste0("cwd_", safe_name, ".Rdata"))

    cat("Saving final file:", out_file, "\n")
    save(cwd_df, file = out_file)
  }

  rm(dt_list, zone_mask); gc()
}

cat("\nALL PROCESSING COMPLETE!\n")
