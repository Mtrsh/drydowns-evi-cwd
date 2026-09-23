# =============================================================================
# DRYDOWNS : LAND-COVER MASK AND STUDY PERIOD, ONE ZONE PER JOB
# Original file name: Prepare_drydowns_LC_masked.R
# Methods 2.5 (land cover) and 2.6 (study period)
# -----------------------------------------------------------------------------
# Last step before computing the peak and low states. Keeps only the years
# 2000 to 2023 and the grid cells with vegetation.
#
# Years: we keep 2000 to 2023 and remove the weeks without EVI or CWD. This
# removes the 1999 weeks from the drydown step.
#
# Land cover: applied for each grid cell and year, with the ESA-CCI land-cover
# map at 0.1 degree (one layer per year, 1992 to 2022). A grid cell is kept in
# a given year if it has a land-cover class in this map and is not irrigated
# cropland (class 20). Urban, bare soil, water, and snow and ice are already
# NA in this map (Methods 2.5), so they are removed too. The map stops in 2022,
# so the 2022 layer is also used for 2023. The map was prepared from the
# ESA-CCI/C3S annual maps before this pipeline.
#
# Coordinates are rounded to one decimal on both sides before matching on
# lat, lon and year.
#
# Input  : <zone>_cwd_drydowns_with_evi.Rdata from 04b_drydowns_add_evi.R
#          dominant_landcover_no_others_1992_2022_0.1deg.tif
# Output : <zone>_cwd_drydowns_with_evi_lc.Rdata (object extracted_drydowns:
#            lat, lon, time, cwd, year, EVI)
#          The column evi is renamed EVI here, and all next scripts use EVI.
#
# Run    : one job per zone, the zone is given as argument:
#            Rscript 04c_drydowns_landcover_mask.R tropical
# =============================================================================

## ---- paths (change these lines to run the script elsewhere) ----------------
source_path <- "/Net/Groups/BGI/scratch/mterristi/PhD/CWD_0_1/New_Climate_Zones/Drydowns/"
output_path <- "/Net/Groups/BGI/scratch/mterristi/PhD/CWD_0_1/New_Climate_Zones_percentiles/Drydowns/"
path_lc     <- "/Net/Groups/BGI/scratch/mterristi/PhD/Data/ESACCI_LC/dominant_landcover_no_others_1992_2022_0.1deg.tif"

## ---- environment -----------------------------------------------------------

.libPaths(c("/Net/Groups/BGI/scratch/mterristi/R_lodomeria/4.4",
            "/Net/Groups/BSI/work_scratch/quincy/model/software/r_packages/r_4.4.x",
            "/opt/ohpc/pub/libs/gnu12/R/4.4.0/lib64/R/library"))

# Libraries needed on our cluster, skipped if they do not exist.
for (lib in c("/opt/ohpc/pub/libs/hwloc/lib/libhwloc.so.15",
              "/opt/ohpc/pub/libs/gnu9/openmpi4/hdf5/1.10.8/lib/libhdf5_hl.so.100",
              "/opt/ohpc/pub/apps/gdal/3.5.1/lib/libgdal.so.31")) {
  if (file.exists(lib)) dyn.load(lib)
}

library(data.table)
library(terra)

## ---- zone to process -------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Provide zone as argument")
zone <- tolower(args[1])

cat("\n===================================\n")
cat("STEP 0 — PREPARING:", toupper(zone), "\n")
cat("===================================\n\n")

dir.create(output_path, showWarnings = FALSE, recursive = TRUE)

input_file  <- file.path(source_path, paste0(zone, "_cwd_drydowns_with_evi.Rdata"))
output_file <- file.path(output_path, paste0(zone, "_cwd_drydowns_with_evi_lc.Rdata"))

if (!file.exists(input_file)) stop("Input not found: ", input_file)

## ---- drydown data and study period -----------------------------------------
cat("Loading:", input_file, "\n")
load(input_file)
setDT(extracted_drydowns)
cat("Source rows:", nrow(extracted_drydowns), "\n")
cat("Columns:", paste(names(extracted_drydowns), collapse = ", "), "\n")

if ("evi" %in% names(extracted_drydowns)) setnames(extracted_drydowns, "evi", "EVI")

# Study period, and weeks with both CWD and EVI (see header).
extracted_drydowns <- extracted_drydowns[year >= 2000 & year <= 2023]
n_before_drop <- nrow(extracted_drydowns)
extracted_drydowns <- extracted_drydowns[!is.na(EVI) & !is.na(cwd)]
cat("After year/NA filter:", nrow(extracted_drydowns),
    "(removed", n_before_drop - nrow(extracted_drydowns), ")\n")

## ---- land-cover map --------------------------------------------------------
cat("\nLoading LC raster...\n")
lc <- rast(path_lc)
years_available <- 1992:(1992 + nlyr(lc) - 1)
cat("LC years available:", min(years_available), "to", max(years_available), "\n")

IRRIGATED_CROPLAND_VALUE <- 20

## ---- grid cells and years to keep ------------------------------------------
# List of (lat, lon, year) to keep. For years after 2022, the 2022 layer is
# used.
cat("\nBuilding per-year LC lookup...\n")
lc_list <- list()
for (yr in 2000:2023) {
  if (yr <= max(years_available)) {
    layer_idx <- which(years_available == yr)
  } else {
    layer_idx <- which(years_available == max(years_available))
    cat("  Year", yr, "-> using", max(years_available), "LC layer\n")
  }
  lc_yr <- lc[[layer_idx]]
  lc_df <- as.data.table(as.data.frame(lc_yr, xy = TRUE))
  setnames(lc_df, c("x", "y", names(lc_yr)), c("lon", "lat", "lc_value"))
  lc_df <- lc_df[!is.na(lc_value) & lc_value != IRRIGATED_CROPLAND_VALUE]
  lc_df[, lat := round(lat, 1)]
  lc_df[, lon := round(lon, 1)]
  lc_df[, year := yr]
  lc_list[[as.character(yr)]] <- lc_df[, .(lat, lon, year)]
}
lc_lookup <- rbindlist(lc_list)
setkey(lc_lookup, lat, lon, year)
cat("LC lookup rows:", nrow(lc_lookup), "\n")
rm(lc_list); gc(verbose = FALSE)

## ---- keep the grid cells and years of the list -----------------------------
# Keep a drydown week only if its grid cell and year are in the list.
cat("\nApplying per-year LC mask...\n")
extracted_drydowns[, lat_r := round(lat, 1)]
extracted_drydowns[, lon_r := round(lon, 1)]

n_before <- nrow(extracted_drydowns)
extracted_drydowns <- extracted_drydowns[
  lc_lookup,
  on = c(lat_r = "lat", lon_r = "lon", year = "year"),
  nomatch = 0
]
extracted_drydowns[, c("lat_r", "lon_r") := NULL]

cat("After LC mask:", nrow(extracted_drydowns),
    "rows (removed", n_before - nrow(extracted_drydowns), ")\n")
cat("Retained:", round(100 * nrow(extracted_drydowns) / n_before, 1), "%\n")

## ---- summary and save ------------------------------------------------------
cat("\n=== SUMMARY ===\n")
cat("Unique pixels:", uniqueN(extracted_drydowns[, .(lat, lon)]), "\n")
cat("Years covered:", paste(sort(unique(extracted_drydowns$year)), collapse = ", "), "\n")

cat("\nSaving:", output_file, "\n")
save(extracted_drydowns, file = output_file)
cat("Done.\n")
