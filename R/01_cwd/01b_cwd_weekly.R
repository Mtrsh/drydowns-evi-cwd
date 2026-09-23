# =============================================================================
# CWD, STEP 2 of 2 : WEEKLY CUMULATIVE WATER DEFICIT
# Original file name: cwd01_array.R
# Methods 2.5
# -----------------------------------------------------------------------------
# Computes the weekly cumulative water deficit (CWD) for one year, on the
# 0.1 degree land grid.
#
# Definition: the weekly water balance is the sum of daily P - ET over the
# week. CWD is the running sum of this balance while it is negative, given as a
# positive value in mm. CWD goes back to zero in two cases:
#   1. the running balance becomes positive (rain has made up for the deficit);
#   2. the week is in the wettest month of the grid cell (from step 1,
#      01a_cwd_wettest_month.R). This stops CWD from growing for years in
#      regions where case 1 rarely happens.
#
# Years are linked: the deficit at the last week of a year is saved in
# final_deficit_<year>.rds and used as the starting value of the next year.
# This way, a drydown going on at the end of December is not cut.
#
# IMPORTANT: for this reason, years must be run in order, starting in 1980.
# If the file of the previous year is missing, the script does not stop: it
# starts the year from zero. Running all years at the same time would
# therefore give wrong results without any error message.
#
# The data are read in blocks of 30 days to save memory, but the weekly sums
# are computed on the whole year at once.
#
# Input  : MSWEP v2.8 daily precipitation for <year>
#          GLEAM v4.2a daily evapotranspiration for <year>
#          GLDAS land-sea mask at 0.25 degree
#          wettest_month.nc from step 1
#          final_deficit_<year-1>.rds from the previous year
# Output : cwd_weekly_<year>.nc, weekly CWD in mm, used by
#          03_cwd_by_climate_zone.R
#          final_deficit_<year>.rds, starting value for the next year
#
# Run    : one job per year, the year is given as argument:
#            Rscript 01b_cwd_weekly.R 1980
#          Progress is written to LOG1-cwd_<year>.out in the output folder.
# =============================================================================

## ---- paths (change these lines to run the script elsewhere) ----------------
path_prec <- "/Net/Groups/BGI/data/DataStructureMDI/DATA/grid/Global/0d10_daily/MSWEP/v280/Data/"
path_et   <- "/Net/Groups/BGI/data/DataStructureMDI/DATA/Incoming/gleam/v4.2a/daily/"
path_lsm  <- "/Net/Groups/BGI/scratch/mterristi/PhD/Data/LSM/GLDASp5_landmask_025d.nc4"
out_path  <- "/Net/Groups/BGI/scratch/mterristi/PhD/CWD_0_1/cwd_computation/"

## ---- environment -----------------------------------------------------------
# Libraries needed on our cluster. They are skipped if they do not exist, so
# the script also runs on other machines.
for (lib in c("/opt/ohpc/pub/libs/hwloc/lib/libhwloc.so.15",
              "/opt/ohpc/pub/libs/gnu9/openmpi4/hdf5/1.10.8/lib/libhdf5_hl.so.100",
              "/opt/ohpc/pub/apps/gdal/3.5.1/lib/libgdal.so.31")) {
  if (file.exists(lib)) dyn.load(lib)
}

library(ncdf4)
library(raster)
library(terra)
library(lubridate)

## ---- year to process -------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
year <- as.numeric(args[1])

setwd(out_path)

options(immediate.write = TRUE)
out <- file(paste0("LOG1-cwd_", year, ".out"), open = "wt")
sink(out, type = "output")
sink(out, type = "message")

print(paste("Processing year:", year))
flush.console()

## ---- land-sea mask ---------------------------------------------------------
# Same grid as in step 1.
print("Loading land-sea mask...")
flush.console()
lsm_025 <- rast(path_lsm)
template_01 <- rast(xmin = -180, xmax = 180, ymin = -60, ymax = 90,
                    resolution = 0.1, crs = "EPSG:4326")
lsm_01 <- resample(lsm_025, template_01, method = "near")
mask_raster <- ifel(lsm_01 == 0, NA, 1)
lsm_global <- crop(mask_raster, ext(-180, 180, -60, 90))

## ---- month of the yearly reset, from step 1 --------------------------------
print("Loading climatological wettest month...")
flush.console()
wettest_month <- rast(file.path(out_path, "wettest_month.nc"))

## ---- deficit from the previous year ----------------------------------------
# Does not exist for 1980, the first year. For later years, if the file is
# missing, the year starts from zero (see header).
previous_deficit <- NULL
if (year > 1980) {
  prev_file <- file.path(out_path, paste0("final_deficit_", year - 1, ".rds"))
  if (file.exists(prev_file)) {
    print("Loading previous year's deficit...")
    flush.console()
    previous_deficit <- readRDS(prev_file)
  }
}

## ---- input files -----------------------------------------------------------
prec_file <- file.path(path_prec, paste0("prec.MSWEPv280.3600.1800.", year, ".nc"))
et_file   <- file.path(path_et, year, paste0("E_", year, "_GLEAM_v4.2a.nc"))

# Here a missing file stops the script: skipping a year would break the link
# between years.
if (!file.exists(prec_file) || !file.exists(et_file)) {
  stop(paste("Files missing for year", year))
}

nc <- nc_open(prec_file)
n_days <- nc$dim$time$len
nc_close(nc)

## ---- daily water balance ---------------------------------------------------
# Read in blocks of 30 days to save memory. The blocks are put back together
# before the weekly sums, so no week is split.
chunk_days <- 30
daily_wb_list <- list()

print("Computing daily water balance in chunks...")
flush.console()
for (start_day in seq(1, n_days, by = chunk_days)) {
  end_day <- min(start_day + chunk_days - 1, n_days)

  print(paste("  Days", start_day, "to", end_day))
  flush.console()

  prec_chunk <- rast(prec_file, lyrs = start_day:end_day)
  et_chunk   <- rast(et_file, lyrs = start_day:end_day)

  prec_chunk <- crop(prec_chunk, ext(lsm_global))
  prec_chunk <- mask(prec_chunk, lsm_global)
  et_chunk   <- crop(et_chunk, ext(lsm_global))
  et_chunk   <- mask(et_chunk, lsm_global)

  wb_chunk <- prec_chunk - et_chunk

  daily_wb_list[[length(daily_wb_list) + 1]] <- wb_chunk

  rm(prec_chunk, et_chunk, wb_chunk)
  gc()
}

print("Combining daily water balance...")
flush.console()
wb_daily <- rast(daily_wb_list)
rm(daily_wb_list)
gc()

## ---- weekly aggregation ----------------------------------------------------
# Done on the whole year. week() numbers the weeks from 1 January, so the
# last week contains the last days of December.
print("Aggregating to weekly...")
flush.console()
dates_daily <- seq(as.Date(paste0(year, "-01-01")), by = "day", length.out = nlyr(wb_daily))
weeks_all <- week(dates_daily)

wb_weekly <- tapp(wb_daily, weeks_all, fun = sum, na.rm = TRUE)
rm(wb_daily)
gc()

n_weeks <- nlyr(wb_weekly)
dates_weekly <- seq(as.Date(paste0(year, "-01-01")), by = "week", length.out = n_weeks)
time(wb_weekly) <- dates_weekly

## ---- CWD -------------------------------------------------------------------
# ongoing_deficit is the running balance (negative when there is a deficit),
# one column per week. CWD is its absolute value, in mm.
print("Computing CWD values...")
flush.console()
cwd <- wb_weekly
values(cwd) <- 0

wb_values <- values(wb_weekly)
wettest_values <- values(wettest_month)

if (!is.null(previous_deficit)) {
  ongoing_deficit <- matrix(0, nrow = nrow(wb_values), ncol = ncol(wb_values))
  ongoing_deficit[, 1] <- previous_deficit + wb_values[, 1]
} else {
  ongoing_deficit <- wb_values
}

for (i in 1:n_weeks) {
  current_month <- month(dates_weekly[i])

  if (i > 1) {
    ongoing_deficit[, i] <- ongoing_deficit[, i - 1] + wb_values[, i]
  }

  # Reset to zero: the balance is positive, or the week is in the wettest
  # month.
  reset_mask <- ongoing_deficit[, i] > 0 | wettest_values == current_month
  ongoing_deficit[reset_mask, i] <- 0

  values(cwd)[, i] <- abs(ongoing_deficit[, i])
}

## ---- outputs ---------------------------------------------------------------
# The deficit of the last week is the starting value of the next year.
print("Saving final deficit state...")
flush.console()
saveRDS(ongoing_deficit[, n_weeks], file.path(out_path, paste0("final_deficit_", year, ".rds")))

print("Saving weekly CWD...")
flush.console()
cwd_file <- sprintf("%s/cwd_weekly_%d.nc", out_path, year)
writeCDF(cwd, filename = cwd_file, varname = "cwd", unit = "mm", overwrite = TRUE)

rm(wb_weekly, cwd, ongoing_deficit)
gc()

print(paste("Successfully processed year:", year))
flush.console()
