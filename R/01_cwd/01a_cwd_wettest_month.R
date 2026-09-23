# =============================================================================
# CWD, STEP 1 of 2 : WETTEST MONTH OF EACH GRID CELL
# Original file name: cwd01_climatology.R
# Methods 2.5
# -----------------------------------------------------------------------------
# The cumulative water deficit (CWD) is the running sum of P - ET. It goes back
# to zero when rain makes up for the deficit. In very dry or very seasonal
# regions this may not happen for several years, and CWD then keeps growing.
# To avoid this, we also reset CWD once a year, in the wettest month of each
# grid cell. This script finds that month.
#
# How: we compute daily P - ET, sum it by week, and add up the weekly values
# for each calendar month over all years. The wettest month is the month with
# the largest total.
#
# Each year is read in blocks of 30 days to save memory, and the weekly sums
# are computed within each block. A week that falls across two blocks is
# therefore counted as two partial weeks. A week that falls across two months
# is counted in both months. (In step 2, the weekly sums are computed on the
# whole year at once.)
#
# Years: 1980 to 2023, all years available in both MSWEP and GLEAM. We use
# more years than the study period (2000 to 2023) to get a more reliable
# seasonal cycle.
#
# Area: land only, from 60 S to 90 N, at 0.1 degree. The land-sea mask is the
# GLDAS mask at 0.25 degree, resampled to 0.1 degree.
#
# Input  : MSWEP v2.8 daily precipitation, one file per year
#          GLEAM v4.2a daily evapotranspiration, one file per year
#          GLDAS land-sea mask at 0.25 degree
# Output : wettest_month.nc, month number (1 to 12) for each land grid cell,
#          used in step 2 (01b_cwd_weekly.R)
#
# Run    : one job, years processed one after the other. Progress is written
#          to LOG1-climatology.out in the output folder.
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

setwd(out_path)

options(immediate.write = TRUE)
out <- file("LOG1-climatology.out", open = "wt")
sink(out, type = "output")
sink(out, type = "message")

## ---- land-sea mask on the 0.1 degree grid ----------------------------------
# All data are cropped and masked to this grid.
print("Creating 0.1 degree land-sea mask...")
lsm_025 <- rast(path_lsm)
template_01 <- rast(xmin = -180, xmax = 180, ymin = -60, ymax = 90,
                    resolution = 0.1, crs = "EPSG:4326")
lsm_01 <- resample(lsm_025, template_01, method = "near")
mask_raster <- ifel(lsm_01 == 0, NA, 1)   # ocean to NA, land to 1
lsm_global <- crop(mask_raster, ext(-180, 180, -60, 90))

## ---- sum the water balance by month ----------------------------------------
# monthly_sums[[m]] is the total of weekly P - ET for month m, summed over
# all years.
print("Calculating climatological wettest month...")
monthly_sums <- vector("list", 12)
for (m in 1:12) {
  monthly_sums[[m]] <- lsm_global * 0
}

chunk_days <- 30   # days read at once; a whole year of daily global fields does
                   # does not fit in memory

for (year in 1980:2023) {
  print(paste("Processing year", year, "for climatology..."))
  flush.console()

  prec_file <- file.path(path_prec, paste0("prec.MSWEPv280.3600.1800.", year, ".nc"))
  et_file   <- file.path(path_et, year, paste0("E_", year, "_GLEAM_v4.2a.nc"))

  # If one of the two files is missing, the whole year is skipped.
  if (!file.exists(prec_file) || !file.exists(et_file)) {
    print(paste("WARNING: Files missing for year", year))
    next
  }

  nc <- nc_open(prec_file)
  n_days <- nc$dim$time$len
  nc_close(nc)

  for (start_day in seq(1, n_days, by = chunk_days)) {
    end_day <- min(start_day + chunk_days - 1, n_days)

    print(paste("  Processing days", start_day, "to", end_day))
    flush.console()

    prec_chunk <- rast(prec_file, lyrs = start_day:end_day)
    et_chunk   <- rast(et_file, lyrs = start_day:end_day)

    prec_chunk <- crop(prec_chunk, ext(lsm_global))
    prec_chunk <- mask(prec_chunk, lsm_global)
    et_chunk   <- crop(et_chunk, ext(lsm_global))
    et_chunk   <- mask(et_chunk, lsm_global)

    # Water balance: positive when P is larger than ET.
    wb_chunk <- prec_chunk - et_chunk

    dates_chunk <- seq(as.Date(paste0(year, "-01-01")) + start_day - 1,
                       length.out = end_day - start_day + 1, by = "day")
    weeks_chunk <- week(dates_chunk)

    # Weekly sums, same time step as CWD in step 2. tapp returns one layer
    # per week present in the block.
    wb_weekly_chunk <- tapp(wb_chunk, weeks_chunk, fun = sum, na.rm = TRUE)

    # The layers returned by tapp follow the sorted week numbers: layer i is
    # unique_weeks[i], not week i. This is why we use unique_weeks below.
    unique_weeks <- sort(unique(weeks_chunk))

    for (m in 1:12) {
      # A week that falls across two months is counted in both.
      weeks_in_month <- unique_weeks[sapply(unique_weeks, function(w) {
        any(month(dates_chunk[weeks_chunk == w]) == m)
      })]

      if (length(weeks_in_month) > 0) {
        layer_indices <- which(unique_weeks %in% weeks_in_month)

        if (length(layer_indices) > 0 && length(layer_indices) <= nlyr(wb_weekly_chunk)) {
          monthly_sums[[m]] <- monthly_sums[[m]] +
            sum(wb_weekly_chunk[[layer_indices]], na.rm = TRUE)
        }
      }
    }

    rm(prec_chunk, et_chunk, wb_chunk, wb_weekly_chunk)
    gc()
  }
}

## ---- wettest month and output ----------------------------------------------
# which.max gives the layer with the largest value, i.e. the month number.
print("Finding wettest month...")
flush.console()
monthly_stack <- rast(monthly_sums)
wettest_month <- which.max(monthly_stack)

print("Saving climatological wettest month...")
flush.console()
writeCDF(wettest_month,
         filename = file.path(out_path, "wettest_month.nc"),
         varname = "wettest_month",
         unit = "month",
         overwrite = TRUE)

print("Climatology calculation complete!")
flush.console()
