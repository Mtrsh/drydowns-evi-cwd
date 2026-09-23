# =============================================================================
# DRYDOWNS : ADD WEEKLY EVI, ALL CLIMATE ZONES
# Original file name: Merging_drydowns_evi_cz.R
# Methods 2.6
# -----------------------------------------------------------------------------
# For each week of the selected drydowns, adds the EVI of the same grid cell
# and the same week. The resulting (CWD, EVI) pairs are used to compute the
# peak and low states.
#
# EVI and CWD are matched by date. The drydown weeks (03_cwd_by_climate_zone.R)
# and the EVI weeks (02_evi_weekly.R) use the same weekly dates starting on
# 1 January, so week n is the same period in both. If no EVI is found for a
# week, a message is written in the log and EVI is left as NA. This should not
# happen, since the dates are built the same way.
#
# EVI is read at the lon, lat of each row, which is the centre of the grid
# cell.
#
# Years: EVI is only available for 2000 to 2023. The 1999 weeks from the
# drydown step are kept with EVI = NA and are removed in
# 04c_drydowns_landcover_mask.R.
#
# Year column: in 04a_drydown_events.R, the year is the start year of the
# drydown. Here it is replaced by the calendar year of each week
# (year := year(time)), and this is the year used in all the next scripts.
# A drydown that goes over 31 December is therefore split between two years
# from this point.
#
# Input  : <zone>_cwd_drydowns.Rdata from 04a_drydown_events.R
#            (object extracted_drydowns: lat, lon, time, cwd, year)
#          EVI_weekly_<year>.nc, 2000 to 2023, from 02_evi_weekly.R
# Output : <zone>_cwd_drydowns_with_evi.Rdata (object extracted_drydowns:
#            lat, lon, time, cwd, year, evi)
#          log_<zone>_evi.txt, progress log with the share of weeks with EVI
#
# Run    : two jobs:
#            Rscript 04b_drydowns_add_evi.R          the seven smaller zones
#            Rscript 04b_drydowns_add_evi.R large    Cold Humid alone
#          Cold Humid is run alone because it is much larger than the other
#          zones. If one zone fails, the script goes on with the next one.
#          The paths are at the top of add_evi_to_drydowns() below.
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

library(terra)
library(data.table)
library(lubridate)
library(dplyr)

write_log <- function(message, log_file) {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  log_message <- paste(timestamp, "-", message)
  cat(log_message, "\n")
  write(log_message, log_file, append = TRUE)
}

## ---- one zone --------------------------------------------------------------
add_evi_to_drydowns <- function(zone_name) {

  # Paths: change these lines to run the script elsewhere.
  drydowns_path <- "/Net/Groups/BGI/scratch/mterristi/PhD/CWD_0_1/New_Climate_Zones/Drydowns"
  evi_path      <- "/Net/Groups/BGI/scratch/mterristi/PhD/CWD_0_1/EVI/Weekly_computation/"
  out_path      <- "/Net/Groups/BGI/scratch/mterristi/PhD/CWD_0_1/New_Climate_Zones/Drydowns/"
  dir.create(out_path, recursive = TRUE, showWarnings = FALSE)

  log_file <- file.path(out_path, paste0("log_", zone_name, "_evi.txt"))
  write_log(paste("Starting EVI extraction for zone:", zone_name), log_file)

  drydown_file <- file.path(drydowns_path, paste0(zone_name, "_cwd_drydowns.Rdata"))
  if (!file.exists(drydown_file)) stop(paste("File not found:", drydown_file))
  write_log(paste("Loading:", drydown_file), log_file)
  load(drydown_file)

  setDT(extracted_drydowns)
  write_log(paste("Rows loaded:", nrow(extracted_drydowns)), log_file)

  # time is a factor with dates as text; convert it to Date.
  if (!inherits(extracted_drydowns$time, "Date")) {
    extracted_drydowns[, time := as.Date(time, origin = "1970-01-01")]
  }

  # Use the calendar year of each week (see header).
  extracted_drydowns[, year := year(time)]
  extracted_drydowns[, evi  := NA_real_]

  all_years <- sort(unique(extracted_drydowns$year))
  years     <- all_years[all_years >= 2000 & all_years <= 2023]

  write_log(paste("Data years:", min(all_years), "to", max(all_years)), log_file)
  write_log(paste("Processing years with EVI:", min(years), "to", max(years)), log_file)

  for (yr in years) {
    write_log(paste("Processing year", yr), log_file)

    evi_file <- file.path(evi_path, paste0("EVI_weekly_", yr, ".nc"))
    if (!file.exists(evi_file)) {
      write_log(paste("WARNING: EVI file not found for year", yr), log_file)
      next
    }

    evi_year  <- rast(evi_file)
    evi_dates <- as.Date(time(evi_year))

    year_data    <- extracted_drydowns[year == yr]
    unique_dates <- sort(unique(year_data$time))

    write_log(paste(" ", length(unique_dates), "unique dates,", nlyr(evi_year), "EVI layers"), log_file)

    # For each drydown week, take the EVI layer with the same date and read it
    # at all grid cells with a drydown that week.
    for (date_obj in unique_dates) {
      date_obj   <- as.Date(date_obj, origin = "1970-01-01")
      date_match <- which(evi_dates == date_obj)

      if (length(date_match) == 0) {
        write_log(paste("  WARNING: No EVI layer for date", date_obj), log_file)
        next
      }

      evi_layer   <- evi_year[[date_match]]
      date_coords <- year_data[time == date_obj, .(lon, lat)]
      if (nrow(date_coords) == 0) next

      evi_values <- extract(evi_layer, date_coords[, .(lon, lat)],
                            method = "simple", ID = FALSE)

      extracted_drydowns[year == yr & time == date_obj, evi := evi_values[[1]]]
    }

    write_log(paste("  Completed year", yr), log_file)
  }

  # Share of weeks with EVI, written to the log.
  n_total    <- nrow(extracted_drydowns)
  n_evi_yrs  <- nrow(extracted_drydowns[year >= 2000 & year <= 2023])
  n_with_evi <- sum(!is.na(extracted_drydowns$evi))
  pct        <- round(100 * n_with_evi / n_evi_yrs, 2)

  write_log(paste("Total rows:", n_total), log_file)
  write_log(paste("Rows in EVI period:", n_evi_yrs), log_file)
  write_log(paste("Rows with EVI:", n_with_evi, "(", pct, "%)"), log_file)

  out_file <- file.path(out_path, paste0(zone_name, "_cwd_drydowns_with_evi.Rdata"))
  save(extracted_drydowns, file = out_file)
  write_log(paste("Saved:", out_file), log_file)

  return(extracted_drydowns)
}

## ---- run -------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
mode <- args[1]

climate_zones <- if (mode == "large") {
  c("cold_humid")
} else {
  c("polar", "tropical", "arid_hot", "arid_cold", "cold_dry", "temperate_dry", "temperate_humid")
}

for (zone in climate_zones) {
  cat("\n================================\n")
  cat("Processing zone:", zone, "\n")
  cat("================================\n\n")

  tryCatch({
    result <- add_evi_to_drydowns(zone)
    cat("Successfully processed:", zone, "\n")
    rm(result); gc()
  }, error = function(e) {
    cat("ERROR processing", zone, ":", conditionMessage(e), "\n")
  })
}

cat("\nAll climate zones processed!\n")
