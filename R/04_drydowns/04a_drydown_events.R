# =============================================================================
# DRYDOWNS : DETECTION AND LARGEST DRYDOWN OF EACH YEAR, ALL CLIMATE ZONES
# Original file name: Drydowns_all_zones.R
# Methods 2.5 (drydown definition) and 2.6 (largest drydown of each year)
# -----------------------------------------------------------------------------
# For each of the eight climate zones from 03_cwd_by_climate_zone.R, the
# script does three things.
#
# 1. Drydowns. For each grid cell, a drydown is a series of consecutive weeks
#    with CWD > 0. It ends at the first week with CWD = 0, or at the end of the
#    data. For each drydown we save its start and end week, its duration in
#    weeks, its maximum CWD and its year (the year of the start week).
#    A drydown that goes over 31 December keeps the year in which it started.
#
# 2. Largest drydown. For each grid cell and year, we keep the drydown with the
#    highest maximum CWD. If two drydowns have exactly the same value, both are
#    kept.
#
# 3. Weekly data. We extract the weekly CWD values of these drydowns, with the
#    year of the drydown. These data are used in the next steps.
#    04b_drydowns_add_evi.R then replaces this year by the calendar year of
#    each week.
#
# Years: we use 1999 to 2023, one year more than the study period (2000 to
# 2023), so that a drydown already going on in January 2000 is captured from
# its start.
#
# The weeks are sorted by the time column. Its dates were created in
# time order in 03_cwd_by_climate_zone.R, so after sorting the weeks
# follow each other in time.
#
# Input  : cwd_<zone>.Rdata (object cwd_df) from 03_cwd_by_climate_zone.R
# Output : for each zone,
#          <zone>_drydown_events.Rdata (object drydown_events):
#            lat, lon, start_date, end_date, duration, max_deficit, year
#          <zone>_cwd_drydowns.Rdata (object extracted_drydowns):
#            lat, lon, time, cwd, year
#          <zone>_drydown_log.txt, progress log
#
# Run    : one job, zones one after the other:
#            Rscript 04a_drydown_events.R
#          Cold Humid has about 930 million weekly rows, so the job needs a
#          node with a lot of memory.
# =============================================================================

## ---- paths (change these lines to run the script elsewhere) ----------------
base_in_path  <- "/Net/Groups/BGI/scratch/mterristi/PhD/CWD_0_1/New_Climate_Zones/CWD_cz"
base_out_path <- "/Net/Groups/BGI/scratch/mterristi/PhD/CWD_0_1/New_Climate_Zones/Drydowns"

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

dir.create(base_out_path, recursive = TRUE, showWarnings = FALSE)

# Zone names as in the file names of 03_cwd_by_climate_zone.R.
zones <- c("tropical", "arid_hot", "arid_cold", "temperate_dry",
           "temperate_humid", "cold_dry", "cold_humid", "polar")

write_log <- function(message) {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  log_message <- paste(timestamp, "-", message)
  cat(log_message, "\n")
  write(log_message, log_file, append = TRUE)
}

## ---- drydown detection -----------------------------------------------------
# Runs on the weekly series of each grid cell (by lat, lon), sorted in time.
# A drydown starts at the first week with CWD > 0 and ends at the last week
# before CWD goes back to zero, or at the last week of the data.
analyze_drought_events <- function(df) {
  df[, {
    in_event = FALSE
    event_start_index = NULL
    events = list()

    for (i in 1:.N) {
      if (cwd[i] > 0) {
        if (!in_event) {
          in_event = TRUE
          event_start_index = i
        }
        if (i == .N || (i < .N && cwd[i + 1] <= 0)) {
          in_event = FALSE
          event_end_index = i
          events[[length(events) + 1]] = .(
            start_date = as.Date(time[event_start_index]),
            end_date = as.Date(time[event_end_index]),
            duration = as.integer(event_end_index - event_start_index + 1),
            max_deficit = max(cwd[event_start_index:event_end_index]),
            year = as.integer(year(time[event_start_index]))
          )
        }
      }
    }

    if (length(events) > 0) {
      rbindlist(events, fill = TRUE)
    } else {
      # Grid cells without any drydown return an empty table, so that the
      # results can still be combined.
      data.table(start_date = as.Date(character()), end_date = as.Date(character()),
                 duration = integer(), max_deficit = numeric(), year = integer())
    }
  }, by = .(lat, lon)]
}

## ---- run, one zone after the other -----------------------------------------
for (target_zone in zones) {

  log_file <- file.path(base_out_path, paste0(target_zone, "_drydown_log.txt"))

  in_file <- file.path(base_in_path, paste0("cwd_", target_zone, ".Rdata"))

  if (!file.exists(in_file)) {
    stop(paste("File not found for", target_zone, "at path:", in_file))
  }

  write_log(paste("Starting CLIMATE drydown analysis for zone:", target_zone))

  load(in_file)
  setDT(cwd_df)

  # One year before the study period (see header).
  write_log("Filtering data for years 1999 to 2023...")
  cwd_df <- cwd_df[year(time) >= 1999 & year(time) <= 2023]
  write_log(paste("Data filtered. Remaining rows:", nrow(cwd_df)))

  setkey(cwd_df, lat, lon, time)

  write_log("Analyzing drought events...")
  drydown_events <- analyze_drought_events(cwd_df)
  write_log(paste("Total drought events identified:", nrow(drydown_events)))

  save(drydown_events, file = file.path(base_out_path, paste0(target_zone, "_drydown_events.Rdata")))

  ## ---- largest drydown per grid cell and year ------------------------------
  write_log("Finding maximum deficit events per year...")
  max_deficits <- drydown_events[drydown_events[, .I[max_deficit == max(max_deficit)], by = .(year, lat, lon)]$V1]

  ## ---- weekly data of the selected drydowns --------------------------------
  write_log("Extracting rows for maximum deficit events...")
  cwd_df[, date_tmp := as.Date(time)]

  # Keep the weeks between the start and end of each selected drydown, for the
  # same grid cell. The year is the year of the drydown (i.year), not the year
  # of the week.
  extracted_drydowns <- cwd_df[
    max_deficits,
    on = .(lat, lon, date_tmp >= start_date, date_tmp <= end_date),
    nomatch = 0L,
    .(lat, lon, time, cwd, year = i.year)
  ]

  write_log(paste("Saving final data for", target_zone))
  save(extracted_drydowns, file = file.path(base_out_path, paste0(target_zone, "_cwd_drydowns.Rdata")))

  write_log(paste("Completed CLIMATE analysis for", target_zone))

  rm(cwd_df, drydown_events, max_deficits, extracted_drydowns); gc()
}
