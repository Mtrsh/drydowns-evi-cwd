# =============================================================================
# POOLED STATES, STEP B of C : PEAK AND LOW STATES
# Original file name: Percentiles_pooled_worker_v2.R
# Methods 2.6
# -----------------------------------------------------------------------------
# For each grid cell of one group, we take together all weekly (CWD, EVI) values of its largest drydowns
# from 2000 to 2023, and compute two states:
#
#   peak state  EVI_E95    median EVI of the 5 % of weeks with the highest EVI
#               CWD_atE95  median CWD of these same weeks
#   low state   EVI_E05    median EVI of the 5 % of weeks with the lowest EVI
#               CWD_atE05  median CWD of these same weeks
#
# 5 % means ceiling(0.05 n) weeks, and at least one week. We need at least five
# weeks of data. With fewer than 20 weeks, each state is a single week.
#
# The low state is taken from all weeks, with no condition on its position
# relative to the peak state. Their order along the CWD axis is saved in
# p05_position and gives the regime:
#   "after"   CWD_atE05 > CWD_atE95: EVI decreases as the deficit grows,
#             water-limited
#   "before"  CWD_atE05 < CWD_atE95: EVI increases as the deficit grows,
#             energy-limited
#   "equal"   both states are at the same CWD
#
# We also save the number of weeks, EVI_max, CWD_max and EVI_range.
#
# Input  : <zone>_cwd_drydowns_with_evi_lc.Rdata
#          <zone>_chunk_assignments.Rdata from step A
# Output : pooled_chunked/<zone>/results/<zone>_chunk_<id>_results.Rdata
#            (object percentiles_summary, one row per grid cell:
#            Latitude, Longitude, n_observations, EVI_E95, CWD_atE95,
#            EVI_E05, CWD_atE05, p05_position, EVI_max, CWD_max, EVI_range,
#            climate_zone)
#
# Run    : one job per zone and group:
#            Rscript 05b_pooled_states.R tropical 3
#          If one grid cell fails, it is skipped and the job goes on.
# =============================================================================

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

## ---- zone and group --------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("Provide: (1) zone (2) chunk_id")
zone     <- tolower(args[1])
CHUNK_ID <- as.integer(args[2])

## ---- paths (change these lines to run the script elsewhere) ----------------
input_path  <- "/Net/Groups/BGI/scratch/mterristi/PhD/CWD_0_1/New_Climate_Zones_percentiles/Drydowns/"
chunk_path  <- "/Net/Groups/BGI/scratch/mterristi/PhD/CWD_0_1/New_Climate_Zones_percentiles/Drydowns/Drydowns_percentiles_pooled_v2/pooled_chunked"
output_path <- file.path(chunk_path, zone, "results")
dir.create(output_path, showWarnings = FALSE, recursive = TRUE)

log_file <- file.path(chunk_path, zone, paste0(zone, "_chunk_", CHUNK_ID, "_log.txt"))
write_log <- function(message) {
  ts <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  msg <- paste(ts, "-", message)
  cat(msg, "\n"); write(msg, log_file, append = TRUE)
}
write_log(paste("STEP 1b POOLED v2 — chunk", CHUNK_ID, "zone:", zone))

## ---- grid cells of this group ----------------------------------------------
chunk_file <- file.path(chunk_path, zone, paste0(zone, "_chunk_assignments.Rdata"))
if (!file.exists(chunk_file)) stop("Chunk file not found: ", chunk_file)
load(chunk_file)
chunk_gridcells <- unique_gridcells[chunk_id == CHUNK_ID, .(lat, lon)]
write_log(paste("Grid cells in chunk:", nrow(chunk_gridcells)))
if (nrow(chunk_gridcells) == 0) stop("No grid cells in this chunk")

## ---- data of these grid cells ----------------------------------------------
data_file <- file.path(input_path, paste0(zone, "_cwd_drydowns_with_evi_lc.Rdata"))
if (!file.exists(data_file)) stop("LC-masked source not found: ", data_file)
load(data_file)
setDT(extracted_drydowns)
drydown_data <- extracted_drydowns[chunk_gridcells, on = .(lat, lon), nomatch = 0]
rm(extracted_drydowns); gc()
write_log(paste("Rows after chunk filter:", nrow(drydown_data)))

## ---- peak and low states (see header) --------------------------------------
process_pixel_percentiles <- function(data) {
  if (nrow(data) < 5) return(NULL)
  n <- nrow(data)
  
  # Peak state: the 5 % of weeks with the highest EVI. With equal values,
  # order() keeps the weeks in time order.
  n_top   <- max(1L, ceiling(n * 0.05))
  top_idx <- order(-data$EVI)[1:n_top]
  EVI_E95   <- median(data$EVI[top_idx], na.rm = TRUE)
  CWD_atE95 <- median(data$cwd[top_idx], na.rm = TRUE)
  if (!is.finite(EVI_E95) || !is.finite(CWD_atE95)) return(NULL)
  
  # Low state: the 5 % of weeks with the lowest EVI, taken from all weeks.
  n_bottom   <- max(1L, ceiling(n * 0.05))
  bottom_idx <- order(data$EVI)[1:n_bottom]
  EVI_E05   <- median(data$EVI[bottom_idx], na.rm = TRUE)
  CWD_atE05 <- median(data$cwd[bottom_idx], na.rm = TRUE)
  if (!is.finite(EVI_E05))   EVI_E05   <- NA_real_
  if (!is.finite(CWD_atE05)) CWD_atE05 <- NA_real_
  
  # Position of the low state compared with the peak state (the regime).
  p05_position <- NA_character_
  if (!is.na(CWD_atE05) && !is.na(CWD_atE95)) {
    p05_position <- ifelse(CWD_atE05 > CWD_atE95, "after",
                    ifelse(CWD_atE05 < CWD_atE95, "before", "equal"))
  }
  
  # Maximum values over the same weeks.
  EVI_max <- max(data$EVI, na.rm = TRUE)
  CWD_max <- max(data$cwd, na.rm = TRUE)
  EVI_range <- EVI_max - min(data$EVI, na.rm = TRUE)
  
  data.table(
    Latitude       = data$lat[1],
    Longitude      = data$lon[1],
    n_observations = n,
    EVI_E95        = EVI_E95,
    CWD_atE95      = CWD_atE95,
    EVI_E05        = EVI_E05,
    CWD_atE05      = CWD_atE05,
    p05_position   = p05_position,
    EVI_max        = EVI_max,
    CWD_max        = CWD_max,
    EVI_range      = EVI_range
  )
}

## ---- run on all grid cells of the group ------------------------------------
# .SD is copied with its keys, so that lat and lon are available as columns.
write_log("Starting pooled processing...")
n_total  <- uniqueN(drydown_data[, .(lat, lon)])
counter  <- 0L
last_log <- Sys.time()

results <- drydown_data[, {
  counter <<- counter + 1L
  if (counter %% 5000L == 0L || difftime(Sys.time(), last_log, units = "mins") > 10) {
    write_log(paste0("Progress: ", counter, " / ", n_total,
                     " (", round(100 * counter / n_total, 1), "%)"))
    last_log <<- Sys.time()
  }
  sd_with_keys <- copy(.SD)
  sd_with_keys[, `:=`(lat = .BY$lat, lon = .BY$lon)]
  list(result = list(tryCatch(process_pixel_percentiles(sd_with_keys), error = function(e) NULL)))
}, by = .(lat, lon)]

results_list <- results$result[!sapply(results$result, is.null)]
write_log(paste("Processed pixels:", length(results_list)))

## ---- save the results of this group ----------------------------------------
percentiles_summary <- rbindlist(results_list, fill = TRUE)
percentiles_summary[, climate_zone := zone]

output_file <- file.path(output_path, paste0(zone, "_chunk_", CHUNK_ID, "_results.Rdata"))
save(percentiles_summary, file = output_file)
write_log(paste("Saved:", output_file))
write_log(paste("Total rows:", nrow(percentiles_summary)))
write_log(paste("p05 after  p95:", sum(percentiles_summary$p05_position == "after",  na.rm = TRUE)))
write_log(paste("p05 before p95:", sum(percentiles_summary$p05_position == "before", na.rm = TRUE)))
write_log(paste("p05 equal  p95:", sum(percentiles_summary$p05_position == "equal",  na.rm = TRUE)))
write_log(paste("Chunk", CHUNK_ID, "done."))
