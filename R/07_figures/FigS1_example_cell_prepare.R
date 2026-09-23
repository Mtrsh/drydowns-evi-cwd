# =============================================================================
# FIGURE S1, STEP 1 of 2 : DATA FOR THE EXAMPLE GRID CELL
# Original file name: Prepare_Fig_drydown.R
# Methods 2.6 (example grid cell in Fig. S1)
# -----------------------------------------------------------------------------
# Extracts the data needed for the example figure for one grid cell and saves
# them in a small .rds file, so that the figure can be drawn on any computer
# without the large zone tables (step 2, FigS1_example_cell_plot.R).
#
# Choice of the grid cell: either set by hand with PIXEL_LAT and PIXEL_LON, or
# chosen automatically. In that case we take the grid cells with data for at
# least 80 % of the weeks from 2000 to 2023 and a maximum CWD of at least
# 30 mm, and keep the one whose maximum CWD is closest to the median of these
# grid cells. If no grid cell meets both conditions, we use the 50 grid cells
# with the most data instead.
#
# Drydowns are found as in 04a_drydown_events.R: consecutive weeks with
# CWD > 0, with the year of the first week, and the largest drydown of each
# year is kept.
#
# Input  : cwd_<zone>.Rdata from 03_cwd_by_climate_zone.R
#          EVI_weekly_<year>.nc, 2000 to 2023
# Output : FigureS1_data_<zone>.rds, a list with the weekly CWD and EVI of the
#          grid cell, the largest drydown of each year (year, start, end, date
#          and value of the CWD peak) and the EVI at each peak
#
# Run    : Rscript FigS1_example_cell_prepare.R temperate_humid
#          Without argument, the zone is temperate_humid, as in the paper.
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
library(data.table); library(terra)
## ---- settings and paths (change these lines to run the script elsewhere) ---
args <- commandArgs(trailingOnly = TRUE)
zone <- if (length(args) >= 1) args[1] else "temperate_humid"
cwd_file    <- file.path("/Net/Groups/BGI/scratch/mterristi/PhD/CWD_0_1/New_Climate_Zones/CWD_cz",
                         paste0("cwd_", zone, ".Rdata"))
evi_dir     <- "/Net/Groups/BGI/scratch/mterristi/PhD/CWD_0_1/EVI/Weekly_computation"
out_dir     <- "/Net/Groups/BGI/scratch/mterristi/PhD/CWD_0_1/New_Climate_Zones_percentiles/Drydowns/Figures"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
out_file    <- file.path(out_dir, paste0("FigureS1_data_", zone, ".rds"))
YEAR0 <- 2000L; YEAR1 <- 2023L
MIN_PEAK  <- 0          # mm; drydowns whose peak CWD is below this are ignored (0: none)
PIXEL_LAT <- NA_real_   # set both to pin a cell; NA = auto-select representative
PIXEL_LON <- NA_real_
## ---- weekly CWD of the zone ------------------------------------------------
cat("Loading", basename(cwd_file), "...\n")
nm <- load(cwd_file)
dt <- NULL; best <- -1L; picked <- NA_character_
for (o in nm) { x <- get(o); if (is.data.frame(x) && nrow(x) > best) { dt <- as.data.table(x); best <- nrow(x); picked <- o } }
if (is.null(dt)) stop("No data.frame found in ", basename(cwd_file))
cat("Loaded object:", picked, "| columns:", paste(names(dt), collapse = ", "), "| rows:", nrow(dt), "\n")
nmv  <- names(dt)
dcol <- nmv[vapply(dt, function(x) inherits(x, c("Date","IDate","POSIXct","POSIXt")), logical(1))]
if (!length(dcol)) dcol <- grep("date|time", nmv, ignore.case = TRUE, value = TRUE)
ccol <- grep("cwd|deficit", nmv, ignore.case = TRUE, value = TRUE); ccol <- ccol[vapply(ccol, function(cc) is.numeric(dt[[cc]]), logical(1))]
latc <- grep("^lat|latitude",  nmv, ignore.case = TRUE, value = TRUE)
lonc <- grep("^lon|longitude", nmv, ignore.case = TRUE, value = TRUE)
if (!length(dcol) || !length(ccol) || !length(latc) || !length(lonc))
  stop("Could not map columns from: ", paste(nmv, collapse = ", "))
setnames(dt, c(dcol[1], ccol[1], latc[1], lonc[1]), c("time", "CWD", "lat", "lon"))
setkey(dt, lat, lon)
## ---- choice of the grid cell (see header) ----------------------------------
if (is.na(PIXEL_LAT) || is.na(PIXEL_LON)) {
  cat("Selecting a representative pixel...\n")
  cov  <- dt[, .(n = .N, maxcwd = max(CWD, na.rm = TRUE)), by = .(lat, lon)]
  cand <- cov[n >= 0.8 * (YEAR1 - YEAR0 + 1) * 52 & maxcwd >= 30]
  if (!nrow(cand)) cand <- cov[order(-n, -maxcwd)][seq_len(min(50, .N))]
  sel <- cand[which.min(abs(maxcwd - median(cand$maxcwd)))]
  PIXEL_LAT <- sel$lat; PIXEL_LON <- sel$lon
  cat(sprintf("Auto-selected pixel: lat %.4f, lon %.4f (n=%d weeks, peakCWD=%.0f mm)\n",
              sel$lat, sel$lon, sel$n, sel$maxcwd))
}
px <- dt[.(PIXEL_LAT, PIXEL_LON), nomatch = 0L]
if (!nrow(px)) {
  u <- unique(dt[, .(lat, lon)]); u[, d2 := (lat-PIXEL_LAT)^2 + (lon-PIXEL_LON)^2]
  s <- u[which.min(d2)]; px <- dt[.(s$lat, s$lon), nomatch = 0L]
  PIXEL_LAT <- s$lat; PIXEL_LON <- s$lon
}
px[, Date := as.Date(time)]
px <- px[Date >= as.Date(paste0(YEAR0,"-01-01")) & Date <= as.Date(paste0(YEAR1,"-12-31"))]
ts_cwd <- unique(px[order(Date), .(Date, CWD)], by = "Date")
## ---- drydowns: consecutive weeks with CWD > 0, largest drydown per year ----
v <- ts_cwd[order(Date)]
v[, indry := CWD > 0]; v[, run := rleid(indry)]
dd <- v[indry == TRUE, .(start = min(Date), end = max(Date),
                         peak_date = Date[which.max(CWD)], peak_cwd = max(CWD),
                         year = year(min(Date))), by = run]
dd <- dd[peak_cwd >= MIN_PEAK]
bands <- dd[dd[, .I[which.max(peak_cwd)], by = year]$V1][order(year)]
## ---- weekly EVI of the grid cell -------------------------------------------
cat("Extracting weekly EVI at the pixel...\n")
evi_px <- rbindlist(lapply(YEAR0:YEAR1, function(yr) {
  f <- file.path(evi_dir, paste0("EVI_weekly_", yr, ".nc"))
  if (!file.exists(f)) { cat("  WARNING: missing", basename(f), "\n"); return(NULL) }
  r <- rast(f)
  vv <- terra::extract(r, data.frame(lon = PIXEL_LON, lat = PIXEL_LAT), method = "simple", ID = FALSE)
  data.table(Date = as.Date(terra::time(r)), EVI = as.numeric(vv[1, ]))
}), fill = TRUE)
evi_px <- evi_px[!is.na(Date)][order(Date)]
## ---- EVI at each CWD peak, and save ----------------------------------------
peaks <- merge(bands[, .(Date = peak_date, CWD = peak_cwd)],
               evi_px[, .(Date, EVI)], by = "Date", all.x = TRUE)
figdat <- list(zone = zone, lat = PIXEL_LAT, lon = PIXEL_LON,
               years = c(YEAR0, YEAR1),
               cwd = ts_cwd, evi = evi_px,
               bands = bands[, .(year, start, end, peak_date, peak_cwd)], peaks = peaks)
saveRDS(figdat, out_file)
cat("\nKept", nrow(bands), "annual-maximum drydowns.\n")
cat("Saved figure data ->", out_file, "\n")
cat("Now run FigureS1_plot.R locally on this .rds.\n")