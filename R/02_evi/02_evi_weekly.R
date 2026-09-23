# =============================================================================
# EVI : WEEKLY MEAN EVI
# Original file name: evi_weekly_computation.R
# Methods 2.3
# -----------------------------------------------------------------------------
# Computes weekly mean EVI from the MOD13C1.061 EVI data, one file per year
# from 2000 to 2023. We need weekly EVI to match it with weekly CWD
# (01b_cwd_weekly.R).
#
# Weeks: we use lubridate::week(), which counts weeks from 1 January (not ISO
# weeks). The CWD script uses the same function and the same dates, so EVI week
# n and CWD week n are the same period. This is needed to match EVI and CWD in
# 04b_drydowns_add_evi.R.
#
# Values below zero are set to NA before the weekly mean.
#
# There are two ways to run the script, with the same result. By default, each
# year is read at once. With use_chunking = TRUE, the grid is split into
# n_chunks bands of rows, to use less memory.
#
# If a year fails, it is skipped and the other years are still processed.
#
# Input  : EVI.3600.1800.<year>.nc, daily MOD13C1.061 EVI at 0.1 degree
# Output : EVI_weekly_<year>.nc, variable evi, weekly mean EVI
#
# Run    : run the whole file. The last line runs main_process(2000, 2023).
#          The input and output paths are at the top of
#          process_global_evi_weekly(). Progress is written to
#          LOG-evi-processing-FIXED.out in the working folder.
# =============================================================================

## ---- environment -----------------------------------------------------------
# Libraries needed on our cluster, skipped if they do not exist.
for (lib in c("/opt/ohpc/pub/libs/hwloc/lib/libhwloc.so.15",
              "/opt/ohpc/pub/libs/gnu9/openmpi4/hdf5/1.10.8/lib/libhdf5_hl.so.100",
              "/opt/ohpc/pub/apps/gdal/3.5.1/lib/libgdal.so.31")) {
  if (file.exists(lib)) dyn.load(lib)
}

required_packages <- c("ncdf4", "raster", "terra", "sf", "dplyr", "lubridate", "remotes")
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

## ---- helpers ---------------------------------------------------------------
setup_logging <- function(filename) {
  options(immediate.write = TRUE)
  log_file <- file(filename, open = "wt")
  sink(log_file, type = "output")
  sink(log_file, type = "message")
  return(log_file)
}

manage_memory <- function() {
  mem_info <- gc()
  mem_used <- sum(mem_info[, 2])
  print(paste("Memory used:", round(mem_used, 2), "MB"))
  flush.console()
}

## ---- loop over the years ---------------------------------------------------
process_global_evi_weekly <- function(start_year, end_year, use_chunking = FALSE, n_chunks = 6) {
  print("Starting global EVI weekly processing...")
  flush.console()

  # Paths: change these lines to run the script elsewhere.
  input_path <- "/Net/Groups/BGI/data/DataStructureMDI/DATA/grid/Global/0d10_daily/MODIS/MOD13C1.061/Data/EVI/"
  output_path <- "/Net/Groups/BGI/scratch/mterristi/PhD/CWD_0_1/EVI/Weekly_computation"

  dir.create(output_path, showWarnings = FALSE, recursive = TRUE)

  for (year in start_year:end_year) {
    print(paste("\n=== Processing year:", year, "==="))
    flush.console()

    evi_file <- file.path(input_path, paste0("EVI.3600.1800.", year, ".nc"))

    if (!file.exists(evi_file)) {
      print(paste("WARNING: File not found:", evi_file))
      next
    }

    # If a year fails, the script goes on with the next one.
    tryCatch({
      print(paste("Loading EVI data for year:", year))
      flush.console()

      if (use_chunking) {
        print(paste("Using spatial chunking with", n_chunks, "chunks"))
        process_evi_chunked(evi_file, year, output_path, n_chunks)

      } else {
        evi_stack <- rast(evi_file)

        print("Processing invalid values...")
        evi_stack[evi_stack < 0] <- NA

        dates <- as.Date(time(evi_stack))
        print(paste("Date range:", min(dates), "to", max(dates)))

        # week() counts weeks from 1 January (not ISO weeks), as in the CWD
        # script, so the EVI and CWD weeks match.
        weeks <- week(dates)

        print(paste("Calculating weekly averages..."))
        flush.console()

        evi_weekly <- tapp(evi_stack, weeks, fun = mean, na.rm = TRUE)

        # Weekly dates starting on 1 January, as in 01b_cwd_weekly.R.
        dates_weekly <- seq(as.Date(paste0(year, "-01-01")), by = "week", length.out = nlyr(evi_weekly))
        time(evi_weekly) <- dates_weekly

        print(paste("Created", nlyr(evi_weekly), "weekly layers"))
        print(paste("Weekly date range:", min(dates_weekly), "to", max(dates_weekly)))
        print(paste("Sample dates:", paste(head(dates_weekly, 3), collapse = ", "), "...",
                    paste(tail(dates_weekly, 3), collapse = ", ")))

        print("Saving weekly EVI...")
        flush.console()
        output_file <- file.path(output_path, paste0("EVI_weekly_", year, ".nc"))
        writeCDF(evi_weekly,
                 filename = output_file,
                 varname = "evi",
                 longname = "Enhanced Vegetation Index (Weekly Mean)",
                 unit = "index",
                 overwrite = TRUE)

        print(paste("Saved:", output_file))

        rm(evi_stack, evi_weekly)
      }

      manage_memory()
      print(paste("✓ Successfully processed year:", year))

    }, error = function(e) {
      print(paste("✗ ERROR processing year", year, ":", e$message))
      flush.console()
    })
  }

  print("\n=== Processing complete ===")
}

## ---- version that uses less memory -----------------------------------------
# Same weekly means as above, computed band by band to use less memory. Each
# band is written into an empty raster of the full size, so the result is the
# same.
process_evi_chunked <- function(evi_file, year, output_path, n_chunks = 6) {

  # Read the grid size and the weeks once.
  evi_template <- rast(evi_file)
  dates <- as.Date(time(evi_template))
  weeks <- week(dates)
  n_weeks <- length(unique(weeks))

  dates_weekly <- seq(as.Date(paste0(year, "-01-01")), by = "week", length.out = n_weeks)

  full_extent <- ext(evi_template)
  nrows <- nrow(evi_template)
  chunk_size <- ceiling(nrows / n_chunks)

  weekly_template <- rast(nrows = nrows, ncols = ncol(evi_template),
                          nlyrs = n_weeks, extent = full_extent,
                          crs = crs(evi_template))
  time(weekly_template) <- dates_weekly

  rm(evi_template)
  gc()

  for (chunk in 1:n_chunks) {
    row_start <- (chunk - 1) * chunk_size + 1
    row_end <- min(chunk * chunk_size, nrows)

    print(paste("  Processing chunk", chunk, "of", n_chunks,
                "(rows", row_start, "to", row_end, ")"))
    flush.console()

    # Only this band is read into memory.
    chunk_ext <- ext(full_extent[1], full_extent[2],
                     full_extent[3] + (row_start - 1) * res(weekly_template)[2],
                     full_extent[3] + row_end * res(weekly_template)[2])

    evi_chunk <- rast(evi_file, win = chunk_ext)
    evi_chunk[evi_chunk < 0] <- NA

    evi_weekly_chunk <- tapp(evi_chunk, weeks, fun = mean, na.rm = TRUE)

    weekly_template[row_start:row_end, ] <- evi_weekly_chunk

    rm(evi_chunk, evi_weekly_chunk)
    gc()
  }

  output_file <- file.path(output_path, paste0("EVI_weekly_", year, ".nc"))
  writeCDF(weekly_template,
           filename = output_file,
           varname = "evi",
           longname = "Enhanced Vegetation Index (Weekly Mean)",
           unit = "index",
           overwrite = TRUE)

  print(paste("Saved chunked output:", output_file))
  print(paste("Date range:", min(dates_weekly), "to", max(dates_weekly)))
  rm(weekly_template)
}

## ---- run -------------------------------------------------------------------
main_process <- function(start_year, end_year, use_chunking = FALSE, n_chunks = 6) {
  log_file <- setup_logging("LOG-evi-processing-FIXED.out")

  tryCatch({
    print(paste("=== EVI Weekly Processing (FIXED DATE ASSIGNMENT) ==="))
    print(paste("Years:", start_year, "to", end_year))
    print(paste("Chunking:", ifelse(use_chunking, paste("YES (", n_chunks, "chunks)"), "NO")))
    print(paste("Start time:", Sys.time()))

    mem_info <- gc()
    print(paste("Initial memory used:", round(sum(mem_info[, 2]), 2), "MB"))
    flush.console()

    process_global_evi_weekly(start_year, end_year, use_chunking, n_chunks)

    print(paste("\nEnd time:", Sys.time()))
    manage_memory()

  }, finally = {
    sink(type = "output")
    sink(type = "message")
    close(log_file)
  })
}

# Study period. Use use_chunking = TRUE if there is not enough memory.
main_process(2000, 2023, use_chunking = FALSE, n_chunks = 6)
