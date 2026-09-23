# =============================================================================
# YEARLY STATES, STEP A of C : SPLIT THE GRID CELLS INTO GROUPS
# Original file name: Percentiles_per_year_chunk_assignments_v2.R
# Methods 2.6
# -----------------------------------------------------------------------------
# Splits the grid cells of each zone into groups of about 5000, so that step B
# can run as many separate jobs. Grid cells are assigned to groups at random,
# with a fixed seed, so the groups are the same every time and each group is a
# mix of the whole zone.
#
# Seven zones are used: Polar is not analysed. Arid Hot is included here and
# removed later in the figure scripts.
#
# Input  : <zone>_cwd_drydowns_with_evi_lc.Rdata
#            from 04c_drydowns_landcover_mask.R
# Output : per_year_chunked/<zone>/<zone>_chunk_assignments.Rdata
#            (object unique_gridcells: lat, lon, chunk_id)
#          per_year_chunked/<zone>/<zone>_chunk_summary.txt
#
# Run    : Rscript 05a_per_year_chunks.R
#          Then run step B once per zone and group, and step C once.
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

library(terra)
library(data.table)

input_path  <- "/Net/Groups/BGI/scratch/mterristi/PhD/CWD_0_1/New_Climate_Zones_percentiles/Drydowns/"
output_path <- "/Net/Groups/BGI/scratch/mterristi/PhD/CWD_0_1/New_Climate_Zones_percentiles/Drydowns/Drydowns_percentiles_per_year_v2/per_year_chunked"
dir.create(output_path, showWarnings = FALSE, recursive = TRUE)

# Group size. Zones have from tens of thousands to several hundred thousand
# grid cells, so this gives between a few and a few dozen jobs per zone.
TARGET_GRIDCELLS_PER_CHUNK <- 5000

zones <- c("tropical", "arid_hot", "arid_cold",
           "cold_dry", "cold_humid", "temperate_dry", "temperate_humid")

for (zone in zones) {
  cat("\n===================================\n")
  cat("STEP 2a — Zone:", toupper(zone), "\n")
  cat("===================================\n")
  
  data_file <- file.path(input_path, paste0(zone, "_cwd_drydowns_with_evi_lc.Rdata"))
  if (!file.exists(data_file)) {
    cat("Skipping - LC-masked file not found. Run Step 0 first.\n"); next
  }
  
  load(data_file)
  setDT(extracted_drydowns)
  cat("Rows:", nrow(extracted_drydowns), "\n")
  
  unique_gridcells <- unique(extracted_drydowns[, .(lat, lon)])
  n_gridcells      <- nrow(unique_gridcells)
  cat("Unique grid cells:", n_gridcells, "\n")
  
  n_chunks <- ceiling(n_gridcells / TARGET_GRIDCELLS_PER_CHUNK)
  cat("Creating", n_chunks, "chunks (~", round(n_gridcells / n_chunks), "gridcells each)\n")
  
  # Fixed seed: the groups are the same every time.
  set.seed(42)
  unique_gridcells[, chunk_id := sample(1:n_chunks, .N, replace = TRUE)]
  
  zone_chunk_dir <- file.path(output_path, zone)
  dir.create(zone_chunk_dir, showWarnings = FALSE, recursive = TRUE)
  
  chunk_file <- file.path(zone_chunk_dir, paste0(zone, "_chunk_assignments.Rdata"))
  save(unique_gridcells, file = chunk_file)
  cat("Saved:", chunk_file, "\n")
  
  summary_file <- file.path(zone_chunk_dir, paste0(zone, "_chunk_summary.txt"))
  writeLines(c(
    paste("Zone:", toupper(zone)),
    paste("Date:", Sys.time()),
    paste("Total grid cells:", n_gridcells),
    paste("Number of chunks:", n_chunks)
  ), summary_file)
  
  rm(extracted_drydowns); gc()
}

cat("\nChunking complete.\n")
