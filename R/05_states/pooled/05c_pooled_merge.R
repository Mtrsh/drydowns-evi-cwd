# =============================================================================
# POOLED STATES, STEP C of C : MERGE THE GROUPS
# Original file name: Percentiles_pooled_merge_v2.R
# Methods 2.6
# -----------------------------------------------------------------------------
# Merges the results of step B into one table per zone and one table for all
# zones. If a group was run twice, the duplicates are removed (by lat, lon).
# Zone names are changed to the form used in the figures ("arid_hot" becomes
# "Arid Hot").
#
# Input  : pooled_chunked/<zone>/results/<zone>_chunk_*_results.Rdata
# Output : <zone>/<zone>_pooled_percentiles_summary_v2.Rdata
#          combined_pooled_percentiles_summary_v2.Rdata (object all_summaries).
#          This combined file is used by the figure and table scripts (states by
#          zone and regime, sensitivity, regime shifts, supplementary tables).
#
# Run    : Rscript 05c_pooled_merge.R, when all step B jobs are finished.
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

chunk_path        <- "/Net/Groups/BGI/scratch/mterristi/PhD/CWD_0_1/New_Climate_Zones_percentiles/Drydowns/Drydowns_percentiles_pooled_v2/pooled_chunked"
final_output_path <- "/Net/Groups/BGI/scratch/mterristi/PhD/CWD_0_1/New_Climate_Zones_percentiles/Drydowns/Drydowns_percentiles_pooled_v2/"

zones <- c("tropical", "arid_hot", "arid_cold",
           "cold_dry", "cold_humid", "temperate_dry", "temperate_humid")

# Zone names as used in the figures.
zone_labels <- list(
  tropical        = "Tropical",
  arid_hot        = "Arid Hot",
  arid_cold       = "Arid Cold",
  cold_dry        = "Cold Dry",
  cold_humid      = "Cold Humid",
  temperate_dry   = "Temperate Dry",
  temperate_humid = "Temperate Humid"
)

per_zone_list <- list()

for (zone in zones) {
  cat("\n=== STEP 1d Merging:", toupper(zone), "===\n")
  
  zone_dir <- file.path(chunk_path, zone, "results")
  if (!dir.exists(zone_dir)) { cat("  no results dir, skipping\n"); next }
  
  result_files <- list.files(zone_dir,
                             pattern = paste0(zone, "_chunk_.*_results\\.Rdata$"),
                             full.names = TRUE)
  if (length(result_files) == 0) { cat("  no files\n"); next }
  cat("  Found", length(result_files), "chunk files\n")
  
  all_results <- lapply(result_files, function(f) {
    tryCatch({
      load(f)
      if (exists("percentiles_summary") && nrow(percentiles_summary) > 0) {
        return(percentiles_summary)
      } else NULL
    }, error = function(e) NULL)
  })
  all_results <- all_results[!sapply(all_results, is.null)]
  if (length(all_results) == 0) { cat("  no valid results\n"); next }
  
  percentiles_summary <- rbindlist(all_results, fill = TRUE)
  
  n_before <- nrow(percentiles_summary)
  percentiles_summary <- unique(percentiles_summary, by = c("Latitude", "Longitude"))
  if (n_before != nrow(percentiles_summary)) {
    cat("  Removed", n_before - nrow(percentiles_summary), "duplicates\n")
  }
  
  percentiles_summary[, climate_zone := zone_labels[[zone]]]
  cat("  Total pixels:", nrow(percentiles_summary), "\n")
  
  zone_output_dir <- file.path(final_output_path, zone)
  dir.create(zone_output_dir, showWarnings = FALSE, recursive = TRUE)
  output_file <- file.path(zone_output_dir, paste0(zone, "_pooled_percentiles_summary_v2.Rdata"))
  save(percentiles_summary, file = output_file)
  cat("  Saved:", output_file, "\n")
  
  n_after  <- sum(percentiles_summary$p05_position == "after",  na.rm = TRUE)
  n_bef    <- sum(percentiles_summary$p05_position == "before", na.rm = TRUE)
  n_eq     <- sum(percentiles_summary$p05_position == "equal",  na.rm = TRUE)
  cat("  after  p95:", n_after, sprintf(" (%.1f%%)\n", 100 * n_after / nrow(percentiles_summary)))
  cat("  before p95:", n_bef,   sprintf(" (%.1f%%)\n", 100 * n_bef   / nrow(percentiles_summary)))
  cat("  equal  p95:", n_eq,    sprintf(" (%.1f%%)\n", 100 * n_eq    / nrow(percentiles_summary)))
  
  per_zone_list[[zone]] <- percentiles_summary
  rm(percentiles_summary, all_results); gc(verbose = FALSE)
}

cat("\n=== Combining all zones ===\n")
all_summaries <- rbindlist(per_zone_list, fill = TRUE)
cat("Combined pixels:", nrow(all_summaries), "\n")

combined_file <- file.path(final_output_path, "combined_pooled_percentiles_summary_v2.Rdata")
save(all_summaries, file = combined_file)
cat("Saved combined:", combined_file, "\n")

cat("\n=== STEP 1d POOLED MERGE v2 DONE ===\n")
