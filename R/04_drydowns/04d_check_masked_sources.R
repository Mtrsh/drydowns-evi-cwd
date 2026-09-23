# =============================================================================
# CHECK : SUMMARY OF THE FILES FROM THE LAND-COVER STEP
# Original file name: Verify_LC_masked_sources.R
# -----------------------------------------------------------------------------
# Not part of the processing. Prints, for each zone, the number of rows and
# grid cells, the years, the EVI range, the maximum CWD and the file size of
# the files from 04c_drydowns_landcover_mask.R, and then the totals. We used it
# to check that all zones were ready before computing the states.
#
# Seven zones are listed: Polar is not used from here on. Arid Hot is still
# included at this stage and is removed later in the figure scripts.
#
# Input  : <zone>_cwd_drydowns_with_evi_lc.Rdata
# Output : printed in the console
# Run    : Rscript 04d_check_masked_sources.R
# =============================================================================

## ---- paths (change these lines to run the script elsewhere) ----------------
input_path <- "/Net/Groups/BGI/scratch/mterristi/PhD/CWD_0_1/New_Climate_Zones_percentiles/Drydowns/"

## ---- environment -----------------------------------------------------------

custom_libs <- c(
  "/Net/Groups/BGI/scratch/mterristi/R_lodomeria/4.4",
  "/Net/Groups/BSI/work_scratch/quincy/model/software/r_packages/r_4.4.x"
)
.libPaths(c(custom_libs[dir.exists(custom_libs)], .libPaths()))

library(data.table)

zones <- c("tropical", "arid_hot", "arid_cold",
           "cold_dry", "cold_humid", "temperate_dry", "temperate_humid")

cat("\n===========================================\n")
cat("LC-MASKED SOURCES — VERIFICATION SUMMARY\n")
cat("===========================================\n\n")

summary_list <- list()

for (zone in zones) {
  f <- file.path(input_path, paste0(zone, "_cwd_drydowns_with_evi_lc.Rdata"))
  
  if (!file.exists(f)) {
    cat(sprintf("%-20s : FILE NOT FOUND\n", zone))
    next
  }
  
  load(f)
  setDT(extracted_drydowns)
  
  summary_list[[zone]] <- data.table(
    zone        = zone,
    rows        = nrow(extracted_drydowns),
    pixels      = uniqueN(extracted_drydowns[, .(lat, lon)]),
    years_n     = uniqueN(extracted_drydowns$year),
    years_range = paste(min(extracted_drydowns$year), max(extracted_drydowns$year), sep = "-"),
    EVI_range   = sprintf("%.2f - %.2f",
                          min(extracted_drydowns$EVI, na.rm = TRUE),
                          max(extracted_drydowns$EVI, na.rm = TRUE)),
    CWD_max     = round(max(extracted_drydowns$cwd, na.rm = TRUE), 1),
    file_MB     = round(file.info(f)$size / 1024^2, 1)
  )
  
  rm(extracted_drydowns); gc(verbose = FALSE)
}

summary_dt <- rbindlist(summary_list, fill = TRUE)
print(summary_dt)

cat("\n=== TOTALS ===\n")
cat("Total rows across all zones:", format(sum(summary_dt$rows), big.mark = ","), "\n")
cat("Total unique pixels:",         format(sum(summary_dt$pixels), big.mark = ","), "\n")
cat("Total disk used (MB):",        round(sum(summary_dt$file_MB), 1), "\n")