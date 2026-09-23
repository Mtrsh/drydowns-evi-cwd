# =============================================================================
# TRENDS OF THE YEARLY STATES : THEIL-SEN SLOPE AND MANN-KENDALL TEST
# Original file name: Percentiles_trends_per_year_V2.R
# Methods 2.8
# -----------------------------------------------------------------------------
# Uses the yearly states of each grid cell (EVI_E95, CWD_atE95, EVI_E05,
# CWD_atE05, one value per grid cell and year) and computes the trend of each
# series with the Theil-Sen slope and the Mann-Kendall test
# (trend::sens.slope). Slopes are per year.
#
# Trends are computed separately for each regime. The years of a grid cell are
# split by the regime of each year ("after": water-limited, "before":
# energy-limited), and one trend is computed for each group of years. A grid
# cell can therefore have two trends, one per regime. We need at least five
# years per grid cell and regime. Years with p05_position "equal" or NA are not
# used.
#
# Before that, values with EVI_E95 > 1 are removed, because EVI cannot be
# larger than 1.
#
# A trend is significant if p < 0.05, and is labelled Significant Increase,
# Significant Decrease or No Significant Trend. The figure scripts use these
# labels and the slopes.
#
# Input  : combined_per_year_percentiles_summary_v2.Rdata (object all_summaries)
# Output : sens_slope_analysis_complete_v2.Rdata (object all_sens_trends), one
#            row per grid cell, regime and variable
#          sens_slope_<metric>_<regime>_v2.Rdata, the same results in eight
#            files, used by the trend figure scripts
#          sens_slope_summary_stats_v2.Rdata and two CSV summaries by variable,
#            regime and zone
#          analysis_metadata_v2.Rdata, settings and R session information
#
# Run    : Rscript 06_trends_theil_sen.R (one job)
# =============================================================================

## ---- paths (change these lines to run the script elsewhere) ----------------
base_path  <- "/Net/Groups/BGI/scratch/mterristi/PhD/CWD_0_1/New_Climate_Zones_percentiles/Drydowns/"
output_dir <- "/Net/Groups/BGI/scratch/mterristi/PhD/CWD_0_1/New_Climate_Zones_percentiles/Trends/Trends_V2/"

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

library(data.table)
library(trend)

## ---- yearly states ---------------------------------------------------------
v2_file    <- file.path(base_path,
                        "Drydowns_percentiles_per_year_v2/combined_per_year_percentiles_summary_v2.Rdata")

cat("Loading V2 per-year data...\n")
load(v2_file)  # -> all_summaries
combined_df <- copy(all_summaries)
setDT(combined_df)
rm(all_summaries); gc(verbose = FALSE)

cat("Rows loaded:", nrow(combined_df), "\n")
cat("Columns:", paste(names(combined_df), collapse = ", "), "\n")

# Rename columns.
if ("Latitude"  %in% names(combined_df)) setnames(combined_df, "Latitude",  "lat")
if ("Longitude" %in% names(combined_df)) setnames(combined_df, "Longitude", "lon")

cat("Year range:", min(combined_df$Year, na.rm = TRUE), "to",
    max(combined_df$Year, na.rm = TRUE), "\n")

## ---- filters (see header) --------------------------------------------------
cat("\nFiltering...\n")

# EVI cannot be larger than 1.
n_before <- nrow(combined_df)
combined_df <- combined_df[is.na(EVI_E95) | EVI_E95 <= 1]
cat("Removed EVI_E95 > 1:", n_before - nrow(combined_df), "rows\n")

# Keep the two regimes only ("equal" and NA have no regime).
combined_df <- combined_df[p05_position %in% c("after", "before")]
cat("Rows with valid regime:", nrow(combined_df), "\n")

cat("\nDistribution by regime:\n")
print(combined_df[, .N, by = p05_position])

cat("\nGrid cells by climate zone:\n")
print(table(combined_df$climate_zone))

# Number of years per grid cell and regime, printed before the five-year
# rule is applied in the trend function.
grid_cell_counts <- combined_df[, .N, by = .(lat, lon, p05_position)]
cat("\n(lat,lon,regime) groups with < 5 years:",  sum(grid_cell_counts$N < 5),  "\n")
cat("(lat,lon,regime) groups with >= 5 years:", sum(grid_cell_counts$N >= 5), "\n")

## ---- trend per grid cell and regime ----------------------------------------
# For each grid cell, zone and regime, the values are sorted by year and
# passed to sens.slope(), which gives the Theil-Sen slope and the
# Mann-Kendall p-value. Groups with fewer than five years get NA and are
# removed.
perform_sens_slope_analysis_v2 <- function(df, metric_col) {
  cat("\nAnalysing Sen's slope for:", metric_col, "(stratified by regime)\n")
  
  df_metric <- df[!is.na(get(metric_col))]
  
  results <- df_metric[, {
    vals <- get(metric_col)[order(Year)]
    
    if (length(vals) >= 5 && !all(is.na(vals))) {
      res <- tryCatch(
        sens.slope(vals),
        error = function(e) list(estimates = NA_real_, p.value = NA_real_)
      )
      mn <- mean(vals, na.rm = TRUE)
      list(
        sens_slope  = as.numeric(res$estimates),
        sens_pvalue = as.numeric(res$p.value),
        mean_value  = mn,
        cv_value    = ifelse(mn != 0, sd(vals, na.rm = TRUE) / mn, NA_real_),
        n_years     = length(vals),
        metric_name = metric_col
      )
    } else {
      list(
        sens_slope  = NA_real_,
        sens_pvalue = NA_real_,
        mean_value  = NA_real_,
        cv_value    = NA_real_,
        n_years     = length(vals),
        metric_name = metric_col
      )
    }
  }, by = .(lon, lat, climate_zone, regime = p05_position)]
  
  results <- results[!is.na(sens_slope) & !is.na(sens_pvalue)]
  
  results[, significance    := ifelse(sens_pvalue < 0.05,
                                      "Significant (p<0.05)", "Not Significant")]
  results[, trend_direction := fcase(
    sens_slope > 0 & sens_pvalue < 0.05, "Significant Increase",
    sens_slope < 0 & sens_pvalue < 0.05, "Significant Decrease",
    default = "No Significant Trend"
  )]
  
  cat("  Valid results:", nrow(results), "\n")
  for (rg in c("after", "before")) {
    sub <- results[regime == rg]
    cat("  Regime", rg, ":\n")
    cat("    Significant increase:", sum(sub$trend_direction == "Significant Increase"), "\n")
    cat("    Significant decrease:", sum(sub$trend_direction == "Significant Decrease"), "\n")
    cat("    No significant trend:", sum(sub$trend_direction == "No Significant Trend"), "\n")
  }
  
  results
}

## ---- the four variables ----------------------------------------------------
metrics <- c("EVI_E95", "CWD_atE95", "EVI_E05", "CWD_atE05")

cat("\nRunning Sen's slope analysis for all V2 thresholds...\n")
all_sens_trends <- rbindlist(lapply(metrics, function(m) {
  perform_sens_slope_analysis_v2(combined_df, m)
}), fill = TRUE)

cat("\nSen's slope analysis completed.\n")
cat("Total trend analyses:", nrow(all_sens_trends), "\n")

## ---- summaries -------------------------------------------------------------
cat("\n=== SUMMARY BY METRIC × REGIME ===\n")
metric_regime_summary <- all_sens_trends[, .(
  total_cells          = .N,
  significant_increase = sum(trend_direction == "Significant Increase"),
  significant_decrease = sum(trend_direction == "Significant Decrease"),
  no_significant_trend = sum(trend_direction == "No Significant Trend"),
  mean_slope           = round(mean(sens_slope, na.rm = TRUE), 5),
  median_slope         = round(median(sens_slope, na.rm = TRUE), 5)
), by = .(metric_name, regime)][order(metric_name, regime)]
print(metric_regime_summary)

cat("\n=== SUMMARY BY METRIC × REGIME × CLIMATE ZONE ===\n")
zone_summary <- all_sens_trends[, .(
  total_cells          = .N,
  significant_increase = sum(trend_direction == "Significant Increase"),
  significant_decrease = sum(trend_direction == "Significant Decrease"),
  no_significant_trend = sum(trend_direction == "No Significant Trend"),
  pct_significant      = round(sum(sens_pvalue < 0.05) / .N * 100, 1),
  median_slope         = round(median(sens_slope, na.rm = TRUE), 5)
), by = .(metric_name, regime, climate_zone)][order(metric_name, regime, climate_zone)]
print(zone_summary)

## ---- outputs ---------------------------------------------------------------
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

save(all_sens_trends,
     file = file.path(output_dir, "sens_slope_analysis_complete_v2.Rdata"))
cat("\nSaved complete results.\n")

# One file per variable and regime, as read by the figure scripts.
for (m in metrics) {
  for (rg in c("after", "before")) {
    metric_result <- all_sens_trends[metric_name == m & regime == rg]
    save(metric_result,
         file = file.path(output_dir,
                          paste0("sens_slope_", m, "_", rg, "_v2.Rdata")))
  }
  cat("Saved", m, "(both regimes)\n")
}

save(metric_regime_summary, zone_summary,
     file = file.path(output_dir, "sens_slope_summary_stats_v2.Rdata"))

fwrite(metric_regime_summary,
       file.path(output_dir, "sens_slope_summary_by_metric_regime.csv"))
fwrite(zone_summary,
       file.path(output_dir, "sens_slope_summary_by_zone.csv"))

analysis_metadata <- list(
  analysis_date          = Sys.time(),
  total_rows_input       = nrow(combined_df),
  metrics_analyzed       = metrics,
  regimes_analyzed       = c("after", "before"),
  min_years_required     = 5,
  significance_threshold = 0.05,
  climate_zones          = unique(combined_df$climate_zone),
  method                 = "V2 unconditional percentiles per pixel-year, stratified by p05_position",
  notes                  = c(
    "Each pixel may have a trend in both regimes if it had >=5 years in each",
    "EVI_E95 and CWD_atE95 are computed independently of the regime, but separate trends are still informative",
    "EVI_E05 and CWD_atE05 are directly tied to the regime classification"
  ),
  session_info           = sessionInfo()
)
save(analysis_metadata,
     file = file.path(output_dir, "analysis_metadata_v2.Rdata"))

cat("\nAll results saved to:", output_dir, "\n")
cat("=== ANALYSIS COMPLETE ===\n")