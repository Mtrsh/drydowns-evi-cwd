# =============================================================================
# SUPPLEMENTARY TABLES ON THE POOLED STATES (TABLES S1 TO S4)
# Original file name: Stats_Part3_1.R
# Methods 2.10
# -----------------------------------------------------------------------------
# Statistics for the variables shown in Fig. 2, by climate zone and regime.
# The script reads the pooled states and writes four CSV tables:
#
#   TableS1_median_spread_skew.csv
#       median, quartiles, IQR and two skewness measures (Bowley and moment
#       skewness) of each variable, by zone and regime
#   TableS2_regime_fractions_areaweighted.csv
#       share of EL and WL grid cells per zone, weighted by cos(latitude);
#       the number of grid cells and the shares without weights are also given
#   TableS3_KW_effectsize.csv, TableS3_cliffs_delta_both_regimes.csv
#       Kruskal-Wallis test between the six zones for each regime (epsilon
#       squared = H / (n - 1)), and Cliff's delta between pairs of zones, for
#       dEVI and dCWD
#   TableS4_ELvsWL_within_zone.csv
#       EL compared with WL within each zone: medians and Cliff's delta
#
# The numbers in the file names come from our working folder and differ from
# the supplement. See the README for the match.
#
# Area fractions are weighted by cos(latitude). All other statistics are
# computed per grid cell, without weights, like the boxplots of Fig. 2. With
# tens to hundreds of thousands of grid cells, all tests are significant, so
# we use effect sizes rather than p-values.
#
# Same filters and variables as in Fig2_states_by_zone_and_regime.R.
#
# Input  : combined_pooled_percentiles_summary_v2.Rdata
# Output : TableS1_median_spread_skew.csv
#          TableS2_regime_fractions_areaweighted.csv
#          TableS3_KW_effectsize.csv
#          TableS3_cliffs_delta_both_regimes.csv
#          TableS4_ELvsWL_within_zone.csv
#
# Run    : Rscript TablesS_state_statistics.R
# =============================================================================

## ---- environment -----------------------------------------------------------

# Our R package folder on the cluster, used only if it exists.
custom_libs <- c(
  "/Net/Groups/BGI/scratch/mterristi/R_lodomeria/4.4",
  "/Net/Groups/BSI/work_scratch/quincy/model/software/r_packages/r_4.4.x"
)
.libPaths(c(custom_libs[dir.exists(custom_libs)], .libPaths()))
library(data.table)

input_path  <- "/Net/Groups/BGI/scratch/mterristi/PhD/CWD_0_1/New_Climate_Zones_percentiles/Drydowns/Drydowns_percentiles_pooled_v2"
output_path <- "/Net/Groups/BGI/scratch/mterristi/PhD/CWD_0_1/New_Climate_Zones_percentiles/Drydowns/Figures"
dir.create(output_path, showWarnings = FALSE, recursive = TRUE)

## ---- pooled states, with the filters of the header -------------------------
cat("Loading V2 pooled data...\n")
load(file.path(input_path, "combined_pooled_percentiles_summary_v2.Rdata"))
setDT(all_summaries)
all_summaries <- all_summaries[!is.na(EVI_E05) & !is.na(CWD_atE05) &
                                 p05_position %in% c("before", "after") &
                                 climate_zone != "Arid Hot"]
all_summaries[, Delta_EVI := 100 * fifelse(
  p05_position == "after",
  (EVI_E05 - EVI_E95) / EVI_E95,
  (EVI_E95 - EVI_E05) / EVI_E95)]
all_summaries[, Delta_CWD := abs(CWD_atE05 - CWD_atE95)]
all_summaries[, w := cos(Latitude * pi / 180)]                        # area weight
zone_levels <- c("Tropical", "Arid Cold", "Temperate Dry",
                 "Temperate Humid", "Cold Dry", "Cold Humid")
all_summaries[, climate_zone := factor(climate_zone, levels = zone_levels)]
all_summaries[, regime := factor(fifelse(p05_position == "before", "EL", "WL"),
                                 levels = c("EL", "WL"))]
cat(sprintf("Grid cells: %d\n", nrow(all_summaries)))
cat(sprintf("Regime split — COUNT:  EL %.1f%%  WL %.1f%%\n",
            100*mean(all_summaries$regime=="EL"), 100*mean(all_summaries$regime=="WL")))
cat(sprintf("Regime split — AREA-WEIGHTED:  EL %.1f%%  WL %.1f%%\n",
            100*all_summaries[regime=="EL", sum(w)]/all_summaries[, sum(w)],
            100*all_summaries[regime=="WL", sum(w)]/all_summaries[, sum(w)]))

## ---- median, spread and skewness, per grid cell ----------------------------
bowley <- function(x){
  q <- quantile(x, c(.25, .5, .75), na.rm = TRUE, names = FALSE)
  if (diff(q[c(1, 3)]) == 0) NA_real_ else (q[3] + q[1] - 2*q[2]) / (q[3] - q[1])
}
mskew <- function(x){
  x <- x[is.finite(x)]; n <- length(x)
  if (n < 3) return(NA_real_)
  m <- mean(x); s <- sd(x); if (s == 0) NA_real_ else (sum((x - m)^3) / n) / s^3
}
vars  <- c(EVI95="EVI_E95", EVI05="EVI_E05", dEVI="Delta_EVI",
           CWD95="CWD_atE95", CWD05="CWD_atE05", dCWD="Delta_CWD")
units <- c(EVI95="-", EVI05="-", dEVI="%", CWD95="mm", CWD05="mm", dCWD="mm")
supp <- rbindlist(lapply(names(vars), function(nm){
  v <- vars[[nm]]
  all_summaries[, .(
    variable = nm, unit = units[[nm]], n = .N,
    median = round(median(get(v), na.rm = TRUE), 3),
    Q25    = round(quantile(get(v), .25, na.rm = TRUE, names = FALSE), 3),
    Q75    = round(quantile(get(v), .75, na.rm = TRUE, names = FALSE), 3),
    IQR    = round(IQR(get(v), na.rm = TRUE), 3),
    skew_q = round(bowley(get(v)), 3),
    skew_m = round(mskew(get(v)),  3)
  ), by = .(climate_zone, regime)]
}))
setcolorder(supp, c("variable","unit","climate_zone","regime","n",
                    "median","Q25","Q75","IQR","skew_q","skew_m"))
setorder(supp, variable, climate_zone, regime)
print(supp)
fwrite(supp, file.path(output_path, "TableS1_median_spread_skew.csv"))

## ---- share of each regime per zone, weighted by area -----------------------
# pct_area is the value cited in the paper; n and pct_count are given for
# information.
frac <- all_summaries[, .(n = .N, wsum = sum(w)), by = .(climate_zone, regime)]
frac[, pct_count := round(100 * n    / sum(n),    1), by = climate_zone]
frac[, pct_area  := round(100 * wsum / sum(wsum), 1), by = climate_zone]
frac[, wsum := NULL]
setorder(frac, climate_zone, regime)
print(frac)
fwrite(frac, file.path(output_path, "TableS2_regime_fractions_areaweighted.csv"))

## ---- effect sizes between zones, per grid cell -----------------------------
# Kruskal-Wallis test and Cliff's delta between zones, for dEVI and dCWD,
# for each regime.
kw_eps2 <- function(x, g){
  ok <- is.finite(x) & !is.na(g); x <- x[ok]; g <- droplevels(as.factor(g[ok]))
  kt <- kruskal.test(x, g); H <- unname(kt$statistic); k <- nlevels(g); n <- length(x)
  data.table(H = round(H,1), k = k, n = n, p = signif(kt$p.value,3),
             epsilon2 = round(H/(n-1), 4))
}
cliffs_delta <- function(a, b){
  a <- a[is.finite(a)]; b <- b[is.finite(b)]
  n1 <- as.numeric(length(a)); n2 <- as.numeric(length(b))
  r <- rank(c(a, b)); U1 <- sum(r[seq_len(length(a))]) - n1*(n1+1)/2
  round(2*U1/(n1*n2) - 1, 3)
}
key <- c(dEVI = "Delta_EVI", dCWD = "Delta_CWD")
kw <- rbindlist(lapply(names(key), function(nm)
  all_summaries[, kw_eps2(get(key[[nm]]), climate_zone), by = regime][, variable := nm][]))
setcolorder(kw, c("variable", "regime")); setorder(kw, variable, regime); print(kw)
fwrite(kw, file.path(output_path, "TableS3_KW_effectsize.csv"))

# Cliff's delta between pairs of zones.
pair_delta <- function(var, reg){
  d <- all_summaries[regime == reg]; zs <- levels(droplevels(d$climate_zone))
  out <- CJ(zone_a = zs, zone_b = zs)[zone_a < zone_b]
  out[, delta := mapply(function(x, y) cliffs_delta(d[climate_zone==x][[var]],
                                                    d[climate_zone==y][[var]]),
                        zone_a, zone_b)]
  out[, `:=`(variable = var, regime = reg)][]
}
cliff <- rbindlist(list(
  pair_delta("Delta_EVI", "WL"), pair_delta("Delta_CWD", "WL"),
  pair_delta("Delta_EVI", "EL"), pair_delta("Delta_CWD", "EL")))
setorder(cliff, variable, regime, delta); print(cliff)
fwrite(cliff, file.path(output_path, "TableS3_cliffs_delta_both_regimes.csv"))

## ---- EL compared with WL in each zone --------------------------------------
# delta_EL_vs_WL > 0 means that EL values are usually higher than WL values.
vars_wz <- c(EVI95="EVI_E95", EVI05="EVI_E05", CWD95="CWD_atE95",
             CWD05="CWD_atE05", dCWD="Delta_CWD")   # dEVI omitted: sign is definitional
elwl <- rbindlist(lapply(names(vars_wz), function(nm){
  v <- vars_wz[[nm]]
  all_summaries[, .(
    variable = nm,
    med_EL = round(median(get(v)[regime=="EL"], na.rm=TRUE), 2),
    med_WL = round(median(get(v)[regime=="WL"], na.rm=TRUE), 2),
    delta_EL_vs_WL = cliffs_delta(get(v)[regime=="EL"], get(v)[regime=="WL"])
  ), by = climate_zone]
}))
setcolorder(elwl, c("variable", "climate_zone")); setorder(elwl, variable, climate_zone)
print(elwl)
fwrite(elwl, file.path(output_path, "TableS4_ELvsWL_within_zone.csv"))

cat("\nDone. Wrote: TableS1_median_spread_skew.csv, TableS2_regime_fractions_areaweighted.csv,\n",
    "TableS3_KW_effectsize.csv, TableS3_cliffs_delta_both_regimes.csv, TableS4_ELvsWL_within_zone.csv\n")