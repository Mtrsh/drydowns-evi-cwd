# =============================================================================
# FIGURE 3 AND TABLES S5a TO S5d : EVI-CWD SENSITIVITY
# Original file name: Stats_Part3_2.R
# Methods 2.7 (sensitivity) and 2.10 (statistics)
# -----------------------------------------------------------------------------
# Sensitivity of EVI to the water deficit, for each grid cell:
#   S = dEVI / dCWD, in % per mm
# dEVI is the relative EVI change with a sign (positive in EL, negative in WL)
# and dCWD = |CWD_EVI05 - CWD_EVI95| in mm, as in
# Fig2_states_by_zone_and_regime.R. S is not computed when dCWD is 5 mm or
# less (CWD_FLOOR), because the two states are then at almost the same CWD
# and the ratio has no meaning.
#
# Figure 3: a) global map of S, with the colour scale limited to +/- 10 % per
# mm; b) boxplots of S by zone and regime.
#
# Tables (the numbers in the file names come from our working folder, see the
# README for the match with the supplement):
#   TableS5a  median, quartiles, IQR, 5th and 95th percentiles and two
#             skewness measures of S, by zone and regime
#   TableS5b  Kruskal-Wallis test between zones for each regime (epsilon
#             squared)
#   TableS5c, TableS5d  Cliff's delta between pairs of zones, for WL and EL
#
# As in TablesS_state_statistics.R, statistics are computed per grid cell,
# without weights. The area fractions and medians cited in Sect. 3.2 are
# weighted by cos(latitude) and printed in the console.
#
# We keep grid cells with both states, EVI95 > 0 and a regime. Arid Hot is
# removed.
#
# The output file name comes from our working folder. This is Fig. 3 of the
# paper.
#
# Input  : combined_pooled_percentiles_summary_v2.Rdata
# Output : Figure4_Sensitivity_v3.png                    Fig. 3
#          TableS5a_Sensitivity_by_zone_regime.csv
#          TableS5b_Sensitivity_KW.csv
#          TableS5c_Sensitivity_cliffs_delta_WL.csv
#          TableS5d_Sensitivity_cliffs_delta_EL.csv
#
# Run    : Rscript Fig3_sensitivity.R
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
library(data.table); library(ggplot2); library(sf); library(rnaturalearth); library(patchwork)
input_path  <- "/Net/Groups/BGI/scratch/mterristi/PhD/CWD_0_1/New_Climate_Zones_percentiles/Drydowns/Drydowns_percentiles_pooled_v2"
output_path <- "/Net/Groups/BGI/scratch/mterristi/PhD/CWD_0_1/New_Climate_Zones_percentiles/Drydowns/Figures"
dir.create(output_path, showWarnings = FALSE, recursive = TRUE)
## ---- pooled states, with the filters of the header -------------------------
load(file.path(input_path, "combined_pooled_percentiles_summary_v2.Rdata"))
setDT(all_summaries)
all_summaries <- all_summaries[!is.na(EVI_E05) & !is.na(CWD_atE05) &
                                 !is.na(EVI_E95) & EVI_E95 > 0 &
                                 p05_position %in% c("before","after") &
                                 climate_zone != "Arid Hot"]
all_summaries[, Delta_EVI := 100 * fifelse(p05_position=="after",
                                           (EVI_E05-EVI_E95)/EVI_E95, (EVI_E95-EVI_E05)/EVI_E95)]  # relative %, signed
all_summaries[, Delta_CWD := abs(CWD_atE05 - CWD_atE95)]                              # mm, >= 0
all_summaries[, regime := factor(fifelse(p05_position=="before","EL","WL"), levels=c("EL","WL"))]
all_summaries[, w := cos(Latitude * pi/180)]                                         # area weight
zone_levels <- c("Tropical","Arid Cold","Temperate Dry","Temperate Humid","Cold Dry","Cold Humid")
all_summaries[, climate_zone := factor(climate_zone, levels=zone_levels)]
## ---- sensitivity S (see header) --------------------------------------------
CWD_FLOOR <- 5  # mm; S undefined where the two states sit at nearly the same CWD
all_summaries[, S := fifelse(Delta_CWD > CWD_FLOOR, Delta_EVI / Delta_CWD, NA_real_)]
n_all <- nrow(all_summaries); n_S <- sum(!is.na(all_summaries$S))
cat(sprintf("S computed for %d / %d cells (%.1f%%); %d dropped (ΔCWD ≤ %d mm)\n",
            n_S, n_all, 100*n_S/n_all, n_all-n_S, CWD_FLOOR))
cat(sprintf("SIGN CHECK — EL median %.3f (+),  WL median %.3f (−)\n",
            median(all_summaries[regime=="EL"]$S, na.rm=TRUE),
            median(all_summaries[regime=="WL"]$S, na.rm=TRUE)))
## ---- numbers cited in Sect. 3.2, weighted by area, printed only ------------
wmed <- function(x, w){ ok<-is.finite(x); x<-x[ok]; w<-w[ok]
o<-order(x); x<-x[o]; w<-w[o]; x[which(cumsum(w)/sum(w) >= 0.5)[1]] }
Wtot <- all_summaries[, sum(w)]
cat("\n[TEXT] Area-weighted regime split:\n")
print(all_summaries[, .(pct = round(100*sum(w)/Wtot,1)), by = regime])
cat(sprintf("[TEXT] |ΔEVI| area-weighted median: %.1f %%\n",
            wmed(abs(all_summaries$Delta_EVI), all_summaries$w)))
cat(sprintf("[TEXT] ΔCWD  area-weighted median: %.1f mm (overall)\n",
            wmed(all_summaries$Delta_CWD, all_summaries$w)))
cat("[TEXT] ΔCWD area-weighted median by regime (mm):\n")
print(all_summaries[, .(cwd_med_aw = round(wmed(Delta_CWD, w),1)), by = regime])
## ---- statistics of S per zone and regime (TableS5a) ------------------------
bowley <- function(x){q<-quantile(x,c(.25,.5,.75),na.rm=TRUE,names=FALSE)
if(diff(q[c(1,3)])==0) NA_real_ else (q[3]+q[1]-2*q[2])/(q[3]-q[1])}
mskew  <- function(x){x<-x[is.finite(x)];n<-length(x);if(n<3)return(NA_real_)
m<-mean(x);s<-sd(x);if(s==0)NA_real_ else (sum((x-m)^3)/n)/s^3}
S_stats <- all_summaries[!is.na(S), .(
  n=.N, median=round(median(S),3),
  Q25=round(quantile(S,.25,names=FALSE),3), Q75=round(quantile(S,.75,names=FALSE),3),
  IQR=round(IQR(S),3),
  P05=round(quantile(S,.05,names=FALSE),3), P95=round(quantile(S,.95,names=FALSE),3),
  skew_q=round(bowley(S),3), skew_m=round(mskew(S),3)
), by=.(climate_zone, regime)][order(regime, median)]
print(S_stats)
fwrite(S_stats, file.path(output_path, "TableS5a_Sensitivity_by_zone_regime.csv"))
## ---- effect sizes (TableS5b to TableS5d) -----------------------------------
# Kruskal-Wallis test between zones, for each regime (S5b).
kw_eps2 <- function(x,g){ok<-is.finite(x)&!is.na(g);x<-x[ok];g<-droplevels(as.factor(g[ok]))
kt<-kruskal.test(x,g);H<-unname(kt$statistic);k<-nlevels(g);n<-length(x)
data.table(H=round(H,1),k=k,n=n,p=signif(kt$p.value,3),epsilon2=round(H/(n-1),4))}
kwS <- all_summaries[!is.na(S), kw_eps2(S, climate_zone), by=regime]
print(kwS); fwrite(kwS, file.path(output_path, "TableS5b_Sensitivity_KW.csv"))
# Cliff's delta between pairs of zones, WL (S5c) and EL (S5d).
cliffs_delta <- function(a,b){ a<-a[is.finite(a)]; b<-b[is.finite(b)]
n1<-as.numeric(length(a)); n2<-as.numeric(length(b))
r<-rank(c(a,b)); U1<-sum(r[seq_len(length(a))]) - n1*(n1+1)/2
round(2*U1/(n1*n2)-1,3) }
cliff_table <- function(reg){
  d  <- all_summaries[regime==reg & !is.na(S)]
  zs <- levels(droplevels(d$climate_zone))
  ct <- CJ(zone_a=zs, zone_b=zs)[zone_a<zone_b]
  ct[, delta := mapply(function(x,y) cliffs_delta(d[climate_zone==x]$S,
                                                  d[climate_zone==y]$S), zone_a, zone_b)]
  ct[order(delta)]
}
cliffWL <- cliff_table("WL"); cliffEL <- cliff_table("EL")
cat("\n== Cliff's δ — WATER-LIMITED ==\n");  print(cliffWL)
cat("\n== Cliff's δ — ENERGY-LIMITED ==\n"); print(cliffEL)
fwrite(cliffWL, file.path(output_path,"TableS5c_Sensitivity_cliffs_delta_WL.csv"))
fwrite(cliffEL, file.path(output_path,"TableS5d_Sensitivity_cliffs_delta_EL.csv"))
## ---- Figure 3: a) map of S, b) boxplots by zone ----------------------------
BASE <- 26; col_EL<-"#2E7D32"; col_WL<-"#8B4513"
robin <- "+proj=robin +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
pts <- st_transform(st_as_sf(all_summaries[!is.na(S)],
                             coords=c("Longitude","Latitude"), crs=4326), crs=robin)
xy  <- st_coordinates(pts)
dmap <- data.table(X=xy[,1], Y=xy[,2], S=all_summaries[!is.na(S)]$S)
world <- st_crop(st_transform(ne_countries(scale="medium", returnclass="sf"), robin),
                 xmin=-17e6, xmax=17e6, ymin=-6e6, ymax=8.5e6)
# ---- graticule: lines and labels from the same lon/lat vectors --------------
grat_lat <- c(-60, -30, 0, 30, 60)                     # parallels (latitude)
grat_lon <- c(-120, -60, 0, 60, 120)                   # meridians (longitude)
LAT_LO   <- -60; LAT_HI <- 84                          # plotted vertical frame
yb <- st_coordinates(st_transform(
  st_as_sf(data.frame(lon = 0, lat = c(LAT_LO, LAT_HI)), coords = c("lon","lat"), crs = 4326),
  robin))[, "Y"]
grat <- st_transform(st_graticule(lon = grat_lon, lat = grat_lat, crs = st_crs(4326)), robin)
grat <- st_crop(grat, xmin = -17e6, xmax = 17e6, ymin = yb[1], ymax = yb[2])
lat_xy <- st_coordinates(st_transform(
  st_as_sf(data.frame(lon = -180, lat = grat_lat), coords = c("lon","lat"), crs = 4326), robin))
lat_lab <- data.table(x = lat_xy[, "X"], y = lat_xy[, "Y"],
                      label = paste0(abs(grat_lat), "°", ifelse(grat_lat > 0, "N", ifelse(grat_lat < 0, "S", ""))))
lon_xy <- st_coordinates(st_transform(
  st_as_sf(data.frame(lon = grat_lon, lat = LAT_LO), coords = c("lon","lat"), crs = 4326), robin))
lon_lab <- data.table(x = lon_xy[, "X"], y = lon_xy[, "Y"],
                      label = paste0(abs(grat_lon), "°", ifelse(grat_lon > 0, "E", ifelse(grat_lon < 0, "W", ""))))
lim <- 10   # display cap, +/- 10 % per mm; values beyond are clamped, not dropped
p_map <- ggplot() +
  geom_tile(data=dmap, aes(X, Y, fill=S), width=11000, height=11000) +
  geom_sf(data=grat, color="grey85", linewidth=0.2) +            # our own grid lines
  geom_sf(data=world, fill=NA, color="gray30", linewidth=.2) +
  geom_text(data=lat_lab, aes(x=x, y=y, label=label), inherit.aes=FALSE,
            hjust=1, nudge_x=-4e5, size=10, colour="grey45") +    # latitude labels (left)
  geom_text(data=lon_lab, aes(x=x, y=y, label=label), inherit.aes=FALSE,
            vjust=1, nudge_y=-4e5, size=10, colour="grey45") +    # longitude labels (bottom)
  scale_fill_gradient2(low=col_WL, mid="grey95", high=col_EL, midpoint=0,
                       limits=c(-lim,lim), oob=scales::squish, name="Sensitivity (% per mm)",
                       breaks=pretty(c(-lim,lim),5)) +
  coord_sf(crs=robin, ylim=yb, expand=FALSE, clip="off") + labs(tag="a)") +
  theme_void(base_size=BASE) +
  theme(legend.position="bottom", legend.key.width=unit(3,"cm"),
        legend.key.height=unit(.5,"cm"), legend.title=element_text(size=24),
        legend.text=element_text(size=22), plot.tag=element_text(size=30),
        plot.margin=margin(6, 8, 28, 36))              # room for bottom/left coord labels
p_box <- ggplot(all_summaries[!is.na(S)], aes(climate_zone, S, fill=regime)) +
  geom_hline(yintercept=0, linetype="dashed", linewidth=.4) +
  geom_boxplot(outlier.shape=NA, position=position_dodge(.75), width=.65) +
  scale_fill_manual(values=c(EL=col_EL, WL=col_WL), name="Regime") +
  coord_cartesian(ylim=c(-10,10)) + labs(x=NULL, y="Sensitivity (% per mm)", tag="b)") +
  theme_bw(base_size=BASE) +
  theme(axis.text.x=element_text(angle=30,hjust=1,size=22),
        axis.text.y=element_text(size=22), panel.grid.major.x=element_blank(),
        plot.tag=element_text(size=30))
fig4 <- p_map / p_box + plot_layout(heights=c(3, 1.4))
ggsave(file.path(output_path,"Figure4_Sensitivity_v3.png"), fig4,
       width=20, height=16, dpi=300, bg="white", limitsize=FALSE)
cat("\nDone. Wrote S5a, S5b, S5c (WL), S5d (EL), and Figure4_Sensitivity_v3.png\n")
