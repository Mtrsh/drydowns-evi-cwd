# =============================================================================
# FIGURE 5 : JOINT EVI AND CWD TRENDS AT THE PEAK AND LOW STATES
# Original file name: Percentiles_KDE_only.R
# Methods 2.9
# -----------------------------------------------------------------------------
# For each climate zone (columns) and regime (rows), each grid cell is placed
# according to its CWD trend (x axis, mm per decade) and its EVI trend (y axis,
# per decade), separately for the peak state (EVI95 and CWD_EVI95) and the low
# state (EVI05 and CWD_EVI05).
#
# For each state we compute a 2D kernel density (MASS::kde2d), shown as
# contour lines, and its mode (the point with the highest density). An arrow
# goes from the mode of the peak state to the mode of the low state. The lines
# at zero split the plot into four parts, named after the sign of the two
# trends. A small inset shows the direction of the arrow when the two modes
# are very close.
#
# The slopes come from 06_trends_theil_sen.R and are multiplied by 10 to get
# values per decade. A grid cell is used for a state only if both its EVI trend
# and its CWD trend are significant (p < 0.05). Arid Hot is removed.
#
# Input  : sens_slope_<metric>_<regime>_v2.Rdata, for the four variables
#          combined_pooled_percentiles_summary_v2.Rdata (climate zones)
# Output : KDE_grid_zone_regime_v2.png                    Fig. 5
#
# Run    : Rscript Fig5_trend_trajectories_kde.R
# =============================================================================

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
library(ggplot2)
library(patchwork)
library(MASS)
library(cowplot)
library(ggnewscale)
## ---- paths and parameters (change these lines to run the script elsewhere) ---
trends_dir <- "/Net/Groups/BGI/scratch/mterristi/PhD/CWD_0_1/New_Climate_Zones_percentiles/Trends/Trends_V2"
v2_dir     <- "/Net/Groups/BGI/scratch/mterristi/PhD/CWD_0_1/New_Climate_Zones_percentiles/Drydowns/Drydowns_percentiles_pooled_v2"
out_dir    <- file.path(trends_dir, "Figures")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
p_threshold <- 0.05
BASE <- 25                   # base font size, as in the other figure scripts
peak_fill     <- "#1b7837"   # green (p95, high-greenness state)
critical_fill <- "#8c510a"   # brown (p05, low-greenness state)
XLIM <- c(-200, 200)
YLIM <- c(-0.30, 0.30)
zone_order <- c("Tropical", "Arid Cold",
                "Temperate Dry", "Temperate Humid",
                "Cold Dry", "Cold Humid")
## ---- climate zone of each grid cell, without Arid Hot ----------------------
cat("Loading V2 lookup...\n")
load(file.path(v2_dir, "combined_pooled_percentiles_summary_v2.Rdata"))
setDT(all_summaries)
zone_lookup <- unique(all_summaries[climate_zone != "Arid Hot",
                                    .(lat = Latitude, lon = Longitude, climate_zone)])
## ---- trends, significant only, per decade ----------------------------------
# regime_lab is the name used in the files ("before" or "after"):
#   before: CWD_atE05 < CWD_atE95, energy-limited (EL)
#   after : CWD_atE05 >= CWD_atE95, water-limited (WL)
load_trend <- function(metric, regime) {
  f <- file.path(trends_dir,
                 paste0("sens_slope_", metric, "_", regime, "_v2.Rdata"))
  env <- new.env(); load(f, envir = env)
  d <- as.data.table(get(ls(env)[1], envir = env))
  d[!is.na(sens_pvalue) & sens_pvalue < p_threshold,
    .(lat, lon, sens_slope)]
}
build_scatter <- function(regime_lab) {
  evi_p <- load_trend("EVI_E95",   regime_lab)[, .(lat, lon, EVI_slope = sens_slope)]
  cwd_p <- load_trend("CWD_atE95", regime_lab)[, .(lat, lon, CWD_slope = sens_slope)]
  peak  <- merge(evi_p, cwd_p, by = c("lat", "lon"))[, threshold := "Threshold p95"]
  
  evi_c <- load_trend("EVI_E05",   regime_lab)[, .(lat, lon, EVI_slope = sens_slope)]
  cwd_c <- load_trend("CWD_atE05", regime_lab)[, .(lat, lon, CWD_slope = sens_slope)]
  crit  <- merge(evi_c, cwd_c, by = c("lat", "lon"))[, threshold := "Threshold p05"]
  
  d <- rbind(peak, crit)
  d <- merge(d, zone_lookup, by = c("lat", "lon"))
  d <- d[!is.na(EVI_slope) & !is.na(CWD_slope)]
  d[, `:=`(EVI_slope = EVI_slope * 10, CWD_slope = CWD_slope * 10)]
  d[, threshold := factor(threshold, levels = c("Threshold p95", "Threshold p05"))]
  d[, regime := c(before = "EL", after = "WL")[[regime_lab]]]   # WL/EL, not before/after
  d
}
## ---- mode of a 2D kernel density -------------------------------------------
kde_mode <- function(x, y, n = 100) {
  if (length(x) < 10) return(c(NA_real_, NA_real_))
  d <- kde2d(x, y, n = n)
  i <- which(d$z == max(d$z), arr.ind = TRUE)[1, ]
  c(d$x[i[1]], d$y[i[2]])
}
# Labels of the four parts: G greening, B browning, D drying, W wetting.
quad_labels <- data.table(
  x     = c( 180, -180,  180, -180),
  y     = c( 0.255, 0.255, -0.255, -0.255),
  label = c("GD", "GW", "BD", "BW")
)

## ---- one panel: zone by regime ---------------------------------------------
make_panel <- function(zone, regime_lab, scatter_dt,
                       show_xlab = TRUE, show_ylab = TRUE,
                       show_xticks = TRUE, show_yticks = TRUE) {
  dt <- scatter_dt[climate_zone == zone & regime == regime_lab]
  if (nrow(dt) == 0) return(plot_spacer())
  
  dt_peak <- dt[threshold == "Threshold p95"]
  dt_crit <- dt[threshold == "Threshold p05"]
  
  # Modes of the two states
  m_peak <- if (nrow(dt_peak) >= 10)
    kde_mode(dt_peak$CWD_slope, dt_peak$EVI_slope) else c(NA, NA)
  m_crit <- if (nrow(dt_crit) >= 10)
    kde_mode(dt_crit$CWD_slope, dt_crit$EVI_slope) else c(NA, NA)
  
  mode_dt <- rbind(
    data.table(CWD_mode = m_peak[1], EVI_mode = m_peak[2],
               threshold = factor("Threshold p95",
                                  levels = c("Threshold p95", "Threshold p05"))),
    data.table(CWD_mode = m_crit[1], EVI_mode = m_crit[2],
               threshold = factor("Threshold p05",
                                  levels = c("Threshold p95", "Threshold p05")))
  )
  
  arrow_dt <- data.table(x = m_peak[1], y = m_peak[2],
                         xend = m_crit[1], yend = m_crit[2])
  
  # --- inset: direction of the arrow from the peak to the low state ----------
  # Drawn with the same axes as the main plot (CWD trend, EVI trend). It shows
  # the direction clearly even when the two modes are almost at the same place.
  inset_box <- list(xmin = -45, xmax = 45, ymin = -0.295, ymax = -0.225)
  inset_layer <- NULL
  if (!is.na(m_peak[1]) && !is.na(m_crit[1])) {
    dx <- m_crit[1] - m_peak[1]
    dy <- m_crit[2] - m_peak[2]
    norm_axes <- sqrt((dx / diff(XLIM))^2 + (dy / diff(YLIM))^2)
    if (norm_axes > 0) {
      ux <- (dx / diff(XLIM)) / norm_axes
      uy <- (dy / diff(YLIM)) / norm_axes
      cx <- mean(c(inset_box$xmin, inset_box$xmax))
      cy <- mean(c(inset_box$ymin, inset_box$ymax))
      half_w <- (inset_box$xmax - inset_box$xmin) * 0.35
      half_h <- (inset_box$ymax - inset_box$ymin) * 0.35
      px <- cx - ux * half_w
      py <- cy - uy * half_h
      qx <- cx + ux * half_w
      qy <- cy + uy * half_h
      inset_pts <- data.frame(
        x  = c(px, qx),
        y  = c(py, qy),
        threshold = factor(c("Threshold p95", "Threshold p05"),
                           levels = c("Threshold p95", "Threshold p05"))
      )
      inset_arrow <- data.frame(x = px, y = py, xend = qx, yend = qy)
      inset_layer <- list(
        annotate("rect",
                 xmin = inset_box$xmin, xmax = inset_box$xmax,
                 ymin = inset_box$ymin, ymax = inset_box$ymax,
                 fill = "white", color = "black", linewidth = 0.5),
        geom_segment(data = inset_arrow,
                     aes(x = x, y = y, xend = xend, yend = yend),
                     inherit.aes = FALSE,
                     color = "grey25", linewidth = 0.7,
                     arrow = arrow(length = unit(0.20, "cm"),
                                   type = "closed", ends = "last")),
        geom_point(data = inset_pts,
                   aes(x = x, y = y, fill = threshold),
                   inherit.aes = FALSE, size = 2.5,
                   shape = 21, color = "grey25", stroke = 0.6)
      )
    }
  }
  
  p <- ggplot() +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.4) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.4) +
    
    # grid cells as light grey points
    geom_point(data = dt_peak,
               aes(x = CWD_slope, y = EVI_slope),
               color = peak_fill, size = 0.55, alpha = 0.40, shape = 16) +
    geom_point(data = dt_crit,
               aes(x = CWD_slope, y = EVI_slope),
               color = critical_fill, size = 0.55, alpha = 0.40, shape = 16) +
    
    # Density contours, with a colour gradient for each state. ggnewscale
    # allows two separate colour scales in the same plot:
    #   - peak state (p95): teal, light outside and dark inside
    #   - low state (p05): brown, light outside and dark inside
    # The filled areas are light grey.
    stat_density_2d(data = dt_peak,
                    aes(x = CWD_slope, y = EVI_slope,
                        fill = after_stat(level)),
                    geom = "polygon", bins = 4,
                    alpha = 0.12, color = NA, show.legend = FALSE) +
    scale_fill_gradient(low = "grey92", high = "grey70", guide = "none") +
    ggnewscale::new_scale_fill() +
    
    # contour lines of the peak state, teal
    geom_density_2d(data = dt_peak,
                    aes(x = CWD_slope, y = EVI_slope,
                        colour = after_stat(level)),
                    bins = 4, linewidth = 0.85,
                    show.legend = FALSE) +
    scale_colour_gradientn(
      colours = c("#C7EAE5", "#5AB4AC", "#1F7C73", "#003C30"),
      guide = "none"
    ) +
    ggnewscale::new_scale_colour() +
    
    # contour lines of the low state, brown
    geom_density_2d(data = dt_crit,
                    aes(x = CWD_slope, y = EVI_slope,
                        colour = after_stat(level)),
                    bins = 4, linewidth = 0.85,
                    show.legend = FALSE) +
    scale_colour_gradientn(
      colours = c("#F5D499", "#D88A21", "#995B0F", "#543005"),
      guide = "none"
    ) +
    
    # arrow from the peak to the low state, dark grey
    geom_segment(data = arrow_dt[!is.na(x) & !is.na(xend)],
                 aes(x = x, y = y, xend = xend, yend = yend),
                 color = "grey25", linewidth = 0.85,
                 arrow = arrow(length = unit(0.25, "cm"),
                               type = "closed", ends = "last")) +
    
    # modes, coloured points with a thin dark grey border
    ggnewscale::new_scale_colour() +
    geom_point(data = mode_dt[!is.na(CWD_mode)],
               aes(x = CWD_mode, y = EVI_mode, fill = threshold),
               shape = 21, color = "grey25", stroke = 0.7, size = 3.8) +
    scale_fill_manual(values = c("Threshold p95" = peak_fill,
                                 "Threshold p05" = critical_fill),
                      name = NULL) +
    
    # inset: direction of the arrow
    inset_layer +
    geom_label(data = quad_labels,
               aes(x = x, y = y, label = label),
               size = 8, fontface = "bold", color = "black",
               label.size = 0, label.padding = unit(0.15, "lines"),
               fill = alpha("white", 0.7)) +
    coord_cartesian(xlim = XLIM, ylim = YLIM) +
    labs(x = NULL, y = NULL) +
    theme_minimal(base_size = BASE) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "grey94"),
      panel.border     = element_rect(color = "black", fill = NA, linewidth = 0.7),
      panel.background = element_rect(fill = "white", color = NA),
      axis.line        = element_blank(),
      axis.text.x      = if (show_xticks) element_text(size = 22) else element_blank(),
      axis.text.y      = if (show_yticks) element_text(size = 22) else element_blank(),
      axis.ticks.x     = if (show_xticks) element_line(color = "black") else element_blank(),
      axis.ticks.y     = if (show_yticks) element_line(color = "black") else element_blank(),
      axis.title       = element_text(size = 25),
      plot.margin      = margin(2, 4, 2, 4),
      legend.position  = "none"
    )
  
  # zone code in the top left corner (left column only)
  p
}
## ---- all twelve panels -----------------------------------------------------
cat("\nBuilding scatter datasets...\n")
scatter_before <- build_scatter("before")   # -> regime "EL"
scatter_after  <- build_scatter("after")    # -> regime "WL"
scatter_all    <- rbind(scatter_before, scatter_after)
scatter_all[, climate_zone := factor(climate_zone, levels = zone_order)]
# The 12 panels: one row per regime, one column per zone.
# Row 1: energy-limited (EL), CWD_atE05 < CWD_atE95
# Row 2: water-limited (WL), CWD_atE05 >= CWD_atE95
# Columns follow zone_order.
panels <- list()
for (regime_lab in c("EL", "WL")) {
  is_bottom_row <- (regime_lab == "WL")
  for (i in seq_along(zone_order)) {
    zone <- zone_order[i]
    is_left_col <- (i == 1)
    panels[[paste0(regime_lab, "_", zone)]] <- make_panel(
      zone, regime_lab, scatter_all,
      show_xlab   = is_bottom_row,
      show_ylab   = is_left_col,
      show_xticks = is_bottom_row,
      show_yticks = is_left_col
    )
  }
}
# Row labels are added with cowplot::draw_label after plot_grid.
## ---- put the panels together, two rows and six columns ---------------------
# Row 1: energy-limited (EL, "before"). Row 2: water-limited (WL, "after").
# The panel list holds row 1 (6 zones) then row 2 (6 zones), so byrow = TRUE.
grid_obj <- plot_grid(
  plotlist = panels,
  ncol = length(zone_order), byrow = TRUE
)
# Zone names above the columns
col_label_grobs <- lapply(zone_order, function(z) {
  ggdraw() + draw_label(z, fontface = "bold", size = 25, hjust = 0.5)
})
col_titles <- plot_grid(plotlist = col_label_grobs,
                        ncol = length(zone_order))
# Regime names on the left of the rows, rotated by 90 degrees
row_label_grobs <- list(
  ggdraw() + draw_label("Energy-limited",
                        fontface = "bold", size = 24, angle = 90, vjust = 0.5),
  ggdraw() + draw_label("Water-limited",
                        fontface = "bold", size = 24, angle = 90, vjust = 0.5)
)
row_label_col <- plot_grid(plotlist = row_label_grobs, ncol = 1, byrow = TRUE)
# One y-axis label for the whole figure, left of the regime names and centred
# on the two rows.
y_axis_label <- ggdraw() +
  draw_label(expression("EVI trend (decade"^{-1}*")"),
             size = 25, angle = 90, vjust = 0.5, hjust = 0.5)
# Put together: regime names, y-axis label, panels
left_and_panels <- plot_grid(
  row_label_col, y_axis_label, grid_obj,
  ncol = 3, rel_widths = c(0.05, 0.03, 1)
)
# Add the zone names on top (with empty space above the row labels)
title_and_grid <- plot_grid(
  plot_grid(
    ggdraw(),       # empty corner above row labels
    ggdraw(),       # empty corner above y-axis label
    col_titles,
    ncol = 3, rel_widths = c(0.05, 0.03, 1)
  ),
  left_and_panels,
  ncol = 1, rel_heights = c(0.07, 1)
)
# One x-axis label for the whole figure, below the panels and centred on
# them.
x_axis_label <- plot_grid(
  ggdraw(),  # empty corner under row labels
  ggdraw(),  # empty corner under y-axis label
  ggdraw() + draw_label(expression("CWD trend (mm decade"^{-1}*")"),
                        size = 25, hjust = 0.5),
  ncol = 3, rel_widths = c(0.05, 0.03, 1)
)
# Legend at the bottom, kept narrow so it stays compact over the full width
# of the figure.
legend_plot <- ggplot() +
  annotate("point", x = 18.0, y = 0, colour = peak_fill,    size = 9) +
  annotate("text",  x = 18.4, y = 0, label = "Threshold p95",
           hjust = 0, size = 8) +
  annotate("point", x = 21.8, y = 0, colour = critical_fill, size = 9) +
  annotate("text",  x = 22.2, y = 0, label = "Threshold p05",
           hjust = 0, size = 8) +
  coord_cartesian(xlim = c(0, 40), ylim = c(-1, 1), clip = "off") +
  theme_void()
final <- plot_grid(
  title_and_grid,
  x_axis_label,
  legend_plot,
  ncol = 1,
  rel_heights = c(1, 0.06, 0.10)
)
## ---- save ------------------------------------------------------------------
out_png <- file.path(out_dir, "KDE_grid_zone_regime_v2.png")
ggsave(out_png, final, width = 30, height = 13, dpi = 350, bg = "white")

cat("\nSaved:\n  ", out_png, "\n")
cat("\n=== DONE ===\n")