# Vegetation dynamics along drydowns across global climate zones: analysis code

This repository contains the R code used in:

> Terristi, M., Sippel, S., Stocker, B. D., and Bastos, A.: Vegetation dynamics along
> drydowns across global climate zones, EGUsphere [preprint],
> https://doi.org/10.5194/egusphere-2026-5509, 2026.

Archived version: https://doi.org/10.5281/zenodo.XXXXXXX

## Overview

We compute a weekly cumulative water deficit (CWD) from precipitation (MSWEP) and
evapotranspiration (GLEAM). A drydown is a period during which CWD stays above zero.
For each grid cell and year, we keep the drydown with the highest CWD and attach
weekly EVI values to it.

From these data we extract two states for each grid cell:
- the peak state: the 5 % of weeks with the highest EVI, and the CWD at those weeks
- the low state: the 5 % of weeks with the lowest EVI, and the CWD at those weeks

If the low state happens at a higher CWD than the peak state, EVI decreases as the
deficit grows and the grid cell is water-limited (WL). If it happens at a lower CWD,
EVI increases as the deficit grows and the grid cell is energy-limited (EL).

We also compute the sensitivity S, the change in EVI (in %) per mm of CWD, and the
trends of the two states over 2000 to 2023. Results are grouped by Koeppen-Geiger
climate zone.

## Input data

The input data are not included here. They can be downloaded from:

| Dataset | Source |
| --- | --- |
| MODIS MOD13C1 v6.1 EVI | https://doi.org/10.5067/MODIS/MOD13C1.061 |
| MSWEP v2.8 precipitation | https://www.gloh2o.org/mswep/ |
| GLEAM v4.2a evapotranspiration | https://www.gleam.eu |
| Koeppen-Geiger climate classification v2 | https://doi.org/10.6084/m9.figshare.21789074 |
| ESA-CCI / C3S land cover | https://doi.org/10.24381/cds.006f2c9a |
| GLDAS p5 land-sea mask | https://ldas.gsfc.nasa.gov/gldas |

Two files used by the code were prepared from these data beforehand:
- `grouped_koppen_1991_2020.nc`: the Koeppen-Geiger map grouped into the eight classes
  listed in stage 03 below
- `dominant_landcover_no_others_1992_2022_0.1deg.tif`: the main ESA-CCI land-cover class
  of each 0.1 degree grid cell and year, without urban, bare soil, water, and snow and
  ice

## How the code is organised

The scripts are numbered in the order in which they are run. Each script starts with
a short description of what it does, its inputs and outputs, and how to run it. It
also gives the name the script had in our working folder.

| Stage | Script | Methods | What it does |
| --- | --- | --- | --- |
| 01 | `01a_cwd_wettest_month.R` | 2.5 | Finds the wettest month of each grid cell, when CWD is reset every year |
| 01 | `01b_cwd_weekly.R` | 2.5 | Computes weekly CWD, one year at a time |
| 02 | `02_evi_weekly.R` | 2.3 | Computes weekly EVI from MODIS, 2000 to 2023 |
| 03 | `03_cwd_by_climate_zone.R` | 2.4 | Splits the weekly CWD by climate zone |
| 04 | `04a_drydown_events.R` | 2.5, 2.6 | Finds the drydowns and keeps the largest one per grid cell and year |
| 04 | `04b_drydowns_add_evi.R` | 2.6 | Adds weekly EVI to these drydowns |
| 04 | `04c_drydowns_landcover_mask.R` | 2.5, 2.6 | Removes non-vegetated and irrigated grid cells, keeps 2000 to 2023 |
| 04 | `04d_check_masked_sources.R` | | Checks the output of 04c (number of rows and grid cells) |
| 05 | `pooled/05a_pooled_chunks.R`, `05b_pooled_states.R`, `05c_pooled_merge.R` | 2.6 | Peak and low states for each grid cell, all years together |
| 05 | `per_year/05a_per_year_chunks.R`, `05b_per_year_states.R`, `05c_per_year_merge.R` | 2.6, 2.8 | Peak and low states for each grid cell and each year |
| 06 | `06_trends_theil_sen.R` | 2.8 | Trends of the yearly states (Theil-Sen slope, Mann-Kendall test) |
| 07 | figure and table scripts | 2.7, 2.9, 2.10 | See the table below |

Stage 05 is run in three steps because of the size of the data: the grid cells are
split into groups (a), each group is processed as a separate job (b), and the results
are merged (c). The two merged files, `combined_pooled_percentiles_summary_v2.Rdata`
and `combined_per_year_percentiles_summary_v2.Rdata`, are used by stages 06 and 07.

## Figures and tables

| Script | Output file | In the paper |
| --- | --- | --- |
| `FigS1_example_cell_prepare.R`, `FigS1_example_cell_plot.R` | `FigureS1_example_timeseries.png` | Fig. S1 |
| `Fig2_states_by_zone_and_regime.R` | `Figure3_V2_all_split_with_map_noAridHot.png` | Fig. 2 |
| `Fig3_sensitivity.R` | `Figure4_Sensitivity_v3.png` | Fig. 3 |
| `Fig4_evi_trend_direction.R` | `Figure_evi_trend_direction_by_regime.png` | Fig. 4 |
| `Fig5_trend_trajectories_kde.R` | `KDE_grid_zone_regime_v2.png` | Fig. 5 |
| `Fig6_regime_shifts_extreme_years.R` | `FigureS_EL_to_WL_flips_map.png` | Fig. 6 |
| `FigS_drydown_timing.R` | `Figure_drydown_timing.png` | Fig. S2 (timing) |
| `FigS_global_maps.R` | `Map_Delta_CWD.png` and other maps | Fig. S2a-b (maps) |
| `Fig4_evi_trend_direction.R` | `FigureS_EVI_trend_magnitude_4maps.png` | Fig. S3 |
| `Fig4_evi_trend_direction.R` | `Figure_EVI_MK_pvalue_maps.png` | Fig. S4 (significance) |
| `Fig6_regime_shifts_extreme_years.R` | `FigureS_Regime_flip_frequency_timeseries.png` | Fig. S4 (time series) |
| `TablesS_state_statistics.R` | `TableS2_regime_fractions_areaweighted.csv` | Table S1 |
| `TablesS_state_statistics.R` | `TableS1_median_spread_skew.csv` | Table S2 |
| `Fig3_sensitivity.R` | `TableS5a` to `TableS5d` | Table S2a to S2d |
| `TablesS_state_statistics.R` | `TableS3_cliffs_delta_both_regimes.csv`, `TableS3_KW_effectsize.csv` | Table S3a, S3b |
| `TablesS_state_statistics.R` | `TableS4_ELvsWL_within_zone.csv` | Table S4 |
| `Fig4_evi_trend_direction.R` | printed in the R console | Table S6a, S6b |
| `Fig5_trend_trajectories_kde.R` | computed within the script, not written to a file | Table S7, S7a, S7b |
| `Fig6_regime_shifts_extreme_years.R` | `TableS_EL_to_WL_flips_by_zone.csv` | Table S8 |
| `FigS_cwd_trend_direction.R`, `FigS_trend_slope_maps.R`, `FigS_drydown_timing.R` (duration, variability) | additional maps | |

The output file names come from our working folder and do not follow the figure
numbers of the paper. The table above gives the match.

## Choices to be aware of

- **Climate zones.** The code processes eight zones: Tropical, Arid Hot, Arid Cold,
  Temperate Dry, Temperate Humid, Cold Dry, Cold Humid and Polar. Polar is not used
  after stage 04. Arid Hot is removed from all figures and tables because the yearly
  CWD reset at the wettest month does not work well there. The results cover six zones.
- **Years.** Drydowns are detected from 1999 so that a drydown already going on in
  January 2000 is captured from its start. The 1999 data are removed in stage 04c.
- **Year of a drydown.** In stage 04a, a drydown gets the year in which it starts. From
  stage 04b onwards, each week gets its own calendar year.
- **Peak and low states.** Each state uses 5 % of the weeks (at least one week). A grid
  cell needs at least five weeks of data.
- **Sensitivity.** S is not computed when the two states are less than 5 mm of CWD apart.
- **Trends.** Trends are computed separately for EL and WL years, with at least five
  years each. Grid cells with EVI above 1 are removed. Slopes are per year in stage 06
  and per decade in the figures.
- **Fig. 5.** A grid cell is shown only if both its EVI trend and its CWD trend are
  significant (p < 0.05).
- **Statistics.** Medians, spreads and effect sizes are computed per grid cell. Area
  fractions are weighted by cos(latitude).

## Running the code

The code was run with R 4.4 on the MPI-BGC computing cluster. To run it elsewhere,
change the paths at the top of each script. Stages 01, 04 and 05 were run as SLURM
jobs (one per year, zone or group of grid cells); each script explains how to call it.

R packages used: `terra`, `ncdf4`, `raster`, `lubridate`, `data.table`, `dplyr`,
`trend`, `ggplot2`, `patchwork`, `cowplot`, `ggnewscale`, `sf`, `rnaturalearth`,
`scales`, `viridis`, `MASS`.

## Licence

MIT, see `LICENSE`. The input datasets keep the licences of their providers.
