# Code for *"Less extreme and earlier outbursts of ice-dammed lakes since 1900"*

## Overview

**This repository contains six scripts to estimate trends in the volume (*V*<sub>0</sub>), peak discharge (*Q*<sub>p</sub>), timing (day of year *doy*), and source elevation (*Z*) of ice-dam failures on regional and local (i.e. lake-) level. In addition, we investigate the consequences of melting glacier dams on the magnitudes of GLOFs.**

- [01_preprocessing.R](#01_lake_area_volumer)
- [02_trends_in_Qp_and_V0.R](#02_trends_in_qp_and_v0r)
- [03_trends_in_doy.R](#03_trends_in_doyr)
- [04_glacier_volumes_and_ice_loss.R](#04_glacier_volumes_and_ice_lossr)
- [05_trends_in_Z.R](#05_trends_in_zr)
- [06_magnitudes_vs_elev_change.R](#06_magnitudes_vs_elev_changer)
- [summary_stats_veh_revision.py](#summary_stats_veh_revisionpy)

The codes are written in the statistical programming language **R** (https://www.r-project.org/), Version 4.2.0, and called within
the Graphical User Interface **RStudio** (https://rstudio.com) under a Microsoft Windows 10 operating system. 
Please install both **R and RStudio** on your machine to successfully run the codes and produce figures and R data objects.

The R codes depend on a number of packages, listed at the beginning of all scripts. Please install those packages before running the scripts. 
The comments within the scripts provide further details on model dependencies and usage of functions. 

Each script will call one or more input data object(s), which are available via ***Zenodo***.  
We also use freely available digital elevation models (DEMs), glaciological data (glacier outlines, estimates of ice thickness and mass loss), and lake outlines. Please download the data from the web sources provided in the scripts.  
Please put all input files into the same folder, and change the folder used in the script to your folder structure. The scripts can be executed one after the other, with the user generating output that is used as input for the next script.
The scripts (and parts thereof) can also be run independent of each other using the input files (in most cases *.RDS* files) from Zenodo.
Each script will produce output in form of a figure (displayed in the associate manuscript and Extended Data figures) or R-objects.

## Scripts

### 01_preprocessing.R

**Script to preprocess a raw OpenOffice table of reported glacier lake outburst floods.**

*Mandatory input data*: 
- "Global_GLOF_database_2022_05_30.ods" (table with all reported GLOFs. Compiliation as of May 30, 2022)
- CRU TS V4.05 temperature data


*Main outputs*: 
- "all_glofs_tibble.RDS" (R-object of all reported GLOFs in the global GLOF database)
- "all_glofs_V0_tibble.RDS" (R-object of all GLOFs that have reported values of *V*<sub>0</sub>)
- "all_glofs_qp_tibble.RDS" (R-object of all GLOFs that have reported values of *Q*<sub>p</sub>)
- "glof_reporting.pdf" (Multi-panel histogram of reported values of all reported GLOFs, reported values of Qp, and reported values of V0 from ice-dammed lakes) 
- "temp_doy_histogram.pdf" (histogram that both shows the number of reported GLOFs and the mean air temperature in a given month)

---

### 02_trends_in_Qp_and_V0.R

**Script to fit quantile regression models (50th and 90th percentile) of peak discharges *Q*<sub>p</sub> and volumes *V*<sub>0</sub> versus time 
from ice-dam failures in six mountain ranges.**

*Mandatory input data*: 
- "all_glofs_V0_tibble.RDS" (R-object of all GLOFs that have reported values of *V*<sub>0</sub>)
- "all_glofs_qp_tibble.RDS" (R-object of all GLOFs that have reported values of *Q*<sub>p</sub>)

*Main outputs*: 
- "qp_models.RDS" (R-object with regional quantile regression models of *Q*<sub>p</sub> versus time for the 50th and 90th for 4 time periods)
- "V0_models.RDS" (R-object with regional quantile regression models of *V*<sub>0</sub> versus time for the 50th and 90th for 4 time periods)
- "fig2.pdf" (PDF figure containing the regional posterior slopes for *Q*<sub>p</sub> and *V*<sub>0</sub> for two different time periods)
- "all_pooled_mods.pdf" (PDF figure containing the pooled trendes of *Q*<sub>p</sub> and *V*<sub>0</sub> for two different time periods)
- "qp_model_median_local.RDS"  (R-object with local quantile regression models of median *Q*<sub>p</sub> versus time)
- "Qp_local.pdf" (PDF figure showing temporal trends of median *Q*<sub>p</sub> for individual glacier lakes) 
- "V0_model_median_local.RDS"  (R-object with regional quantile regression models of median *V*<sub>0</sub> versus time)
- "V0_local.pdf" (PDF figure showing temporal trends of median *V*<sub>0</sub> for individual glacier lakes) 
- "local_posterior_trends.pdf" (PDF figure showing posterior distributions of the trends in local *Q*<sub>p</sub> and V*<sub>0</sub>)

---

### 03_trends_in_doy.R

**Script to estimate trends in the annual timing (*doy*, i.e. day in a given year) of ice-dam failures on regional and local scale.**

*Mandatory input data*: 
- "all_glofs_tibble.RDS" (R-object with a preprocessed table of all reported GLOFs)

*Output*: 
- "doy_trends_per_region.RDS" (R-object with regression models of *doy* versus time for all dated GLOFs in the six regions)
- "doy_change.pdf" (Plot of the temporal trends in *doy* for each region, including the posterior differences in *doy* between 2021 and 1900)
- "doy_trends_per_glacier.RDS"  (R-object with regression models of *doy* versus time for lakes with repeat GLOFs)
- "doy_local.pdf" (Plot of local changes in *doy* versus time)
- "post_trend_doy_per_lake.pdf"  (Plot of local  posterior differences in *doy* for each lake)

---

### 04_glacier_volumes_and_ice_loss.R

**Script to obtain the total volumes of glaciers and their volume loss between 2000 and 2019 in 100-m elevation bins.**

*Mandatory input data (Data sources from external repositories are provided in the script)*: 
- Folder "Region_extents" (Contains the ESRI shapefile *Extent_pol.shp* to display the extent of the study regions)
- Glacier outlines from the Randolph Glacier Inventory (RGI)
- Glacier surface DEMs from Farinotti et al. (2019)
- Glacier volume DEMs from Farinotti et al. (2019)
- Glacier elevation change data from Hugonnet et al. (2021)

*Output*: 
- "Regional_glacier_and_melt_volumes.rds" (R-object containing the total volume of glacier volume and volume change between 2000 and 2019 in 100-m elevation bins)

---

### 05_trends_in_Z.R

**Script to estimate regional trends in the source elevation (*Z*) of ice-dammed failures.**

*Mandatory input data*: 
- Digital Elevation models from ALOS World 3D - 30m (AW3D30)
- Files from the "GDL_database" (We created a merged lake inventory in ESRI shapefile format from regional lake databases. This lake database is available upon request)   
- "Regional_glacier_and_melt_volumes.rds" (R-object containing the total volume of glacier volume and volume change between 2000 and 2019 in 100-m elevation bins)

*Major outputs*: 
- "gdl_database_centroid.RDS" (R-object of glacier lake centroids in the six study regions)
- "glofs_ice_with_z.RDS" (R-object of first reported GLOF from a given lake and its elevation)
- "Z_trends_per_region.RDS" (R-object with a hierarchical regression models of *Z* versus time for dated GLOFs in the six regions between 1900 and 2021)
- "elev_trend.pdf" (Plot of the change in GLOF source elevation for six regions between 1900 and 2021, including the posterior regression slope)
- "Lake_GLOF_elevation.pdf" / "Lake_GLOF_elevation.png" (Plot of the elevation distribution of historic burst ice-dammed lakes and present-day ice-dammed lakes for six regions between 1900 and 2021)

---

### 06_magnitudes_vs_elev_change.R

**Script to estimate local trends of  V<sub>0</sub> and  Q<sub>p</sub> with elevation change of the glacier dam.**

*Mandatory input data*: 

- Folder "dh_pergla_cut" (Tables of cumulative elevation change (in m) for glaciers with repeat GLOFs between 2000 and 2019)
- "all_glofs_tibble.RDS" (R-object with a preprocessed table of all reported GLOFs)
- "all_glofs_V0_tibble.RDS" (Table of lakes with repeat GLOFs and reported V<sub>0</sub>)
- "all_glofs_qp_tibble.RDS" (Table of lakes with repeat GLOFs and reported Q<sub>p</sub>)
- Folder "Region_extents" (Contains the ESRI shapefile *Extent_pol.shp* to display the extent of the study regions)

*Output*: 

- "local_Qp_vs_dhdt_model.RDS" (R-Object containing a hierarchical model of local changes in Q<sub>p</sub> versus glacier elevation change)
- "local_V0_vs_dhdt_model.RDS" (R-Object containing a hierarchical model of local changes in V<sub>0</sub> versus glacier elevation change)
- "map_and_trends.pdf" (Map of lakes with repeat GLOFs between 2000 and 2019; local trends of V<sub>0</sub> and Q<sub>p</sub> with cumulative changes in glacier dam elevation)
- "dam_thinning_rats.shp" (ESRI shapefile showing mean annual elevation change of glacier dams with repeat outbursts between 2000 and 2019)
- "elev_change_per_glacier.pdf" / "elev_change_per_glacier.png" (Plot of cumulative elevation change for each glacier that produced repeated GLOFs between 2000 and 2019)

---

### summary_stats_veh_revision.py

** Script by Romain Hugonnet to obtain elevation changes from glacier dams. Please contact R. Hugonnet, if you have further questions.**


## Input data

Please visit repository on Zenodo to obtain the input files.


## References

Georg Veh, Natalie Lützow, Jenny Tamm, Lisa V. Luna, Romain Hugonnet, Kristin Vogel, Marten Geertsema, John J Clague, and Oliver Korup: *Less extreme and earlier outbursts from ice-dammed lakes since 1900*.

## See also

http://glofs.geoecology.uni-potsdam.de

## Contact

**Georg Veh**  
Postdoctoral researcher in the working group on natural hazards  
Institute of Environmental Sciences and Geography  
University of Potsdam  
georg.veh@uni-potsdam.de  
https://www.uni-potsdam.de/de/umwelt/forschung/ag-naturgefahren.html
