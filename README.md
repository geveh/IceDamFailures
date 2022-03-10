# Changes in ice-dam failures with ongoing glacier decay

## Overview

**This repository contains seven scripts to identify trends in the volume (*V*<sub>0</sub>), peak discharge (*Q*<sub>p</sub>), timing (day of year *doy*), and source elevation (*Z*) of ice-dam failures on regional and local (i.e. lake-) level. Furthermore, we investigate the consequences of melting glacier dams on the magnitudes of GLOFs.**

- [01_lake_area_volume.R](#01_lake_area_volumer)
- [02_quantile_regression.R](#02_quantile_regressionr)
- [03_local_trends_in_V0_and_Qp.R](#03_local_trends_in_v0_and_qpr)
- [04_trends_in_doy.R](#04_trends_in_doyr)
- [05_glacier_volumes_and_ice_loss.R](#05_glacier_volumes_and_ice_lossr)
- [06_trends_in_Z.R](06_trends_in_zr)
- [07_magnitudes_vs_elev_change.R](07_magnitudes_vs_elev_changer)

The codes are written in the statistical programming language **R** (https://www.r-project.org/), Version 4.0.0, and called within
the Graphical User Interface **RStudio** (https://rstudio.com) under a Microsoft Windows 10 operating system. 
Please install both **R and RStudio** on your machine to successfully run the codes and produce figures and R data objects.

The R codes depend on a number of packages, listed at the beginning of all scripts. Please install those packages before running the scripts. 
The comments within the scripts provide further details on model dependencies and usage of functions. 

Each script will call one or more input data object(s), which will be made available ***via Zenodo soon***.  
We also use freely available digital elevation models (DEMs) and glaciological data (glacier outlines, estimates of ice thickness and mass loss). Please download the data from the web sources provided in the scripts.  
Please put all input files into the same folder. The scripts can be executed one after the other, with the user generating output that is used as input for the next script.
The scripts (and parts thereof) can also be run independent of each other using the input files from Zenodo.
Each script will produce output in form of a figure (displayed in the associate manuscript and supplementary figures) or R-objects.

## Scripts

### 01_lake_area_volume.R

**Script to predict glacier lake volume *V* from glacier lake area *A*. We use a Bayesian piece-wise
regression model that objectively learns the location of a breakpoint in the empirical relationship between lake volume and lake area.**

*Mandatory input data*: 
- "Global_GLOF_database_2021_12_08.ods" (table with all reported GLOFs according to the Global GLOF database)
- "lake_area_volume_compiliaton.xlsx" (table with previously reported data pairs of volume and area of glacier lakes)

*Output*: 
- "va_breakpoint.pdf" / "va_breakpoint.png" (Plot of the piece-wise *V-A* regression model)
- "va_model.RDS" (R-object of the fitted volume-area-relationship using the package mcp)
---

### 02_quantile_regression.R

**Script to fit quantile regression models (50th and 90th percentile) of peak discharges *Q*<sub>p</sub> and volumes *V*<sub>0</sub> versus time 
from ice-dam failures in six mountain ranges.**

*Mandatory input data*: 
- "Global_GLOF_database_2021_12_08.ods" (table with all reported GLOFs according to the Global GLOF database)

*Output*: 
- "all_glofs_tibble.RDS" (R-object with a preprocessed table of all reported GLOFs)
- "qp_regional_quantreg_fits.RDS" (R-object with a quantile regression model of V<sub>0</sub> and Q<sub>p</sub> versus time for the 50th and 90th percentile for each region)
- "qp_regional_quantreg_posteriors.RDS" (R-object of the posterior predictive distribution of V<sub>0</sub> and Q<sub>p</sub> for each region and both percentiles for each year in the period 1900-2021)
- "qp_regional_quantreg_trends.RDS" (R-object of posterior trends in V<sub>0</sub> and Q<sub>p</sub> for each region and both percentiles)
- "qp_regional_medline.RDS" (R-object of median temporal change in V<sub>0</sub> and Q<sub>p</sub> for each region and both percentiles)
---

### 03_local_trends_in_V0_and_Qp.R

**Script to estimate trends in V<sub>0</sub> and Q<sub>p</sub> for individual ice-dammed lakes with repeat outburst.** 

*Mandatory input data*: 
- "all_glofs_tibble.RDS" (R-object with a preprocessed table of all reported GLOFs)
- "va_model.RDS" (R-object of the fitted volume-area-relationship. Will be used to make new predictions of lake volume from mapped lake areas)

*Output*: 

GLOF volume V<sub>0</sub>
- "all_glofs_V0_tibble.RDS" (Table of lakes with repeat GLOFs and reported V<sub>0</sub>)
- "local_V0_model.RDS" (Hierarchical regression model of V<sub>0</sub> versus time for individual glacier lakes)
- "regression_V0_per_lake.pdf" / "regression_V0_per_lake.png" (Plot of the posterior trends in V<sub>0</sub> for each lake)
- "post_trend_V0_per_lake.pdf" / "post_trend_V0_per_lake.png" (Plot of the posterior regression slope of V<sub>0</sub> versus time)
- "all_glofs_V0_tibble.RDS" (Table of lakes with repeat GLOFs and reported V<sub>0</sub>)

GLOF peak discharge Q<sub>p</sub>
- "all_glofs_qp_tibble.RDS" (Table of lakes with repeat GLOFs and reported Q<sub>p</sub>)
- "local_qp_model.RDS" (Hierarchical regression model of Q<sub>0</sub> versus time for individual glacier lakes)
- "regression_Qp_per_lake.pdf" / "regression_Qp_per_lake.png" (Plot of the posterior trends in Q<sub>p</sub> for each lake)
- "post_trend_Qp_per_lake.pdf" / "post_trend_Qp_per_lake.png" (Plot of the posterior regression slope of Q<sub>p</sub> versus time)
- "all_glofs_qp_tibble.RDS" (Table of lakes with repeat GLOFs and reported Q<sub>p</sub>)
---

### 04_trends_in_doy.R

**Script to estimate trends in the annual timing (i.e. day in a given year) of ice-dam failures on regional and local scale.**

*Mandatory input data*: 
- "all_glofs_tibble.RDS" (R-object with a preprocessed table of all reported GLOFs)

*Output*: 
- "doy_trends_per_region.RDS" (R-object with regression models of *doy* versus time for all dated GLOFs in the six regions)
- "doy_change.pdf" / "doy_change.png" (Plot of the temporal trends in *doy* for each region, including the posterior differences in *doy* between 2021 and 1900)
- "doy_trends_per_glacier.RDS"  (R-object with regression models of *doy* versus time for lakes with repeat GLOFs)
- "regression_doy_per_lake.pdf" / "regression_doy_per_lake.png" (Plot of local changes in *doy* versus time)
- "post_trend_doy_per_lake.pdf" / "post_trend_doy_per_lake.png" (Plot of local  posterior differences in *doy* for each lake)
---

### 05_glacier_volumes_and_ice_loss.R

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

### 06_trends_in_Z.R

**Script to estimate regional trends in the source elevation (*Z*) of ice-dammed failures**

*Mandatory input data*: 
- Digital Elevation models from ALOS World 3D - 30m (AW3D30)
- Folder "GDL_database" (Regional inventories of ice-dammed lakes in ESRI shapefile format)
- "Regional_glacier_and_melt_volumes.rds" (R-object containing the total volume of glacier volume and volume change between 2000 and 2019 in 100-m elevation bins)

*Output*: 
- "Z_trends_per_region.RDS" (R-object with a hierarchical regression models of *Z* versus time for dated GLOFs in the six regions between 1900 and 2021)
- "elev_trend.pdf" / "elev_trend.png" (Plot of the change in GLOF source elevation for six regions between 1900 and 2021, including the posterior regression slope)
- "Lake_GLOF_elevation.pdf" / "Lake_GLOF_elevation.png" (Plot of the elevation distribution of historic burst ice-dammed lakes and present-day ice-dammed lakes for six regions between 1900 and 2021)

---

### 07_magnitudes_vs_elev_change.R

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

## Input data

Please visit repository on Zenodo (***link will appear here soon***) to obtain the input files.


## References

Georg Veh, Natalie Lützow, Jenny Tamm, Romain Hugonnet, Marten Geertsema, John J Clague, and Oliver Korup: *Smaller and earlier outbursts from ice-dammed lakes with ongoing glacier decay* (submitted).

## See also

http://glofs.geoecology.uni-potsdam.de

## Contact

**Georg Veh**  
Postdoctoral researcher in the working group on natural hazards  
Institute of Environmental Sciences and Geography  
University of Potsdam  
georg.veh@uni-potsdam.de  
https://www.uni-potsdam.de/de/umwelt/forschung/ag-naturgefahren.html
