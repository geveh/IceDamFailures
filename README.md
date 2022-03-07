# Ice-Dam Failures

## Overview

This repository contains **XXX scripts** to analyse trends in the GLOF volume V<sub>0</sub>, Q<sub>p</sub>, timing, and source elevation of ice-dam failures on regional and local (i.e. lake-level) scale. 

- 01_lake_area_volume.R
- 02_quantile_regression.R
- 03_local_trends_in_V0_and_Qp.R
- 04_trends_in_doy.R
- 05_XXX.R
- 06_XXX.R
- 07_magnitudes_vs_elev_change.R

The codes are written in the statistical programming language **R** (https://www.r-project.org/), Version 4.0.0, and called within
the Graphical User Interface **RStudio** (https://rstudio.com) under a Microsoft Windows 10 operating system. 
Please install both **R and RStudio** on your machine to successfully run the codes and produce figures and R data objects.

The R codes also depend on a number of packages, listed at the beginning of both scripts. Please install those packages before running the codes. 
The comments within the scripts provide further details on model dependencies and usage of functions. 

Each script will call one or more input data object(s), which will be made available via Zenodo soon. 
Please put all Input files into the same folder.
The Zenodo folder also contains the results of each script, which you can use to compare with the results of your own analysis.

## Scripts

### 01_lake_area_volume.R

Script to predict glacier lake volume from glacier lake area. We use a Bayesian piece-wise
regression model that objectively learns the location of breakpoint in the empirical relationship between lake volume and lake area.

*Mandatory input data*: 
- "Global_GLOF_database_2021_12_08.ods" (table with all reported GLOFs according to the Global GLOF database)
- "lake_area_volume_compiliaton.xlsx" (table with previously reported data pairs of volume and area of glacier lakes)

*Output*: 
- "va_breakpoint.pdf" / "va_breakpoint.png" (Plot of the piece-wise regression model)
- "va_model.RDS" (R-object of the fitted volume-area-relationship using the package mcp)


### 02_quantile_regression.R

Script to fit quantile regression models (50th and 90th percentile) of peak discharges and volumes versus time 
from ice-dam failures in six mountain ranges. 

*Mandatory input data*: 
- "Global_GLOF_database_2021_12_08.ods" (table with all reported GLOFs according to the Global GLOF database)

*Output*: 
- "all_glofs_tibble.RDS" (R-object with a preprocessed table of all reported GLOFs)
- "qp_regional_quantreg_fits.RDS" (R-object with a quantile regression model of V<sub>0</sub> and Q<sub>p</sub> versus time for the 50th and 90th percentile for each region)
- "qp_regional_quantreg_posteriors.RDS" (R-object of the posterior predictive distribution of V<sub>0</sub> and Q<sub>p</sub> for each region and both percentiles for each year in the period 1900-2021)
- "qp_regional_quantreg_trends.RDS" (R-object of posterior trends in V<sub>0</sub> and Q<sub>p</sub> for each region and both percentiles)
- "qp_regional_medline.RDS" (R-object of median temporal change in V<sub>0</sub> and Q<sub>p</sub> for each region and both percentiles)


### 03_local_trends_in_V0_and_Qp.R

Script to estimate trends in V<sub>0</sub> and Q<sub>p</sub> for individual ice-dammed lakes with repeat outburst. 

*Mandatory input data*: 
- "all_glofs_tibble.RDS" (R-object with a preprocessed table of all reported GLOFs)
- "va_model.RDS" (R-object of the fitted volume-area-relationship. Will be used to make new predictions of lake volume from mapped lake area)

*Output*: 

GLOF volume  V<sub>0</sub>
- "all_glofs_V0_tibble.RDS" (Table of lakes with repeat GLOFs and reported V<sub>0</sub>)
- "local_V0_model.RDS" (Hierarchical regression model of V<sub>0</sub> versus time for individual glacier lakes)
- "regression_V0_per_lake.pdf" / "regression_V0_per_lake.png" (Plot of the posterior trends in V<sub>0</sub> for each lake)
- "post_trend_V0_per_lake.pdf" / "post_trend_V0_per_lake.png" (Plot of the posterior regression slope of V<sub>0</sub> versus time)
- "all_glofs_qp_tibble.RDS" (Table of lakes with repeat GLOFs and reported Q<sub>p</sub>)

GLOF peak discharge Q<sub>p</sub>
- "all_glofs_qp_tibble.RDS" (Table of lakes with repeat GLOFs and reported Q<sub>p</sub>)
- "local_qp_model.RDS" (Hierarchical regression model of Q<sub>0</sub> versus time for individual glacier lakes)
- "regression_Qp_per_lake.pdf" / "regression_Qp_per_lake.png" (Plot of the posterior trends in Q<sub>p</sub> for each lake)
- "post_trend_Qp_per_lake.pdf" / "post_trend_Qp_per_lake.png" (Plot of the posterior regression slope of Q<sub>p</sub> versus time)
- "all_glofs_qp_tibble.RDS" (Table of lakes with repeat GLOFs and reported Q<sub>p</sub>)


### 04_trends_in_doy.R

Script to estimate trends in the annual timing (i.e. day in a given year) of ice-dam failures on regional and local scale.

*Mandatory input data*: 
*Output*: 


### 05_XXX.R


*Description TBA*
*Mandatory input data*: TBA
*Output*: TBA


### 06_XXX.R

*Description TBA*
*Mandatory input data*: TBA
*Output*: TBA


### 07_magnitudes_vs_elev_change.R

*Mandatory input data*: 
*Output*: 



## Input data

### Global_GLOF_database_2021_06_09.ods

Open-Office spreadsheet as of 09 June 2021 with seven sheets named after the regions, for which we obtained historical GLOF occurrences. 
Each sheet has 32 columns containing the attributes that we were able to collect for each GLOF. Empty cells mean 'No Data'. 
The first row is the column name, followed by two rows with further description of the content and the data structure.
The content of the columns 'Major_RGI_Region', 'Mountain_range_Region', 'Glacier',	'RGI_Glacier_Id', and	'RGI_Glacier_Area' is from the
Randolph Glacier Inventory, V6.0 (https://www.glims.org/RGI/rgi60_dl.html).


### Region_extents.zip

Extents of study regions in a WGS 84 / World Mercator projection


### regional_glof_stats.rds

R-Data object (a list with 8 entries) containing regional annual statistics of GLOF occurrences, temperatures, and research activity.
Description of the column names:
- 'year': Year;
- 'freq': Total number of reported GLOFs per year, including GLOFs from volcanic eruptions;
- 'moraine': Number of moraine-dam failures per year;
- 'ice': Number of ice-dam failures per year;
- 'other': Number of GLOFs from other (bedrock, water pockets, supraglacial) or unknown sources;
- 'volc': Number of GLOFs from subglacial lakes beneath ice-covered volcanoes;
- 'mb_meas': Annual number of glacier surveys measuring in-situ mass balances from the WGMS database;
- 'front_meas': Annual number of glacier surveys measuring in-situ front variations;
- 'dch_meas': Annual number of glacier surveys measuring geodetic mass balances (includes also remote sensing studies);
- 'all_meas': Annual sum of mb_meas, front_meas, and dch_meas;
- 'mb_and_front': Annual sum of mb_meas and front_meas;
- 'region': Name of the study region;
- 'year_scale': Standardised years (zero mean and unit standard deviation);
- 'temp_mean': Mean annual air temperature extracted from the CRU TS 4.05 dataset from all lakes that produced at least one GLOF in a given region;
- 'temp_q25': 25th percentile of annual air temperatures in a given region;
- 'temp_q75': 75th percentile of annual air temperatures in a given region;
- 'pre_sum': total amount of precipitation in a given region.



## References

Georg Veh, Natalie Lützow, Jenny Tamm, Romain Hugonnet, Marten Geertsema, John J Clague, and Oliver Korup: *Smaller and earlier outbursts from ice-dammed lakes with ongoing glacier decay* (submitted).

## See also

http://glofs.geoecology.uni-potsdam.de

## Contact

Georg Veh

Working group on natural hazards

University of Potsdam

georg.veh@uni-potsdam.de

https://www.uni-potsdam.de/de/umwelt/forschung/ag-naturgefahren.html
