# Ice-Dam Failures

This repository contains XXX scripts to analyse trends in the magnitude, timing, and source elevation of ice-dam failures on regional and local (i.e. lake-level) scale. 

- 01_lake_area_volume.R
- 02_quantile_regression.R
- 03_local_trends_in_V0_and_Qp.R
- 04_trends_in_doy.R
- XXX.R
- XXX.R
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

Script to predict glacier lake volume from glacier lake area, aka the volume-area relationship

- GLOF counts by dam type;
- glaciological research activity;
- temperature; and
- precipitation

from the original Open-Office spreadsheet 'Global_GLOF_database_2021_06_09.ods', available at this page (see a detailed description in the section on [Input data](#global_glof_database_2021_06_09ods)).
The script produces the R-Data object *regional_glof_stats.rds*, which is already available on this page.

### assessing_bias.R

Main script to 

- find change points in time series of reported GLOFs, air temperatures, and glacier surveys;
- estimate trends in GLOF reporting for each study region and dam type;
- predict annual GLOF counts from air temperatures and the number of glacier surveys;
- to hind- and forecast the number of GLOFs before and after the global break in GLOF reporting.
- reproduce all figures in the manuscript.




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
