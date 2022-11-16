################################################################################
#####   Obtain the total volumes of glaciers and their volume loss      ########
#####       between 2000 and 2019 in 100-m elevation bins               ########
#####                                                                   ########
#####                         by Georg Veh                              ########
#####                         12 Oct, 2021                              ########
#####               comments added on 09 March, 2022                    ########   
#####             more comments added on 15 Nov, 2022                   ########
################################################################################

# This workflow uses
# - Digitized glacier outlines as reported by the RGI Consortium (2017): 
#   Randolph Glacier Inventory - A Dataset of Global Glacier Outlines: 
#   Version 6.0: Technical Report, Global Land Ice Measurements from Space, 
#   Colorado, USA. Digital Media. DOI: https://doi.org/10.7265/N5-RGI-60
# - Gridded data on glacier thickness as reported by Farinotti et al. (2019): 
#   A consensus estimate for the ice  thickness distribution of all glaciers 
#   on Earth. Nature Geoscience, 12(3), 168-173.
# - Gridded data on glacier elevation change between 2000 and 2019 as reported 
#   by Hugonnet et al. (2021): Accelerated global glacier mass loss in the early 
#   twenty-first century. Nature, 592(7856), 726-731.
#
# You will be prompted to download the files from these publications in the 
# course of this script. Make sure you have some tens of gigabytes of free 
# local hard drive storage on your machine. Furthermore, note that these 
# scripts were executed on a Windows Server environment with 52 core and 
# 256 GB RAM. If necessary, adapt the number of cores to reduce RAM consumption, 
# although this might increase computation time.

# Load packages. Use install.packages("packagename") if one or more of these 
# packages are not yet available on your machine.

library(s2)
library(sf)
library(tidyverse)
library(parallel)
library(stringr)
library(pbapply)
library(foreach)
library(raster)
library(exactextractr)

# Load the shapefile containing the extent of the study region.

sr <-  read_sf("D:/nrc_user/veh/LW_F/Shapefiles/Region_extents/Extent_pol.shp") %>% 
  filter(Region != "Other") %>%
  st_make_valid()

# Reproject this shapefile to WGS84 (a geographic coordinate system).

sr.wgs <- st_transform(sr, 4326)

# Dissolve boundaries of study region and reproject.

sr.uni <- st_union(sr) %>%
  st_transform(4326) %>%
  st_make_valid()

# List the names of the outlines of all glaciers in the Randolph Glacier Inventory (RGI) 6.0.
# To this end, download the entire RGI from https://www.glims.org/RGI/rgi60_files/00_rgi60.zip
# Unzip the folder and change the path to your local folder.

# List the files by changing the folder to the folder structure on your computer.

gl.list <- list.files(path = "D:/nrc_user/veh/LW_F/Shapefiles/00_rgi60", 
                      recursive = T, 
                      include.dirs = T, 
                      full.names = T) %>%
  as_tibble() %>% 
  filter(str_detect(value, ".shp$")) %>%
  filter(str_detect(value, ".Regions", negate = T))

# Identify all glaciers that are covered by the six study regions.  We use 
# parallel processing,and export the 19 RGI regions to the clusters. 
# For each RGI region, we would like to obtain the 
# outlines of all glaciers that intersect with the extent of the study regions.

cl <- makeCluster(20) # Adapt to the number of available cores on YOUR machine.

# Export R-packages to the cluster, which will be used in the apply loop. Make sure 
# you installed these packages before.

clusterEvalQ(cl = cl, c(library("sf"),
                        library("stringr"),
                        library("tidyverse")))

# Export the file list of the glacier shapefiles and the extent of the study 
# regions to the clusters.

clusterExport(cl = cl, list("gl.list", "sr.wgs"))

# Iterate over the list of glaciers, read the glaciers to memory, and find the 
# glaciers in each RGI region that overlap with the study regions.

all.gl <- pblapply(1: nrow(gl.list), cl = cl, function (m) {
  
  reg.gl <- st_read(gl.list[m, ]) %>%
    st_make_valid()
  
  int <- st_intersects(reg.gl, sr.wgs, sparse = F)
  
  row.true <- apply(int, 1, any)
  
  if(any(row.true)) {
    
    sub.reg.gl <- reg.gl[row.true, ] %>%
      mutate(Region = sr.wgs$Region[apply(int[row.true, ], 1, which)])
    
    return (sub.reg.gl) }
  
})

# Close the cluster environment.

stopCluster(cl)
gc()

# Bind all glaciers in the study regions to one large multi-polygon.

all.gl.pol <- bind_rows(all.gl)
all.gl.rgiid <- all.gl.pol$RGIId

# We would now like to know at which altitudes in our study regions glacier ice 
# is stored. We subtract the gridded data of the glacier volumes from gridded
# data of the glacier surface. We first obtain a list of the surface DEMs of 
# all glaciers worldwide. Data are available under "surface_DEMs_RGI60.zip" at 
# https://www.research-collection.ethz.ch/handle/20.500.11850/315707. 
# Download the files and unzip them to one folder.

# List the files by changing the folder to the folder structure on your computer.

gl.dem.list <- list.files("D:/nrc_user/veh/LW_F/Farinotti_Ice_thickness_Raw_DEMs", 
                          pattern = ".tif$",
                          recursive = T, include.dirs = T, full.names = T) %>%
  as_tibble() %>%
  filter(str_detect(value, ".tif$")) 

# Find all glacier surface DEMs that exist for the study regions.

gl.dem.list <- gl.dem.list[str_sub(gl.dem.list$value, -18, -5) %in% all.gl.rgiid , ]

# Then obtain a list of ice thickness data from the same data source are available
# under "composite_thickness_RGI60-all_regions.zip" at
# https://www.research-collection.ethz.ch/handle/20.500.11850/315707.
# Download the files and unzip them to one folder.

# List the files by changing the folder to the folder structure on your computer.

gl.thick.list <- list.files("D:/nrc_user/veh/LW_F/Farinotti_Ice_thickness_consensus", 
                            recursive = T, include.dirs = T, full.names = T) %>%
  as_tibble() %>%
  filter(str_detect(value, ".tif$")) 

# Find all glacier thickness DEMs that exist for the study regions.

gl.thick.list <- gl.thick.list[str_sub(gl.thick.list$value, -28, -15) %in% all.gl.rgiid, ]

################################################################################
#######  Estimate the regional glacier volume in 100-m elevation bins   ########

# Set up the cluster environment. We use 50 cores - adjust to your system 
# specifications.

cl <- makeCluster(50)
registerDoParallel(cl)

# We iterate over all glacier polygons. We generate text files as output that  
# have the glacier volume in 100-m bins between -1,000 and 9,000 m a.s.l. 
# (thus 100 values). 

# Within the foreach-loop, adjust the folders for I/O operations to your 
# folder structure.

vol.per.bin <- foreach(i = 1: nrow(all.gl.pol),  
                       .packages =  c('sf', 'raster', 'exactextractr', 'tidyverse'),
                       .combine = cbind) %dopar% {
                         
                         # Select the shapefile for the glacier.
                         
                         gl.shp <- all.gl.pol[i, ]
                         
                         # Obtain the location of the file containing the 
                         # associate ice thickness. If there is no file,
                         # we report No Data for each elevation bin.
                         
                         thick.file <- gl.thick.list$value[grep(gl.shp$RGIId, gl.thick.list$value)]
                         if (length(thick.file) == 0) {return(rep(NA, 100))}
                         
                         # Read the ice thickness data to memory.
                         
                         th <- raster(thick.file)
                         
                         # Obtain the spatial resolution of the ice thickness data. 
                         # Note that resolution of the grids varies between the files.
                         
                         res.th <- res(th)
                         
                         # Obtain the location of the file containing the glacier surface DEM.
                         # If there is no file, we report No Data for each elevation bin.
                         
                         surf.file <- gl.dem.list$value[grep(gl.shp$RGIId, gl.dem.list$value)]
                         if (length(surf.file) == 0) {return(rep(NA, 100))}
                         
                         # Read the glacier surface DEM to memory.
                         
                         surf <- raster(surf.file)
                         
                         # Calculate the glacier bed topography.
                         
                         bed  <- surf-th
                         
                         # Generate a stack of rasters with three layers: the surface DEM, 
                         # the ice thickness, and the glacier bed topography.
                         
                         gl.stack <- stack(surf, th , bed)
                         names(gl.stack) <- c("surf", "thick", "bed")
                         
                         # Reproject the glacier polygon using the projection of 
                         # the gridded glacier data.
                         
                         if(st_crs(gl.shp) != st_crs(gl.stack)) { 
                           gl.shp <- st_transform(gl.shp, st_crs(gl.stack)) %>%
                                     st_make_valid()}
                         
                         # Obtain the cell values from these three layers. 
                         # Do some rough filtering of outliers, i.e. 
                         # glacier thicknesses that are 10 times higher or lower 
                         # than the standard deviation of observed thicknesses.
                         
                         df <- exact_extract(gl.stack, gl.shp) %>% 
                           bind_rows() %>%
                           as_tibble() %>%
                           drop_na() %>%
                           filter((thick < (median(thick) + (10* sd(thick)))) &  
                                    (thick > 0) &
                                    (surf > -20) &
                                    (surf < (median(surf) + (10* sd(surf))))) 
                         
                         # If the data are pure noise and no data remain after 
                         # filtering, write NA to disk.
                         
                         if (nrow(df) <= 1) {
                           
                           # Define a location where to store the output. 
                           
                           write.table(rep(NA, 100), 
                                       file = paste0("D:/nrc_user/veh/LW_F/Ice_vol_per_bin/", 
                                                     all.gl.pol$RGIId[i], ".txt"), 
                                       quote = F, 
                                       row.names = F, 
                                       col.names =  all.gl.pol$RGIId[i])
                           
                           return(rep(NA, 100))}
                         
                         # Otherwise, iterate in 100-m steps from -1000 m a.s.l. 
                         # (ice below sea level at outlet glaciers) 
                         # to the highest peaks on Earth (9000 m a.s.l.). 
                         # Count the number of cells in each bin and multiply
                         # with the cell resolution.
                         
                         mat <- mapply(FUN = function(x, y, z) {
                           
                           hist(seq(x, y, by = 1), 
                                breaks = seq(-1000, 9000, by = 100), 
                                plot = F)$count * z * res.th[1] * res.th[2]},
                           
                           x = df$bed,
                           y = df$surf,
                           z = df$coverage_fraction)
                         
                         rs <- rowSums(mat) 
                         
                         # Define a location where to store the output, 
                         # that is the glacier volume in each 100-m bin. 
                         # Use the same folder location as above. 
                         
                         write.table(rs, 
                                     file = paste0("D:/nrc_user/veh/LW_F/Ice_vol_per_bin/", 
                                                   all.gl.pol$RGIId[i], ".txt"), 
                                     quote = F, 
                                     row.names = F, 
                                     col.names =  all.gl.pol$RGIId[i])
                         
                         return(rs) }

# Stop the cluster environment.

stopCluster(cl)

# List the glacier volumes that were stored as text files to disk.
# Change the folder structure to your needs.

vol.txt <- list.files(path = "D:/nrc_user/veh/LW_F/Ice_vol_per_bin",
                pattern = ".txt", 
                include.dirs = T, 
                full.names = T)

# Read the binned glacier volumes to memory. 
# We use parallel computing, and adjust the number of cores, if necessary.

cl <- makeCluster(40)
clusterExport(cl = cl, list("vol.txt"))
vol.tab <- pbsapply(vol.txt, cl = cl, read.csv2, dec = ".")
parallel::stopCluster(cl)

# Bind all binned glacier volume to one large tibble and shorten the column names 
# to the RGI ID.

t.bind <- bind_cols(vol.tab)
colnames(t.bind) <- str_remove(basename(vol.txt), ".txt")

# Transpose the tibble.

t.bind2 <- as_tibble(t(t.bind)) %>%
  mutate(RGIId = names(t.bind))

# Add the binned glacier volume to the associate shapefile for each glacier.

glacier.elev <- all.gl.pol %>% 
  transmute(RGIId, Region) %>%
  left_join(., t.bind2, by = "RGIId")

# Summary statistics: Calculate the glacier volume in each region in 100-m 
# elevation bins. 

elev.stat <- glacier.elev %>%
  st_drop_geometry() %>%
  group_by(Region) %>%
  summarise(across(starts_with("V"), ~ sum(.x, na.rm = TRUE))) 

# The transposed version of the binned statistics.

elev.stat.2 <- elev.stat %>%
  dplyr::select(!Region) %>%
  t() %>%
  as_tibble() 

colnames(elev.stat.2) <- elev.stat$Region

elev.stat.glaciers <- elev.stat.2 %>% 
  mutate(mids = seq(-950, 8950, by = 100)) %>%
  pivot_longer(cols = !mids, 
               names_to = "Region",
               values_to = "Glacier_vol") 

################################################################################
#######       Estimate the regional glacier volume loss                 ########
#######     between 2000 and 2019 in 100-m elevation bins               ######## 

# Use the dataset by Hugonnet et al. (2021) to calculate the change in glacier 
# volume in each study region. Download the elevation change maps for the 
# period Jan 2000 to Dec 2019 from 
# https://www.sedoo.fr/theia-publication-products/?uuid=c428c5b9-df8f-4f86-9b75-e04c778e29b9.
# Put all GeoTiffs into the same folder. 

# Generate a list of elevation change maps. The elevation change maps are provided as rates 
# (meters per year) during the 20-year study period. 

# Change the folder structure to your folder.

dhdts <- list.files("D:/nrc_user/veh/LW_F/Hugonnet_glacier_elevation_change",
                    pattern = "_dhdt.tif$", 
                    full.names = T, 
                    include.dirs = T, recursive = T)

# We iterate over the list of elevation change maps to identify the maps that 
# intersect with  the study regions.

cl <- makeCluster(50) # Change the number of cores, if necessary.

# Export packages and the list of elevation change maps (dhdts) to the cluster

clusterEvalQ(cl = cl, c(library("sf"),
                        library("raster"),
                        library("tidyverse")))

clusterExport(cl = cl, list("dhdts"))

dem.ex <- pblapply(dhdts, cl = cl, function (y) {
  
  x <- raster(y) 
  
  # We generate a polygon from the dhdt map.
  
  p <- st_as_sfc(st_bbox(x)) %>%
    st_sf() %>%
    st_transform(4326) %>%
    mutate(scene = y)

  return(p)})

# Stop the cluster environment.

stopCluster(cl)
gc()

# Bind the polygons containing the extents of the dhdt maps to one large 
# multi-polygon.

dem.ex <- bind_rows(dem.ex)

# Again, we iterate over the glacier polygons in each study region. We would like to 
# obtain the loss in glacier surface elevation in 100-m elevation bins, which we
# convert to volumes using the spatial resolution of the dataset.

library(doParallel)
cl <- makeCluster(50)
registerDoParallel(cl)

vol.loss.per.bin <- foreach(i = 1:nrow(all.gl.pol), 
                            .packages =  c('sf', 'raster', 'exactextractr', 'tidyverse'),
                            .combine = cbind) %dopar%   {
                              
                              # Select the glacier polygon.
                              
                              gl.shp <- all.gl.pol[i, ]
                              
                              # Select the corresponding surface DEM file. We will subtract the
                              # elevation change rate from this DEM.
                              
                              surf.file <- gl.dem.list$value[grep(gl.shp$RGIId, gl.dem.list$value)]
                              
                              # If there is no surface DEM, write a text file containing only NA
                              # for the elevation bins to disk.
                              
                              if (length(surf.file) == 0) {
                                
                                return(rep(NA, 100))
                                
                                # Adjust the path to YOUR output folder.
                                
                                write.table(rep(NA, 100), 
                                            file = paste0("D:/nrc_user/veh/LW_F/Ice_vol_per_bin/", 
                                                          all.gl.pol$RGIId[i], ".txt"), 
                                            quote = F, 
                                            row.names = F, 
                                            col.names =  all.gl.pol$RGIId[i])
                              }
                              
                              # Read the surface DEM to disk.
                              
                              surf <- raster(surf.file)
                              
                              # Obtain the resolution of the surface DEM.
                              
                              res.surf <- res(surf)
                              
                              # The spatial resolution of the surface DEMs are 
                              # not the same for all files. However,
                              # the dhdt maps always have a resolution of 
                              # 100 x 100 m throughout. Therefore,
                              # reproject the surface DEMs to a spatial resolution 
                              # of 100 m, if necessary.
                              
                              if (res.surf[1] < 100) {
                                
                                surf.agg <- aggregate(surf, 
                                                      fact = 100/ res.surf[1], 
                                                      fun = median)
                                
                              } else if (res.surf[1] == 100) { 
                                
                                surf.agg <- surf
                                
                              } else { surf.agg <- disaggregate(surf, 
                                                                fact = res.surf[1] / 100, 
                                                                method = 'bilinear') }
                              
                              # Identify the tile of the elevation change map that intersect with the 
                              # glacier polygon.
                              
                              loss.dem.ex <- dem.ex[st_intersects(dem.ex, gl.shp, sparse = F)[, 1], ]
                              
                              # The polygon can intersect with more than one dhdt tile. If so,
                              # aggregate them to one large tile.
                              
                              if(nrow(loss.dem.ex) > 1) {
                                
                                # Read all dhdt maps as rasters to memory.
                                
                                loss.dems <- lapply(loss.dem.ex$scene, raster)
                                
                                # Reproject to the same UTM zone as the surface 
                                # DEM, and use a consistent spatial resolution
                                #  of 100 x 100 m.
                                
                                loss.dems <- lapply(loss.dems, function (x) {
                                  
                                  if(st_crs(x) != st_crs(surf.agg)) {
                                    
                                    reproj <- projectRaster(from = x,
                                                            to = surf.agg,
                                                            res = c(100, 100)) 
                                    
                                    if(all(is.na(getValues(reproj)))) { return() } else { return(reproj) }
                                    
                                  } else { 
                                    
                                    ex.dem <- st_as_sfc(st_bbox(x)) 
                                    ex.surf.agg <- st_as_sfc(st_bbox(surf.agg)) 
                                    
                                    if (st_intersects(ex.dem, ex.surf.agg, sparse = F)[, 1])  {
                                      
                                      resamp <- resample(x = x, y = surf.agg)
                                      
                                      if(all(is.na(getValues(resamp)))) { return() } else { return(resamp) }
                                      
                                    } else { return() } 
                                    
                                  } 
                                  
                                })
                                
                                loss.dems <- loss.dems[sapply(loss.dems, function (x) !is.null(x))]
                                
                                # Check if list is empty
                                
                                if (length(loss.dems) == 0) { 
                                  
                                  # If there is only No Data in the dhdt tile, 
                                  # write a text file to disk that contains 
                                  # only NA for each elevation bin.
                                  # Change the folder location, if necessary.
                                  
                                  write.table(rep(NA, 100), 
                                              file = paste0("D:/nrc_user/veh/LW_F/Vol_loss_per_bin/", 
                                                            all.gl.pol$RGIId[i], ".txt"), 
                                              quote = F, 
                                              row.names = F, 
                                              col.names =  all.gl.pol$RGIId[i])
                                  
                                  return(rep(NA, 100)) 
                                  
                                } else if (length(loss.dems) == 1)  { 
                                  
                                  loss.dem <-  loss.dems[[1]]
                                  
                                  } else {
                                  
                                   # Finally, merge all dhdt tiles to a single tile. 
                                    
                                   loss.dem <-  do.call(merge, loss.dems) }
                                
                              } else {  loss.dem <- raster(loss.dem.ex$scene) # case when there is only one dhdt tile.
                              
                              if(st_crs(loss.dem) != st_crs(surf.agg)) {
                                
                                loss.dem <- projectRaster(from = loss.dem, 
                                                          to = surf.agg,
                                                          res = c(100, 100))}
                              
                              }
                              
                              crop.loss <- resample(loss.dem, surf.agg)
                              
                              # Loss data are average rates of mass loss in a
                              # 20-year period (2000-2019). Therefore,
                              #  multiply the elevation change maps with 20.
                              # We subtract the elevation loss from the surface 
                              # and call this elevation the 'bed', although this 
                              # is not the bed in reality, given that ice may 
                              # still remain in a given glacier cell.
                              
                              bed <- surf.agg - (crop.loss*20)
                              
                              # Generate a stack of all three layers.
                              
                              gl.stack <- stack( surf.agg, crop.loss, bed)
                              names(gl.stack) <- c("surf", "loss", "bed")
                              
                              # Reproject the shapefile to the projection of the stack.
                              
                              gl.shp.repro <- st_transform(gl.shp, st_crs(gl.stack))
                              
                              # Extract cell values from the stack, and do some raw 
                              # filtering of outliers in the mass loss and elevation data.
                              
                              df <- exact_extract(gl.stack, gl.shp.repro) %>%
                                bind_rows() %>%
                                as_tibble() %>%
                                drop_na() %>%
                                filter((loss < (median(loss) + (10* sd(loss)))) & (loss > (median(loss) - (10* sd(loss)))) &
                                         (surf > -20 ) & (surf < (median(surf) + (10* sd(surf))))) 
                              
                              # If there is only noise (i.e. no cells remain 
                              # after filtering), write NA to disk.
                              
                              if (nrow(df) <= 1) {
                                
                                write.table(rep(NA, 100), 
                                            file = paste0("D:/nrc_user/veh/LW_F/Vol_loss_per_bin/", 
                                                          all.gl.pol$RGIId[i], ".txt"), 
                                            quote = F, 
                                            row.names = F, 
                                            col.names =  all.gl.pol$RGIId[i])
                                
                                return(rep(NA, 100))}
                              
                              # Otherwise, calculate the elevation loss 
                              # (negative values) or elevation gain (positive values)
                              # in 100-m elevation bins
                              
                              mat <- mapply(FUN = function(x, y, z) {
                                
                                if( x > y ) { # i.e. mass gain
                                  
                                  hist(seq(round(y), round(x), by = 1),
                                       breaks = seq(-1000, 9000, by = 100), 
                                       plot = F)$count * z * 100 * 100
                                  
                                }  else { # i.e. mass loss
                                  
                                  -1 * (hist(seq(round(x), round(y),  by = 1),
                                             breaks = seq(-1000, 9000, by = 100), 
                                             plot = F)$count * z * 100 * 100)
                                }
                              },
                              
                              x = df$surf,  y = df$bed,  z = df$coverage_fraction)
                              
                              # Calculate the sum of elevation change in each bin
                              
                              rs <- rowSums(mat) 
                              
                              # Write the text file to disk. Change the output 
                              # location to YOUR folder.
                              
                              write.table(rs, 
                                          file = paste0("D:/nrc_user/veh/LW_F/Vol_loss_per_bin/", 
                                                        all.gl.pol$RGIId[i], ".txt"), 
                                          quote = F, 
                                          row.names = F, 
                                          col.names =  all.gl.pol$RGIId[i])
                              
                              return(rs)
                              
                            } 

# Stop the cluster environment.

stopCluster(cl)

# Similar to workflow on glacier volumes, read the ice loss data from 
# disk using parallel computing.

# Change the folder accordingly.

ice.loss.txt <- list.files(path = "D:/nrc_user/veh/LW_F/Vol_loss_per_bin",
                           pattern = ".txt", 
                           include.dirs = T, 
                           full.names = T)

# Set up cluster, export the list of text files, and read the elevation change 
# data in parallel to memory.

cl <- makeCluster(50)
clusterExport(cl = cl, list("ice.loss.txt"))
tab.ice.loss <- pbsapply(ice.loss.txt, cl = cl, read.csv2, dec = ".")
stopCluster(cl)

# Bind the binned data on glacier elevation change to one large tibble.

t.ice.loss.bind <- bind_cols(tab.ice.loss)
colnames(t.ice.loss.bind) <- str_remove(basename(ice.loss.txt), ".txt")

# Add the RGI ID to the table.

t.ice.loss.bind2 <- as_tibble(t(t.ice.loss.bind)) %>%
  mutate(RGIId = names(t.ice.loss.bind))

# Add the elevation change per bin to the glacier shapefile.

vol.loss.elev <- all.gl.pol %>% 
  transmute(RGIId, Region) %>%
  left_join(., t.ice.loss.bind2, by = "RGIId")

# Calculate regional statistics of glacier volume change.

vol.loss.stat <- vol.loss.elev %>%
  st_drop_geometry() %>%
  group_by(Region) %>%
  summarise(across(starts_with("V"), ~ sum(.x, na.rm = TRUE))) 

# Generate a table with volume change for each 100-m elevation bin per region.

vol.loss.stat.2 <- vol.loss.stat %>%
  dplyr::select(!Region) %>%
  t() %>%
  as_tibble() 

colnames(vol.loss.stat.2) <- vol.loss.stat$Region

# Convert to one long tibble.

vol.loss.stat.glaciers <- vol.loss.stat.2 %>% 
  mutate(mids = seq(-950, 8950, by = 100)) %>%
  pivot_longer(cols = !mids, 
               names_to = "Region",
               values_to = "Vol_loss")


# Join both the statistics of total glacier volume and glacier volume loss 
# between 2000 and 2019.

full.data <- left_join(elev.stat.glaciers, vol.loss.stat.glaciers) 

# Change to your location.

saveRDS(full.data, 
        file = "U:/R-Scripts_Full_Himalaya/Regional_glacier_and_melt_volumes.rds" )
