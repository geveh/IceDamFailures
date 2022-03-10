####################################################################################################
#######     Estimate regional trends in the source elevation of ice-dammed failures    #############
#######                                                                                #############
#######                                  by Georg Veh                                  #############
#######                                   12 Oct, 2021                                 #############
#######   added code to extract elevations of GLOF source locations on 06 Dec, 2021    #############
#######                      more comments added on March 10, 2021                     #############
####################################################################################################

require(sf)
require(tidyverse)
require(pbapply)
require(modelr)
require(parallel)
require(exactextractr)
require(brms)
require(tidybayes)
require(ggpubr)
require(scales)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

#######################

# Set YOUR working directory

setwd("D:/nrc_user/veh/LW_F/Shapefiles")

# Load input data

# Shapefile containing the extent of the six study regions.

sr <-  read_sf("Region_extents/Extent_pol.shp") %>% 
  st_transform(4326)

# Database of present-day outlines of ice-dammed lakes. The database is merge of several regional inventories.
# 
# - Norway:             Nagy & Andreassen (2019): Glacier lake mapping with Sentinel-2 imagery in Norway., 54p. 
# - NW North America:   Rick, B. et al. (2022): Dam type and lake location characterize ice-marginal 
#                       lake area change in Alaska and NW Canada between 1984 and 2019. The Cryosphere 
#                       16, 297-314 (2022).
# - High Mountain Asia: Chen et al. (2021): Annual 30 m dataset for glacial lakes in High Mountain Asia from 
#                       2008 to 2017. ESSD 13, 741-766.
# - Andes:              Wilson et al. (2018): Glacial lakes of the Central and Patagonian Andes. Global and
#                       Planetary Change 162, 275-291.
# - European Alps:      Buckel et al. (2018): Glacial lakes in Austria-Distribution and formation since the 
#                       Little Ice Age. Global and Planetary Change, 164, 39-51.
#                       Mölg et al. (2021): Inventory and evolution of glacial lakes since the Little Ice Age:
#                       Lessons from the case of Switzerland. Earth Surface Processes and Landforms, 46(13), 2551-2564.
# - Iceland:            Own mapping by Natalie Luetzow, University of Potsdam.

gdl.database.updated <- st_read("gdl_database/gdl_database.shp")

# Load the tibble with all reported ice-dammed GLOFs between 1900 and 2021.

all.glofs <- readRDS("all_glofs_tibble.RDS")

# Useful functions

scale_this <- function(x){
  (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
}

################################################################################################
#######    Add elevation (Z) to present-day ice-dammed lakes     ###############################

# Register to get access to the ALOS Global Digital Surface Model "ALOS World 3D - 30m (AW3D30)"
# at https://www.eorc.jaxa.jp/ALOS/en/dataset/aw3d30/aw3d30_e.htm
# Download the DEMs intersecting with the study regions, unzip, and store them in one folder.

# Generate a list of all available ALOS DEMs.

alos <- list.files("D:/nrc_user/veh/LW_D/DEMs/ALOS_3_2",
                   pattern = "_DSM.tif$", 
                   full.names = T, 
                   include.dirs = T, 
                   recursive = T)

# Derive the extent of each DEM tile, and convert to a polgyon.

dem.ex <- pblapply(alos,  function (y) {
  
  x <- raster::raster(y) 
  
  p <- st_as_sfc(st_bbox(x)) %>%
    st_sf() %>%
    st_transform(4326) %>%
    mutate(scene = y)
  
  return(p)})

# Bind all extents to one large multi-polygon.

dem.ex <- bind_rows(dem.ex)

# Write the polygon with the DEM extents to disk.

st_write(dem.ex , "dem_tiles.shp", delete_layer =  T)
# dem.ex <- st_read("dem_tiles.shp")

# We iterate over all ice-dammed lakes and extract the median elevation of the 
# lake surface from the ALOS DEM. We use parallel processing with 15 cores.

cl <- makeCluster(15) # Change the number of cores, if necessary.

# Export the packages, the lake outlines, and the DEM extents to the clusters.

clusterEvalQ(cl = cl, c(library("sf"),
                        library("stringr"),
                        library("tidyverse"),
                        library("exactextractr"),
                        library("raster")))

clusterExport(cl = cl, list("gdl.database.updated", "dem.ex"))

# Identify the median elevation of each ice-dammed lake.

all.z <- pblapply(1: nrow(gdl.database.updated), cl = cl, function (m) {
  
  # Select the lake of interest.
  
  reg.gl <- gdl.database.updated[m, ] %>%
    st_make_valid()
  
  # Identify the DEM tiles that intersect with the lake.
  
  int.dem.ex <- dem.ex[st_intersects(dem.ex, reg.gl, sparse = F)[, 1], ]
  
  # Return NA, if there is no DEM intersecting with the lake (which should not be the case, if you
  # have downloaded all ALOS DEMs.)
  
  if(nrow(int.dem.ex) == 0) {return(NA)}
  
  # If more than one DEM intersect with the DEM, merge the individual tiles to one large DEM.
  
  if(nrow(int.dem.ex) > 1 ) {
    
    int.dems <- lapply(int.dem.ex$scene, raster)
    
    int.dem <-  do.call(merge, int.dems)
    
  } else {  int.dem <- raster(int.dem.ex$scene) }
  
  # Extract the median elevation.
  
  z <- exact_extract(int.dem, reg.gl, fun = 'median') 
  
  return(z)
  
})

# Stop the cluster environment.

stopCluster(cl)

# Add the elevation to the lake polygons and select only the lakes with the highest level
# of confidence of being ice-dammed according to visual image interpretation. 
# In this case, the column "ice.dammed" is flagged "yes" (otherwise "no" or "unsure").

gdl.database.updated <- gdl.database.updated %>%
  mutate(Elevation_m = unlist(all.z)) %>% 
  filter(ice.dammed == "yes")

# From the regional (in many cases multi-temporal) inventories of ice-dammed lakes,
# we select only the lakes with the most recent time stamps, and convert the polygons to points.

gdl.database.centroid <- gdl.database.updated %>%
  filter(!is.na(Elevation_m)) %>%
  filter((Region == "Andes" & Year >= 2015) |
         (Region == "European Alps" & Year >= 2015) |
         (Region == "Pacific NW" & Year == 2018) |
         (Region == "High Mountain Asia" & Year == 2017) |
         (Region == "Scandinavia" & Year == 2018) |
         (Region == "Iceland" & Year > 2018)) %>%
  st_centroid() %>%
  rename(region = Region) %>%
  mutate(region = replace(region, region == "Pacific NW", "NW North America"))

################################################################################################
#######    Add elevation (Z) to the locations of historic ice-dam failures     #################

# Convert the table of reported GLOFs to a point feature geometry.

all.glofs.pts <- all.glofs %>%
  filter(!is.na(Latitude))  %>% 
  st_as_sf(coords = c( "Longitude", "Latitude"), crs = 4326) %>%
  mutate(Latitude = all.glofs %>% filter(!is.na(Latitude)) %>%.$Latitude,
         Longitude = all.glofs %>% filter(!is.na(Latitude)) %>%.$Longitude) 

# We iterate over each observed GLOFs in parallel. Setup the clusters, and export the necessary
# packages and files to the clusters.

cl <- makeCluster(15)
clusterEvalQ(cl = cl, c(library("sf"),
                        library("tidyverse"),
                        library("raster")))

clusterExport(cl = cl, list("all.glofs.pts", "dem.ex"))

# Obtain the elevation for each GLOF source location.

all.z.glofs <- pblapply(1: nrow(all.glofs.pts), cl = cl, function (m) {
  
  # Select the GLOF of interest.
  
  reg.gl <- all.glofs.pts[m, ] %>%
    st_make_valid()
  
  # Find the DEM tiles that intersect with the GLOF.
  
  int.dem.ex <- dem.ex[st_intersects(dem.ex, reg.gl, sparse = F)[, 1], ]
  
  if(nrow(int.dem.ex) == 0) {return(NA)}
  
  if(nrow(int.dem.ex) > 1 ) {
    
    int.dems <- lapply(int.dem.ex$scene, raster)
    
    int.dem <-  do.call(merge, int.dems)
    
  } else {  int.dem <- raster(int.dem.ex$scene) }
  
  # Extract the elevation of the GLOF from the DEM.
  
  z <- raster::extract(int.dem, reg.gl) 
  
  return(z)
  
})

# Stop the cluster enviroment.

stopCluster(cl)

# Add the elevation to the point shapefile containing all historical GLOFs.

all.glofs.pts <- all.glofs.pts %>%
  mutate(Elevation_m = unlist(all.z.glofs))

# Select only GLOFs that happened between 1900 and 2021. For any lake with repeat outbursts, 
# we use only the first reported GLOF that happened during this period

glofs.with.z <- all.glofs.pts %>%
  mutate(region = replace(region, region == "Pacific NW", "NW North America")) %>%
  filter(!is.na(rounded_year),
         !is.na(Elevation_m),
         rounded_year >= 1900) %>% 
  group_by(RGI_Glacier_Id, Glacier, Lake, Lake_type, region) %>%
  summarise(RGI_Glacier_Id = unique(RGI_Glacier_Id), 
            Lake = unique(Lake), 
            Lake_type = unique(Lake_type), 
            region = unique(region),
            first_glof = min(rounded_year, na.rm = T),
            Z = median(Elevation_m),
            Latitude = mean(Latitude),
            Longitude = mean(Longitude))

# Only select ice-dammed lakes. Scale both the predictor (year) and the outcome (Z) using
# a mean of zero and unit standard deviation.

glofs.ice.with.z <- glofs.with.z %>%
  ungroup() %>%
  filter(Lake_type == "ice") %>%
  filter(!((region == "Andes") & (Latitude > -40 ))) %>%
  mutate(Z_scale = scale_this(Z),
         year_scale = scale_this(first_glof)) %>%
  st_drop_geometry()

# Run the Bayesian hierarchical model of GLOF source elevation Z versus time.
# Set t-distributed priors on the intercept and slope, maintain all other default priors
# set in brms.

bprior <- prior(student_t(3, 0, 5), class = "Intercept") + 
          prior(student_t(3, 0, 2.5), class = "b") 

z.year.brm <- brm(Z_scale ~ year_scale + (year_scale | region),
                  data    = glofs.ice.with.z,
                  family  = student(),
                  prior   = bprior,
                  cores   = 4,
                  chains  = 4,
                  warmup  = 2000,
                  iter    = 6000,
                  control = list(adapt_delta = 0.9,
                                 max_treedepth = 14))


# Save the regional models of Z versus year to disk.

saveRDS(z.year.brm, "Z_trends_per_region.RDS")
# z.year.brm <- readRDS("Z_trends_per_region.RDS")

# Assess model fit 

summary(z.year.brm)
plot(z.year.brm)

# Obtain the posterior trend (i.e. the regression slope) in GLOF source elevation with time.
# Re-transform the standardized posterior distributions to original scale, and plot them
# for each study region.

elev.trend <- z.year.brm %>%
  spread_draws(b_year_scale, `r_region`[region, param]) %>%
  filter(param == "year_scale") %>%
  mutate(Z_mean = (b_year_scale + r_region) * sd(glofs.ice.with.z$Z) / sd(glofs.ice.with.z$first_glof) ) %>%
  ungroup() %>%
  mutate(region = str_replace_all(region, "[.]", " "),
         region = reorder(region, Z_mean, median)) %>%
  ggplot(aes(x = Z_mean,
             y = region,
             fill = stat(x > 0))) +
  stat_halfeye(.width = 0.95) +
  scale_fill_manual("Posterior trend > 0", 
                    values = RColorBrewer::brewer.pal(3, "Greens")[1:2]) +
  theme_bw() +
  labs(x = expression(paste("Elevation change of GLOF source lakes [m ", yr^{-1}, ']')),
       y = "Region") +
  geom_vline(xintercept = 0)  + 
  theme( axis.text =   element_text(size = 8.5),
         axis.text.x = element_text(size = 8.5),
         axis.title =  element_text(size = 9),
         strip.text =  element_text(size = 8.5)) + 
  theme(legend.position = "none")

# Predict the change in GLOF source elevation for new data, 
# i.e. for the entire period 1900-2021.

preds <-  glofs.ice.with.z %>%
  data_grid(year_scale = seq_range(year_scale, n = 51),
            region  = unique(region)) %>%
  add_epred_draws(
    object = z.year.brm, 
    newdata = .,
    value = "Z", 
    ndraws = 1000) %>%
  mutate(first_glof = year_scale * sd(glofs.ice.with.z$first_glof) + mean(glofs.ice.with.z$first_glof),
         Z = Z * sd(glofs.ice.with.z$Z) + mean(glofs.ice.with.z$Z)) 

# Load the regional glacier volumes and ice losses between 1900 and 2021. Those will be 
# added at the right margin of the plot showing the trend in GLOF elevation with time. 

gl.melt <- readRDS(file = "Regional_glacier_and_melt_volumes.rds") %>%
  rename(region = Region)

gl.reg <- unique(gl.melt$region)

# Convert the regional glacier volumes and melt rates to density-like objects that
# are plotted on the right margin of the plot (in a value range between 2025 and 2040).
# The values of the glacier melt rates are rescaled such that they are relative to 
# the total glacier volume.

reg.vol.melt <- lapply(gl.reg, function (x) {
  
  gl.vol <- gl.melt %>%
    filter(region == x) %>%
    filter(Glacier_vol > 0)
  
  csum <- cumsum(gl.vol$Glacier_vol) / sum( gl.vol$Glacier_vol)
  idx <- min(which(csum >= 0.005)) : max(which(csum < 0.995 ))
  
  gl.vol <- gl.vol[idx, ]  
  gl.vol <- gl.vol %>%
    add_row(mids = min(gl.vol$mids),
            Glacier_vol = 0,
            region = x,
            Vol_loss = 0,
            .before = 1) %>%
    add_row(mids = max(gl.vol$mids),
            Glacier_vol = 0,
            region = x,
            Vol_loss = 0) %>% 
    mutate(Vol_loss = replace(Vol_loss, Vol_loss >= 0, 0)) %>%
    mutate(Glacier_vol_rescale = rescale(Glacier_vol, 
                                         from = c(0, max(Glacier_vol)), 
                                         to = c(2040, 2025)), 
           Vol_loss_rescale = rescale(-Vol_loss, 
                                      from = c(0, max(Glacier_vol)), 
                                      to   = c(2040, 2025))) 
  
  return(gl.vol)
  
})

reg.vol.melt <- bind_rows(reg.vol.melt) %>%
  mutate(region = str_replace(region, "Pacific NW", "NW North America"))

# Count the number of individual GLOF source locations per region. This information
# will be plotted in the lower left corner (ymin).

n.glofs <- glofs.ice.with.z %>% 
  group_by(region) %>%
  summarise(ymin = min(Z, na.rm = T),
            n_glofs = n()) %>%
  mutate(x = 1900,
         col = "navy")

# Count the number of present-day lakes in each region. This information
# will be plotted in the lower right corner (ymin).

n.lakes <- gdl.database.centroid %>% 
  group_by(region) %>% 
  summarise(n_lakes = n()) %>%
  mutate(x = 2040,
         col = "darkorange") %>%
  st_drop_geometry() %>%
  left_join(., n.glofs %>% select(region, ymin), by = 'region')

# Generate the plot showing the change in GLOF elevation with time (Figure 3).
# Add density-like polygons on the right margin showing the total ice volume in grey,
# and ice loss between 2000 and 2019 in black.

elev.ch.plot <- preds %>%
  ungroup() %>%
  mutate(region = replace(region, region == "Pacific NW", "NW North America")) %>%
  ggplot(aes(x = first_glof, 
             y = Z)) +
  facet_wrap(~region, scales = "free_y", ncol = 2) +
  stat_lineribbon(aes(y = Z), 
                  .width = 0.95,
                  point_interval = mean_qi,
                  fill = RColorBrewer::brewer.pal(3, "Greens")[2]) +
  scale_fill_discrete(name = "Posterior\ninterval") +
  geom_point(data = glofs.ice.with.z, 
             mapping = aes(x = first_glof,
                           y = Z), 
             shape = 16, 
             size = 2.2,
             alpha = 0.7,
             color = "navy") + 
  geom_polygon(data = reg.vol.melt,
               mapping = aes(x = Glacier_vol_rescale,
                             y = mids),
               color = NA,
               fill = "grey80") +
  geom_polygon(data = reg.vol.melt,
               mapping = aes(x = Vol_loss_rescale,
                             y = mids),
               fill= "grey20" ) +
  geom_rug(data = gdl.database.centroid,
           mapping = aes(y = Elevation_m),
           inherit.aes = F,
           sides = "r",
           color = "darkorange") +
  geom_text(aes(x = x, 
                y = ymin, 
                label = n_lakes),
            color = "darkorange",
            size = 3,
            data = n.lakes) +
  geom_text(aes(x = x, 
                y = ymin, 
                label = n_glofs),
            color = "navy",
            size = 3,
            data = n.glofs) + 
  labs(x = "First documented year of an outburst from a given ice-dammed lake",
       y = "Elevation of glacier lakes [m a.s.l.]") +
  theme_bw()  +
  theme( axis.text =   element_text(size = 8.5),
         axis.text.x = element_text(size = 8.5),
         axis.title =  element_text(size = 9),
         strip.text =  element_text(size = 8.5)) +
  scale_x_continuous(breaks = seq(1900, 2025, by = 25),
                     labels = seq(1900, 2025, by = 25),
                     limits = c(1895, 2045),
                     expand = expansion())


# Combine the plots showing the change over time and the posterior slope.

arr.elev.trend <- ggarrange(plotlist = list(elev.ch.plot , elev.trend), 
                            labels = "auto",
                            ncol = 2, 
                            align = "v", 
                            legend = "none",
                            common.legend = T, 
                            widths = c(2,1)) +
  theme(plot.margin = margin(0.2, 1.3, 0.2, 0.2, "cm")) 

# Write the plots showing the change in Z to disk (Figure 3).

ggsave(
  filename = "elev_trend.pdf",
  plot = arr.elev.trend ,
  width = 250,
  height = 150,
  units = "mm"
)

ggsave(
  filename = "elev_trend.png",
  plot = arr.elev.trend ,
  width = 250,
  height = 150,
  units = "mm"
)

# Show differences in elevation (Z) between historic GLOF source locations and 
# present-day lakes.

glofs.and.present <- bind_rows(gdl.database.centroid %>% 
                                 dplyr::select(region, Elevation_m) %>%
                                 rename("Z" = Elevation_m) %>%
                                 st_drop_geometry() %>%
                                 mutate(Lakes = "present-day"),
                               glofs.ice.with.z %>%
                                 dplyr::select(region, Z) %>%
                                 mutate(Lakes = "with GLOFs\n(1900-2021)"))

# Obtain median elevation for the two groups.

med.elev <- glofs.and.present %>%
  group_by(region, Lakes) %>%
  summarise(med_elev = round(median(Z)),
            max_elev = max(Z))

# Show both distributions of Z grouped by study region (Extended Data Figure 10).

lake.elev.plot <- glofs.and.present %>%
  mutate(region = reorder(region, Z, median)) %>%
  ggplot(aes(x = Z,
             y = region,
             fill = Lakes)) +
  stat_halfeye(.width = 0.95,
               slab_size = 1,
               interval_size = 2,
               shape = 21,
               point_color = "black",
               point_fill = "white",
               position = position_dodge(width = -0.5),
               interval_color = "grey50",
               slab_alpha = 0.8) +
  scale_fill_manual(values = c("darkorange", "navy")) +
  geom_text(data = med.elev,
            aes(-150,
                region, 
                label = med_elev, 
                group = Lakes), 
            position = position_dodge(width = -0.5),
            size  = 3, 
            vjust = "center",
            hjust = "center") +
  theme_bw() +
  xlab("Elevation [m a.s.l.]") +
  ylab("Region") +
  theme(legend.position = c(0.75, 0.25),
        legend.background = element_rect(fill = "white", color = "black"))

# Write the plot showing the distribution of Z to disk.

ggsave(filename = "Lake_GLOF_elevation.pdf",
       plot = lake.elev.plot,
       width = 130,
       height = 100, 
       unit = "mm")  

ggsave(filename = "Lake_GLOF_elevation.png",
       plot = lake.elev.plot,
       width = 130,
       height = 100, 
       unit = "mm")  

###############################################################################################
