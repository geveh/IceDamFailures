####################################################################################################
#######               Estimate trends of local GLOF volume and peak discharge          #############
#######                         with elevation change of the glacier dam               #############
#######                                                                                #############
#######                                  by Georg Veh                                  #############
#######                                 03 March, 2022                                 #############
#######                      more comments added on 08 March, 2022                     #############
####################################################################################################

# Load the following packages, or use install.packages("nameofpackage"), if some of them
# are not pre-installed. In some cases you need to restart your R session.

require(tidyverse)
require(lubridate)
require(brms)
require(modelr)
require(tidybayes)
require(patchwork)
require(rnaturalearth)
require(sf)
require(ggthemes)

# Navigate to YOUR working directory where you stored the input files.

setwd("C:/Users/local-admin/Desktop/Lake_area_volume/")

# Useful functions

scale_this <- function(x){
  (x - mean(x, na.rm=TRUE)) / sd(x, na.rm=TRUE)
}

# Load the tables of local cumulative elevation change of glaciers impounding lakes between 2000 and 2019.
# The elevation changes are aggregated by RGI region such that one file may contain data on elevation
# changes for more than one glacier dam. The folder of glacier elevation changes (dhdt) is called dh_pergla_cut. 
# Please contant Romain Hugonnet (romain.hugonnet@legos.obs-mip.fr), if you have any questions on data
# processing. We use the cumulative mass change, in which the first observation on Jan 01, 2000, is set to zero.

dhdt <- list.files(pattern = "cumul.csv$",
                   path = "dh_pergla_cut",
                   full.names = T, 
                   include.dirs = T)

# Load the tables with all GLOFs, and those tables with dated flood volumes V0 or peak discharge Qp.

all.glofs    <- readRDS("all_glofs_tibble.RDS")
all.glofs.V0 <- readRDS("all_glofs_V0_tibble.RDS")
all.glofs.qp <- readRDS("all_glofs_qp_tibble.RDS")

################################################################################################
#######     Correlation between glacier elevation change and GLOF peak discharge Qp    #########

# Select only dated ice-dammed GLOFs with reported Qp that happened between 2000 and 2019.

glofs.after2000.qp <- all.glofs.qp %>% 
  group_by(Glacier_and_lake) %>% 
  filter(nchar(Date) == 10) %>%
  filter((rounded_year >= 2000) & (rounded_year < 2020)) %>%
  filter(n() >= 5) %>%
  ungroup() %>%
  mutate(RGI_Glacier_Id = str_replace_all(RGI_Glacier_Id,  "_", "."),
         year_dec_glofs = decimal_date(ymd(Date)),
         dh = NA,
         dh_err = NA,
         time_dh = NA) 

# We now iterate over the tables containing the glacier elevation changes for a given RGI region.

for (i in 1: length(dhdt)) {
  
  # Read the file with cumulative elevation changes into memory.
  
  f <- read_csv(dhdt[i], 
                col_types = cols_only(rgiid = col_character(), 
                                      time = col_character(),
                                      dh = col_double(),
                                      err_dh = col_double())) 
  
  # Select all glaciers (using their RGI ID) that have reported values of Qp.
  
  f2 <- f %>% 
    filter(rgiid %in% unique(glofs.after2000.qp$RGI_Glacier_Id)) %>%
    mutate(year_dec_dh = decimal_date(ymd(time)),
           rounded_year = as.numeric(str_sub(time, 1, 4))) %>% 
    rename(RGI_Glacier_Id = rgiid) 
  
  # Iterate over table containing the dated GLOFs with reported values of Qp.
  # Add the observed cumulative glacier elevation change (+ error) immediately before 
  # the outburst to the table.
  
  for (j in 1: nrow(glofs.after2000.qp)) {
    
    if(!(any(glofs.after2000.qp$RGI_Glacier_Id[j] == f2$RGI_Glacier_Id))) next
    
    f3 <- f2 %>% 
      filter(RGI_Glacier_Id == glofs.after2000.qp$RGI_Glacier_Id[j])
    
    idx <- max(which( (f3$year_dec_dh - glofs.after2000.qp$year_dec_glofs[j]) < 0))
    
    glofs.after2000.qp$dh[j] <- f3$dh[idx]
    glofs.after2000.qp$dh_err[j] <- f3$err_dh[idx]
    glofs.after2000.qp$time_dh[j] <- f3$time[idx]
    
  }
  
}

# Some post-processing of the table, which now contains reported values of Qp and glacier elevation changes.

dhdt.list.qp <- glofs.after2000.qp  %>%
  filter(!is.na(dh)) %>%
  mutate(dhdt = dh,
         dhdt_err_low  = dhdt - dh_err,
         dhdt_err_up   = dhdt + dh_err,
         dhdt_scale    = scale_this(dhdt),
         qp_scale = scale_this(Peak_discharge_Qp)) %>%
  mutate(Lake = str_sub(Glacier_and_lake, 16),
         Lake = str_replace_all(Lake, "_", " "),
         Glacier_and_lake = str_replace_all(Glacier_and_lake, " ", "_")) 

# Obtain the mean and error of cumulative glacier elevation change.

dhdt.scaling.qp <- dhdt.list.qp %>%
  dplyr::select(dhdt, dhdt_err_low, dhdt_err_up) %>%
  pivot_longer(cols = c(dhdt, dhdt_err_low, dhdt_err_up), 
               names_to = "var", values_to = "DHDT") %>%
  summarise(dhdt_mean = mean(DHDT),
            dhdt_sd = sd(DHDT))

# Use these values to scale the cumulative elevation change and error to mean of zero and
# unit standard deviation.

scaled.dhdt.qp.df <- dhdt.list.qp %>%
  mutate(mean_dhdt_scaled = (dhdt - dhdt.scaling.qp$dhdt_mean) / dhdt.scaling.qp$dhdt_sd,
         se_dhdt_scaled = ((dhdt_err_up - dhdt.scaling.qp$dhdt_mean) / dhdt.scaling.qp$dhdt_sd) - mean_dhdt_scaled)



bprior <- prior(student_t(3, 0, 5), class = "Intercept") +
          prior(student_t(3, 0, 2.5), class = "b") 

# Run the Bayesian hierarchical model of Qp versus cumulative glacier elevation change.

dhdt.qp.brm <- brm(qp_scale ~  me(mean_dhdt_scaled, se_dhdt_scaled) + (me(mean_dhdt_scaled, se_dhdt_scaled) | Glacier_and_lake),
                   data    = scaled.dhdt.qp.df,
                   family  = student(),
                   prior   = bprior,
                   chains  = 4,
                   warmup  = 2000,
                   iter    = 20000,
                   cores = 4,
                   save_pars = save_pars(latent = TRUE),
                   control = list(adapt_delta = 0.9,
                                  max_treedepth = 14))

# Write the hierarchical model of local changes in Qp versus glacier elevation change to disk.

saveRDS(dhdt.qp.brm , "local_Qp_vs_dhdt_model.RDS")
# dhdt.qp.brm  <- readRDS("local_Qp_vs_dhdt_model.RDS")

# Obtain posterior summary and assess convergence of chains.

summary(dhdt.qp.brm)
plot(dhdt.qp.brm)
pp_check(dhdt.qp.brm)

# Extract the range of cumulative elevation change for each glacier.

conds <- scaled.dhdt.qp.df %>% 
  group_by(Glacier_and_lake) %>% 
  summarise(min_r = min(dhdt_scale, na.rm = T) - 0.25, 
            max_r = max(dhdt_scale, na.rm = T) + 0.25)

# Obtain the posterior distribution of Qp with changes in glacier elevation.
# Keep the error in dhdt fixed.

preds <- scaled.dhdt.qp.df %>%
  data_grid(mean_dhdt_scaled = seq_range(dhdt_scale, n = 101, expand = 0.2),
            se_dhdt_scaled = 0.1,
            Glacier_and_lake  = unique(Glacier_and_lake),
            region = unique(region)) %>%
  add_epred_draws(object = dhdt.qp.brm, 
                  value = "qp", 
                  ndraws = 1000)

preds.sub <- list()

for (i in 1:nrow(conds)) {
  
  preds.sub[[i]] <- filter(preds, 
                           (Glacier_and_lake == conds$Glacier_and_lake[i]) & 
                             (mean_dhdt_scaled >= conds$min_r[i]) &  
                             (mean_dhdt_scaled <= conds$max_r[i]))
}

preds.sub <- bind_rows(preds.sub) %>% 
             ungroup()

# Plot the posterior trend in peak discharge Qp with glacier elevation change.

plot_trend_dhdt.qp <- preds.sub %>%
  mutate(dhdt = (mean_dhdt_scaled * dhdt.scaling.qp$dhdt_sd) + dhdt.scaling.qp$dhdt_mean,
         qp = (qp * sd(dhdt.list.qp$Peak_discharge_Qp)) + mean(dhdt.list.qp$Peak_discharge_Qp)) %>%
  ggplot(aes(x = dhdt, 
             y = qp)) +
  facet_wrap(~Glacier_and_lake, 
             scales = "free",
             ncol = 3) +
  scale_fill_manual(name = "Posterior rate", 
                    values = c("#52c8c8c8")) +
  stat_lineribbon(aes(y = qp), 
                  .width = 0.95,
                  point_interval = mean_qi) +
  geom_point(data = scaled.dhdt.qp.df, 
             aes(x = dhdt,
                 y = Peak_discharge_Qp,
                 color = rounded_year), 
             shape = 16,
             size = 1.5) +
  scale_color_viridis_c(name = "Year") +
  geom_linerange(data = scaled.dhdt.qp.df, 
                 aes(y = Peak_discharge_Qp,
                     x = dhdt,
                     xmin = dhdt_err_low,
                     xmax = dhdt_err_up, 
                     color = rounded_year)) +
  labs(x = expression(paste("Cumulative glacier elevation change 2000-2019 [m]")) ,
       y = expression(paste("Peak discharge ", Q[p], " [",  m^{3}, ']'))) +
  theme_bw() +
  theme( axis.text   = element_text(size = 8.5),
         axis.text.x = element_text(size = 8.5),
         axis.title  = element_text(size = 9),
         strip.text.x  = element_text(size = 6),
         legend.position = "bottom") 
 
# Obtain the posterior distribution of the regression slope, i.e. the trend of 
# peak discharge with the cumulative change in glacier elevation.

slopes.dhdt.qp <- dhdt.qp.brm %>%
  spread_draws(bsp_memean_dhdt_scaledse_dhdt_scaled, `r_Glacier_and_lake`[Glacier, param]) %>%
  filter(param == "memean_dhdt_scaledse_dhdt_scaled") %>%
  mutate(mean_qp_scaled = bsp_memean_dhdt_scaledse_dhdt_scaled + r_Glacier_and_lake,
         mean_qp = (mean_qp_scaled * sd(scaled.dhdt.qp.df$Peak_discharge_Qp)) / dhdt.scaling.qp$dhdt_sd ) %>%
  ungroup() %>%
  mutate(Region = str_sub(Glacier, 1, 8)) %>%
  rename(Glacier_and_lake = Glacier) %>%
  mutate(Glacier_and_lake = str_replace_all(Glacier_and_lake, "[.]", "_"),
         RGI = str_sub(Glacier_and_lake, 7, 14),
         Lake = str_sub(Glacier_and_lake, 16),
         Lake = str_replace_all(Lake, "[.]", " "),
         Lake = str_replace_all(Lake, "_", " "),
         type = "Peak discharge") 


################################################################################################
#######     Correlation between glacier elevation change and GLOF volume  V0    ################

# Select only dated ice-dammed GLOFs with reported V0 that happened between 2000 and 2019.

glofs.after2000.V0 <- all.glofs.V0 %>% 
  group_by(Glacier_and_lake) %>% 
  filter(nchar(Date) == 10) %>%
  filter((rounded_year >= 2000) & (rounded_year < 2020)) %>%
  filter(n() >= 5) %>%
  ungroup() %>%
  mutate(RGI_Glacier_Id = str_replace_all(RGI_Glacier_Id,  "_", "."),
         year_dec_glofs = decimal_date(ymd(Date)),
         dh = NA,
         dh_err = NA,
         time_dh = NA) 

# We now iterate over the tables containing the glacier elevation changes for a given RGI region.

for (i in 1: length(dhdt)) {
  
  # Read the cumulative elevation change for each glacier, aggregated per RGI region, to memory. 
  
  f <- read_csv(dhdt[i], 
                col_types = cols_only(rgiid = col_character(), 
                                      time = col_character(),
                                      dh = col_double(),
                                      err_dh = col_double())) 
  
  # Select all glaciers (using their RGI ID) that have reported values of V0.
  
  f2 <- f %>% 
    filter(rgiid %in% unique(glofs.after2000.V0$RGI_Glacier_Id)) %>%
    mutate(year_dec_dh = decimal_date(ymd(time)),
           year = as.numeric(str_sub(time, 1,4))) %>% 
    rename(RGI_Glacier_Id = rgiid,
           rounded_year = year) 
  
  # Iterate over table containing the dated GLOFs with reported values of V0.
  # Add the observed cumulative glacier elevation change (+ error) immediately before 
  # the outburst to the table.
  
  for(j in 1:nrow(glofs.after2000.V0)) {
    
    if(!(any(glofs.after2000.V0$RGI_Glacier_Id[j] == f2$RGI_Glacier_Id))) next
    
    f3 <- f2 %>% 
      filter(RGI_Glacier_Id == glofs.after2000.V0$RGI_Glacier_Id[j])
    
    idx <- max(which( (f3$year_dec_dh - glofs.after2000.V0$year_dec_glofs[j]) < 0))
    
    glofs.after2000.V0$dh[j] <- f3$dh[idx]
    glofs.after2000.V0$dh_err[j] <- f3$err_dh[idx]
    glofs.after2000.V0$time_dh[j] <- f3$time[idx]
    
  }
  
}

# Some post-processing of the table containing values of V0 and glacier elevation change.

dhdt.list.V0 <- glofs.after2000.V0  %>%
  filter(!is.na(dh)) %>%
  mutate(dhdt = dh,
         dhdt_err_low  = dhdt - dh_err,
         dhdt_err_up   = dhdt + dh_err,
         dhdt_scale    = scale_this(dhdt),
         V0_scale = scale_this(Mean_Flood_Volume_V0)) %>%
  mutate(Lake = str_sub(Glacier_and_lake, 16),
         Lake = str_replace_all(Lake, "_", " "),
         Glacier_and_lake = str_replace_all(Glacier_and_lake, " ", "_")) 

# Obtain the mean and error of cumulative glacier elevation change.

dhdt.scaling.V0 <- dhdt.list.V0 %>%
  dplyr::select(dhdt, dhdt_err_low, dhdt_err_up) %>%
  pivot_longer(cols = c(dhdt, dhdt_err_low, dhdt_err_up), 
               names_to = "var", values_to = "DHDT") %>%
  summarise(dhdt_mean = mean(DHDT),
            dhdt_sd = sd(DHDT))

# Use these values to scale the cumulative elevation change and error to mean of zero and
# unit standard deviation.

scaled.dhdt.V0.df <- dhdt.list.V0 %>%
  mutate(mean_dhdt_scaled = (dhdt - dhdt.scaling.V0$dhdt_mean) / dhdt.scaling.V0$dhdt_sd,
         se_dhdt_scaled = ((dhdt_err_up - dhdt.scaling.V0$dhdt_mean) / dhdt.scaling.V0$dhdt_sd) - mean_dhdt_scaled)

# Set weakly informed priors on the slope and intercept. Maintain all other default priors set in brms.

bprior <- prior(student_t(3, 0, 5), class = "Intercept") +
          prior(student_t(3, 0, 2.5), class = "b") 

# Run the Bayesian hierarchical model of Qp versus cumulative glacier elevation change.

dhdt.V0.brm <- brm(V0_scale ~  me(mean_dhdt_scaled, se_dhdt_scaled) + (me(mean_dhdt_scaled, se_dhdt_scaled) | Glacier_and_lake),
                   data    = scaled.dhdt.V0.df,
                   family  = student(),
                   prior   = bprior,
                   chains  = 4,
                   warmup  = 2000,
                   iter    = 20000,
                   cores = 4,
                   save_pars = save_pars(latent = TRUE),
                   control = list(adapt_delta = 0.9,
                                  max_treedepth = 14))


# Write the hierarchical model of local changes in V0 versus glacier elevation change to disk.

saveRDS(dhdt.V0.brm , "local_V0_vs_dhdt_model.RDS")
# dhdt.V0.brm  <- readRDS("local_V0_vs_dhdt_model.RDS")

# Obtain posterior summary and assess convergence of chains.

summary(dhdt.V0.brm)
plot(dhdt.V0.brm)
pp_check(dhdt.V0.brm)

# Extract the range of cumulative elevation change for each glacier.

conds <- scaled.dhdt.V0.df %>% 
  group_by(Glacier_and_lake) %>% 
  summarise(min_r = min(dhdt_scale, na.rm = T) - 0.25, 
            max_r = max(dhdt_scale, na.rm = T) + 0.25)

# Obtain the posterior distribution of V0 with changes in glacier elevation.
# Keep the error in dhdt fixed.

preds <- scaled.dhdt.V0.df %>%
  data_grid(mean_dhdt_scaled = seq_range(dhdt_scale, n = 101, expand = 0.2),
            se_dhdt_scaled = 0.1,
            Glacier_and_lake  = unique(Glacier_and_lake),
            region = unique(region)) %>%
  add_epred_draws(object = dhdt.V0.brm, 
                  value = "V0", 
                  ndraws = 1000)

# Truncate the posterior to the range of observed glacier elevation changes.

preds.sub <- list()

for (i in 1:nrow(conds)) {
  
  preds.sub[[i]] <- filter(preds, 
                           (Glacier_and_lake == conds$Glacier_and_lake[i]) & 
                             (mean_dhdt_scaled >= conds$min_r[i]) &  
                             (mean_dhdt_scaled <= conds$max_r[i]))
}

preds.sub <- bind_rows(preds.sub) %>% 
             ungroup()

# Plot the posterior trend in peak discharge Qp with glacier elevation change.

plot_trend_dhdt.V0 <- preds.sub %>%
  mutate(dhdt = (mean_dhdt_scaled * dhdt.scaling.V0$dhdt_sd) + dhdt.scaling.V0$dhdt_mean,
         V0 = (V0 * sd(dhdt.list.V0$Mean_Flood_Volume_V0)) + mean(dhdt.list.V0$Mean_Flood_Volume_V0)) %>%
  ggplot(aes(x = dhdt, 
             y = V0)) +
  facet_wrap(~Glacier_and_lake, 
             scales = "free",
             ncol = 4) +
  scale_fill_manual(name = "Posterior rate", 
                    values = c("#52c8c8c8")) +
  stat_lineribbon(aes(y = V0), 
                  .width = 0.95,
                  point_interval = mean_qi) +
  geom_point(data = scaled.dhdt.V0.df, 
             aes(x = dhdt,
                 y = Mean_Flood_Volume_V0,
                 color = rounded_year), 
             shape = 16,
             size = 1.5) +
  scale_color_viridis_c(name = "Year") +
  geom_linerange(data = scaled.dhdt.V0.df, 
                 aes(y = Mean_Flood_Volume_V0,
                     x = dhdt,
                     xmin = dhdt_err_low,
                     xmax = dhdt_err_up, 
                     color = rounded_year)) +
  labs(x = expression(paste("Cumulative glacier elevation change 2000-2019 [m]")) ,
       y = expression(paste("Flood volume ", V[0], " [", 10^{6}, " ", m^{3}, ']'))) +
  theme_bw() +
  theme( axis.text   = element_text(size = 6.5),
         axis.text.x = element_text(size = 8.5),
         axis.title  = element_text(size = 9),
         strip.text  = element_text(size = 8.5),
         legend.position = "bottom")  

# Obtain the posterior trend of V0 versus glacier elevation change.

slopes.dhdt.V0 <- dhdt.V0.brm %>%
  spread_draws(bsp_memean_dhdt_scaledse_dhdt_scaled, `r_Glacier_and_lake`[Glacier, param]) %>%
  filter(param == "memean_dhdt_scaledse_dhdt_scaled") %>%
  mutate(mean_V0_scaled = bsp_memean_dhdt_scaledse_dhdt_scaled + r_Glacier_and_lake,
         mean_V0 = (mean_V0_scaled * sd(scaled.dhdt.V0.df$Mean_Flood_Volume_V0)) / dhdt.scaling.V0$dhdt_sd ) %>%
  ungroup() %>%
  mutate(Region = str_sub(Glacier, 1, 8)) %>%
  rename(Glacier_and_lake = Glacier) %>%
  mutate(Glacier_and_lake = str_replace_all(Glacier_and_lake, "[.]", "_"),
         RGI = str_sub(Glacier_and_lake, 7, 14),
         Lake = str_sub(Glacier_and_lake, 16),
         Lake = str_replace_all(Lake, "[.]", " "),
         Lake = str_replace_all(Lake, "_", " "),
         type = "Flood volume") 


################################################################################################
#######     Plot the location of the investigated glacier dams and the local trends     ########
#######                    in V0 and Qp with glacier elevation change                   ########

# Read the polygon containing the outline of the study regions to memory.

sr <- st_read("Region_extents/Extent_pol.shp")

# Obtain an outline of the world continents from the NaturalEarth data repository.
# Reproject the line to a Robinson projection and convert to a polygon.

world <- ne_coastline(scale = "medium", returnclass = "sf") %>% 
  st_make_valid() %>%
  st_transform(., "+proj=robin") %>%
  st_cast("MULTIPOLYGON")

# Combine the locations of glacier dams with repeat GLOFs to one geometry (point features).

glofs.sf <- bind_rows(glofs.after2000.V0,  glofs.after2000.qp)  %>% 
  distinct(Glacier_and_lake, .keep_all = T) %>% 
  st_as_sf(
    coords = c("Longitude", "Latitude" ),
    crs = st_crs(4326)) %>%
  mutate(Lake = str_replace(Lake, " lake", ""),
         Lake = str_replace(Lake, " Lake", ""),
         Lake = str_replace(Lake, "Lake No", "No Lake"))

# Plot a gray world map in the background. Add the locations of the glacier dams on top.
# Add labels (using ggrepel) to each point feature.

world.map <- ggplot() +
  geom_sf(data = world,  
          fill = "gray80", 
          color = "gray80")  +
  theme_map()  +
  ggrepel::geom_text_repel(
    data = glofs.sf,
    aes(label = Lake, 
        geometry = geometry),
    stat = "sf_coordinates",
    max.overlaps = 21, 
    nudge_x = .15,
    size = 3,
    box.padding = 0.5,
    nudge_y = 1,
    segment.curvature = -0.1,
    segment.ncp = 3,
    segment.angle = 20,
    segment.size = 0.6,
    segment.color = "saddlebrown") +
  geom_sf(data = glofs.sf,
          mapping = aes(fill = region),
          color= "white", 
          size = 2,
          pch  = 21) +
  scale_fill_viridis_d(guide = guide_legend(reverse = TRUE, nrow = 2) ) +
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"),
        legend.position = "none")

# Extract the median and 95% HDI of the trends in V0 and Qp with glacier elevation change.
# To this end, combine both tibbles containing the posterior distributions of the trends.
# Change the labels of the RGI regions to a more intuitive description of the region. 

trends.V0.qp <- bind_rows(
  slopes.dhdt.qp %>% rename("val" = mean_qp),
  slopes.dhdt.V0 %>% rename("val" = mean_V0)) %>%
  group_by(Region, type, Lake) %>%
  filter((val >= quantile(val, 0.01)) & (val <= quantile(val, 0.99))) %>%
  ungroup() %>%
  mutate(Region = replace(Region, Region == "RGI60-13", "High Mountain Asia"),
         Region = replace(Region, Region == "RGI60-14", "High Mountain Asia"),
         Region = replace(Region, Region == "RGI60-15", "High Mountain Asia"),
         Region = replace(Region, Region == "RGI60-01", "NW North America"),
         Region = replace(Region, Region == "RGI60-02", "NW North America"),
         Region = replace(Region, Region == "RGI60-17", "Andes"),
         Region = replace(Region, Region == "RGI60-11", "European Alps"),
         Region = replace(Region, Region == "RGI60-08", "Scandinavia"),
         Region = replace(Region, Region == "RGI60-06", "Iceland")) %>%
  group_by(Region, type) %>%
  mutate(Lake = str_replace(Lake, " lake", ""),
         Lake = str_replace(Lake, " Lake", ""),
         Lake = str_replace(Lake, "Lake No", "No Lake"),
         Lake = str_replace(Lake, "Øvre Messingmalmvatn", "Øvre\nMessingm."),
         Lake = str_replace(Lake, "Lac de Faverges", "Lac de\nFaverges"),
         Lake = reorder(Lake,  val, median)) %>%
  ggplot(aes(x = val,
             y = Lake,
             order = Region)) +
  geom_vline(xintercept = 0)  + 
  stat_pointinterval( mapping = aes(color = Region),
                      .width = c(.95), 
                      point_size = 1.5,
                      stroke = 1.5) + 
  scale_color_viridis_d(guide = guide_legend(reverse = TRUE, nrow = 2) ) +
  theme_bw() +
  facet_wrap(~type, scales = "free_x") +
  labs(x = expression(paste("Change in ", V[0] , " [" , 10^{6}, " ", m^{3},  "]",
                            " and ", Q[p], " [" ,  m^{3}, " ", s^{-1}, "]")),
       y = "") +
  theme(axis.text   = element_text(size = 8),
        axis.text.x = element_text(size = 8),
        axis.title  = element_text(size = 9),
        strip.text  = element_text(size = 8),
        plot.margin = unit(c(0, 0.5, 0, 0), "cm"),
        legend.position = "none") + 
  guides(colour = guide_legend(nrow = 3))

# Combine the map and the posterior distributions of the regional trends.

map_and_trends <- world.map + 
  wrap_elements(panel = trends.V0.qp)+
  plot_layout(nrow = 2, byrow = FALSE) +
  plot_annotation(tag_levels = c('a', 'b')) 

# Save the plot to disk (Figure 4 a, b)

ggsave(plot = map_and_trends, 
       filename = "map_and_trends.pdf",
       width = 95,
       height = 150,
       unit = "mm")

################################################################################################
#######                  Obtain mean annual thinning rates of the dams                  ########

# Load the tables of local mean annual elevation change of glaciers impounding lakes between 2000 
# and 2019. The elevation changes are aggregated by RGI region such that one file may contain data
# on elevation changes from more than one glacier dam. The folder of glacier elevation changes  
# is called dh_pergla_cut. Please contant Romain Hugonnet (romain.hugonnet@legos.obs-mip.fr), 
# if you have any questions on data processing.

rates <- list.files(pattern = "rates.csv$",
                    path = "C:/Users/local-admin/Desktop/Lake_area_volume/dh_pergla_cut",
                    full.names = T)

# For each table, extract only the average calculated for the period from Jan 01, 2000, to
# Jan 01, 2020.

dhdt.pts <- list()

for (i in 1: length(rates)) {
  
  f1 <- read_csv(rates[i]) %>% filter(period == "2000-01-01_2020-01-01")
  
  f2 <- read_csv(dhdt[i]) %>%
    dplyr::select(c(rgiid, lat, lon)) %>%
    distinct()
  
  # Convert to a point feature using the columns 'lat' and 'lon', the centroid of the glacier.
  
  glacier.pts <- left_join(f1, f2, by = "rgiid")  %>% 
    sf::st_as_sf(coords = c("lon","lat")) %>% 
    sf::st_set_crs(4326)
  
  dhdt.pts[[i]] <- glacier.pts
  
}

# Combine all points to one geometry.

mean.gl.ch <- bind_rows(dhdt.pts)

# Write the point geometry of mean glacier elevation changes as an ESRI shapefile to disk.

st_write(mean.gl.ch, "dam_thinning_rats.shp")

################################################################################################
#######      Generate a plot of the cumulative elevation change for each glacier      ########## 
#######      that produced more than five GLOFs between 2000 and 2019.                ##########

# Obtain the RGI ID from all glaciers with reported values of Qp and V0.

gl.id <- unique(c(dhdt.list.qp$RGI_Glacier_Id, dhdt.list.V0$RGI_Glacier_Id))

# Extract the cumulative elevation change for these glaciers from the tables provided
# by Romain Hugonnet.

all.gl <- list()

for (i in 1: length(dhdt)) {
  
  f <- read_csv(dhdt[i], 
                col_types = cols_only(rgiid = col_character(), 
                                      time = col_character(),
                                      dh = col_double(),
                                      err_dh = col_double())) 
  
  f2 <- f %>% 
    filter(rgiid %in% unique(gl.id)) %>%
    mutate(year_dec_dh = decimal_date(ymd(time)),
           year = as.numeric(str_sub(time, 1,4))) %>% 
    rename(RGI_Glacier_Id = rgiid,
           rounded_year = year) 
  
  
  all.gl[[i]] <- f2
  
}

# Obtain the dates for all GLOFs, for which we have either reported values of Qp or V0.

uni.glofs <- bind_rows(dhdt.list.qp %>% select(RGI_Glacier_Id, year_dec_glofs, dh),
                       dhdt.list.V0 %>% select(RGI_Glacier_Id, year_dec_glofs, dh)) %>% 
  distinct()

# Generate the plot of mean glacier elevation, including a shade for the standard error.
# Add ticks at the bottom showing the dates when the GLOFs have happened at these glaciers.

elev.ch.plot <- bind_rows(all.gl) %>%
  ggplot(aes(x = year_dec_dh, 
             y = dh)) +
  geom_ribbon(mapping = aes(x = year_dec_dh, 
                            ymax = dh-err_dh , 
                            ymin = dh+err_dh),  fill = "skyblue") +
  geom_line() + 
  theme_bw() +
  geom_rug(data = uni.glofs,
           mapping = aes(x = year_dec_glofs),
           sides = "b") +
  facet_wrap(~RGI_Glacier_Id, scales = "free_y") + 
  theme(axis.text   = element_text(size = 8),
        axis.text.x = element_text(size = 8),
        axis.title  = element_text(size = 9),
        strip.text  = element_text(size = 8),
        plot.margin = unit(c(0, 0.5, 0, 0), "cm"),
        legend.position = "none") +
  xlab("Year") +
  ylab("Mean elevation change of the ice dam [m]")

# Write the plot of glacier elevation change to disk.

ggsave(plot = elev.ch.plot, 
       filename = "elev_change_per_glacier.pdf",
       width = 180,
       height = 180,
       unit = "mm")

ggsave(plot = elev.ch.plot, 
       filename = "elev_change_per_glacier.png",
       width = 180,
       height = 180,
       unit = "mm")

# What is the absolute elevation change of these glaciers between 2000 and 2019?

cumul.elev.ch <- bind_rows(all.gl) %>% 
  group_by(RGI_Glacier_Id) %>% 
  filter(row_number()==n()) %>% 
  ungroup() 

cumul.elev.ch %>% summarise(dh_median = median(dh),
                            dh_q025 = quantile(dh, 0.025),
                            dh_q975 = quantile(dh, 0.975))
