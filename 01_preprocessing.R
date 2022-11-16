################################################################################
#######     Preprocess data to calculate trends in GLOF volume,           ######
#######           peak discharge, timing, and elevation                   ######
#######                                                                   ######
#######                            by Georg Veh                           ######
#######                           03 March, 2022                          ######
#######                      comments added 15 Nov, 2022                  ######
################################################################################

# Load the following packages, or use install.packages("nameofpackage"), if some 
# of them are not pre-installed. In some cases you need to restart your R session.

require(tidyverse)
require(tidybayes)
require(modelr)
require(scales)
require(readODS)
require(patchwork)
require(ncdf4)
require(raster)
require(sf)
require(lubridate)

# Set YOUR working directory folder where to find all files, necessary to run 
# this script. Change the location appropriately.

setwd("D:/data/Lake_area_volume/")

# Useful functions

scale_this <- function(x){
  (x - mean(x, na.rm=TRUE)) / sd(x, na.rm=TRUE)
}

# Open-office spreadsheet with GLOFs per region in separate sheets.

glof.file <- "glofdatabase_2022_05_30.ods"

# Get names of the sheets in the Open Office document. 
# Exclude 'Global', 'Other', and 'Greenland'.

sheetnames <- list_ods_sheets(glof.file)

sheetnames <- sheetnames[!(sheetnames %in% "Global")]
sheetnames <- sheetnames[!(sheetnames %in% "Other")]
sheetnames <- sheetnames[!(sheetnames %in% "Greenland")]

region.list <- list()

# Load the table of regional GLOF reports into memory. 
# We iterate over the names of the spreadsheet.

for(r in sheetnames) {
  
  data <- as_tibble(read_ods(glof.file, sheet = r), .name_repair = "unique")
  
  data <- data[-c(1:2), 1:59] 
  data$region <- r 
  
  y <- as.numeric(str_sub(data$Date, 1, 4))
  
  # Extract the years of reported GLOF occurrences. 
  
  for (i in 1 : length(y)) {
    
    # Some years are NA because these GLOFs have no fixed date of occurrence, 
    # but a range of possible dates. For example, some GLOFs were detected from 
    # satellite images and offer only the last image before 
    # and the next image after the GLOF.
    # If there is NA, we first check whether there is a given range of dates. 
    # If so, we then randomly sample for the range of plausible years. 
    # Finally we increase the observed GLOF count for that dam type in that year 
    # by +1.
    
    if (is.na(y[i])) {
      
      min.date <- as.numeric(str_sub(data$Date_Min[i], 1, 4))
      max.date <- as.numeric(str_sub(data$Date_Max[i], 1, 4))
      
      if((!is.na(min.date)) & (!is.na(max.date))) {
        
        obs.period <- min.date:max.date 
        
        if (length(obs.period) == 1 ) {
          
          random.year <- obs.period
          
        } else {   random.year <- sample(obs.period, size = 1)}
        
        y[i] <- random.year
        
      } 
      
    } 
  }
  
  # Append the year column to the data.frame
  
  data <- data %>% 
    mutate(rounded_year = y) %>%
    mutate(
      Longitude = as.numeric(Longitude),
      Latitude  = as.numeric(Latitude),
      Mean_Lake_Volume_VL = as.numeric(Mean_Lake_Volume_VL),
      Min_VL    = as.numeric(Min_VL), 
      Max_VL    = as.numeric(Max_VL), 
      Mean_Flood_Volume_V0 = as.numeric(Mean_Flood_Volume_V0), 
      Min_V0    = as.numeric(Min_V0), 
      Max_V0    = as.numeric(Max_V0),
      Peak_discharge_Qp = as.numeric(Peak_discharge_Qp), 
      Min_Qp    = as.numeric(Min_Qp), 
      Max_Qp    = as.numeric(Max_Qp),
      Reference = as.character(Reference),
      D_buildings = as.character(D_buildings),
      reported_impacts = as.character(reported_impacts),
      economic_losses	 = as.character(economic_losses),
      D_buildings	     = as.character(D_buildings),
      D_bridges	       = as.character(D_bridges),
      D_roads_paths	   = as.character(D_roads_paths),
      D_railroads	  	 = as.character(D_railroads),
      D_utilities	     = as.character(D_utilities),
      D_flood_protection = as.character(D_flood_protection),
      D_environmental	 = as.character(D_environmental),
      resettlement	   = as.character(resettlement),
      reported_fatalities = as.character(reported_fatalities),
      Image_date_after    = as.character(Image_date_after),
      Image_date_before   = as.character(Image_date_before),
      Lake_area_before    = as.numeric(Lake_area_before), 
      Certainty_level_before = as.numeric(Certainty_level_before),
      Lake_area_after = as.numeric(Lake_area_after), 
      Certainty_level_after = as.numeric(Certainty_level_after)
    )
  
  region.list[[r]] <- data
  
}

# Concatenate all tibbles of regional GLOF occurrences to one long tibble

all.glofs <- bind_rows(region.list) %>%
  mutate(Lake_area_after = as.numeric(Lake_area_after)) %>%
  mutate(Lake_area_after = if_else(Lake_area_after == 0, 1, Lake_area_after)) %>%
  mutate(la_after_log = log10(Lake_area_after/10^6),
         la_before_log = log10(as.numeric(Lake_area_before)/10^6),
         region = str_replace(region, "Pacific NW", "NW North America"))

# Write this table to disk.

saveRDS(all.glofs, "all_glofs_tibble.RDS")
# all.glofs <- readRDS("all_glofs_tibble.RDS")

# Number of all ice dam failures between 1900 and 2021. 

all.glofs %>% 
  filter(Lake_type == "ice", 
         rounded_year >= 1900,
         rounded_year <= 2021) %>%
  summarise(n())

################################################################################
#########    FLOOD VOLUMES PER REGION   ########################################

# Select all ice-dammed with reported flood volume V0 between 1900 and 2021.

all.glofs.V0 <- all.glofs %>%
  filter(!is.na(RGI_Glacier_Id)) %>%
  filter(!is.na(Mean_Flood_Volume_V0)) %>%
  dplyr::select(RGI_Glacier_Id,
                rounded_year,
                Date,
                Mean_Flood_Volume_V0,
                Min_V0,
                Max_V0,
                Major_RGI_Region,
                Longitude,
                Latitude,
                Lake,
                region,
                Lake_type) %>%
  filter(rounded_year >= 1900) %>%
  group_by(RGI_Glacier_Id, Lake,
           .add = TRUE) %>%
  filter(Lake_type == "ice") %>%
  mutate(Glacier_and_lake = paste0(RGI_Glacier_Id, "_", Lake)) %>%
  ungroup() %>%
  mutate(V0_scale = scale_this(log10(Mean_Flood_Volume_V0)),
         year_scale = scale_this(rounded_year)) %>%
  mutate(RGI_Glacier_Id = str_replace_all(RGI_Glacier_Id, "[.]", "_"))

# Plot V0 per region

ggplot(data = all.glofs.V0,
       mapping = aes(x = year_scale,
                     y = V0_scale)) +
  geom_point() +
  facet_wrap(~region) +
  theme_bw() +
  xlab("Standardised year") + 
  ylab("log10-transformed and scaled flood volume V0")

# Write this table to disk.

saveRDS(all.glofs.V0, "all_glofs_V0_tibble.RDS")
# all.glofs.V0 <- readRDS("all_glofs_V0_tibble.RDS") 

# What is the median flood volume according to our database?

q.V0 <- quantile(all.glofs.V0$Mean_Flood_Volume_V0, c(0.025, 0.5, 0.975))
c(q.V0[2], q.V0[2]-q.V0[1], q.V0[3]- q.V0[2])

################################################################################
#########    PEAK DISCHARGES PER REGION   ######################################

# Select all ice-dammed with reported peak discharge Qp between 1900 and 2021.

all.glofs.qp <- all.glofs %>%
  filter(!is.na(RGI_Glacier_Id)) %>%
  filter(!is.na(Peak_discharge_Qp)) %>%
  dplyr::select(RGI_Glacier_Id,
                rounded_year,
                Date,
                Peak_discharge_Qp,
                Major_RGI_Region,
                Longitude,
                Latitude,
                Lake,
                region,
                Lake_type) %>%
  filter(rounded_year >= 1900) %>%
  group_by(RGI_Glacier_Id, Lake,
           .add = TRUE) %>%
  filter(Lake_type == "ice") %>%
  mutate(Glacier_and_lake = paste0(RGI_Glacier_Id, "_", Lake)) %>%
  ungroup() %>%
  mutate(qp_scale = scale_this(log10(Peak_discharge_Qp)),
         year_scale = scale_this(rounded_year)) %>%
  mutate(RGI_Glacier_Id = str_replace_all(RGI_Glacier_Id, "[.]", "_"))

# Plot Qp per region

ggplot(data = all.glofs.qp,
       mapping = aes(x = year_scale,
                     y = qp_scale)) +
  geom_point() +
  facet_wrap(~region) +
  theme_bw() +
  xlab("Standardised year") + 
  ylab("log10-transformed and scaled Peak discharge Qp")


# Write this table to disk.

saveRDS(all.glofs.qp, "all_glofs_qp_tibble.RDS")
# all.glofs.qp <- readRDS("all_glofs_qp_tibble.RDS") 

# What is the median peak discharge according to our database?

q.qp <- quantile(all.glofs.qp$Peak_discharge_Qp, c(0.025, 0.5, 0.975))
c(q.qp[2], q.qp[2]-q.qp[1], q.qp[3]- q.qp[2])

################## Generate a plot that summarizes GLOF reporting ##############

# Multi-panel histogram of reported values of all reported GLOFs,
# reported values of Qp, and reported values of V0 from ice-dammed lakes. 
# Data are aggregated in 30-year bins

a <- rbind(all.glofs %>%
  filter(Lake_type == "ice",
          rounded_year >= 1900) %>%
  transmute(rounded_year, region, type = "all"),
  all.glofs.qp %>% transmute(rounded_year, region, type = "Qp"),
  all.glofs.V0 %>% transmute(rounded_year, region, type = "V0")) %>%
  ggplot(mapping = aes(x = rounded_year, 
                       fill = type)) +
  geom_histogram(breaks   = c(1900, 1930, 1960, 1990, 2021),
                 position = "identity" ) +
  scale_fill_manual(values  = c("black", "blue", "darkorange"))+
  scale_x_continuous(breaks = c(1900, 1930, 1960, 1990, 2021))+
  stat_bin(breaks = c(1900, 1930, 1960, 1990, 2021), 
           aes(label = ..count..),
           size = 2, vjust= -0.5, geom = "text") +
  facet_grid(vars(type), vars(region)) +
  theme_bw() +
  xlab("Year (aggregated to 3 decades each)") + 
  ylab("Count") +
  ylim(c(0, 400)) +
  theme( axis.text   = element_text(size = 6),
         axis.text.x = element_text(size = 6),
         axis.text.y = element_text(size = 6),
         axis.title  = element_text(size = 6),
         strip.text  = element_text(size = 6),
         legend.position = "none")

ice.dammed <- all.glofs %>%
  filter(!is.na(RGI_Glacier_Id)) %>%
  filter(rounded_year >= 1900) %>%
  filter(Lake_type == "ice") %>%
  mutate(First_reference_found = as.numeric(First_reference_found)) 

# Year of GLOF occurrence versus year of first report.

b <- ice.dammed %>%
  ggplot(mapping = aes(y = First_reference_found,
                       x = rounded_year)) +
  geom_bin2d(binwidth = 5)  +
  scale_fill_viridis_b("Count in\n5-year bins",
                       # type = "viridis", 
                       option = "C",
                       breaks = c(1, 2, 5, 10, 25, 50)) +
  theme_bw() +
  geom_point(mapping = aes(y = First_reference_found,
                           x = rounded_year),
             data = ice.dammed,
             shape = 1, 
             size = 0.8,
             color = "white", 
             alpha = 0.1) +
  xlab("Year of GLOF occurrence") +
  ylab("Year of first publication / reporting") +
  theme( axis.text   = element_text(size = 6),
         axis.text.x = element_text(size = 6),
         axis.text.y = element_text(size = 6),
         axis.title  = element_text(size = 6),
         strip.text  = element_text(size = 6),
         legend.title = element_text(size = 6), 
         legend.text = element_text(size = 6),
         legend.box  = "horizontal", 
         legend.position = c(0.8, 0.3))

# Histogram of GLOF reporting.

c <- ice.dammed %>%
  ggplot(mapping = aes(x = First_reference_found)) +
  geom_histogram(breaks = seq(1900, 2020, by = 10)) +
  theme_bw() + 
  xlab("Year of first available report")  +
  ylab("Decadal count of reports")+
  theme( axis.text = element_text(size = 6),
         axis.text.x = element_text(size = 6),
         axis.text.y = element_text(size = 6),
         axis.title = element_text(size = 6),
         strip.text = element_text(size = 6))

# Combine all panels.

abc <- a / (b | c)  + 
  plot_annotation(tag_levels = 'a') & 
  theme(plot.tag = element_text(size = 8, face = "bold"))

# Export and save to disc.

ggsave(
  filename = "glof_reporting.pdf",
  plot = abc,
  width = 180,
  height = 180,
  units = "mm"
)

# Get statistics of GLOFs that had impacts.

imp.dam <- ice.dammed %>% 
  filter(!is.na(Impact_and_destruction))  %>%
  filter(Impact_and_destruction != "none",
         Impact_and_destruction != "none reported",
         Impact_and_destruction != "no damage",
         Impact_and_destruction != "No impacts reported") %>%
  transmute(RGI_Glacier_Id, Lake, rounded_year, Impact_and_destruction ) %>% 
  View()

############## HISTOGRAM OF GLOF COUNTS AND TEMPERATURE PER MONTH ##############

# Extract the date of each GLOF

doy.data <- all.glofs %>% 
  filter(!is.na(Date)) %>%
  filter(nchar(Date) == 10) %>%
  filter(Lake_type == "ice") %>%
  mutate(Glacier_and_lake = paste0(RGI_Glacier_Id, "_", Lake)) %>%
  ungroup() %>%
  mutate(doy = yday(Date),
         year = year(Date),
         month_num = month(Date),
         month_char = month(month_num, label = T, abbr = T)) %>%
  filter(year >= 1900) %>%
  mutate(year_scale = scale_this(year),
         doy_rescale = rescale(doy, to = c(-pi, pi), from = c(0,366)))

# Download the CRU time series from the following source and deposit it in 
# your working directory: 
# https://crudata.uea.ac.uk/cru/data/hrg/cru_ts_4.05/cruts.2103051243.v4.05/tmp/cru_ts4.05.1901.2020.tmp.dat.gz

cru.tmp <- stack("cru_ts4.05.1901.2020.tmp.dat.nc")

glof.doy.points <- doy.data %>% 
  filter(!is.na(Latitude)) %>%
  mutate(coord_paste = paste0(Latitude, Longitude)) %>%
  distinct(coord_paste, .keep_all = TRUE) %>%
  st_as_sf(coords = c( "Longitude", "Latitude"), crs = 4326)

# Extract temperatures from all grid cells in the stacked layers of temperature. 
# Each year has 12 layers (i.e. one per month) for the period 1901 to 2020.

temp.doy.ts <- raster::extract(cru.tmp, glof.doy.points) %>%
  as_tibble() %>%
  mutate(region = glof.doy.points$region) %>%
  pivot_longer( cols = starts_with("X"), names_to = "Month", values_to = "Temp") %>%
  mutate(month_date = as.Date(Month, format = "X%Y.%m.%d"),
         month_num = month(month_date, label = F, abbr = T),
         month_char = month(month_date, label = T, abbr = T))

# Generate the monthly mean annual air temperature

mean.temp <- temp.doy.ts %>%
  group_by(region, month_num) %>%
  summarise(Monthly_mean_temp = mean(Temp))

# Obtain the number of GLOFs reported in each month.

n.doy <- doy.data %>%
  group_by(region, month_num) %>%
  summarise(doy_count = n())

# Generate a histogram that both shows the number of reported GLOFs 
# and the mean air temperature in a given month

temp.doy.plot <- left_join(mean.temp, 
                           n.doy, 
                           by = c("region", "month_num")) %>% 
  pivot_longer(cols      = c("Monthly_mean_temp", "doy_count"), 
               names_to  = "class", 
               values_to = "val") %>%
  ggplot(mapping = aes(x = as_factor(month_num), 
                       y = val, 
                       fill = class)) +
  geom_bar(position = "dodge", stat = "identity") +
  facet_wrap(~region, scales = "free_y") + 
  scale_fill_manual(values = c("navy", "red"),
                    name = "Variable", 
                    labels = c("Number of GLOFs", 
                               "Mean monthly air temperature")) +
  labs(x = "Month",
       y = "Number of GLOFs / Mean monthly air temperature [°C]")  +
  theme_bw() + 
  theme( axis.text   = element_text(size = 7),
         axis.text.x = element_text(size = 7),
         axis.text.y = element_text(size = 7),
         axis.title  = element_text(size = 7),
         strip.text  = element_text(size = 7),
         legend.position = "bottom",
         legend.title = element_text(size = 7), 
         legend.text = element_text(size = 7))


ggsave(
  filename = "temp_doy_histogram.pdf",
  plot = temp.doy.plot,
  width = 180,
  height = 130,
  units = "mm"
)

ggsave(
  filename = "temp_doy_histogram.tiff",
  plot = temp.doy.plot,
  width = 180,
  height = 130,
  units = "mm",
  dpi = 300
)

################ NUMBER OF REPORTED GLOF PER LAKE (USED IN OVERVIEW FIGURE) ####

glof_n <- all.glofs %>%
  filter(rounded_year >= 1900) %>%
  group_by(RGI_Glacier_Id, Lake,
           .add = TRUE) %>%
  filter(Lake_type == "ice") %>%
  summarise(nGLOFs = n(),
            Lon = mean(Longitude),
            Lat = mean(Latitude)) %>%
  mutate(Lake = if_else(is.na(Lake), "unknown", Lake)) %>%
  ungroup() %>%
  drop_na() %>%
  st_as_sf(
    coords = c("Lon", "Lat"),
    crs = st_crs(4326)
  )

###### FIN! ####################################################################

