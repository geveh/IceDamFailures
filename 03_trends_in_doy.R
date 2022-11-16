################################################################################
#######     Estimate trends in annual timing of ice-dam failures on     ########
#######                   regional and local scale                      ########
#######                                                                 ########
#######                        by Georg Veh                             ########
#######                       03 March, 2022                            ########
#######                 comments added 15 Nov, 2022                     ########
################################################################################


# Load the following packages, or use install.packages("nameofpackage"), if some 
# of them are not pre-installed. In some cases you need to restart your R session.

require(brms)
require(tidybayes)
require(tidyverse)
require(lubridate)
require(scales)
require(ggpubr)
require(modelr)

# Set YOUR working directory folder where to find all files, necessary to 
# run this script. Change the location appropriately.

setwd("D:/data/Lake_area_volume/")

# Useful functions

scale_this <- function(x){
  (x - mean(x, na.rm=TRUE)) / sd(x, na.rm=TRUE)
}

################################################################################
#######     REGIONAL TRENDS IN DOY WITH TIME     ###############################

# Extract all dated GLOFs. Dates are given as YYYY-MM-DD, so the date vector  
# needs to have 10 characters. We only consider GLOFs that occurred in the
# period 1900-2021. We also rescale the dates to range between -Pi and + Pi to 
# match the parametrization of the von Mises  distribution in brms. 

# Load the tibble with all reported ice-dammed GLOFs between 1900 and 2021.

all.glofs <- readRDS("all_glofs_tibble.RDS")

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
         doy_rescale = rescale(doy, to = c(-pi, pi), from = c(0,366)),
         region = str_replace(region, "Pacific NW", "NW North America"))

# Specify normally distributed priors.

bprior <-  prior(normal(0, 2.5), class = "Intercept") +
  prior(normal(0, 0.75), class = "b") 

# Set inits for better sampling.

inits      <- list(Intercept = -2)
inits_list <- list(inits, inits, inits)

# This model is NOT hierarchical. We run it separately for each region.

doy.list.brm <- list()

uni.reg <- unique(doy.data$region)

# Iterate over the six study regions and fit a regression model of doy versus
# year with a von Mises likelihood. Add all model fits to a list.

for(reg in 1: length(uni.reg)) {
  
  trend.doy <- brm(doy_rescale ~ year_scale, 
                   data = doy.data %>% filter(region == uni.reg[reg]),
                   family = von_mises(),
                   prior = bprior,
                   warmup = 2000,
                   iter = 6000,
                   cores = 3,
                   chains = 3,
                   inits = inits_list,
                   control = list(adapt_delta = 0.92,
                                  max_treedepth = 15),
                   backend = "cmdstanr",
                   threads = threading(3)) 

  doy.list.brm[[reg]] <- trend.doy
  
}

# Save the regional models to disk.

saveRDS(doy.list.brm, "doy_trends_per_region.RDS")
# doy.list.brm <- readRDS("doy_trends_per_region.RDS")

# Obtain the posterior change in GLOF timing between 1900 and 2021.
# Note that both the predictor and the response were transformed (to
# a mean of zero and unit standard deviation, and -pi to pi, respectively), 
# so we need to re-transform the predictions to the original scale.

doy_post_trends.list <- list()

for(reg in 1: length(uni.reg)) {
  
  doy_post_trends <- doy.list.brm[[reg]] %>%
    add_epred_draws(newdata = data.frame(
      year_scale = (seq(1900, 2021, length.out = 51) - mean(doy.data$year)) / 
                    sd(doy.data$year)),
      value = "doy_rescale",
      ndraws = 1000) %>%
    mutate(doy = rescale(doy_rescale, to = c(0,366), from = c(-pi, pi)),
           year = year_scale * sd(doy.data$year) + mean(doy.data$year),
           region = uni.reg[reg])
  
  doy_post_trends.list[[reg]] <- doy_post_trends
}

# Obtain the number of dated GLOFs per region.

nobs <- doy.data %>%
  group_by(region) %>%
  summarise(n = n(),
            n_eq_n = paste0("n = ", n()))

# Plot the change in GLOF timing between 1900 and 2021.

doy_trend_plot <- bind_rows(doy_post_trends.list) %>%
  ggplot() +
  geom_point(data = doy.data,
             mapping = aes(x = year,
                           y = doy),
             alpha = 0.4,
             size = 0.8,
             color = "blue2") +
  stat_lineribbon(aes(x = year, 
                      y = doy), 
                  .width = c(.95)) + 
  scale_fill_manual("Posterior\ninterval", 
                    values =  c( "lightblue")) +
  scale_color_manual("Posterior\ninterval", 
                     values = c( "navy")) + 
  
  facet_wrap(~region) +
  labs(x = "Year" ,
       y =  "GLOF timing [doy]") +
  theme_bw()   +
  geom_text(data  = nobs, 
            aes(x = 1900, y = 10, label = n), 
            size = 2.3,
            colour = "gray10", 
            inherit.aes = FALSE, 
            parse = FALSE,
            hjust = "left") +
  xlim(c(1895, 2025))  +
  theme( axis.text = element_text(size = 7),
         axis.text.x = element_text(size = 7),
         axis.title = element_text(size = 7),
         strip.text = element_text(size = 7))

# Calculate the change in annual GLOF timing, that is the shift in days within  
# the year between 1900 and 2021 for all regions. Extract the regional posterior 
# distribution of  doy from 1900, and subtract it from the posterior
# distribution of doy in 2021.

doy_post_change.list <- list()

for(reg in 1: length(uni.reg)) {
  
  doy.1900 <- bind_rows(doy_post_trends.list) %>%
    filter(region == uni.reg[reg],
           (year == 1900) )
  
  doy.2021 <- bind_rows(doy_post_trends.list) %>%
    filter(region == uni.reg[reg],
           (year == 2021) )
  
  doy_change_region <- tibble(doy_change = doy.2021$doy - doy.1900$doy  ,
                              region =  uni.reg[reg])
  
  doy_post_change.list[[reg]] <- doy_change_region
  
}

# Bind the output per region into one tibble and obtain summary statistics
# as reported in the manuscript.

bind_rows(doy_post_change.list) %>% 
  group_by(region) %>% 
  summarise(round(quantile(doy_change, c(0.025, 0.5, 0.975))) )

# Plot the change in doy for each region.

doy_change_plot <- bind_rows(doy_post_change.list) %>%
  group_by(region) %>%
  filter(doy_change > quantile(doy_change, 0.005) & 
           doy_change < quantile(doy_change, 0.997)) %>%
  ungroup() %>%
  mutate(region = reorder(region, doy_change, median)) %>%
  ggplot(aes(x = doy_change,
             y = region,
             fill = stat(x > 0))) +
  stat_halfeye(.width = 0.95,
               slab_size = 1,
               interval_size = 2,
               interval_color = "black") + 
  scale_fill_manual("Posterior trend > 0", 
                    values = rev(RColorBrewer::brewer.pal(3, "Blues")[1:2])) +
  theme_bw() +
  labs(y = "Region",
       x = "Change in doy\nbetween 1900 and 2021") +
  geom_vline(xintercept = 0) +
  theme( axis.text = element_text(size = 7),
         axis.text.x = element_text(size = 7),
         axis.title = element_text(size = 7),
         strip.text = element_text(size = 7),
         legend.position = "none") 

# Combine both plots.

doy.arr <- ggarrange(plotlist = list(doy_trend_plot, 
                                     doy_change_plot), 
                     labels = c("a", "b"),
                     ncol = 2, 
                     font.label = list(size = 8),
                     widths = c(2,1),
                     align = "h", 
                     legend = "none",
                     common.legend = T) +
  theme(plot.margin = margin(0.2, 0.7, 0.2, 0.2, "cm")) 

# Write the combined plot to disk.

# Export figure 3.

ggsave(
  filename = "doy_change.pdf",
  plot = doy.arr,
  width = 180,
  height = 90,
  units = "mm"
)

ggsave(
  filename = "doy_change.png",
  plot = doy.arr,
  width = 180,
  height = 90,
  units = "mm"
)

# Show the regional changes in GLOF timing

bind_rows(doy_post_change.list) %>% 
  group_by(region) %>%
  filter(doy_change > quantile(doy_change, 0.005) & 
           doy_change < quantile(doy_change, 0.995)) %>%
  summarise(med_day = median(doy_change),
            q025 = quantile(doy_change, 0.025),
            q975 = quantile(doy_change, 0.975))


################################################################################
#######   LOCAL TRENDS IN DOY WITH TIME FOR LAKES WITH REPEAT OUTBURSTS   ######

# Extract all dated lakes with repeat (>5) dated GLOFs. 

doy.data.per.glacier <- all.glofs %>% 
  filter(!is.na(Date)) %>%
  filter(nchar(Date) == 10) %>%
  filter(Lake_type == "ice") %>%
  mutate(doy = yday(Date),
         year = year(Date)) %>%
  filter(year >= 1900) %>%
  mutate(Glacier_and_lake = paste0(RGI_Glacier_Id, "_", Lake)) %>%
  group_by(Glacier_and_lake) %>%
  filter(n() >= 5) %>%
  filter(length(unique(rounded_year)) > 2) %>%
  ungroup() %>%
  mutate(year_scale = scale_this(year),
         doy_rescale = rescale(doy, to = c(-pi, pi), from = c(0, 366)))

# Approach to fit Bayesian models of GLOF timing versus time for individual
# glacier lakes is the same as above for regions. Comments are largely omitted.

uni.gl <- unique(doy.data.per.glacier$Glacier_and_lake)

inits_list <- list(inits, inits)

# Set the `inits`
inits      <- list(Intercept = -2)
inits_list <- list(inits, inits, inits)

bprior <-  prior(normal(0, 2), class = "Intercept") +
  prior(normal(0, 0.75), class = "b") 

gl.doy.list <- list()

for(gl in 1:length(uni.gl)) {
  
  trend.doy <- brm(doy_rescale ~ year_scale,
                   data = doy.data.per.glacier %>% 
                      filter(Glacier_and_lake == uni.gl[gl]) ,
                   family = von_mises(),
                   prior = bprior,
                   warmup = 2000,
                   iter = 6000,
                   cores = 3,
                   chains = 3,
                   inits = inits_list,
                   control = list(adapt_delta = 0.95,
                                  max_treedepth = 15),
                   backend = "cmdstanr",
                   threads = threading(3)) 
  
  
  plot(conditional_effects(trend.doy), points = T)[[1]] +
    ggtitle(uni.gl[gl]) +
    xlab("DOY (scaled)") +
    ylab("Year (scaled)")
  
  gl.doy.list[[gl]] <- trend.doy
  
  
}

# Save the local models of doy versus time to disk.

saveRDS(object = gl.doy.list,  file = "doy_trends_per_glacier.RDS")
# gl.doy.list <- readRDS("doy_trends_per_glacier.RDS")

# Obtain the posterior distribution of changes in annual GLOF timing 
# for the period spanning the first and last observed GLOF from a given lake.

preds.list <- list()

for (i in 1: length(gl.doy.list)) {
  
  preds <- gl.doy.list[[i]]$data %>%
    data_grid(year_scale = seq_range(year_scale, n = 51),
              Glacier_and_lake  = uni.gl[i]) %>%
    add_epred_draws(object = gl.doy.list[[i]], 
                    value = "doy", 
                    ndraws = 1000) %>%
    ungroup() %>%
    dplyr::select(year_scale, Glacier_and_lake, doy)
  
  preds.list[[i]] <- preds
  
}

# Combine all posterior distributions to one large tibble and
# re-transform the data to the original scale.

preds <- bind_rows(preds.list) %>%
  mutate(rounded_year = (year_scale * sd(doy.data.per.glacier$rounded_year)) + 
                         mean(doy.data.per.glacier$rounded_year),
         doy_orig = rescale(doy, from = c(-pi, pi), to = c(0, 366)))

# Plot the posterior changes in annual GLOF timing for all glacier lakes that
# had at least 5 outbursts.

plot_trend_year.doy <- preds %>%
  mutate(Glacier_and_lake = str_replace(Glacier_and_lake, "_", " ")) %>%
  mutate(Glacier_and_lake = str_replace_all(Glacier_and_lake, 
                                           "RGI60-08.03416 Leirskardsvatnet/Kalvtjønna", 
                                           "RGI60-08.03416 Leirskardsvatnet"),
         Glacier_and_lake = str_replace_all(Glacier_and_lake, 
                                            "RGI60-06.00483 Lake at Breiðarmerkurfjall", 
                                            "RGI60-08.03416 NA")) %>%
  drop_na() %>%
  ggplot(aes(x = rounded_year, 
             y = doy_orig)) +
  facet_wrap(~Glacier_and_lake, 
             scales = "free",
             ncol = 5) +
  scale_fill_manual(name = "Posterior rate", 
                    values = c("#52c8c8c8")) +
  stat_lineribbon(aes(y = doy_orig), 
                  .width = 0.95,
                  point_interval = mean_qi) +
  geom_point(data = doy.data.per.glacier %>%
               mutate(Glacier_and_lake = str_replace(Glacier_and_lake, "_", " "))  %>%
               mutate(Glacier_and_lake = str_replace_all(Glacier_and_lake, 
                                                         "RGI60-08.03416 Leirskardsvatnet/Kalvtjønna", 
                                                         "RGI60-08.03416 Leirskardsvatnet"),
                      Glacier_and_lake = str_replace_all(Glacier_and_lake, 
                                                         "RGI60-06.00483 Lake at Breiðarmerkurfjall", 
                                                         "RGI60-08.03416 NA")), 
             aes(x = rounded_year,
                 y = doy),
             shape = 16 ,
             color = "blue3", 
             size = 1,
             alpha = 0.5) +
  labs(x = "Year",
       y = "Date of year of outburst [doy]") +
  theme_bw() +
  theme( axis.text = element_text(size = 6),
         axis.text.x = element_text(size = 6),
         axis.title = element_text(size = 6),
         strip.text = element_text(size = 5),
         legend.position = "none") 


# Obtain the posterior 95% highest density interval (HDI) for the posterior 
# change in GLOF timing for all glaciers between the first and the last 
# reported GLOF from a given lake. 

# First, obtain the GLOF timing in the first reported GLOF from a given lake.

preds.min.y <- preds %>% 
  group_by(Glacier_and_lake) %>%
  filter( rounded_year == min(rounded_year)) %>%
  rename(doy_early = doy_orig,
         year_early = rounded_year)

# Then, obtain the posterior distribution of doy in the year of the last 
# reported GLOF from a given lake.

preds.max.y <- preds %>% 
  group_by(Glacier_and_lake) %>%
  filter( rounded_year == max(rounded_year)) %>%
  rename(doy_late = doy_orig,
         year_late = rounded_year)

# Calculate the difference between the two distributions, i.e. the shift in days 
# of GLOF occurrence in the observation period.

change.doy <- bind_cols(preds.min.y, preds.max.y ) %>%
  drop_na() %>%
  mutate(diff_doy = doy_late - doy_early) %>%
  dplyr::select(Glacier_and_lake...2, diff_doy) %>%
  rename(Glacier_and_lake = Glacier_and_lake...2) %>%
  mutate(Glacier_and_lake = str_replace(Glacier_and_lake, "_", " ")) %>%
  mutate(Glacier_and_lake = str_replace_all(Glacier_and_lake, 
                                             "RGI60-08.03416 Leirskardsvatnet/Kalvtjønna", 
                                                     "RGI60-08.03416 Leirskardsvatnet"),
          Glacier_and_lake = str_replace_all(Glacier_and_lake, 
                                             "RGI60-06.00483 Lake at Breiðarmerkurfjall", 
                                              "RGI60-08.03416 NA")) %>%
  mutate(Glacier_and_lake = reorder(Glacier_and_lake, diff_doy, median)) %>%
  ggplot(aes(x = diff_doy,
             y = Glacier_and_lake)) +
  stat_pointinterval( .width = c( .95), 
                      color = "navy",
                      interval_size = 0.9,
                      point_size = 1.5) +
  theme_bw() +
  labs(x = "Change in doy between\nfirst and last reported GLOF from a given lake",
       y = "Glacier and lake") +
  geom_vline(xintercept = 0) +
  theme( axis.text   = element_text(size = 7),
         axis.text.x = element_text(size = 7),
         axis.title  = element_text(size = 7),
         strip.text  = element_text(size = 7),
         legend.position = "none") 

# Show the statistics from the posterior distribution of local
# changes in GLOF timing.

bind_cols(preds.min.y, preds.max.y ) %>%
  mutate(diff_doy = doy_late - doy_early) %>%
  dplyr::select(Glacier_and_lake...2, diff_doy) %>% 
  rename(Glacier_and_lake = Glacier_and_lake...2) %>% 
  group_by(Glacier_and_lake) %>% 
  summarise(q025 = quantile(diff_doy, 0.025), 
            q50 = quantile(diff_doy, 0.5), 
            q0975 = quantile(diff_doy, 0.975)) %>%
  View()

# Write the plots of local changes in doy to disk. 

# For the supplementary information.

ggsave(filename = "doy_local.pdf",
       plot = plot_trend_year.doy,
       width = 183,
       height = 210,
       units = "mm")

ggsave(filename = "doy_local.png",
       plot = plot_trend_year.doy,
       width = 183,
       height = 210,
       units = "mm")

# Extended Data Figure 7

ggsave(filename = "post_trend_doy_per_lake.pdf",
       plot = change.doy,
       width = 180,
       height = 200,
       units = "mm")

ggsave(filename = "post_trend_doy_per_lake.png",
       plot = change.doy,
       width = 180,
       height = 200,
       units = "mm")
