####################################################################################################
#######        Quantile regression models of peak discharges and volumes               #############
#######                         from ice-dam failures                                  #############
#######                                  by Georg Veh                                  #############
#######                                 03 March, 2022                                 #############
####################################################################################################

# Load the following packages, or use install.packages("nameofpackage"), if some of them
# are not pre-installed. In some cases you need to restart your R session.

require(readODS)
require(brms)
require(tidybayes)
require(scales)
require(ggpubr)
require(tidyverse)

# Set YOUR working directory folder where to find all files, necessary to run this script.
# Change the location appropriately.

setwd("C:/Users/local-admin/Desktop/Lake_area_volume/")

# Useful functions

scale_this <- function(x){
  (x - mean(x, na.rm=TRUE)) / sd(x, na.rm=TRUE)
}

# Open-office spreadsheet with GLOFs per region in separate sheets.

glof.file <- "Global_GLOF_database_2021_12_08.ods"

# Get names of the sheets in the Open Office document. Exclude 'Global' or 'Other'.

sheetnames <- list_ods_sheets(glof.file)

sheetnames <- sheetnames[!(sheetnames %in% "Global")]
sheetnames <- sheetnames[!(sheetnames %in% "Other")]

region.list <- list()

# Load the table of regional GLOF reports into memory. 
# We iterate over the names of the spreadsheet.

for(r in sheetnames) {
  
  data <- as_tibble(read_ods(glof.file, sheet = r), .name_repair = "unique")
  
  data <- data[-c(1:2), 1:44] 
  data$region <- r 
  
  y <- as.numeric(str_sub(data$Date, 1, 4))
  
  # Extract the years of reported GLOF occurrences. 
  
  for (i in 1 : length(y)) {
    
    # Some years are NA because these GLOFs have no fixed date of occurrence, but a range of possible dates.
    # For example, some GLOFs were detected from satellite images and offer only the last image before 
    # and the next image after the GLOF.
    # If there is NA, we first check whether there is a given range of dates. 
    # If so, we then randomly sample for the range of plausible years. 
    # Finally we increase the observed GLOF count for that dam type in that year by +1.
    
    if (is.na(y[i])) {
      
      min.date <- as.numeric(str_sub(data$Date_Min[i], 1, 4))
      max.date <- as.numeric(str_sub(data$Date_Max[i], 1, 4))
      
      if((!is.na(min.date)) & (!is.na(max.date))) {
        
        obs.period <- min.date:max.date
        
        random.year <- sample(obs.period, size = 1)
        
        y[i] <- random.year
        
      } 
      
    } 
  }
  
  # Append the year column to the data.frame
  
  data <- data %>% 
    mutate(rounded_year = y) %>%
    mutate(
      Longitude = as.numeric(Longitude),
      Latitude = as.numeric(Latitude),
      Mean_Lake_Volume_VL = as.numeric(Mean_Lake_Volume_VL),
      Min_VL = as.numeric(Min_VL), 
      Max_VL = as.numeric(Max_VL), 
      Mean_Flood_Volume_V0 = as.numeric(Mean_Flood_Volume_V0), 
      Min_V0 = as.numeric(Min_V0), 
      Max_V0 = as.numeric(Max_V0),
      Peak_discharge_Qp = as.numeric(Peak_discharge_Qp), 
      Min_Qp = as.numeric(Min_Qp), 
      Max_Qp = as.numeric(Max_Qp),
      Reference = as.character(Reference),
      image_date_after = as.character(image_date_after),
      image_date_before = as.character(image_date_before),
      Lake_area_before  = as.numeric(Lake_area_before), 
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
         la_before_log = log10(as.numeric(Lake_area_before)/10^6))

# Write this table to disk.

saveRDS(all.glofs, "all_glofs_tibble.RDS")
# all.glofs <- readRDS("all_glofs_tibble.RDS")

# Number of all ice dam failures between 1900 and 2021. 

all.glofs %>% 
  filter(Lake_type == "ice", 
         rounded_year >= 1900) %>%
  summarise(n())

################################################################################################
#######     QUANTILE REGRESSION for PEAK DISCHARGE Qp     ######################################

#### Select values of Qp per region

all.glofs.qp.quant <- all.glofs %>%
  filter(!is.na(RGI_Glacier_Id)) %>%
  filter(!is.na(Peak_discharge_Qp)) %>%
  dplyr::select(RGI_Glacier_Id,
                rounded_year,
                Peak_discharge_Qp,
                Min_Qp,
                Max_Qp,
                Major_RGI_Region,
                Longitude,
                Latitude,
                Lake,
                region,
                Lake_type)  %>%
  filter(Lake_type == "ice",
         rounded_year >= 1900) %>%
  mutate(region = str_replace(region, "Pacific NW", "NW North America"))

# We calculate quantile regressions of QP versus time for the 50th and 90th percentile,
# for all ice-dammed lakes in a given region.

quants <- c(0.5, 0.9)
regions <- unique(all.glofs.qp.quant$region)

# We iterate over each quantile and compute the quantile regression models. 
# Each entry contains three models (all, ice, and moraine-dammed lakes). 

qp.region.fits.list <- list()
qp.region.epreds.list <- list()
qp.region.trends.list <- list()
qp.region.med.line.list  <- list()

for (i in 1 :length(regions)) {
  
  # (Could set a clock to monitor the progress of the loop.)
  
  start <- Sys.time()
  
  selected.region <- all.glofs.qp.quant %>% 
    filter(region == regions[i]) %>%
    mutate(qp_log = log10(Peak_discharge_Qp),
           qp_scale = scale_this(qp_log),
           year_scale = scale_this(rounded_year))
  
  # Recover the original parameters of year and QP to plot the data on original scale.
  
  sd_qp     <- sd(selected.region$qp_log)
  mean_qp   <- mean(selected.region$qp_log) 
  sd_year   <- sd(selected.region$rounded_year)
  mean_year <- mean(selected.region$rounded_year)
  
  # Define the range, for which new predictions will be made.
  # ... on original scale
  
  seq.years.orig <- 1895:2025
  
  # ... and on scaled range
  
  seq.years.scaled <- (seq.years.orig - mean_year) / sd_year
  
  # Set weakly informative priors on the regression intercept, slope, and the noise.
  
  bprior <- c(prior(prior = "student_t(3, 0, 5)", class = "Intercept"),
              prior(prior = "student_t(3, 0, 2.5)", class = "b"),
              prior(normal(0, 2.5), class = "sigma"))
  
  # Run a loop over both the 50th and 90th percentile. We generate lists, which we fill 
  # with the fitted models, the posterior distribtutions, the posterior predictive
  # distributions, and the median prediction.

  qp.region.quant.fits  <- list()
  qp.region.post.samps  <- list()
  qp.region.post.epreds <- list()
  qp.region.post.medline <- list()
  
  for (j in 1: length(quants)) {
    
    if (j == 1) {
      
      # Note that there will be warnings that R-hat is NA, and ESS is too low
      # This is because of the quantile being a fixed parameter, for which we do not set 
      # a prior and, hence fail to obtain a posterior distribution.
      # Note that we use 4 cores - decrease or increase, if necessary.
      
      fit.quant <- brm(bf(qp_scale ~  year_scale, 
                          quantile = quants[j]),
                       data    = selected.region,
                       cores   = 4,
                       prior   = bprior,
                       chains  = 4,
                       warmup  = 2000,
                       iter    = 6000,
                       family  = asym_laplace(),
                       control = list(adapt_delta = 0.95,
                                      max_treedepth = 14)) 
      
    } else {
      
      fit.quant <- update(fit.quant, 
                          cores   = 4,
                          chains  = 4,
                          formula = bf(qp_scale ~  year_scale,
                                       quantile = quants[j]))
    }
    
    
    # Obtain the standardized posterior of the regression slope
    # and convert to original scale.
    
    post_samps <- as_draws_df(fit.quant, 
                              variable = "b_year_scale") %>% 
      mutate(region = regions[i]) %>%
      mutate(b_orig = b_year_scale * (sd_qp / sd_year),
             quantile = quants[j])
    
    # Obtain the standardized posterior predictive distribution for new observations
    # and convert the predictions to original scale.
    
    post_epred <- posterior_epred(fit.quant,
                                  dpar = "mu",
                                  ndraws = 1500,
                                  newdata = data.frame(
                                            year_scale = seq.years.scaled)) %>%
      t() %>%
      as_tibble( .name_repair = "unique") %>%
      mutate(region = regions[i],
             quantile = quants[j],
             year = seq.years.orig) %>%
      pivot_longer( cols = starts_with("..."), 
                    values_to = "pred") %>%
      mutate(pred = (pred*sd_qp) + mean_qp, 
             pred = 10^pred) %>%
      dplyr::select(-name)
    
    # Obtain the median posterior intercept and slope of the regression model.
    
    med_int_b <- as_draws_df(fit.quant, 
                             variable = c("b_Intercept", "b_year_scale")) %>% 
      summarise(m_int = median(b_Intercept),
                m_b = median(b_year_scale)) 
    
    # Use these two parameteres to obtain the median prediction of Qp with time.
    
    med_line <- tibble(year = seq.years.orig, 
                       med_pred = med_int_b$m_int + med_int_b$m_b * seq.years.scaled,
                       pred = (med_pred * sd_qp) + mean_qp) %>%
      mutate(pred = 10^pred) %>%
      mutate(quantile = quants[j],
             region = regions[i])
    
    # Write each output to a separate list item.
    
    qp.region.quant.fits[[j]]  <- fit.quant
    qp.region.post.samps[[j]]  <- post_samps
    qp.region.post.epreds[[j]] <- post_epred
    qp.region.post.medline[[j]] <- med_line
    
  }
  
  # Add name for the quantile to the list element.
  
  names(qp.region.quant.fits) <- paste0("Q_", quants)
  
  # Finally, write the regression models and predictions for both quantiles to the top list element.
  
  qp.region.fits.list[[i]] <- qp.region.quant.fits
  qp.region.trends.list[[i]] <- bind_rows(qp.region.post.samps)
  qp.region.epreds.list[[i]] <- bind_rows(qp.region.post.epreds) 
  qp.region.med.line.list[[i]]  <- bind_rows(qp.region.post.medline)
  
  # (Print message of runtime length).
  
  stopp <- Sys.time()
  message(stopp-start)
  
}

# Save these objects to disk.

saveRDS(object = qp.region.fits.list,     file = "qp_regional_quantreg_fits.RDS")
saveRDS(object = qp.region.epreds.list,   file = "qp_regional_quantreg_posteriors.RDS")
saveRDS(object = qp.region.trends.list,   file = "qp_regional_quantreg_trends.RDS")
saveRDS(object = qp.region.med.line.list, file = "qp_regional_medline.RDS")
# qp.region.fits.list <-    readRDS("qp_regional_quantreg_fits.RDS")
# qp.region.epreds.list <-   readRDS("qp_regional_quantreg_posteriors.RDS")
# qp.region.trends.list <-   readRDS( "qp_regional_quantreg_trends.RDS")
# qp.region.med.line.list <- readRDS("qp_regional_medline.RDS")

all_regional_preds <- bind_rows(qp.region.epreds.list)

# Obtain the number of observed values of Qp for each region.

nobs <- all.glofs.qp.quant %>%
  group_by(region) %>%
  summarise(n = paste0("n = ", n()))

# Reproduce Figure 1a: Quantile regression (50th and 90th percentile) of Peak Discharge
# versus time for each region

quants.05.09.qp <- all_regional_preds %>%
  mutate(region = str_replace(region, "Pacific NW", "NW North America")) %>%
  ggplot() +
  stat_lineribbon(aes(x = year, 
                      y = pred, 
                      color = ordered(quantile),
                      fill = ordered(quantile)), 
                  .width = c(.95), 
                  alpha = 0.5) + 
  scale_fill_manual("Posterior\ninterval", 
                    values =  c( "navy", "darkorange")) +
  scale_color_manual("Posterior\ninterval", 
                     values = c( "navy", "darkorange")) + 
  geom_point(data = all.glofs.qp.quant,
             mapping = aes(x = rounded_year,
                           y = Peak_discharge_Qp),
             alpha = 0.75) +
  facet_wrap(~region, scales = "free_y") +
  labs(x = "Year" ,
       y = expression(paste("Peak discharge ", Q[p], " [",  m^{3}, ' ', s^{-1}, ']'))) +
  theme_bw()   +
  scale_y_continuous(trans  = log10_trans(),
                     breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x)))  + 
  geom_text(data  = nobs, 
            aes(x = 1900, y = 100, label = n), 
            colour = "gray10", 
            inherit.aes = FALSE, 
            parse = FALSE,
            hjust = "left")

# Reproduce Figure 1b: Posterior trends of the 50th and 90th percentile of Peak Discharge
# versus time for each region

trends.quants.qp <- bind_rows(qp.region.trends.list) %>%
  mutate(region = str_replace_all(region, "[.]", "\n"),
         region = str_replace(region, "Pacific NW", "NW North America")) %>%
  mutate(region = as_factor(region) %>% fct_rev(),
         quantile = as_factor(quantile) %>% fct_rev()) %>%
  ggplot(aes(x = b_orig,
             y = region,
             color = quantile)) +
  coord_cartesian(xlim = c(-0.032, 0.032)) + 
  stat_pointinterval( .width = c(.95), 
                      size = 1,
                      point_size = 1.25,
                      position = position_dodge(width = -0.4)) + 
  scale_color_manual(values = c("darkorange", "navy")) +
  theme_bw() +
  labs(y = "",
       x = expression(paste("Trend in ", log[10], " ", Q[p], " [",  m^{3}, ' ', s^{-1}, " ", yr^{-1}, ']'))) +
  geom_vline(xintercept = 0) +
  theme( axis.text = element_text(size = 8.5),
         axis.text.x = element_text(size = 8.5),
         axis.title = element_text(size = 9, face = "bold"),
         strip.text = element_text(size = 8.5))

# Combine both plots (change over time and posterior trend) into one plot.

arr.quant.qp <- ggarrange(plotlist = list(quants.05.09.qp, trends.quants.qp), 
                          labels = "auto",
                          ncol = 2, 
                          align = "h", 
                          common.legend = T, 
                          widths = c(2, 1))

# Save plots to disk.

ggsave(
  filename = "qp_quants05_09.pdf",
  plot = arr.quant.qp,
  width = 250,
  height = 170,
  units = "mm"
)

ggsave(
  filename = "qp_quants05_09.png",
  plot = arr.quant.qp,
  width = 250,
  height = 170,
  units = "mm"
)

ggsave(
  filename = "qp_quants05_09.pdf",
  plot = quants.05.09.qp,
  width = 180,
  height = 100,
  units = "mm"
)


gc()
closeAllConnections()


################################################################################################
#######     QUANTILE REGRESSION for FLOOD VOLUME V0     ########################################

#### Select values of Qp per region

all.glofs.v0.quant <- all.glofs %>%
  filter(!is.na(RGI_Glacier_Id)) %>%
  filter(!is.na(Mean_Flood_Volume_V0)) %>%
  dplyr::select(RGI_Glacier_Id,
                rounded_year,
                Mean_Flood_Volume_V0,
                Min_V0,
                Max_V0,
                Major_RGI_Region,
                Longitude,
                Latitude,
                Lake,
                region,
                Lake_type)  %>%
  filter(Lake_type == "ice",
         rounded_year >= 1900) %>%
  mutate(region = str_replace(region, "Pacific NW", "NW North America"))


# We calculate quantile regressions of V0 versus time for the 50th and 90th percentile in the 
# six study regions.

quants <- c(0.5, 0.9)
regions <- unique(all.glofs.v0.quant$region)

# We iterate over the two percentile and compute the quantile regression models. 

v0.region.fits.list <- list()
v0.region.epreds.list <- list()
v0.region.trends.list <- list()
v0.region.med.line.list  <- list()

for (i in 1 :length(regions)) {
  
  # Set a clock to monitor the progress of the loop.
  
  start <- Sys.time()
  
  selected.region <- all.glofs.v0.quant %>% 
    filter(region == regions[i]) %>%
    mutate(v0_log = log10(Mean_Flood_Volume_V0),
           v0_scale = scale_this(v0_log),
           year_scale = scale_this(rounded_year))
  
  # Recover the original parameters of year and V0 to plot the data on original scale.
  
  sd_v0     <- sd(selected.region$v0_log)
  mean_v0   <- mean(selected.region$v0_log) 
  sd_year   <- sd(selected.region$rounded_year)
  mean_year <- mean(selected.region$rounded_year)
  
  # Define the range, for which the data will be plotted.
  # ... on original scale
  
  seq.years.orig <- 1895:2025
  
  # ... and on scaled range
  
  seq.years.scaled <- (seq.years.orig - mean_year) / sd_year
  
  # Set weakly informative priors on the regression intercept, slope, and the noise.
  
  bprior <- c(prior(prior = "student_t(3, 0, 5)", class = "Intercept"),
              prior(prior = "student_t(3, 0, 2.5)", class = "b"),
              prior(normal(0, 2.5), class = "sigma"))
  
  # Model for all lake types.

  v0.region.quant.fits  <- list()
  v0.region.post.samps  <- list()
  v0.region.post.epreds <- list()
  v0.region.post.medline <- list()
  
  for (j in 1: length(quants)) {
    
    if (j == 1) {
      
      # Note that there brms will throw warnings because of the quantile parameter
      # being fixed and thus missing posterior draws.
      
      fit.quant <- brm(bf(v0_scale ~  year_scale, 
                          quantile = quants[j]),
                       data    = selected.region,
                       cores   = 4,
                       prior   = bprior,
                       chains  = 4,
                       warmup  = 2000,
                       iter    = 6000,
                       family  = asym_laplace(),
                       control = list(adapt_delta = 0.95,
                                      max_treedepth = 14)) 
      
    } else {
      
      fit.quant <- update(fit.quant, 
                          cores   = 4,
                          chains  = 4,
                          formula = bf(v0_scale ~  year_scale,
                                       quantile = quants[j]))
    }
    
    # Obtain the standardized posterior of the regression intercept and slope, and 
    # transform to original scale.
    
    post_samps <- as_draws_df(fit.quant, 
                              variable = "b_year_scale") %>% 
      mutate(region = regions[i]) %>%
      mutate(b_orig = b_year_scale * (sd_v0 / sd_year),
             quantile = quants[j])
    
    # Generate predictions of V0 for new data and transform to original scale. 
    
    post_epred <- posterior_epred(fit.quant,
                                  dpar = "mu",
                                  ndraws = 1500,
                                  newdata = data.frame(year_scale = seq.years.scaled)) %>%
      t() %>%
      as_tibble(.name_repair = "unique") %>%
      mutate(region = regions[i],
             quantile = quants[j],
             year = seq.years.orig) %>%
      pivot_longer( cols = starts_with("..."), 
                    values_to = "pred") %>%
      mutate(pred = (pred*sd_v0) + mean_v0, 
             pred = 10^pred) %>%
      dplyr::select(-name)
    
    # Obtain the median prediction of the trend in V0 over time.
    
    med_int_b <- as_draws_df(fit.quant, 
                             variable = c("b_Intercept", "b_year_scale")) %>% 
      summarise(m_int = median(b_Intercept),
                m_b = median(b_year_scale)) 
    
    med_line <- tibble(year = seq.years.orig, 
                       med_pred = med_int_b$m_int + med_int_b$m_b * seq.years.scaled,
                       pred = (med_pred * sd_v0) + mean_v0) %>%
      mutate(pred = 10^pred) %>%
      mutate(quantile = quants[j],
             region = regions[i])
    
    # Write the model outputs to a list with two elements, one for each percentile.
    
    v0.region.quant.fits[[j]]   <- fit.quant
    v0.region.post.samps[[j]]   <- post_samps
    v0.region.post.epreds[[j]]  <- post_epred
    v0.region.post.medline[[j]] <- med_line
    
  }
  
  # Add names to the fits.
  
  names(v0.region.quant.fits) <- paste0("Q_", quants)
  
  v0.region.fits.list[[i]]     <- v0.region.quant.fits
  v0.region.trends.list[[i]]   <- bind_rows(v0.region.post.samps)
  v0.region.epreds.list[[i]]   <- bind_rows(v0.region.post.epreds) 
  v0.region.med.line.list[[i]] <- bind_rows(v0.region.post.medline)
  
  # (Print message of runtime length).
  
  stopp <- Sys.time()
  message(stopp-start)
  
}

# Write all objects to disk.

saveRDS(object = v0.region.fits.list,     file = "v0_regional_quantreg_fits.RDS")
saveRDS(object = v0.region.epreds.list,   file = "v0_regional_quantreg_posteriors.RDS")
saveRDS(object = v0.region.trends.list,   file = "v0_regional_quantreg_trends.RDS")
saveRDS(object = v0.region.med.line.list, file = "v0_regional_medline.RDS")

# v0.region.fits.list <- readRDS( file = "v0_regional_quantreg_fits.RDS")
# v0.region.epreds.list <- readRDS(  file = "v0_regional_quantreg_posteriors.RDS")
# v0.region.trends.list <- readRDS(file = "v0_regional_quantreg_trends.RDS")
# v0.region.med.line.list <- readRDS(file = "v0_regional_medline.RDS")

all_regional_preds <- bind_rows(v0.region.epreds.list)

# Obtain the number of observed values of V0 in each region.

nobs <- all.glofs.v0.quant %>%
  group_by(region) %>%
  summarise(n = paste0("n = ", n()),
            minv0 = min(Mean_Flood_Volume_V0, na.rm = T))

# Reproduce Extended Data Figure 2.

quants.05.09.v0 <- all_regional_preds %>%
  mutate(region = str_replace(region, "Pacific NW", "NW North America")) %>%
  ggplot() +
  stat_lineribbon(aes(x = year, 
                      y = pred, 
                      color = ordered(quantile),
                      fill  = ordered(quantile)), 
                  .width = c(.95), 
                  alpha = 0.5) + 
  scale_fill_manual("Posterior\ninterval", 
                    values =  c( "navy", "darkorange")) +
  scale_color_manual("Posterior\ninterval", 
                     values = c( "navy", "darkorange")) + 
  geom_point(data = all.glofs.v0.quant,
             mapping = aes(x = rounded_year,
                           y = Mean_Flood_Volume_V0),
             alpha = 0.7) +
  facet_wrap(~region, scales = "free_y") +
  labs(x = "Year" ,
       y = expression(paste("Flood volume ", V[0], " [", 10^{6}, " ", m^{3}, ']'))) +
  theme_bw()   +
  scale_y_continuous(trans = log10_trans(),
                     breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x))) + 
  geom_text(data  = nobs, 
            aes(x = 1900, y = minv0, label = n), 
            colour = "gray10", 
            inherit.aes = FALSE, 
            parse = FALSE,
            hjust = "left")

# Reproduce Figure 1 c.

trends.quants.v0 <- bind_rows(v0.region.trends.list) %>%
  mutate(region = str_replace_all(region, "[.]", "\n"),
         region = str_replace(region, "Pacific NW", "NW North America")) %>%
  mutate(region = as_factor(region) %>% fct_rev(),
         quantile = as_factor(quantile) %>% fct_rev()) %>%
  ggplot(aes(x = b_orig,
             y = region,
             color = quantile)) +
  coord_cartesian(xlim = c(-0.032, 0.032)) + 
  stat_pointinterval( .width = c(.95), 
                      size = 1,
                      point_size = 1.25,
                      position = position_dodge(width = -0.4)) + 
  scale_color_manual(values = c("darkorange", "navy")) +
  theme_bw() +
  labs(y = "",
       x = expression(paste("Trend in ", log[10], " ", V[0], " [", 10^{6}, " ", m^{3}, ' ', yr^{-1}, ']'))) +
  geom_vline(xintercept = 0) +
  theme( axis.text = element_text(size = 8.5),
         axis.text.x = element_text(size = 8.5),
         axis.title = element_text(size = 9),
         strip.text = element_text(size = 8.5),
         axis.text.y=element_blank(),
         axis.title.y=element_blank())

# Combine the plot of the temporal trend in V0 and the posterior regression slope.

arr.quant.v0 <- ggarrange(plotlist = list(quants.05.09.v0, 
                                          trends.quants.v0), 
                          labels = "auto",
                          ncol = 2, 
                          align = "h", 
                          common.legend = T, 
                          widths = c(2,1))

# Combine the posterior trends in one plot, labelled b and c in Figure 1.

arr.quant.trends.qp.v0 <- ggarrange(plotlist = list(trends.quants.qp, 
                                                    trends.quants.v0), 
                                    labels = c("b", "c"),
                                    ncol = 2, 
                                    align = "hv", 
                                    legend = "none",
                                    common.legend = T) +
  theme(plot.margin = margin(0.2, 0.7, 0.2, 0.2, "cm")) 

# Write the output to disk.

ggsave(
  filename = "v0_quants05_09.pdf",
  plot = quants.05.09.v0 ,
  width = 180,
  height = 100,
  units = "mm"
)

ggsave(
  filename = "v0_quants05_09.png",
  plot = quants.05.09.v0 ,
  width = 180,
  height = 100,
  units = "mm"
)

ggsave(
  filename = "v0_quants05_09.png",
  plot = arr.quant.v0,
  width = 250,
  height = 170,
  units = "mm"
)

ggsave(
  filename = "trends_quants05_09_qp_V0.pdf",
  plot = arr.quant.trends.qp.v0 ,
  width = 130,
  height = 60,
  units = "mm"
)

##########################################################################################
