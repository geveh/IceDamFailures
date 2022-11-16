################################################################################
######                 Calculate regional and local trends in             ######
######                        GLOF volume and peak discharge              ######
######                                                                    ######
######                                  by Georg Veh                      ######
######                                 03 March, 2022                     ######
######                 multi-level structure revised 26 May, 2022         ######
######                    comments added 11 and 15 Nov, 2022              ######
################################################################################

# Load the following packages, or use install.packages("nameofpackage"), 
# if some of them are not pre-installed. In some cases you need to restart 
# your R session.

require(tidyverse)
require(brms)
require(tidybayes)
require(modelr)
require(scales)
require(readODS)
require(patchwork)

# Pay special attention on installing the cmdstanr package, which
# we use to speed up the Bayesian models.
# https://mc-stan.org/cmdstanr/


# Set YOUR working directory folder where to find all files to run this script.
# Change the location appropriately.

setwd("D:/data/Lake_area_volume/")

# Useful functions

scale_this <- function(x){
  (x - mean(x, na.rm=TRUE)) / sd(x, na.rm=TRUE)
}

# Write a hierarchical Bayesian quantile regression model that uses an Asymmetric
# LaPlace likelihood.

# Create a function that can eat different 
# - percentiles, 
# - responses and predictors, 
# - conditioning variables such as regions or individual lakes, and 
# - time periods.

qr_mod_pred <- function(dat, resp, pred, cond, q, tstart, tend, levels) {

  # Fit a hierarchical QR model in brms for a specified response variable, 
  # one predictor, a quantile, conditional variable, and time period.
  
  # Variable names:
  # dat: tibble containing the data.
  # resp: character name of the response variable in that tibble.
  # pred: character name of the predictor variable in that tibble.
  # cond: character conditioning variable.
  # q: numeric conditional quantile.
  # tstart: integer start year, for which new predictions are desired.
  # tend: integer end year, for which new predictions are desired.
  # levels: vector of group levels for which new predictions are desired.
  
  # Extract predictor and response
  
  x <- dat[, pred] %>% as_vector() %>% unname()
  y <- dat[, resp] %>% as_vector() %>% unname()
  
  # Scale predictor and response variables.
  
  x_scale <- scale(x)[ ,1]
  y_scale <- scale(log10(y))[ ,1]
  
  sd_y     <- sd(log10(y))
  mean_y   <- mean(log10(y)) 
  sd_x   <- sd(x)
  mean_x <- mean(x)
  
  dat.new <- dat %>%
    mutate(x_scale = x_scale,
           y_scale = y_scale)
  
  # Specify priors
  
  bprior <- prior(normal(0, 2.5), class = "Intercept") +
            prior(normal(0, 2.5), class = "b") +
            prior(normal(0, 2.5), class = "sd") +
            prior(normal(0, 2.5), class = "sigma")
  
  # We model the trend in a given quantile for all data in a given 
  # time interval.
  
  mod <- brm(bf(as.formula(paste0("y_scale ~ x_scale + (x_scale |", cond, ")")),  
                          quantile = q),
                       family = asym_laplace(),
                       data = dat.new %>% filter(
                         rounded_year >= tstart,
                         rounded_year <= tend),
                       prior = bprior,
                       cores  = 3,
                       chains = 3,
                       warmup = 2000,
                       iter   = 6000,
                       control = list(adapt_delta = 0.92,
                                      max_treedepth = 15),
                       backend = "cmdstanr",
                       threads = threading(3))
  
  # Define the range, for which new predictions will be made.
  # ... on original scale
  
  seq.x.orig <- seq(from = (tstart-5),
                    to   = (tend+5),
                    length.out = 132)
  
  # ... and on scaled range
  
  seq.x.scaled <- (seq.x.orig - mean_x) / sd_x
  
  # Obtain the standardized posterior predictive distribution for new observations
  # and convert the predictions to original scale.
  # Create an a tibble that we will use to specify the predictions.
  
  pred.tibble <- tibble(x_scale := (rep(seq.x.scaled,  times = length(levels))),
                        !!pred  := rep(seq.x.orig,     times = length(unique(levels))),
                        !!cond  := rep(unique(levels), each  = length(seq.x.scaled)),
                        quantile = q)
  
  # Make 500 predictions for each observation.
  # Convert predictions to original scale.
  
  post_epred  <- posterior_epred(mod,
                                 dpar   = "mu",
                                 ndraws = 500,
                                 re_formula = NULL,
                                 newdata    = pred.tibble) %>%
    t() %>%
    as_tibble(.name_repair = "unique") %>%  
    mutate(x_scale := pred.tibble[, 3] %>% as_vector(),
           !!pred  := pred.tibble[, 2] %>% as_vector(),
           quantile = q) %>%
    pivot_longer( cols = starts_with("..."), 
                  values_to = "prediction") %>%
    mutate(prediction = (prediction*sd_y) + mean_y, 
           prediction = 10^prediction,
           first_year = tstart,
           last_year  = tend,
           response   = resp) %>%
    dplyr::select(-name) 
  
  # Extract posterior trends for each group. 
  
  slopes.cond <- tidy_draws(mod) %>% 
    pivot_longer(starts_with(paste0("r_", cond)),
                 names_to = "param",
                 values_to = "post") %>% 
    filter(str_detect(param, 'x_scale')) %>%
    transmute(b_x_scale, param, post) %>%
    mutate(param = str_replace_all(param, "[.]", " "),
           !!cond := str_extract_all(param,"(?<=\\[).+(?=,)") %>% unlist(),
           post_orig = b_x_scale + post,
           post_orig = post_orig  * sd_y / sd_x) 
  
  # Extract posterior trends for the pooled model. 
    
  slopes.all <- mod %>%
    gather_draws(b_x_scale) %>% 
    mutate(!!cond := "Global", param = "x_scale") %>% 
    mutate(post_orig = .value  * sd_y / sd_x) %>%
    ungroup() %>%
    select(!starts_with("."))
  
  # Combine slopes from pooled estimate and from group levels.
  
  slopes <- bind_rows(slopes.cond ,
                      slopes.all) %>%
    mutate(first_year = tstart,
           last_year  = tend,
           response   = resp,
           quantile = q) 
  
  # The output is a list that has the model, predictions for specific 
  # observations, and the posterior slopes.
  
  mod_and_pred <- list(model = mod,
                       prediction = post_epred,
                       posterior_slopes = slopes)
  
  return(mod_and_pred)

}


###################### REGIONAL MODELS     #####################################

###################### Peak discharge      #####################################
###################### Group-level model   #####################################

# Load Qp data from your working directory.

all.glofs.qp <- readRDS("all_glofs_qp_tibble.RDS")

# We calculate quantile regression models for 4 different time periods,
# starting in 1900, 1930, 1960, and 1990, and ending in 2021

cut_years <- c(1900, 1930, 1960, 1990)

# We calculate the 50th (median) and 90th (10% highest) percentile.

qp_tib <- crossing(q = c(0.5, 0.9), tstart = cut_years)

qp.mods <- list()

# Iterate over combinations of periods and quantiles.

for (i in 1: nrow(qp_tib)) {
  
  
  qp.quant.period <-  qr_mod_pred(dat = all.glofs.qp , 
                              resp = "Peak_discharge_Qp", 
                              pred = "rounded_year",
                              cond = "region", 
                              q = qp_tib$q[i],
                              tstart =  qp_tib$tstart[i],
                              tend = 2021,
                              levels = unique(all.glofs.qp$region))
  
  
  qp.mods[[i]] <- qp.quant.period
  
}

# Save the output to disc.

saveRDS(qp.mods, "qp_models.RDS")
# qp.mods <- readRDS("qp_models.RDS")

# In the following, we only focus on two periods, 1900 - 2021, and 1990 - 2021.
# All other periods remain untouched here.

all.preds.qp <- lapply(qp.mods, function (x) x$prediction) %>% 
             bind_rows() %>%
             mutate(period = paste0(first_year, " - ", last_year)) %>%
             filter(first_year == 1900 | first_year == 1990)

nobs.qp <- list()

# Period 1900- 2021 is the first, 1990 - 2021 the fourth object in the list
# generated above, according to the structure of qp_tib, which we iterated over.

# Obtain number of observations per period.

for(i in c(1,4)) {
  
  qp.orig.n <- all.glofs.qp %>%
               group_by(region) %>%
               summarise(n_orig = n(),
                         min_y = min(Peak_discharge_Qp))
  
  qp.sub.n <- all.glofs.qp %>%
    group_by(region) %>%
    filter(rounded_year >= cut_years[i]) %>%
    summarise(n_sub = n())
  
  label.tibble <- left_join(qp.sub.n, qp.orig.n, by = "region") %>%
    mutate(perc = round(n_sub/n_orig*100, digits = 0),
           period = paste0(cut_years[i], " - ", 2021),
           n_perc = paste0("n = ", n_sub, " (", perc, "%)" ),
           rounded_year = 1895) %>%
    rename(x_scale = region)
    
  nobs.qp[[i]] <- label.tibble
  
}

# Bind output per period together.

nobs.qp <- bind_rows(nobs.qp)

# Add colors to points within and outside of the study period.

qp.points.list <- list()

for(i in c(1,4)) {
  
  new.tib <- all.glofs.qp %>% 
    rename(x_scale = region) %>%
    mutate(pt_color = if_else(rounded_year < cut_years[i], "gray75", "black"),
            period =  paste0(cut_years[i], " - ", 2021))
  
  qp.points.list[[i]] <- new.tib
}

qp.points <- do.call(rbind, qp.points.list)

# Generate a plot showing the trend in Qp for the 50th and 90th percentile with
# time for the period 1900 - 2021. 

quants.qp.q.50.q90 <- all.preds.qp %>%
  ggplot() +
  geom_point(data = qp.points,
             mapping = aes(x = rounded_year,
                           y = Peak_discharge_Qp),
             color = qp.points$pt_color,
             alpha = 0.66,
             size  = 1) +
  stat_lineribbon(aes(x = rounded_year,
                      y = prediction,
                      color = ordered(quantile),
                      fill = ordered(quantile)),
                  .width = c(.95),
                  alpha = 0.6) +
  scale_fill_manual("Quantile", 
                     values =  c("navy", "darkorange")) +
  scale_color_manual("Quantile", 
                      values = c("navy", "darkorange")) + 
  facet_grid(rows = vars(x_scale),
             cols = vars(period)) +
  labs(x = "Year" ,
       y = expression(paste("Peak discharge ", Q[p], " [", m^{3}, " ", 
                            s^{-1}, ']'))) +
  theme_bw()   +
  scale_y_continuous(trans  = log10_trans(),
                     breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x))) +
  geom_text(data  = nobs.qp,
            aes(x = rounded_year, y = 1, label = n_perc),
            size = 2.5,
            colour = "gray10",
            inherit.aes = FALSE,
            parse = FALSE,
            hjust = "left") +
  theme( axis.text   = element_text(size = 7),
         axis.text.x = element_text(size = 7),
         axis.title  = element_text(size = 7),
         strip.text  = element_text(size = 7),
         legend.position = "bottom")

# Extract posterior slopes for the two periods.

all.slopes.qp <- lapply(qp.mods, function (x) x$posterior_slopes) %>% 
  bind_rows() %>%
  mutate(period = paste0(first_year, " - ", last_year),
         posterior_slope = "Posterior slope") %>%
  filter(first_year == 1900 | first_year == 1990)

# Calculate how much of the posterior trend is smaller or larger than zero.

qp.perc <- all.slopes.qp %>%
  filter(region != "Global") %>%
  mutate(period = paste0(first_year, " - ", last_year))%>%
  group_by(region, period, quantile) %>%
  summarise(perc_sm_0 = (sum((post_orig) <0) / n())*100,
            perc_gt_0 = 100-perc_sm_0,
            perc_sm_0 = round(perc_sm_0, digits = 1),
            perc_gt_0 = round(perc_gt_0, digits = 1),
            perc_sm_0 = paste0(perc_sm_0, "%"),
            perc_gt_0 = paste0(perc_gt_0, "%")) 

# Draw the posterior slopes for each regions as HDI bars, and add labels
# for the posterior masses being greater or smaller than zero.

mod_param.qp <- all.slopes.qp %>%
  filter(region != "Global") %>%
  mutate(decadal_trend = post_orig*10) %>%
  ggplot(aes(x = decadal_trend,
             y = as.character(quantile), 
             color = period)) +
  stat_pointinterval( .width = c(0.95), 
                      point_size = 1,
                      size = 0.8,
                      position = position_dodge(width = .6)) + 
  geom_text(data = qp.perc,
            mapping = aes(x = -0.7, 
                          y = as.character(quantile),
                          group = period,
                          label = perc_sm_0), 
            position = position_dodge(width = .6),
            color = "black",
            size = 2) +
  geom_text(data = qp.perc,
            mapping = aes(x = 0.6, 
                          y = as.character(quantile),
                          group = period,
                          label = perc_gt_0), 
            position = position_dodge(width = .6),
            color = "black",
            size = 2) +
  scale_color_manual("Period", values = c( "darkgreen", "darkolivegreen3")) +
  facet_grid(rows = vars(region),
             cols = vars(posterior_slope)) +
  theme_bw() +
  xlim(c(-0.8, 0.8)) +
  labs(x = expression(paste("Trend in ", log[10], "-", Q[p], " [",  m^{3}, ' ', 
                            s^{-1}, ' ', dec^{-1}, "]")),
       y = "Quantile")  +
  geom_vline(xintercept = 0) +
  theme( axis.text   = element_text(size = 7),
         axis.text.x = element_text(size = 7),
         axis.title  = element_text(size = 7),
         strip.text  = element_text(size = 7),
         legend.position = "bottom") +
  guides(color = guide_legend(reverse=TRUE)) 


# Output figure for Qp for Extended data

qp.fin <- quants.qp.q.50.q90 + mod_param.qp + 
          plot_layout(widths = c(3, 1)) & theme(legend.position = "bottom")

qp.plot <- qp.fin + plot_layout(guides = "collect") + 
  plot_annotation(tag_levels = 'a') & 
  theme(plot.tag = element_text(size = 8, face = "bold"),
        legend.title = element_text(size = 7), 
        legend.text = element_text(size = 7))

ggsave(filename = "qp_sensitivity.pdf", qp.plot,
       width = 183, height = 200,
       units = "mm")

###################### Pooled model   ##########################################
###################### Peak discharge ##########################################

# Same plotting approach as above, left largely uncommented.

# Obtain standard deviation and peak discharge on original scale.

sd_year <- sd(all.glofs.qp$rounded_year)
m_year  <- mean(all.glofs.qp$rounded_year)
sd_qp   <- sd(log10(all.glofs.qp$Peak_discharge_Qp))
m_qp    <- mean(log10(all.glofs.qp$Peak_discharge_Qp))

# Generate range of predicted values

pred_range <- seq((1900 - m_year) / sd_year, 
                  (2021 - m_year) / sd_year, 
                  by = 0.02)

# Draw from the posterior predictive distribution.

# 50th percentile of Qp since 1900

qp.pool.50.1900 <- posterior_epred(qp.mods[[1]]$model, 
                              newdata = data.frame(x_scale = pred_range ), 
                              re_formula = NA,
                              dpar = "mu",
                              ndraws = 1000) %>%
  t() %>%
  as_tibble(.name_repair = "unique") %>%  
  mutate(x_scale = pred_range ) %>%
  pivot_longer( cols = starts_with("..."), 
                values_to = "prediction") %>%
  mutate(q = 50,
         period = "1900 - 2021")

# 90th percentile of Qp since 1900

qp.pool.90.1900 <- posterior_epred(qp.mods[[5]]$model, 
                              newdata = data.frame(x_scale = pred_range ), 
                              re_formula = NA,
                              dpar = "mu",
                              ndraws = 1000) %>%
  t() %>%
  as_tibble(.name_repair = "unique") %>%  
  mutate(x_scale = pred_range) %>%
  pivot_longer( cols = starts_with("..."), 
                values_to = "prediction") %>%
  mutate(q = 90,
         period = "1900 - 2021")

# As above, draw from the posterior predictive distribution, but this time
# from 50th percentile of Qp since 1990.

pred_range <- seq((1990 - m_year) / sd_year, 
                  (2021 - m_year) / sd_year, 
                  by = 0.02)

qp.pool.50.1990 <- posterior_epred(qp.mods[[4]]$model, 
                                   newdata = data.frame(x_scale = pred_range ), 
                                   re_formula = NA,
                                   dpar = "mu",
                                   ndraws = 1000) %>%
  t() %>%
  as_tibble(.name_repair = "unique") %>%  
  mutate(x_scale = pred_range ) %>%
  pivot_longer( cols = starts_with("..."), 
                values_to = "prediction") %>%
  mutate(q = 50,
         period = "1990 - 2021")

# And from 90th percentile of Qp since 1990.

qp.pool.90.1990 <- posterior_epred(qp.mods[[8]]$model, 
                                   newdata = data.frame(x_scale = pred_range ), 
                                   re_formula = NA,
                                   dpar = "mu",
                                   ndraws = 1000) %>%
  t() %>%
  as_tibble(.name_repair = "unique") %>%  
  mutate(x_scale = pred_range) %>%
  pivot_longer( cols = starts_with("..."), 
                values_to = "prediction") %>%
  mutate(q = 90,
         period = "1990 - 2021")

# Generate two similar plots: 

# - one that shows the trends of 50th and 90th percentile of Qp since 1900.

qp.pooled.mod.1900 <- rbind(qp.pool.50.1900,
                            qp.pool.90.1900 ) %>%
  mutate(x_orig = x_scale *sd_year + m_year,
         prediction = prediction * 
           sd_qp + 
           m_qp,
         y_orig = 10^prediction) %>%
  ggplot() +
  geom_point(data = all.glofs.qp,
             mapping = aes(x = rounded_year,
                           y = Peak_discharge_Qp),
             alpha = 0.66,
             size  = 1) +
  stat_lineribbon(aes(x = x_orig,
                      y = y_orig ,
                      color = ordered(q),
                      fill = ordered(q)),
                  .width = c(.95),
                  alpha = 0.6) +
  scale_fill_manual("Quantile", 
                    values =  c("navy", "darkorange")) +
  scale_color_manual("Quantile", 
                     values =  c("navy", "darkorange")) +
  labs(x = "Year" ,
       y = expression(paste("Peak discharge ", Q[p], " [", m^{3}, " ", 
                            s^{-1}, ']'))) +
  theme_bw()   +
  scale_y_continuous(trans  = log10_trans(),
                     breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x))) +
  theme( axis.text   = element_text(size = 7),
         axis.text.x = element_text(size = 7),
         axis.title  = element_text(size = 7),
         strip.text  = element_text(size = 7),
         legend.position = "none")

# - and one that shows the trends of 50th and 90th percentile of Qp since 1990.

qp.pooled.mod.1990 <- rbind(qp.pool.50.1990,
                            qp.pool.90.1990 ) %>%
  mutate(x_orig = x_scale *sd_year + m_year,
         prediction = prediction * 
           sd_qp + 
           m_qp,
         y_orig = 10^prediction) %>%
  ggplot() +
  geom_point(data = all.glofs.qp %>% filter(rounded_year < 1990),
             mapping = aes(x = rounded_year,
                           y = Peak_discharge_Qp),
             color = "gray75",
             alpha = 0.66,
             size  = 1) +
  geom_point(data = all.glofs.qp %>% filter(rounded_year >= 1990),
             mapping = aes(x = rounded_year,
                           y = Peak_discharge_Qp),
             color = "black",
             alpha = 0.66,
             size  = 1) +
  stat_lineribbon(aes(x = x_orig,
                      y = y_orig ,
                      color = ordered(q),
                      fill = ordered(q)),
                  .width = c(.95),
                  alpha = 0.6) +
  scale_fill_manual("Quantile", 
                    values =  c("navy", "darkorange")) +
  scale_color_manual("Quantile", 
                     values =  c("navy", "darkorange")) +
  labs(x = "Year" ,
       y = expression(paste("Peak discharge ", Q[p], " [", m^{3}, " ", 
                            s^{-1}, ']'))) +
  theme_bw()   +
  scale_y_continuous(trans  = log10_trans(),
                     breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x))) +
  theme( axis.text   = element_text(size = 7),
         axis.text.x = element_text(size = 7),
         axis.title  = element_text(size = 7),
         strip.text  = element_text(size = 7),
         legend.position = "none")

# Obtain for each model (50th and 90th percentile for the periods 1900-2021
# and 1990-2021) the posterior decadal trend in Qp.

qp.pooled.slopes <- rbind(as_draws_df(qp.mods[[1]]$model, "b_x_scale") %>% 
                            mutate(q = 50, 
                                   period = "1900-2021"), 
                          as_draws_df(qp.mods[[5]]$model, "b_x_scale") %>% 
                            mutate(q = 90, 
                                   period = "1900-2021"), 
                          as_draws_df(qp.mods[[4]]$model, "b_x_scale") %>% 
                            mutate(q = 50, 
                                   period = "1990-2021"), 
                          as_draws_df(qp.mods[[8]]$model, "b_x_scale") %>% 
                            mutate(q = 90, 
                                   period = "1990-2021")) %>%
  as_tibble() %>%
  mutate(b_orig = (b_x_scale * 
                     sd_qp / 
                     sd_year)* 10 ) %>%
  mutate(q = as_factor(q),
         period = as_factor(period)) 

# Calculate how much of the posterior trend is smaller and larger than zero,
# respectively.

qp.pooled.perc <- qp.pooled.slopes  %>%
  group_by(period, q) %>%
  summarise(perc_sm_0 = (sum((b_orig) <0) / n())*100,
            perc_gt_0 = 100-perc_sm_0,
            perc_sm_0 = round(perc_sm_0, digits = 1),
            perc_gt_0 = round(perc_gt_0, digits = 1),
            perc_sm_0 = paste0(perc_sm_0, "%"),
            perc_gt_0 = paste0(perc_gt_0, "%")) 

# Draw the posterior slopes for the 50th and 90th percentile for both periods.

qp.pooled.post <- qp.pooled.slopes  %>%
  ggplot(aes(x = b_orig,
             y = q, 
             fill = period)) +
  stat_halfeye( mapping = aes(fill = period),
                .width = c( .95),
                interval_size = 0.8,
                shape = 21,
                point_color = "black",
                point_fill = "gray80",
                point_size = 1,
                position = position_dodge(width = 0.6))+ 
  scale_fill_manual("Period", values = c( "darkgreen", "darkolivegreen3")) +
  geom_text(data = qp.pooled.perc,
            mapping = aes(x = -0.4,
                          y = as.factor(q),
                          group = period,
                          label = perc_sm_0),
            position = position_dodge(width = .6),
            color = "black",
            size = 2)  +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_bw() +
  labs(x = expression(paste("Pooled decadal trend in ", 
                            log[10], "-", Q[p],
                            " [",  m^{3}, ' ', s^{-1}, ' ', dec^{-1}, "]")),
       y = "Quantile") +
  xlim(c(-0.5, 0.5)) +
  theme( axis.text   = element_text(size = 7),
         axis.text.x = element_text(size = 7),
         axis.title  = element_text(size = 7),
         strip.text  = element_text(size = 7),
         legend.position = "none")

###################### Flood volume ############################################   

# Load data on flood volumes V0

all.glofs.V0 <- readRDS("all_glofs_V0_tibble.RDS")

# As above, iterate over different periods, starting in 1900, 1930, 1960, and 1990
# and ending in 2021.

cut_years <- c(1900, 1930, 1960, 1990)

V0_tib <- crossing(q = c(0.5, 0.9), tstart = cut_years)

# Generate the same model outputs as above, including the Bayesian hierachical 
# model and the posterior trends, and predictive distributions.

V0.mods <- list()

for (i in 1: nrow(V0_tib)) {
  
  
  V0.quant.period <-  qr_mod_pred(dat = all.glofs.V0 , 
                                  resp = "Mean_Flood_Volume_V0", 
                                  pred = "rounded_year",
                                  cond = "region", 
                                  q = V0_tib$q[i],
                                  tstart =  V0_tib$tstart[i],
                                  tend = 2021,
                                  levels = unique(all.glofs.V0$region))
  
  
  V0.mods[[i]] <- V0.quant.period
  
}

# Save output to disc.

saveRDS(V0.mods, "V0_models.RDS")
# V0.mods <- readRDS("V0_models.RDS")

all.preds.V0 <- lapply(V0.mods, function (x) x$prediction) %>% 
  bind_rows() %>%
  mutate(period = paste0(first_year, " - ", last_year)) %>%
  filter(first_year == 1900 | first_year == 1990)

# Focus only on the period 1900 - 2021 (first element) and 1990 - 2021 (fourth
# element in the list).

# Obtain the number of observations in each period and region,
# which we will add as labels to the plot.

nobs.V0 <- list()

for(i in c(1, 4)) {
  
  V0.orig.n <- all.glofs.V0 %>%
    group_by(region) %>%
    summarise(n_orig = n(),
              min_y = min(Mean_Flood_Volume_V0))
  
  V0.sub.n <- all.glofs.V0 %>%
    group_by(region) %>%
    filter(rounded_year >= cut_years[i]) %>%
    summarise(n_sub = n())
  
  label.tibble <- left_join(V0.sub.n, V0.orig.n, by = "region") %>%
    mutate(perc = round(n_sub/n_orig*100, digits = 0),
           period = paste0(cut_years[i], " - ", 2021),
           n_perc = paste0("n = ", n_sub, " (", perc, "%)" ),
           rounded_year = 1895) %>%
    rename(x_scale = region)
  
  nobs.V0[[i]] <- label.tibble
  
}

nobs.V0 <- bind_rows(nobs.V0)


V0.points.list <- list()

for(i in c(1,4)) {
  
  new.tib <- all.glofs.V0 %>% 
    rename(x_scale = region) %>%
    mutate(pt_color = if_else(rounded_year < cut_years[i], "gray75", "black"),
           period =  paste0(cut_years[i], " - ", 2021))
  
  V0.points.list[[i]] <- new.tib
}

V0.points <- do.call(rbind, V0.points.list) 

# Plot the trends in the quantile.

quants.V0.q.50.q90 <- all.preds.V0 %>%
  ggplot() +
  geom_point(data = all.glofs.V0 %>% rename(x_scale = region),
             mapping = aes(x = rounded_year,
                           y = Mean_Flood_Volume_V0),
             color = V0.points$pt_color,
             alpha = 0.66,
             size  = 1) +
  stat_lineribbon(aes(x = rounded_year,
                      y = prediction,
                      color = ordered(quantile),
                      fill = ordered(quantile)),
                  .width = c(0.95),
                  alpha = 0.6) +
  scale_fill_manual("Quantile", 
                    values =  c("navy", "darkorange")) +
  scale_color_manual("Quantile", 
                     values = c("navy", "darkorange")) + 
  facet_grid(rows = vars(x_scale),
             cols = vars(period)) +
  labs(x = "Year" ,
       y = expression(paste("Flood volume ", V[0], " [", 10^{6}, " ", m^{3}, ']'))) +
  theme_bw()   +
  scale_y_continuous(trans  = log10_trans(),
                     breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x))) +
  geom_text(data  = nobs.V0,
            aes(x = rounded_year, y = 10^-2, label = n_perc),
            size = 2.5,
            colour = "gray10",
            inherit.aes = FALSE,
            parse = FALSE,
            hjust = "left") +
  theme( axis.text   = element_text(size = 7),
         axis.text.x = element_text(size = 7),
         axis.title  = element_text(size = 7),
         strip.text  = element_text(size = 7),
         legend.position = "bottom")

# Obtain the posterior slopes in V0.

all.slopes.V0 <- lapply(V0.mods, function (x) x$posterior_slopes) %>% 
  bind_rows() %>%
  mutate(period = paste0(first_year, " - ", last_year),
         posterior_slope = "Posterior slope") %>%
  filter(first_year == 1900 | first_year == 1990)

# Find out how much of the posterior trends of V0 is smaller or larger than 
# zero, respectively.

V0.perc <- all.slopes.V0 %>%
  filter(region != "Global") %>%
  mutate(period = paste0(first_year, " - ", last_year))%>%
  group_by(region, period, quantile) %>%
  summarise(perc_sm_0 = (sum((post_orig) <0) / n())*100,
            perc_gt_0 = 100-perc_sm_0,
            perc_sm_0 = round(perc_sm_0, digits = 1),
            perc_gt_0 = round(perc_gt_0, digits = 1),
            perc_sm_0 = paste0(perc_sm_0, "%"),
            perc_gt_0 = paste0(perc_gt_0, "%")) 

# Plot the posterior slopes of V0.

mod_param.V0 <- all.slopes.V0 %>%
  filter(region != "Global") %>%
  mutate(decadal_trend = post_orig*10) %>%
  ggplot(aes(x = decadal_trend,
             y = as.character(quantile), 
             color = period)) +
  stat_pointinterval( .width = 0.95, 
                      point_size = 1,
                      size = 0.8,
                      position = position_dodge(width = .6)) + 
  scale_color_manual("Period", values = c( "darkgreen", "darkolivegreen3")) +
  geom_text(data = V0.perc,
            mapping = aes(x = -0.65, 
                          y = as.character(quantile),
                          group = period,
                          label = perc_sm_0), 
            position = position_dodge(width = .6),
            color = "black",
            size = 2) +
  geom_text(data = V0.perc,
            mapping = aes(x = 0.65, 
                          y = as.character(quantile),
                          group = period,
                          label = perc_gt_0), 
            position = position_dodge(width = .6),
            color = "black",
            size = 2) +
  facet_grid(rows = vars(region),
             cols = vars(posterior_slope)) +
  theme_bw() +
  xlim(c(-0.8, 0.8)) +
  labs(x = expression(paste("Trend in ", log[10], "-", V[0], " [", 10^{6}, ' ', m^{3}, ' ', dec^{-1}, "]")),
       y = "Quantile")  +
  geom_vline(xintercept = 0) +
  theme( axis.text   = element_text(size = 7),
         axis.text.x = element_text(size = 7),
         axis.title  = element_text(size = 7),
         strip.text  = element_text(size = 7),
         legend.position = "bottom") +
  guides(color = guide_legend(reverse=TRUE)) 

# Plot trends in flood volumes for the Extended Data.

V0.fin <- quants.V0.q.50.q90 + mod_param.V0 + 
  plot_layout(widths = c(3, 1)) & theme(legend.position = "bottom")

V0.plot <- V0.fin + plot_layout(guides = "collect") + 
  plot_annotation(tag_levels = 'a') & 
  theme(plot.tag = element_text(size = 8, face = "bold"),
        legend.title = element_text(size = 7), 
        legend.text = element_text(size = 7))

ggsave(filename = "V0_sensitivity.pdf", V0.plot,
       width = 183, height = 200,
       units = "mm")


# Generate a figure containing the posterior slopes
# for Qp and V0 for two different time periods.

# Obtain number of observations (n)

nobs.tot <- rbind(nobs.qp %>% 
                    mutate(response = "Peak_discharge_Qp"),
                  nobs.V0 %>% mutate(response = "Mean_Flood_Volume_V0")) %>%
  filter(period == "1900 - 2021") %>%
  rename("region" = "x_scale")

# Extract slopes for the period 1990 - 2021

slopes.since.1990 <- rbind(all.slopes.V0 %>% 
                             filter(first_year == 1990),
                           all.slopes.qp %>%
                             filter(first_year == 1990)) %>%
  filter(region != "Global") %>%
  mutate(decadal_trend = post_orig*10) 

fig2 <- rbind(all.slopes.V0 %>% 
                filter(first_year == 1900) %>%
                mutate(col = if_else(quantile == 0.5,
                                     "coral1", "coral4")),
              all.slopes.qp %>%
                filter(first_year == 1900)%>%
                mutate(col = if_else(quantile == 0.5,
                                     "deepskyblue1", "deepskyblue4"))) %>%
  filter(region != "Global") %>%
  mutate(decadal_trend = post_orig*10) %>%
  ggplot(aes(x = decadal_trend,
             y = as.character(response))) +
  stat_slab(data = slopes.since.1990,
            aes(x = decadal_trend,
                y = as.character(response),
                group = quantile),
            size = 0.5, 
            fill = "gray82",
            slab_color = NA,
            interval_size = 0, 
            position = position_dodge(width = 0.6)) + 
  stat_halfeye( mapping = aes(fill = col),
                .width = c( .95),
                interval_size = 0.8,
                shape = 21,
                point_color = "black",
                point_fill = "gray80",
                point_size = 1,
                position = position_dodge(width = 0.6)) +
  scale_fill_manual("Quantile", values = c("coral4",
                                           "coral1",
                                           "deepskyblue4",
                                           "deepskyblue1")) +
  facet_wrap(~region) +
  theme_bw()+
  xlim(-0.71, 0.71) +
  labs(x = expression(paste("Trend in ", log[10], "-", italic(V)[0], 
                            " [", 10^{6}, ' ', m^{3}, ' ', dec^{-1}, "]",
                            " and ", log[10], "-", italic(Q)[p], 
                            " [",  m^{3}, ' ', s^{-1}, ' ', dec^{-1}, "]")),
       y = "GLOF magnitude")  +
  geom_vline(xintercept = 0, linetype = "longdash") +
  scale_y_discrete(labels = c(expression(paste(italic(V)[0])),
                              expression(paste(italic(Q)[p])))) +
  geom_text(data  = nobs.tot,
            aes(x = -0.7, y = response, label = n_orig),
            size = 2.5,
            colour = "gray10",
            inherit.aes = FALSE,
            parse = FALSE,
            hjust = "left") +
  theme( axis.text = element_text(size = 7),
         axis.text.x = element_text(size = 7),
         axis.title = element_text(size = 7),
         strip.text = element_text(size = 7),
         legend.position = "none") 
         
ggsave(filename = "fig2.pdf", 
       fig2,
       width = 170, height = 90,
       units = "mm")

###################### Pooled model   ##########################################
###################### Flood volume   ##########################################

# Obtain the mean and standard deviation of the model parameters.

sd_year <- sd(all.glofs.V0$rounded_year)
m_year  <- mean(all.glofs.V0$rounded_year)
sd_V0   <- sd(log10(all.glofs.V0$Mean_Flood_Volume_V0))
m_V0    <- mean(log10(all.glofs.V0$Mean_Flood_Volume_V0))

# As above generate predictions of the periods 1900-2021 and 1990-2021 of 
# the trend in the 50th and 90th percentile of flood volumes V0.

pred_range <- seq((1900 - m_year) / sd_year, 
                  (2021 - m_year) / sd_year, 
                  by = 0.02)

V0.pool.50.1900 <- posterior_epred(V0.mods[[1]]$model, 
                                   newdata = data.frame(x_scale = pred_range ), 
                                   re_formula = NA,
                                   dpar = "mu",
                                   ndraws = 1000) %>%
  t() %>%
  as_tibble(.name_repair = "unique") %>%  
  mutate(x_scale = pred_range ) %>%
  pivot_longer( cols = starts_with("..."), 
                values_to = "prediction") %>%
  mutate(q = 50,
         period = "1900 - 2021")

V0.pool.90.1900 <- posterior_epred(V0.mods[[5]]$model, 
                                   newdata = data.frame(x_scale = pred_range ), 
                                   re_formula = NA,
                                   dpar = "mu",
                                   ndraws = 1000) %>%
  t() %>%
  as_tibble(.name_repair = "unique") %>%  
  mutate(x_scale = pred_range) %>%
  pivot_longer( cols = starts_with("..."), 
                values_to = "prediction") %>%
  mutate(q = 90,
         period = "1900 - 2021")


pred_range <- seq((1990 - m_year) / sd_year, 
                  (2021 - m_year) / sd_year, 
                  by = 0.02)

V0.pool.50.1990 <- posterior_epred(V0.mods[[4]]$model, 
                                   newdata = data.frame(x_scale = pred_range ), 
                                   re_formula = NA,
                                   dpar = "mu",
                                   ndraws = 1000) %>%
  t() %>%
  as_tibble(.name_repair = "unique") %>%  
  mutate(x_scale = pred_range ) %>%
  pivot_longer( cols = starts_with("..."), 
                values_to = "prediction") %>%
  mutate(q = 50,
         period = "1990 - 2021")

V0.pool.90.1990 <- posterior_epred(V0.mods[[8]]$model, 
                                   newdata = data.frame(x_scale = pred_range ), 
                                   re_formula = NA,
                                   dpar = "mu",
                                   ndraws = 1000) %>%
  t() %>%
  as_tibble(.name_repair = "unique") %>%  
  mutate(x_scale = pred_range) %>%
  pivot_longer( cols = starts_with("..."), 
                values_to = "prediction") %>%
  mutate(q = 90,
         period = "1990 - 2021")

# Plotting of the trends with time is exactly the same as for peak discharges.

V0.pooled.mod.1900 <- rbind(V0.pool.50.1900,
                            V0.pool.90.1900 ) %>%
  mutate(x_orig = x_scale *sd_year + m_year,
         prediction = prediction * 
           sd_V0 + 
           m_V0,
         y_orig = 10^prediction) %>%
  ggplot() +
  geom_point(data = all.glofs.V0,
             mapping = aes(x = rounded_year,
                           y = Mean_Flood_Volume_V0),
             alpha = 0.66,
             size  = 1) +
  stat_lineribbon(aes(x = x_orig,
                      y = y_orig ,
                      color = ordered(q),
                      fill = ordered(q)),
                  .width = c(.95),
                  alpha = 0.6) +
  scale_fill_manual("Quantile", 
                    values =  c("navy", "darkorange")) +
  scale_color_manual("Quantile", 
                     values =  c("navy", "darkorange")) +
  labs(x = "Year" ,
       y = expression(paste("Flood volume ", 
                            V[0], " [", 10^{6}, " ", m^{3}, ']'))) +
  theme_bw()   +
  scale_y_continuous(trans  = log10_trans(),
                     breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x))) +
  theme( axis.text   = element_text(size = 7),
         axis.text.x = element_text(size = 7),
         axis.title  = element_text(size = 7),
         strip.text  = element_text(size = 7),
         legend.position = "none")

# Plotting of the posterior trends is exactly the same as for peak discharges.

V0.pooled.mod.1990 <- rbind(V0.pool.50.1990,
                            V0.pool.90.1990 ) %>%
  mutate(x_orig = x_scale *sd_year + m_year,
         prediction = prediction * 
           sd_V0 + 
           m_V0,
         y_orig = 10^prediction) %>%
  ggplot() +
  geom_point(data = all.glofs.V0 %>% filter(rounded_year < 1990),
             mapping = aes(x = rounded_year,
                           y = Mean_Flood_Volume_V0),
             color = "gray75",
             alpha = 0.66,
             size  = 1) +
  geom_point(data = all.glofs.V0 %>% filter(rounded_year >= 1990),
             mapping = aes(x = rounded_year,
                           y = Mean_Flood_Volume_V0),
             color = "black",
             alpha = 0.66,
             size  = 1) +
  stat_lineribbon(aes(x = x_orig,
                      y = y_orig ,
                      color = ordered(q),
                      fill = ordered(q)),
                  .width = c(.95),
                  alpha = 0.6) +
  scale_fill_manual("Quantile", 
                    values =  c("navy", "darkorange")) +
  scale_color_manual("Quantile", 
                     values =  c("navy", "darkorange")) +
  labs(x = "Year" ,
       y = expression(paste("Flood volume ", V[0], 
                            " [", 10^{6}, " ", m^{3}, ']'))) +
  theme_bw()   +
  scale_y_continuous(trans  = log10_trans(),
                     breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x))) +
  theme( axis.text   = element_text(size = 7),
         axis.text.x = element_text(size = 7),
         axis.title  = element_text(size = 7),
         strip.text  = element_text(size = 7),
         legend.position = "none")

# Obtain the posterior slopes for the 2 quantiles and 2 periods.

V0.pooled.slopes <- rbind(as_draws_df(V0.mods[[1]]$model, "b_x_scale") %>% 
                            mutate(q = 50, 
                                   period = "1900-2021"), 
                          as_draws_df(V0.mods[[5]]$model, "b_x_scale") %>% 
                            mutate(q = 90, 
                                   period = "1900-2021"), 
                          as_draws_df(V0.mods[[4]]$model, "b_x_scale") %>% 
                            mutate(q = 50, 
                                   period = "1990-2021"), 
                          as_draws_df(V0.mods[[8]]$model, "b_x_scale") %>% 
                            mutate(q = 90, 
                                   period = "1990-2021")) %>%
  as_tibble() %>%
  mutate(b_orig = (b_x_scale * sd_V0 / sd_year)* 10 ) %>%
  mutate(q = as_factor(q),
         period = as_factor(period)) 

V0.pooled.perc <- V0.pooled.slopes  %>%
  group_by(period, q) %>%
  summarise(perc_sm_0 = (sum((b_orig) <0) / n())*100,
            perc_gt_0 = 100-perc_sm_0,
            perc_sm_0 = round(perc_sm_0, digits = 1),
            perc_gt_0 = round(perc_gt_0, digits = 1),
            perc_sm_0 = paste0(perc_sm_0, "%"),
            perc_gt_0 = paste0(perc_gt_0, "%")) 

# Generate the plot of the posterior trends.

V0.pooled.post <- V0.pooled.slopes  %>%
  ggplot(aes(x = b_orig,
             y = q, 
             fill = period)) +
  stat_halfeye( mapping = aes(fill = period),
                .width = c( .95),
                interval_size = 0.8,
                shape = 21,
                point_color = "black",
                point_fill = "gray80",
                point_size = 1,
                position = position_dodge(width = 0.6))+ 
  scale_fill_manual("Period", values = c( "darkgreen", "darkolivegreen3")) +
  geom_text(data = V0.pooled.perc,
            mapping = aes(x = -0.4,
                          y = as.factor(q),
                          group = period,
                          label = perc_sm_0),
            position = position_dodge(width = .6),
            color = "black",
            size = 2)  +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_bw() +
  labs(x = expression(paste("Pooled decadal trend in ", log[10], "-", V[0], 
                            " [", 10^{6}, " ", m^{3}, ' ',  dec^{-1}, "]")),
       y = "Quantile") +
  xlim(c(-0.5, 0.5)) +
  theme( axis.text   = element_text(size = 7),
         axis.text.x = element_text(size = 7),
         axis.title  = element_text(size = 7),
         strip.text  = element_text(size = 7),
         legend.position = "none")


all.pooled.mods <- (V0.pooled.mod.1900 | V0.pooled.mod.1990 | V0.pooled.post) /
  (qp.pooled.mod.1900 | qp.pooled.mod.1990 | qp.pooled.post)  + 
  plot_annotation(tag_levels = 'a') & 
  theme(plot.tag = element_text(size = 8, face = "bold"))

ggsave(filename = "all_pooled_mods.pdf", 
       all.pooled.mods,
       width = 183, height = 140,
       units = "mm")


###################### LOCAL MODELS     ########################################

# Now e focus on local trends, i.e. changes of flood magnitudes for lakes that
# had repeated outbursts

###################### Peak discharge   ########################################

# Select only lakes that had more than five outbursts and reported values of Qp.

all.glofs.qp.local <- all.glofs.qp %>%
  group_by(Glacier_and_lake) %>%
  mutate(Glacier_and_lake = str_replace_all(Glacier_and_lake, " ", "_")) %>%
  filter(n() >= 5) %>%
  ungroup()

# Now condition the model on individual glaciers and lakes.

# We are only interested in the median trend (Q50).

qp.median.local <-  qr_mod_pred(dat    = all.glofs.qp.local, 
                                resp   = "Peak_discharge_Qp", 
                                pred   = "rounded_year",
                                cond   = "Glacier_and_lake", 
                                q      = 0.5,
                                tstart = 1900,
                                tend   = 2021,
                                levels = unique(all.glofs.qp.local$Glacier_and_lake))

# Save the output to disc.

saveRDS(qp.median.local, "qp_model_median_local.RDS")
# qp.median.local <- readRDS("qp_model_median_local.RDS")

# Extract the posterior prediction.

all.preds.qp.local <- qp.median.local$prediction 

conds <- all.glofs.qp.local %>% 
  group_by(Glacier_and_lake) %>% 
  summarise(min_r = min(rounded_year)-5, 
            max_r = max(rounded_year)+5)

# Truncate the predicted values to period of observed V0.

all.preds.qp.local.sub <- list()

for (i in 1:nrow(conds)) {
  
  all.preds.qp.local.sub[[i]] <- filter(all.preds.qp.local, 
                                        (x_scale == conds$Glacier_and_lake[i]) & 
                                          (rounded_year >= conds$min_r[i]) &  
                                          (rounded_year <= conds$max_r[i]))
}

all.preds.qp.local.sub <- bind_rows(all.preds.qp.local.sub) %>%
  mutate(x_scale = str_replace_all(x_scale, "_", " "))

# Generate a plot that shows the trend in median Qp for every lake with time.

qp.local.timeseries <- all.preds.qp.local.sub %>%
  ggplot() +
  stat_lineribbon(aes(x = rounded_year, 
                      y = prediction, 
                      color = ordered(quantile),
                      fill = ordered(quantile)), 
                  .width = c(.95), 
                  alpha = 0.6) + 
  scale_fill_manual( "Quantile", values =  "coral1") +
  scale_color_manual("Quantile", values =  "coral1") + 
  geom_point(data = all.glofs.qp.local %>% 
               rename(x_scale = Glacier_and_lake) %>%
               mutate(x_scale = str_replace_all(x_scale, "_", " ")),
             mapping = aes(x = rounded_year,
                           y = Peak_discharge_Qp),
             alpha = 0.66) +
  facet_wrap(~x_scale, scales = "free", ncol = 4) +
  labs(x = "Year" ,
       y = expression(paste("Peak discharge ", 
                            Q[p], " [", m^{3}, " ", s^{-1}, ']'))) +
  theme_bw()  +
  theme( axis.text   = element_text(size = 6),
         axis.text.x = element_text(size = 6),
         axis.title  = element_text(size = 6),
         strip.text  = element_text(size = 5),
         legend.position = "none")

# Write plots of local Qp to disc.

ggsave(filename = "Qp_local.pdf", 
       qp.local.timeseries,
       width = 183, 
       height = 183,
       units = "mm")

# Generate a plot that shows the posterior decadal the trend in the
# 50th percentile of local Qp.

all.slopes.qp.local <- qp.median.local$posterior_slopes %>% 
  mutate(posterior_slope = "Posterior slope") %>%
  filter(Glacier_and_lake != "Global") %>%
  mutate(decadal_trend = post_orig*10) %>%
  mutate(Glacier_and_lake = str_replace_all(Glacier_and_lake, " ", "."),
         Glacier_and_lake = str_replace_all(Glacier_and_lake, "_", " "))

# Find out how much of the posterior mass is smaller or larger than zero, 
# respectively.

qp.post.perc <- all.slopes.qp.local %>%
  group_by(Glacier_and_lake) %>%
  summarise(perc_sm_0 = (sum((decadal_trend) <0) / n())*100,
            perc_gt_0 = 100-perc_sm_0,
            perc_sm_0 = round(perc_sm_0, digits = 1),
            perc_gt_0 = round(perc_gt_0, digits = 1),
            perc_sm_0 = paste0(perc_sm_0, "%"),
            perc_gt_0 = paste0(perc_gt_0, "%")) 

# Plot the posterior distributions of local median peak discharge.

qp.local.posteriors <- all.slopes.qp.local %>%
  mutate(Glacier_and_lake = reorder(Glacier_and_lake, decadal_trend, median)) %>%
  ggplot(aes(x = decadal_trend,
             y = Glacier_and_lake,
             fill = stat(x > 0))) +
  stat_halfeye( .width = c( .95), 
                interval_size = 0.8, 
                shape = 21,
                point_color = "black",
                point_size = 1,
                #alpha = 0.6,
                position = position_dodge(width = 0.5)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_fill_manual(values = c("coral4", "coral1")) +
  theme_bw() +
  labs(x = expression(paste("Decadal trend in ", log[10], "-", italic(Q)[p], 
                            " [",  m^{3}, ' ', s^{-1}, ' ', dec^{-1}, "]")),
       y = "")  +
  xlim(c(-0.8, 0.8)) + 
  geom_text(data  = qp.post.perc,
            aes(x = -0.75, y = Glacier_and_lake, label = perc_sm_0),
            size = 2.2,
            colour = "gray10",
            inherit.aes = FALSE,
            parse = FALSE,
            hjust = "left")  +
  geom_text(data  = qp.post.perc,
            aes(x = 0.75, y = Glacier_and_lake, label = perc_gt_0),
            size = 2.2,
            colour = "gray10",
            inherit.aes = FALSE,
            parse = FALSE,
            hjust = "left") +
  theme( axis.text   = element_text(size = 7),
         axis.text.x = element_text(size = 7),
         axis.title  = element_text(size = 7),
         strip.text  = element_text(size = 7),
         legend.position = "none")

###################### Flood volume   ##########################################

# Select only lakes that had more than five outbursts and reported values of V0.

all.glofs.V0.local <- all.glofs.V0 %>%
  group_by(Glacier_and_lake) %>%
  mutate(Glacier_and_lake = str_replace_all(Glacier_and_lake, " ", "_")) %>%
  filter(n() >= 5) %>%
  ungroup()

# Fit a model of local median flood volumes V0.

V0.median.local <-  qr_mod_pred(dat    = all.glofs.V0.local, 
                                resp   = "Mean_Flood_Volume_V0", 
                                pred   = "rounded_year",
                                cond   = "Glacier_and_lake", 
                                q      = 0.5,
                                tstart = 1900,
                                tend   = 2021,
                                levels = unique(all.glofs.V0.local$Glacier_and_lake))

# Save the output to disc.

saveRDS(V0.median.local, "V0_model_median_local.RDS")
# V0.median.local <- readRDS("V0_model_median_local.RDS")

# Approach is exactly the same as above for Qp, so left largely uncommented.

all.preds.V0.local <- V0.median.local$prediction 

conds <- all.glofs.V0.local %>% 
  group_by(Glacier_and_lake) %>% 
  summarise(min_r = min(rounded_year)-5, 
            max_r = max(rounded_year)+5)


# Now truncate the predicted values to period of observed V0.

all.preds.V0.local.sub <- list()

for (i in 1:nrow(conds)) {
  
  all.preds.V0.local.sub[[i]] <- filter(all.preds.V0.local, 
                                        (x_scale == conds$Glacier_and_lake[i]) & 
                                          (rounded_year >= conds$min_r[i]) &  
                                          (rounded_year <= conds$max_r[i]))
}

all.preds.V0.local.sub <- bind_rows(all.preds.V0.local.sub) %>%
  mutate(x_scale = str_replace_all(x_scale, "_", " "))

# Generate a plot of local flood volumes with time.

V0.local.timeseries <- all.preds.V0.local.sub %>%
  ggplot() +
  stat_lineribbon(aes(x = rounded_year, 
                      y = prediction, 
                      color = ordered(quantile),
                      fill = ordered(quantile)), 
                  .width = c(0.95), 
                  alpha = 0.6) + 
  scale_fill_manual( "Quantile", values =  "deepskyblue1") +
  scale_color_manual("Quantile", values =  "deepskyblue1") + 
  geom_point(data = all.glofs.V0.local %>% 
               rename(x_scale = Glacier_and_lake) %>%
               mutate(x_scale = str_replace_all(x_scale, "_", " ")),
             mapping = aes(x = rounded_year,
                           y = Mean_Flood_Volume_V0),
             alpha = 0.66) +
  facet_wrap(~x_scale, scales = "free", ncol = 4) +
  labs(x = "Year" ,
       y = expression(paste("Flood volume ", 
                            italic(V)[0], " [", 10^{6}, " ", m^{3}, ']'))) +
  theme_bw()  +
  theme( axis.text = element_text(size = 6),
         axis.text.x = element_text(size = 6),
         axis.title = element_text(size = 6),
         strip.text = element_text(size = 5),
         legend.position = "none")

# Save plot of local V0 vs. time to disc.

ggsave(filename = "V0_local.pdf", 
       V0.local.timeseries,
       width = 183, 
       height = 200,
       units = "mm")

# Extract the posterior of the trend in median local V0.

all.slopes.V0.local <- V0.median.local$posterior_slopes %>% 
  mutate(posterior_slope = "Posterior slope") %>%
  filter(Glacier_and_lake != "Global") %>%
  mutate(decadal_trend = post_orig*10) %>%
  mutate(Glacier_and_lake = str_replace_all(Glacier_and_lake, " ", "."),
         Glacier_and_lake = str_replace_all(Glacier_and_lake, "_", " "))

# Identify how much of the posterior is greater or smaller than zero, 
# respecitvely.

V0.post.perc <- all.slopes.V0.local %>%
  group_by(Glacier_and_lake) %>%
  summarise(perc_sm_0 = (sum((decadal_trend) <0) / n())*100,
            perc_gt_0 = 100-perc_sm_0,
            perc_sm_0 = round(perc_sm_0, digits = 1),
            perc_gt_0 = round(perc_gt_0, digits = 1),
            perc_sm_0 = paste0(perc_sm_0, "%"),
            perc_gt_0 = paste0(perc_gt_0, "%")) 

# Draw the posterior distribution of the trend for each lake.

V0.local.posteriors <- all.slopes.V0.local %>%
  mutate(Glacier_and_lake = reorder(Glacier_and_lake, decadal_trend, median)) %>%
  ggplot(aes(x = decadal_trend,
             y = Glacier_and_lake,
             fill = stat(x > 0))) +
  stat_halfeye( .width = c( .95), 
                interval_size = 0.8, 
                shape = 21,
                point_color = "black",
                point_size = 1,
                #alpha = 0.6,
                position = position_dodge(width = 0.5)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_fill_manual(values = c("deepskyblue4", "deepskyblue1")) +
  theme_bw() +
  labs(x = expression(paste("Decadal trend in ", log[10], "-", italic(V)[0], " [", 10^{6}, 
                            " ", m^{3}, ' ',  dec^{-1}, "]")),
       y = "")  +
  xlim(c(-0.8, 0.8)) +
  geom_text(data  = V0.post.perc,
            aes(x = -0.75, y = Glacier_and_lake, label = perc_sm_0),
            size = 2.2,
            colour = "gray10",
            inherit.aes = FALSE,
            parse = FALSE,
            hjust = "left")  +
  geom_text(data  = V0.post.perc,
            aes(x = 0.75, y = Glacier_and_lake, label = perc_gt_0),
            size = 2.2,
            colour = "gray10",
            inherit.aes = FALSE,
            parse = FALSE,
            hjust = "left") +
  theme( axis.text   = element_text(size = 7),
         axis.text.x = element_text(size = 7),
         axis.title  = element_text(size = 7),
         strip.text  = element_text(size = 7),
         legend.position = "none")

# Combine the posterior distributions of the trends in local Qp and V0 in one 
# plot. This will go into the Extended data.

local.posteriors <- V0.local.posteriors + qp.local.posteriors + 
  plot_layout(heights = c(1, 1)) + 
  plot_annotation(tag_levels = 'a') & 
  theme(plot.tag = element_text(size = 8, face = "bold"))


ggsave(filename = "local_posterior_trends.pdf", 
       local.posteriors,
       width = 150, 
       height = 200,
       units = "mm")

ggsave(filename = "local_posterior_trends.tiff", 
       local.posteriors,
       width = 150, 
       height = 200,
       units = "mm",
       dpi = 260)


##### FIN! #####################################################################