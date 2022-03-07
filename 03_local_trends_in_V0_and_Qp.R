####################################################################################################
#######               Estimate local trends in flood volume and peak discharge         #############
#######                         of ice-dammed lakes with repeat outburst               #############
#######                                                                                #############
#######                                  by Georg Veh                                  #############
#######                                 03 March, 2022                                 #############
####################################################################################################

# Load the following packages, or use install.packages("nameofpackage"), if some of them
# are not pre-installed. In some cases you need to restart your R session.

require(brms)
require(tidybayes)
require(modelr)
require(tidyverse)

# Set YOUR working directory folder where to find all files, necessary to run this script.
# Change the location appropriately.

setwd("C:/Users/local-admin/Desktop/Lake_area_volume/")

# Useful functions

scale_this <- function(x){
  (x - mean(x, na.rm=TRUE)) / sd(x, na.rm=TRUE)
}

# Load the tibble with all reported ice-dammed GLOFs between 1900 and 2021.

all.glofs <- readRDS("all_glofs_tibble.RDS")

################################################################################################
#######     LOCAL TRENDS IN FLOOD VOLUME V0 WITH TIME     ######################################

# Load the piece-wise regression model to predict lake volumes for lakes with known lake area.

fit.va <- readRDS("va_model.RDS")

# Find GLOFs that do not have observed flood volumes, but have a mapped lake area before 
# and after the outburst. The table includes a "Certainty_level" that indicates how accurately 
# we were able to map the lake outline. We use only Level 2 data, which represent lake 
# outlines with the highest degree of reliability.

newdata.v0 <- all.glofs %>%
  filter(!is.na(la_after_log),
         !is.na(la_before_log),
         is.na(Mean_Flood_Volume_V0)) %>%
  filter(la_after_log < la_before_log) %>%
  filter((Certainty_level_before == 2) & (Certainty_level_after == 2)) %>%
  pivot_longer(cols = c(la_before_log, la_after_log),
               names_to = "before_or_after",
               values_to = "area_log")

# Predict flood volumes V0 for these lakes from their mapped surface areas.

predicted.V0 <- mcp:::predict.mcpfit(object = fit.va,
                                    newdata = newdata.v0,
                                    summary = F,
                                    samples_format = "matrix")

# Derive the flood volume for each lake: Subtract the predicted volume
# after the GLOF from the predicted volume before the GLOF had happened.
# Then derive percentiles from the posterior distributions.

prctl.predicted.V0 <- sapply(seq(1, ncol(predicted.V0), by = 2), function (i) {
  
  a <- (sort(10^predicted.V0[, i], decreasing = T) - 
          sort(10^predicted.V0[, i+1], decreasing = T) )*10^6

  return(quantile(a, c(0.025, 0.165, 0.25, 0.5, 0.75, 0.835, 0.975)))
  
})

# Add the predictions to the original tibble.

post.pred.V0 <- t(prctl.predicted.V0) %>%
  as_tibble() %>%
  bind_cols(newdata.v0[ seq(1, nrow(newdata.v0), by = 2), ], .)

# We interpret the data between 16.5% and 83.5% of the predicted value range as 
# 1 sigma standard error, given that our response is normally distributed 
# for log-transformed data.

post.pred.V0$Mean_Flood_Volume_V0 <- post.pred.V0$`50%`
post.pred.V0$Min_V0 <- post.pred.V0$`16.5%`
post.pred.V0$Max_V0 <- post.pred.V0$`83.5%`

# Select only specific columns and add a further column expressing that V0 was estimated
# for the V-A scaling relationship

post.pred.V0.sel <- post.pred.V0 %>%
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
  mutate(Mean_Flood_Volume_V0 = Mean_Flood_Volume_V0 / 10^6,
         Min_V0 = Min_V0 / 10^6,
         Max_V0 = Max_V0 / 10^6) %>%
  mutate(method = "V-A")

# Select the same data from the original table, and add a column expressing that 
# V0 is from previous literature.

all.glofs.sel <- all.glofs %>%
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
  mutate(method = "gaged/lit")

# Now combine the tables with estimated V0 from the V-A scaling relationship
# and with values reported in the literature.
# Select only lakes that are ice-dammed and that have at least 5 outbursts between 1900 and 2021.

all.glofs.V0 <- bind_rows(post.pred.V0.sel, all.glofs.sel) %>%
  filter(rounded_year >= 1900) %>%
  group_by(RGI_Glacier_Id, Lake,
           .add = TRUE) %>%
  filter(Lake_type == "ice") %>%
  filter(n() >= 5) %>%
  filter(length(unique(rounded_year)) > 2) %>%
  mutate(Glacier_and_lake = paste0(RGI_Glacier_Id, "_", Lake)) %>%
  ungroup() %>%
  filter(!(Lake == "Lago Argentino")) %>%
  mutate(V0_scale = scale_this(Mean_Flood_Volume_V0),
         year_scale = scale_this(rounded_year)) %>%
  mutate(RGI_Glacier_Id = str_replace_all(RGI_Glacier_Id, "[.]", "_"))

# Write this table to disk.

saveRDS(all.glofs.V0, "all_glofs_V0_tibble.RDS")
# all.glofs.V0 <- readRDS("all_glofs_V0_tibble.RDS") 

# What is the median flood volume according to our database?

q.V0 <- quantile(all.glofs.V0$Mean_Flood_Volume_V0, c(0.025, 0.5, 0.975))
c(q.V0[2], q.V0[2]-q.V0[1], q.V0[3]- q.V0[2])

# Run the hierarchical model of V0 versus time. We set weakly informed priors on the slope
# and intercept and maintain the default priors set in brms.

bprior <- prior(student_t(3, 0, 5), class = "Intercept") + 
  prior(student_t(3, 0, 2.5), class = "b") 

# Note that we use 4 cores - decrease or increase, if necessary.

trend.v0  <- brm(bf(V0_scale ~ year_scale + ( year_scale | Glacier_and_lake)),
                 family = student(),
                 data = all.glofs.V0,
                 prior = bprior,
                 cores = 4,
                 chains = 4,
                 warmup = 2000,
                 iter = 6000,
                 control = list(adapt_delta = 0.94,
                                max_treedepth = 15))

# Write the hierarchical model V0 model to disk.

saveRDS(trend.v0, "local_V0_model.RDS")
# trend.V0 <- readRDS("local_V0_model.RDS")

# Assess model parameters, check for divergences, and do posterior predictive checks.

summary(trend.v0)
plot(trend.v0)
brms::pp_check(trend.v0)

# Obtain the change in GLOF volume with time for each lake.
# We define the range of years for which we want to obtain posterior draws. This period
# is different for each lake.

conds <- all.glofs.V0 %>% 
  group_by(Glacier_and_lake) %>% 
  summarise(min_r = min(year_scale)-0.25, 
            max_r = max(year_scale)+0.25)

# Obtain the posterior distribution of GLOF volumes for each lake and its parent glacier.

preds <- all.glofs.V0 %>%
  data_grid(year_scale = seq_range(year_scale, n = 101, expand = 0.5),
            Glacier_and_lake  = unique(Glacier_and_lake)) %>%
  add_epred_draws(object = trend.v0, 
                  value = "V0", 
                  ndraws = 1000)

# Now truncate the predicted values to period of observed V0.

preds.sub <- list()

for (i in 1:nrow(conds)) {
  
  preds.sub[[i]] <- filter(preds, 
                           (Glacier_and_lake == conds$Glacier_and_lake[i]) & 
                             (year_scale >= conds$min_r[i]) &  
                             (year_scale <= conds$max_r[i]))
}

preds.sub <- bind_rows(preds.sub)

# Plot the posterior trends in V0 for each lake.

plot_trend_year.V0 <- preds.sub %>%
  mutate(rounded_year = year_scale * sd(all.glofs.V0$rounded_year) + mean(all.glofs.V0$rounded_year),
         V0 = V0 * sd(all.glofs.V0$Mean_Flood_Volume_V0) + mean(all.glofs.V0$Mean_Flood_Volume_V0),
         Glacier_and_lake = str_replace_all(Glacier_and_lake, "_", " ")) %>%
  ggplot(aes(x = rounded_year, 
             y = V0)) +
  facet_wrap(~Glacier_and_lake, 
             scales = "free",
             ncol = 4) +
  scale_color_manual(name = "Method",
                     values = c("navy", "orange")) +
  scale_fill_manual(name = "Posterior rate", 
                    values = c("#52c8c8c8")) +
  stat_lineribbon(aes(y = V0), 
                  .width = 0.95,
                  point_interval = mean_qi) +
  geom_pointrange(data = all.glofs.V0 %>%
                    mutate(method = str_replace(method, 
                                                "gaged/lit", 
                                                "Gaged (obtained from literature)"),
                           method = str_replace(method, 
                                                "V-A", 
                                                "Volume-Area relationship"),
                           Glacier_and_lake = str_replace_all(Glacier_and_lake, "_", " ")), 
                  aes(x = rounded_year,
                      y = Mean_Flood_Volume_V0,
                      ymin = Min_V0, 
                      ymax = Max_V0, 
                      color = method, size = 0.5), shape = 16, size = 0.5) +
  labs(x = "Year",
       y = expression(paste("Flood volume ", V[0], " [", 10^{6}, " ", m^{3}, ']'))) +
  theme_bw()  +
  theme( axis.text = element_text(size = 8.5),
         axis.text.x = element_text(size = 8.5),
         axis.title = element_text(size = 9),
         strip.text = element_text(size = 6.5),
         legend.position = "bottom")

# Obtain the posterior 95% highest density interval (HDI) for the trends in
# V0 for all lakes. Note that these refer to the standardised predictor year.
# Therefore, convert the posterior slope to the original scale.

slopes.V0 <- trend.v0 %>%
  spread_draws(b_year_scale, `r_Glacier_and_lake`[Glacier, param]) %>%
  filter(param == "year_scale") %>%
  mutate(mean_V0_scaled = b_year_scale + r_Glacier_and_lake,
         mean_V0 = mean_V0_scaled * sd(all.glofs.V0$Mean_Flood_Volume_V0) / sd(all.glofs.V0$rounded_year) ) %>%
  ungroup() %>%
  mutate(Region = str_sub(Glacier, 1, 8)) %>%
  rename(Glacier_and_lake = Glacier) %>%
  mutate(Glacier_and_lake = str_replace_all(Glacier_and_lake, "[.]", " "),
         Glacier_and_lake = str_replace_all(Glacier_and_lake, "_", " "),
         RGI = str_sub(Glacier_and_lake, 7, 14),
         Lake = str_sub(Glacier_and_lake, 16),
         Lake = str_replace_all(Lake, "[.]", " "),
         type = "V0 trend") 

# Show the posterior the median and 95% HDI of the posterior regression slope.

slopes.V0 %>%
  group_by(Glacier_and_lake) %>%
  summarise(q025 = quantile(mean_V0, 0.025),
            q50 = quantile(mean_V0, 0.5),
            q975 = quantile(mean_V0, 0.975)) %>%
  View()

# Plot the posterior regression slope of V0 versus time.

mod_param.V0 <- slopes.V0 %>%
  mutate(Glacier_and_lake = reorder(Glacier_and_lake, mean_V0, median)) %>%
  ggplot(aes(x = mean_V0,
             y = Glacier_and_lake)) +
  stat_pointinterval( .width = c( .95), 
                      color = "navy",
                      point_size = 1.5) + 
  scale_color_brewer("Posterior\ninterval") + 
  theme_bw() +
  labs(x = expression(paste("Regression slope of ", V[0], " [", 10^{6}, ' ', m^{3}, ' ', yr^{-1}, "] versus year")),
       y = "Glacier and lake")  +
  geom_vline(xintercept = 0) +
  theme( axis.text = element_text(size = 8.5),
         axis.text.x = element_text(size = 8.5),
         axis.title = element_text(size = 9),
         strip.text = element_text(size = 8.5))

# Write the plots to disk, as pngs and PDFs.

# Extended Data Figure 5.

ggsave(
  filename = "regression_V0_per_lake.pdf",
  plot = plot_trend_year.V0,
  width = 180,
  height = 240,
  units = "mm"
)

ggsave(
  filename = "regression_V0_per_lake.png",
  plot = plot_trend_year.V0,
  width = 180,
  height = 240,
  units = "mm"
)

# Extended Data Figure 6.

ggsave(
  filename = "post_trend_V0_per_lake.pdf",
  plot = mod_param.V0,
  width = 180,
  height = 200,
  units = "mm"
)

ggsave(
  filename = "post_trend_V0_per_lake.png",
  plot = mod_param.V0,
  width = 180,
  height = 200,
  units = "mm"
)


################################################################################################
#######     LOCAL TRENDS IN PEAK DISCHARGE Qp WITH TIME     ####################################

# Select only lakes that are ice-dammed and that have at least 5 outbursts between 1900 and 2021.

all.glofs.qp <- all.glofs %>%
  filter(!is.na(RGI_Glacier_Id)) %>%
  filter(!is.na(Peak_discharge_Qp)) %>%
  dplyr::select(RGI_Glacier_Id,
                Date,
                rounded_year,
                Peak_discharge_Qp,
                Min_Qp,
                Max_Qp,
                Major_RGI_Region,
                Longitude,
                Latitude,
                Lake,
                Lake_type,
                region) %>%
  group_by(RGI_Glacier_Id, Lake,
           .add = TRUE) %>%
  filter(Lake_type == "ice") %>%
  filter(n() >= 5) %>%
  filter(length(unique(rounded_year)) > 2) %>%
  mutate(Glacier_and_lake = paste0(RGI_Glacier_Id, "_", Lake)) %>%
  ungroup() %>%
  filter(!(Lake == "Gornersee" & Peak_discharge_Qp > 100) ) %>%
  mutate(qp_scale = scale_this(Peak_discharge_Qp),
         year_scale = scale_this(rounded_year)) %>%
  mutate(RGI_Glacier_Id = str_replace_all(RGI_Glacier_Id, "[.]", "_"))

# Write the table of repeat GLOFs with reported Qp to disk.

saveRDS(all.glofs.qp, "all_glofs_qp_tibble.RDS")
# all.glofs.qp <- readRDS("all_glofs_qp_tibble.RDS")

# What is the median reported peak discharge by GLOFs?

q.qp <- quantile(all.glofs.qp$Peak_discharge_Qp, c(0.025, 0.5, 0.975))
c(q.qp[2], q.qp[2]-q.qp[1], q.qp[3]- q.qp[2])

# Set weakly informed priors on the slope and intercept. Maintain the
# default priors in brms for the other parameter

bprior <- prior(student_t(3, 0, 5), class = "Intercept") + 
          prior(student_t(3, 0, 2.5), class = "b") 

# Run the hierarchical regression model to estimate the change in peak discharge
# with time for all ice-dammed lakes with repeat GLOFs.

trend.qp  <- brm(bf(qp_scale ~ year_scale + ( year_scale | Glacier_and_lake)),
                 family = student(),
                 data = all.glofs.qp,
                 prior = bprior,
                 cores = 4,
                 chains = 4,
                 warmup = 2000,
                 iter = 6000,
                 control = list(adapt_delta = 0.94,
                                max_treedepth = 12))

# Write the hierarchical model Qp model to disk.

saveRDS(trend.qp, "local_qp_model.RDS")

# Assess model parameters, check for divergences, and do posterior predictive checks.

summary(trend.qp)
plot(trend.qp)
brms::pp_check(trend.qp)

# Obtain the change in GLOF peak discharge with time for each lake.
# We define the range of years for which we want to obtain posterior draws. This period
# is different for each lake.

conds <- all.glofs.qp %>% 
  group_by(Glacier_and_lake) %>% 
  summarise(min_r = min(year_scale)-0.25, 
            max_r = max(year_scale)+0.25)

# Obtain the posterior distribution of the change in Qp with time

preds <- all.glofs.qp %>%
  data_grid(year_scale = seq_range(year_scale, n = 101, expand = 0.5),
            Glacier_and_lake  = unique(Glacier_and_lake)) %>%
  add_epred_draws(object = trend.qp, 
                  value = "Qp", 
                  ndraws =  1000)

# Truncate the range of predicted posterior values to the period, for which
# we have observed values of Qp for a given lake.

preds.sub <- list()

for (i in 1:nrow(conds)) {
  
  preds.sub[[i]] <- filter(preds, 
                           (Glacier_and_lake == conds$Glacier_and_lake[i]) & 
                             (year_scale >= conds$min_r[i]) &  
                             (year_scale <= conds$max_r[i]))
}

preds.sub <- bind_rows(preds.sub) %>% 
              ungroup()

# Plot the posterior trend in peak discharge for each lake that produced at least 5 outbursts.

plot_trend_year.qp <- preds.sub %>%
  mutate(rounded_year = (year_scale * sd(all.glofs.qp$rounded_year)) + mean(all.glofs.qp$rounded_year),
         Qp = (Qp * sd(all.glofs.qp$Peak_discharge_Qp)) + mean(all.glofs.qp$Peak_discharge_Qp),
         Glacier_and_lake = str_replace_all(Glacier_and_lake, "_", " ")) %>%
  ggplot(aes(x = rounded_year, 
             y = Qp)) +
  facet_wrap(~Glacier_and_lake, 
             scales = "free",
             ncol = 4) +
  scale_fill_manual(name = "Posterior rate", 
                    values = c("#52c8c8c8")) +
  stat_lineribbon(aes(y = Qp), 
                  .width = 0.95,
                  point_interval = mean_qi) +
  geom_pointrange(data = all.glofs.qp %>%
                    mutate(Glacier_and_lake = str_replace_all(Glacier_and_lake, "_", " ")), 
                  aes(x = rounded_year,
                      y = Peak_discharge_Qp,
                      ymin = Min_Qp, 
                      ymax = Max_Qp,
                      size = 0.5), shape = 16, size = 0.5) +
  labs(x = "Year",
       y = expression(paste("Peak discharge ", Q[p], " [",  m^{3}, " ", s^{-1}, ' ', yr^{-1},  ']'))) +
  theme_bw()  +
  theme( axis.text = element_text(size = 8.5),
         axis.text.x = element_text(size = 8.5),
         axis.title = element_text(size = 9),
         strip.text = element_text(size = 6.5),
         legend.position = "bottom")

# Obtain the posterior 95% highest density interval (HDI) for the trends in
# Qp for all regions. Note that these refer to the standardised predictor year,
# so we need to re-transform this parameter to the original scale.

slopes.qp <- trend.qp %>%
  spread_draws(b_year_scale, `r_Glacier_and_lake`[Glacier, param]) %>%
  filter(param == "year_scale") %>%
  mutate(mean_qp_scaled = b_year_scale + r_Glacier_and_lake,
         mean_qp = mean_qp_scaled * sd(all.glofs.qp$Peak_discharge_Qp) / sd(all.glofs.qp$rounded_year) ) %>%
  ungroup() %>%
  mutate(Region = str_sub(Glacier, 1, 8)) %>%
  rename(Glacier_and_lake = Glacier) %>%
  mutate(Glacier_and_lake = str_replace_all(Glacier_and_lake, "[.]", " "),
         Glacier_and_lake = str_replace_all(Glacier_and_lake, "_", " "),
         RGI = str_sub(Glacier_and_lake, 7, 14),
         Lake = str_sub(Glacier_and_lake, 16),
         Lake = str_replace_all(Lake, "[.]", " "),
         type = "QP trend") 

# Show the posterior regression slope of Qp versus time.

slopes.qp %>%
  group_by(Glacier_and_lake) %>%
  summarise(q025 = quantile(mean_qp, 0.025),
            q50 = quantile(mean_qp, 0.5),
            q975 = quantile(mean_qp, 0.975)) %>%
  View()

# Plot the posterior regression slope.

mod_param.qp <- slopes.qp %>%
  mutate(Glacier_and_lake = reorder(Glacier_and_lake, mean_qp, median)) %>%
  ggplot(aes(x = mean_qp,
             y = Glacier_and_lake)) +
  stat_pointinterval( .width = .95, 
                      color = "navy",
                      point_size = 1.5) +
  geom_point(data = med_trends.qp, 
             mapping = aes(x = med_trend, y = Glacier_and_lake)) +
  theme_bw() +
  labs(x = expression(paste("Regression slope of ", Q[p], " [",  m^{3}, ' ', s^{-1}, ' ', yr^{-1}, '] versus year')),
       y = "Glacier") +
  geom_vline(xintercept = 0) +
  theme( axis.text = element_text(size = 8.5),
         axis.text.x = element_text(size = 8.5),
         axis.title = element_text(size = 9),
         strip.text = element_text(size = 8.5))

# Write the plots to disk, as PDF and png.

# Extended Data Figure 3.

ggsave(
  filename = "regression_Qp_per_lake.pdf",
  plot = plot_trend_year.qp,
  width = 180,
  height = 240,
  units = "mm"
)

ggsave(
  filename = "regression_Qp_per_lake.png",
  plot = plot_trend_year.qp,
  width = 180,
  height = 240,
  units = "mm"
)

# Extended Data Figure 4.

ggsave(
  filename = "post_trend_Qp_per_lake.pdf",
  plot = mod_param.qp,
  width = 180,
  height = 200,
  units = "mm"
)

ggsave(
  filename = "post_trend_Qp_per_lake.png",
  plot = mod_param.qp,
  width = 180,
  height = 200,
  units = "mm"
)

##############################################################################################