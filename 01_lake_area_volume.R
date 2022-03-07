####################################################################################################
#######               Predict glacier lake volume from glacier lake area               #############
#######                   also known as the volume-area relationship                   #############
#######                                                                                #############
#######                                  by Georg Veh                                  #############
#######                                 02 March, 2022                                 #############
#######                       more comments added on 07 March, 2022                    #############                        
####################################################################################################

# Load the following packages, or use install.packages("nameofpackage"), if some of them
# are not pre-installed. In some cases you need to restart your R session.

require(openxlsx)
require(mcp)
require(tidybayes)
require(scales)
require(ggplot2)
require(tidyverse)

# Set YOUR working directory folder where to find the input files, necessary to successfully run 
# this script. Change the location appropriately.

setwd("C:/Users/local-admin/Desktop/Lake_area_volume/")

# Useful functions

scale_this <- function(x){
  (x - mean(x, na.rm=TRUE)) / sd(x, na.rm=TRUE)
}

# Open-office spreadsheet with GLOFs per region in separate sheets.

glof.file <- "Global_GLOF_database_2021_12_08.ods"

# Input data for the empirical scaling relationship lake volume V versus lake area A.

va.file <- read.xlsx("lake_area_volume_compiliaton.xlsx") %>%
  as_tibble()

# log10 transform of lake area and lake volume.

va.file <- va.file %>% 
  mutate(vol_log = log10(Bathymetric_volume_mil_m3),
         area_log = log10(Area_m2/10^6))

# Formulate a piecewise regression model, here with two linear segments using the package 'mcp'.
# To learn more about the model syntax, please see: 
# https://lindeloev.github.io/mcp/articles/formulas.html

model <- list(
  vol_log ~  area_log ,  # intercept +  slope of segment 1
  ~ 0  + area_log        # intercept +  slope of segment 2
)

# Setting priors

# The change point is centred on 0.5 km² as reported by Shugar et al. (2020). 
# We set robust priors on the regression intercepts and the slopes. 

priors = list(
  int_1         = "dt( 2, 3, 3)",  # Intercept of the model
  area_log_1    = "dt( 0, 2, 3)",  # Slope of the first segment
  area_log_2    = "dt( 0, 2, 3)",  # Slope of the second segment
  cp_1          = "dnorm( -0.3, 1.5) T(-3, 2.5)", # -0.3 = log10(0.5 km²) 
  sigma_1       = "dnorm(0, 1) T(0, )"
)

# Run the piecewise regression model. We choose a Poisson likelihood, which
# is the default distribution for count data.

fit.va <- mcp(model  = model, 
              data   = va.file, 
              prior  = priors, 
              family = gaussian(), 
              chains = 4, 
              iter   = 20000, 
              cores  = 4, 
              adapt  = 8000)

# Show a summary of the posterior model fit

summary(fit.va)

# Assess model performance using leave-one-out cross-validation.

loo1 <- loo(fit.va)

# We compare the model with one breakpoint against a simpler model without a change point.

model <- list(
  vol_log ~  area_log # Only with intercept and slope, no breakpoint
)

# Setting priors.

priors <- list(
  int_1           = "dt( 2, 3, 3)", # Intercept
  area_log_1    = "dt( 0, 2, 3)",  # slope
  sigma_1       = "dnorm(0, 1) T(0, )"
)

# Running the model without a breakpoint.

fit.va.no.bp <- mcp(model =  model,
                    data =   va.file,
                    prior =  priors,
                    family = gaussian(), 
                    chains = 4, 
                    iter  =  20000, 
                    cores =  4, 
                    adapt =  8000)

# Print the model summary.

summary(fit.va.no.bp)

# Assess ELPD for this model.

loo2 <- loo(fit.va.no.bp)

# Compare the two models.

m1 <- data.frame(print(loo::loo_compare(loo1, loo2), 
                       simplify = F))


#### Plotting

# Obtain the prior distribution of the change point.

prior.cp1.x.freq <- seq(from = -2.5, to = 2.5, length.out = 150)

prior.xy.cp1.freq <- dnorm(x = prior.cp1.x.freq , mean = 0.5, sd = 1.5)

# Rescale the prior on the change point location to be located at the bottom of the plot.

prior.line <- tibble(x = prior.cp1.x.freq ,
                     y = rescale(x = prior.xy.cp1.freq, 
                                 to = c(-6, -5)))

# Obtain the posterior location of the change point.

post_cp_freq <- do.call(rbind, fit.va$mcmc_post)

post.x.orig.freq <- post_cp_freq[, "cp_1"]

density.post.orig.freq <- density(post.x.orig.freq)

# Rescale the posterior on the change point locationto be located at the bottom of the plot.

density.post.rescaled.freq <- rescale(x = density.post.orig.freq$y, to = c(-6, -5))

post.cp.polygon <- tibble(x = density.post.orig.freq$x,
                          y = density.post.rescaled.freq )

q.bp <- unname(round(quantile(post.x.orig.freq, c(0.025, 0.5, 0.975)), digits = 2))
q.bp.up <- paste0("+", q.bp[3] - q.bp[2])
q.bp.low <- paste0("-", q.bp[2] - q.bp[1])

# Obtain the posterior predictive distribution for all data points.

predict.freq <- predict(fit.va, 
                        summary = F, 
                        nsamples = 2000, 
                        samples_format = "tidy")


# Plot the posterior predictive distribution of the V-A model (Extended Data Figure 12).

plot_breakpoint <- predict.freq %>%
  ggplot(aes(x = area_log, y = predict)) +
  stat_lineribbon(aes(y = predict), 
                  point_interval = mean_qi, 
                  color = "navy",
                  fill = "lightblue",
                  .width  = 0.95 ) +
  geom_point(data   = va.file,
             aes(x  = area_log, y = vol_log),
             size   = 1.5, 
             shape  = 21, 
             fill   = "black", 
             color  = "white", 
             stroke = 0.5) +
  theme_bw() +
  labs(x = expression(paste(log[10], "-transformed lake area [", km^{2}, ']')),
       y = expression(paste(log[10], "-transformed lake volume [", 10^{6}, " ", m^{3}, ']'))) +
  geom_line(data = prior.line, 
            aes(x = x, y = y),
            color = "grey60") +
  annotate("text", 
           label =  "Prior CP",
           x = -2.5, y = -5.5, hjust = 0, size = 3, color = "grey60", fontface = 'italic') + 
  geom_polygon(data = post.cp.polygon, aes(x = x, y = y)) +
  annotate("text", 
           label =  "Posterior CP",
           x = 0.7, y = -5.5, hjust = 0, size = 3,  fontface = 'italic') + 
  annotate("text", 
           label =  substitute(paste("Change point: "*  x^y * "/"[z]*" log10 lake area"), 
                               list(x = q.bp[2],
                                    y = q.bp.up, 
                                    z = q.bp.low)),
           x = -4.5, y = 4.5, hjust = 0, size = 3) +
  theme(legend.position = 'none')

# Write plot to disk.

ggsave(
  filename = "va_breakpoint.pdf",
  plot = plot_breakpoint,
  width = 150,
  height = 110,
  units = "mm"
)

ggsave(
  filename = "va_breakpoint.png",
  plot = plot_breakpoint,
  width = 150,
  height = 110,
  units = "mm"
)

# Write the V-A model to disk (for predictions  of lake volume later in the workflow).

saveRDS(fit.va, "va_model.RDS")
