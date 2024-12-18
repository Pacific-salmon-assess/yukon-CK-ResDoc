library(here)
library(tidyverse)
library(rstan)
library(bayesplot)
library(shinystan)

harvest <- read.csv(here("analysis/data/raw/harvest-data.csv")) |>
  dplyr::rename(stock = population, 
         harv_cv = cv)

sp_har <- read.csv(here("analysis/data/raw/esc-data.csv")) |>
  dplyr::rename(spwn = mean, 
         spwn_cv = cv) |>
  select(-obs, - se) |>
  left_join(harvest, by = c("stock", "year")) |>
  dplyr::rename(cu = stock) |>
  mutate(N = spwn+harv)

ages <- read.csv(here("analysis/data/raw/age-data-aggregate.csv"))

A_obs <- ages |>
  select(a4:a7) |>
  as.matrix()

rm(harvest)

# fit to 1 cu ----------------------------------------------------------------------------
sp_har1 <- filter(sp_har, cu == unique(sp_har$cu)[1]) 

a_min <- 4
a_max <- 7 
nyrs <- nrow(sp_har1) #number of spawning years
A <- a_max - a_min + 1 #total age classes
nRyrs <- nyrs + A - 1 #number of recruitment years, i.e. add in unobserved age classes at start to predict 1st year of spawners

stan.data <- list("nyrs" = nyrs,
                  "a_min" = a_min,
                  "a_max" = a_max,
                  "A" = A,
                  "nRyrs" = nyrs + A - 1,
                  "A_obs" = A_obs,
                  "S_obs" = sp_har1$spwn,
                  "H_obs" = sp_har1$harv,
                  "S_cv" = sp_har1$spwn_cv,
                  "H_cv" = sp_har1$harv_cv)

AR1.fit <- stan(file = here("analysis/Stan/SS-SR_AR1.stan"), 
                 data = stan.data)

TV.fit <- stan(file = here("analysis/Stan/SS-SR_TVA.stan"), 
               data = stan.data)

#diagnose/play with priors to squish ESS warnings, generate prior & posterior predictive check 
# basic diagnostics ----------------------------------------------------------------------
AR1.model.summary <- as.data.frame(rstan::summary(AR1.fit)$summary)
model.pars.AR1 <- rstan::extract(AR1.fit)

#plot PPC
R <- (sp_har1$harv+sp_har1$spwn) #recruit data
R_rep <- model.pars.AR1$H_rep[1:500,] + model.pars.AR1$S_rep[1:500,]

ppc_dens_overlay(R, R_rep) +
  #xlim(quantile(R_rep, .01), quantile(R_rep, .99)) + #why are some <0? Seems ok based on vignette, not sure why xmin default =0 tho
  theme(legend.position = "none") +
  labs(y = "density", x = "y_est")


# Ideally n_eff for individual parmeters is >200
min(AR1.model.summary$n_eff, na.rm = TRUE)
ggplot(AR1.model.summary, aes(n_eff)) +
  geom_histogram() +
  geom_vline(xintercept=200, lty=2) +
  labs(y = "frequency",
       x = "ESS") +
  theme_minimal()

# If chains have not mixed well (i.e., the between- and within-chain estimates don't agree),
# R-hat is > 1. Only use the sample if R-hat is less than 1.05.
ggplot(AR1.model.summary, aes(Rhat)) +
  geom_histogram() +
  labs(y = "frequency",
       x =  expression(hat("R"))) +
  theme_minimal()

max(AR1.model.summary$Rhat, na.rm = T)

# check the chains directly
#leading pars
mcmc_combo(AR1.fit, pars = c("beta", "lnalpha", "sigma_R", "phi"),
           combo = c("dens_overlay", "trace"),
           gg_theme = legend_none())

#age pars that were having a hard time
mcmc_combo(AR1.fit, pars = c("D_scale", "D_sum"),
           combo = c("dens_overlay", "trace"),
           gg_theme = legend_none())

mcmc_combo(AR1.fit, pars = c("Dir_alpha[1]", "Dir_alpha[2]", 
                             "Dir_alpha[3]", "Dir_alpha[4]"),
           combo = c("dens_overlay", "trace"),
           gg_theme = legend_none())


# how do correlations in leading parameters look?
pairs(AR1.fit, pars = c("beta", "lnalpha", "sigma_R", "phi"))
