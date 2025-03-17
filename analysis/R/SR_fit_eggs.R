# fit egg based SR models for all CUs, run diagnostics, make initial plots
library(here)
library(tidyverse)
library(ggplotify) #for as.ggplot() to help mcmc_combo() plotting
library(rstan)
library(bayesplot)
library(shinystan)
#library(gsl) #for lambertw0() to calc U_MSY
source(here("analysis/R/data_functions.R"))

#refit <- FALSE #toggle T/F if you want to refit models
refit <- TRUE

# fit AR1 models with eggs as measure of reproductive output --------------------------------------
for(i in unique(sp_har$cu)){
    
    sp_har1 <- filter(sp_har, cu == i) 
    
    a_min <- 4
    a_max <- 7 
    nyrs <- nrow(sp_har1) #number of years of observations
    A <- a_max - a_min + 1 #total age classes
    nRyrs <- nyrs + A - 1 #number of recruitment years: unobserved age classes at ? to predict last year of spawners
    
    stan.data <- list("nyrs" = nyrs,
                      "a_min" = a_min,
                      "a_max" = a_max,
                      "A" = A,
                      "nRyrs" = nRyrs,
                      "A_obs" = A_obs,
                      "S_obs" = sp_har1$spwn,
                      "H_obs" = sp_har1$harv,
                      "S_cv" = sp_har1$spwn_cv,
                      "H_cv" = sp_har1$harv_cv, 
                      "Smax_p" = 0.75*max(sp_har1$spwn), #data for priors in semi_inform models, can tinker based on what assumed Smax is 
                      "Smax_p_sig" = 0.75*max(sp_har1$spwn),
                      "fem_S_comps"=fem_S_comps, 
                      "fem_fec_comps"=fem_fec_comps,
                      "fem_egg_mass_comps"=fem_egg_mass_comps)
    
    #AR1.eggs.fit <- stan(file = here("analysis/Stan/SS-SR_AR1_eggs.stan"), 
    #                data = stan.data,
    #                iter = 2000)
    
    AR1.eggs.fit <- stan(file = here("analysis/Stan/SS-SR_AR1_egg_mass.stan"), 
                         data = stan.data,
                         iter = 2000)
    
    #saveRDS(AR1.eggs.fit, here("analysis/data/generated/model_fits/AR1_eggs/", 
    #                      paste0(i, "_AR1_eggs.rds")))

    saveRDS(AR1.eggs.fit, here("analysis/data/generated/model_fits/AR1_egg_mass/", 
                          paste0(i, "_AR1_eggs.rds")))    
}
