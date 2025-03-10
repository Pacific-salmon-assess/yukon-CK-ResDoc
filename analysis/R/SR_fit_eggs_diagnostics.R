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
if(refit == TRUE){
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
}else{
  if(!dir.exists(here("analysis/data/generated/model_fits/AR1_egg_mass"))){
    stop(print("Make sure you have fit  models at least once!"))
  }
  AR1.eggs.fits <- lapply(list.files(here("analysis/data/generated/model_fits/AR1_egg_mass"),
                                full.names = T), 
                     readRDS)
  names(AR1.eggs.fits) <- unique(sp_har$CU)[order(unique(sp_har$CU))]
  
}

# describe diagnostics in loop -----------------------------------------------------------
min_ESS_AR1_eggs <- NULL
max_Rhat_AR1_eggs <- NULL

for(i in unique(sp_har$cu)){
  sub_dat <- filter(sp_har, cu==i)
  
  #diagnostics associated w/ AR1 model---
  sub_AR1_eggs_summary <- as.data.frame(rstan::summary(AR1.eggs.fits[[i]])$summary)
  
  sub_AR1_eggs_pars <- rstan::extract(AR1.eggs.fits[[i]])
  
  min_ESS_AR1_eggs <- rbind(min_ESS_AR1_eggs, data.frame(CU=i, 
                                               ESS=round(min(sub_AR1_eggs_summary$n_eff, na.rm=T))))
  
  max_Rhat_AR1_eggs <- rbind(max_Rhat_AR1_eggs, data.frame(CU=i, 
                                                 Rhat=round(max(sub_AR1_eggs_summary$Rhat, na.rm=T), 3)))
  
  R <- (sub_dat$harv+sub_dat$spwn)
  R_rep <- sub_AR1_eggs_pars$H_rep[1:500,] + sub_AR1_eggs_pars$S_rep[1:500,]
  
  ppc_dens_overlay(R, R_rep) +
    xlim(NA, quantile(R_rep, 0.99)) +
    theme(legend.position = "none") +
    labs(y = "density", x = "y_est", title = paste(i, "posterior predictive check"))
  
  my.ggsave(here("analysis/plots/diagnostics/AR1_egg_mass", paste0("PPC_", i, ".PNG")))
  
  p <- mcmc_combo(AR1.eggs.fits[[i]], pars = c("beta", "lnalpha", "sigma_R", "phi"), #for AR1
                  combo = c("dens_overlay", "trace"),
                  gg_theme = legend_none()) |> 
    as.ggplot() +
    labs(title = paste(i, "leading parameters"))
  my.ggsave(here("analysis/plots/diagnostics/AR1_egg_mass", paste0("mcmc_leading_", i, ".PNG")), plot = as.ggplot(p))
  
  p <- mcmc_combo(AR1.eggs.fits[[i]], pars = c("D_scale", "D_sum"),
                  combo = c("dens_overlay", "trace"),
                  gg_theme = legend_none())|> 
    as.ggplot() +
    labs(title = paste(i, "age pars"))
  
  my.ggsave(here("analysis/plots/diagnostics/AR1_egg_mass", paste0("mcmc_age_par_", i, ".PNG")), p)
  
  p <- mcmc_combo(AR1.eggs.fits[[i]], pars = c("Dir_alpha[1]", "Dir_alpha[2]", 
                                          "Dir_alpha[3]", "Dir_alpha[4]"),
                  combo = c("dens_overlay", "trace"),
                  gg_theme = legend_none())|> 
    as.ggplot() +
    labs(title = paste(i, "age probs"))
  
  my.ggsave(here("analysis/plots/diagnostics/AR1_egg_mass", paste0("mcmc_ages_", i, ".PNG")), p)
  
}

