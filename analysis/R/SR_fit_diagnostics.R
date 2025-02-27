# fit SR models for all CUs, run diagnostics, make initial plots
library(here)
library(tidyverse)
library(ggplotify) #for as.ggplot() to help mcmc_combo() plotting
library(rstan)
library(bayesplot)
library(shinystan)
library(gsl) #for lambertw0() to calc U_MSY
source(here("analysis/R/data_functions.R"))

#refit <- FALSE #toggle T/F if you want to refit models
refit <- TRUE

# fit AR1 and time varying productivity (TVA) models--------------------------------------
if(refit == TRUE){
  for(i in unique(sp_har$CU)){
    
    sp_har1 <- filter(sp_har, CU == i) 
    
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
                      "Smax_p_sig" = 0.75*max(sp_har1$spwn))
    
    AR1.fit <- stan(file = here("analysis/Stan/SS-SR_AR1.stan"), 
                    data = stan.data,
                    iter = 4000)
    
    saveRDS(AR1.fit, here("analysis/data/generated/model_fits/AR1/", 
                          paste0(i, "_AR1.rds")))
    
    TV.fit <- stan(file = here("analysis/Stan/SS-SR_TVA.stan"), 
                    data = stan.data,
                    iter = 4000)
    #launch_shinystan(TV.fit)
    saveRDS(TV.fit, here("analysis/data/generated/model_fits/TVA/", 
                          paste0(i, "_TVA.rds")))
  }
}else{
  if(!dir.exists(here("analysis/data/generated/model_fits/AR1")) &
     !dir.exists(here("analysis/data/generated/model_fits/TVA"))){
    stop(print("Make sure you have fit both models at least once!"))
  }
  AR1.fits <- lapply(list.files(here("analysis/data/generated/model_fits/AR1"),
                                full.names = T), 
                     readRDS)
  names(AR1.fits) <- unique(sp_har$CU)[order(unique(sp_har$cu))]
  
  TVA.fits <- lapply(list.files(here("analysis/data/generated/model_fits/TVA"), 
                                full.names = T), 
                     readRDS)
  names(TVA.fits) <- unique(sp_har$CU)[order(unique(sp_har$CU))]
}

# describe diagnostics in loop -----------------------------------------------------------
min_ESS_AR1 <- NULL
max_Rhat_AR1 <- NULL
min_ESS_TVA <- NULL
max_Rhat_TVA <- NULL
for(i in unique(sp_har$CU)){
  sub_dat <- filter(sp_har, CU==i)
  
  #diagnostics associated w/ AR1 model---
  sub_AR1_summary <- as.data.frame(rstan::summary(AR1.fits[[i]])$summary)
  
  sub_AR1_pars <- rstan::extract(AR1.fits[[i]])
  
  min_ESS_AR1 <- rbind(min_ESS_AR1, data.frame(CU=i, 
                                       ESS=round(min(sub_AR1_summary$n_eff, na.rm=T))))
  
  max_Rhat_AR1 <- rbind(max_Rhat_AR1, data.frame(CU=i, 
                                         Rhat=round(max(sub_AR1_summary$Rhat, na.rm=T), 3)))
  
  R <- (sub_dat$harv+sub_dat$spwn)
  R_rep <- sub_AR1_pars$H_rep[1:500,] + sub_AR1_pars$S_rep[1:500,]
  
  ppc_dens_overlay(R, R_rep) +
    xlim(NA, quantile(R_rep, 0.99)) +
    theme(legend.position = "none") +
    labs(y = "density", x = "y_est", title = paste(i, "posterior predictive check"))
  
  my.ggsave(here("analysis/plots/diagnostics/AR1", paste0("PPC_", i, ".PNG")))
  
  p <- mcmc_combo(AR1.fits[[i]], pars = c("beta", "lnalpha", "sigma_R", "phi"), #for AR1
             combo = c("dens_overlay", "trace"),
             gg_theme = legend_none()) |> 
    as.ggplot() +
    labs(title = paste(i, "leading parameters"))
  my.ggsave(here("analysis/plots/diagnostics/AR1", paste0("mcmc_leading_", i, ".PNG")), plot = as.ggplot(p))
  
  p <- mcmc_combo(AR1.fits[[i]], pars = c("D_scale", "D_sum"),
             combo = c("dens_overlay", "trace"),
             gg_theme = legend_none())|> 
    as.ggplot() +
    labs(title = paste(i, "age pars"))
  
  my.ggsave(here("analysis/plots/diagnostics/AR1", paste0("mcmc_age_par_", i, ".PNG")), p)
  
  p <- mcmc_combo(AR1.fits[[i]], pars = c("Dir_alpha[1]", "Dir_alpha[2]", 
                               "Dir_alpha[3]", "Dir_alpha[4]"),
             combo = c("dens_overlay", "trace"),
             gg_theme = legend_none())|> 
    as.ggplot() +
    labs(title = paste(i, "age probs"))
  
  my.ggsave(here("analysis/plots/diagnostics/AR1", paste0("mcmc_ages_", i, ".PNG")), p)
  
  #diagnostics associated w/ TVA model---
  sub_TVA_summary <- as.data.frame(rstan::summary(TVA.fits[[i]])$summary)
  
  sub_TVA_pars <- rstan::extract(TVA.fits[[i]])
  
  min_ESS_TVA <- rbind(min_ESS_TVA, data.frame(CU=i, 
                                               ESS=round(min(sub_TVA_summary$n_eff, na.rm=T))))
  
  max_Rhat_TVA <- rbind(max_Rhat_TVA, data.frame(CU=i, 
                                                 Rhat=round(max(sub_TVA_summary$Rhat, na.rm=T), 3)))
  
  R <- (sub_dat$harv+sub_dat$spwn)
  R_rep <- sub_TVA_pars$H_rep[1:500,] + sub_TVA_pars$S_rep[1:500,]
  
  ppc_dens_overlay(R, R_rep) +
    xlim(NA, quantile(R_rep, 0.99)) +
    theme(legend.position = "none") +
    labs(y = "density", x = "y_est", title = paste(i, "posterior predictive check"))
  
  my.ggsave(here("analysis/plots/diagnostics/TVA", paste0("PPC_", i, ".PNG")))
  
  p <- mcmc_combo(TVA.fits[[i]], pars = c("beta", "sigma_alpha", "sigma_R", "sigma_tot", "F_rw", 
                                          "mean_ln_R0"), 
                  combo = c("dens_overlay", "trace"),
                  gg_theme = legend_none()) |> 
    as.ggplot() +
    labs(title = paste(i, "leading parameters"))
  my.ggsave(here("analysis/plots/diagnostics/TVA", paste0("mcmc_leading_", i, ".PNG")), plot = as.ggplot(p))
  
  p <- mcmc_combo(TVA.fits[[i]], pars = c("D_scale", "D_sum"),
                  combo = c("dens_overlay", "trace"),
                  gg_theme = legend_none())|> 
    as.ggplot() +
    labs(title = paste(i, "age pars"))
  
  my.ggsave(here("analysis/plots/diagnostics/TVA", paste0("mcmc_age_par_", i, ".PNG")), p)
  
  p <- mcmc_combo(TVA.fits[[i]], pars = c("Dir_alpha[1]", "Dir_alpha[2]", 
                                          "Dir_alpha[3]", "Dir_alpha[4]"),
                  combo = c("dens_overlay", "trace"),
                  gg_theme = legend_none())|> 
    as.ggplot() +
    labs(title = paste(i, "age probs"))
  
  my.ggsave(here("analysis/plots/diagnostics/TVA", paste0("mcmc_ages_", i, ".PNG")), p)
}
