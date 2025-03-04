# fit SR models for all CUs, run diagnostics, make initial plots
library(here)
library(tidyverse)
library(rstan)
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
                    iter = 8000)
    
    saveRDS(AR1.fit, here("analysis/data/generated/model_fits/AR1/", 
                          paste0(i, "_AR1.rds")))
    
    TV.fit <- stan(file = here("analysis/Stan/SS-SR_TVA.stan"), 
                    data = stan.data,
                    iter = 8000)
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
  names(AR1.fits) <- unique(sp_har$CU)[order(unique(sp_har$CU))]
  
  TVA.fits <- lapply(list.files(here("analysis/data/generated/model_fits/TVA"), 
                                full.names = T), 
                     readRDS)
  names(TVA.fits) <- unique(sp_har$CU)[order(unique(sp_har$CU))]
}
