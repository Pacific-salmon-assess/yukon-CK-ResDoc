# functional programming to loop through CU fits and diagnostics
library(here)
library(tidyverse)
library(rstan)
library(bayesplot)
library(shinystan)
source(here("analysis/R/functions.R"))

# load data ------------------------------------------------------------------------------
refit <- FALSE #toggle T/F if you want to refit models

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

stan_data <- unique(sp_har$cu) |> 
  purrr::set_names() |> 
  map(make_stan_data)

# fit AR1 --------------------------------------------------------------------------------
if(refit == TRUE){
  for(i in unique(sp_har$cu)){ #just checking to see
    sp_har1 <- filter(sp_har, cu == i) 
    
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
    
    saveRDS(AR1.fit, here("analysis/data/generated/model_fits/AR1/", 
                          paste0(i, "_AR1.rds")))
  }
}else{
  #make big summary list 
  AR1.fits <- lapply(list.files(here("analysis/data/generated/model_fits/AR1"), 
                                full.names = T), 
                     readRDS)
  names(AR1.fits) <- unique(sp_har$cu)
}

AR1.summary <- AR1.fits |> 
  map(rstan::summary) |> 
  map(\(x) x$summary) |> 
  map(as.data.frame) |> 
  mutate(par = rownames()) |>
  list_rbind(names_to = "cu")

AR1.pars <- AR1.fits |> 
  map(\(x) x |> 
      rstan::extract() |> 
      map(as.data.frame) |> 
      enframe()) |> 
  list_rbind(names_to = "cu")
# diagnostics


# fit TV-alpha (fix in in the single CU one) ---------------------------------------------

# diagnostics
lead_pars <- c("beta", "lnalpha", "sigma_R", "lnresid_0", "phi")

AR1.fits |> 
  map(\(x) mcmc_combo(x, 
      pars = lead_pars,
      combo = c("dens_overlay", "trace"),
      gg_theme = legend_none())) |>
  ggsave()


# age pars
AR1_mods |> 
  map(
    \(x) mcmc_combo(
      x, 
      pars = c("D_scale", "D_sum"),
      combo = c("dens_overlay", "trace"),
      gg_theme = legend_none()
    )
  )

AR1_mods |> 
  map(
    \(x) mcmc_combo(
      x, 
      pars = paste0("Dir_alpha[", 1:4, "]"),
      combo = c("dens_overlay", "trace"),
      gg_theme = legend_none()
    )
  )
# wrangle important outputs and write them