library(here)
library(tidyverse)
library(rstan)
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

stan.fit <- stan(file = here("analysis/Stan/SS-SR_AR1.stan"), 
                 data = stan.data)



#loop to make stan data and fit? ---------------------------------------------------------