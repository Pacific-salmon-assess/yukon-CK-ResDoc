# script to calculate reference points for demographic models
library(here)
library(tidyverse)
library(rstan)

source(here("analysis/R/data_functions.R"))

male_S_comps <- read.csv(here("analysis/data/raw/male_age_comps.csv"))[,3:6] # male spawner age composition
fem_S_comps <- read.csv(here("analysis/data/raw/female_age_comps.csv"))[,3:6] # female spawner age composition
fem_fec_comps <- read.csv(here("analysis/data/raw/female_length_comps_eggs.csv")) # female spawner eggs by age
fem_egg_mass_comps <- read.csv(here("analysis/data/raw/female_length_comps_egg_mass.csv")) # female spawner eggs by age

AR1.eggs.fits <- lapply(list.files(here("analysis/data/generated/model_fits/AR1_egg_mass"),
                                   full.names = T), 
                        readRDS)
names(AR1.eggs.fits) <- unique(sp_har$CU)[order(unique(sp_har$CU))]

bench_AR1_eggs <- NULL

for(i in unique(sp_har$CU)){
  sub_dat <- filter(sp_har, CU==i)
  
  #load AR1 model with eggs parameters ---
  sub_AR1_eggs_pars <- rstan::extract(AR1.eggs.fits[[i]])
  
  #benchamrks based on early demographics
  bench_early <- matrix(NA,length(sub_AR1_eggs_pars$lnalpha),2,
                  dimnames = list(seq(1:length(sub_AR1_eggs_pars$lnalpha)), c("Smsy","Smsr")))
  
  for(j in 1:length(sub_AR1_eggs_pars$lnalpha)){ 
    ln_a <- sub_AR1_eggs_pars$lnalpha[j]
    beta <- sub_AR1_eggs_pars$beta[j]
    year_index = c(1:5)  
    #year_index = c(1:40)  
    #year_index = c(30:40)  
    eta_as_F <- colMeans(fem_S_comps[year_index,]) # probability of returning by age for females
    eta_as_M <- colMeans(male_S_comps[year_index,]) # probability of returning by age for males
    eta_as <- c(eta_as_F,eta_as_M) # probability of returning by age and sex
    z_as <- c(colMeans(fem_egg_mass_comps[year_index,]),c(0,0,0,0)) #  age/sex class-specific reproductive output
    opt <- optimize(calcCeq, c(0,5), maximum=TRUE)
    bench_early[j,1] <- sum(Seq_as)
    opt <- optimize(calcNeq, c(0,5), maximum=TRUE)
    bench_early[j,2] <- sum(Seq_as)
  }
  bench_early <- as.data.frame(bench_early)
  bench_early$Smsy[bench_early$Smsy <0] <- NA
  bench_early$Smsr[bench_early$Smsr <0] <- NA
  bench_early$period <- "early"
  bench_early$CU <- i

  #benchamrks based on long-term average demographics
  bench_avg <- matrix(NA,length(sub_AR1_eggs_pars$lnalpha),2,
                        dimnames = list(seq(1:length(sub_AR1_eggs_pars$lnalpha)), c("Smsy","Smsr")))
  
  for(j in 1:length(sub_AR1_eggs_pars$lnalpha)){ 
    ln_a <- sub_AR1_eggs_pars$lnalpha[j]
    beta <- sub_AR1_eggs_pars$beta[j]
    #year_index = c(1:10)  
    year_index = c(1:40)  
    #year_index = c(35:40)  
    eta_as_F <- colMeans(fem_S_comps[year_index,]) # probability of returning by age for females
    eta_as_M <- colMeans(male_S_comps[year_index,]) # probability of returning by age for males
    eta_as <- c(eta_as_F,eta_as_M) # probability of returning by age and sex
    z_as <- c(colMeans(fem_egg_mass_comps[year_index,]),c(0,0,0,0)) #  age/sex class-specific reproductive output
    opt <- optimize(calcCeq, c(0,5), maximum=TRUE)
    bench_avg[j,1] <- sum(Seq_as)
    opt <- optimize(calcNeq, c(0,5), maximum=TRUE)
    bench_avg[j,2] <- sum(Seq_as)
  }
  bench_avg <- as.data.frame(bench_avg)
  bench_avg$Smsy[bench_avg$Smsy <0] <- NA
  bench_avg$Smsr[bench_avg$Smsr <0] <- NA
  bench_avg$period <- "avg"
  bench_avg$CU <- i

  #benchamrks based on recent demographics
  bench_recent <- matrix(NA,length(sub_AR1_eggs_pars$lnalpha),2,
                      dimnames = list(seq(1:length(sub_AR1_eggs_pars$lnalpha)), c("Smsy","Smsr")))
  
  for(j in 1:length(sub_AR1_eggs_pars$lnalpha)){ 
    ln_a <- sub_AR1_eggs_pars$lnalpha[j]
    beta <- sub_AR1_eggs_pars$beta[j]
    #year_index = c(1:5)  
    #year_index = c(1:40)  
    year_index = c(35:40)  
    eta_as_F <- colMeans(fem_S_comps[year_index,]) # probability of returning by age for females
    eta_as_M <- colMeans(male_S_comps[year_index,]) # probability of returning by age for males
    eta_as <- c(eta_as_F,eta_as_M) # probability of returning by age and sex
    z_as <- c(colMeans(fem_egg_mass_comps[year_index,]),c(0,0,0,0)) #  age/sex class-specific reproductive output
    opt <- optimize(calcCeq, c(0,5), maximum=TRUE)
    bench_recent[j,1] <- sum(Seq_as)
    opt <- optimize(calcNeq, c(0,5), maximum=TRUE)
    bench_recent[j,2] <- sum(Seq_as)
  }
  bench_recent <- as.data.frame(bench_recent)
  bench_recent$Smsy[bench_recent$Smsy <0] <- NA
  bench_recent$Smsr[bench_recent$Smsr <0] <- NA
  bench_recent$period <- "recent"
  bench_recent$CU <- i
  
  bench_AR1_eggs <- rbind(bench_AR1_eggs, bench_early, bench_avg, bench_recent)
}

  summary_bench_AR1_eggs <- bench_AR1_eggs |>
    group_by(CU,period) |>
    summarize(upper.BM = median(0.8*Smsy, na.rm=T),
              upper.alt = median(Smsr, na.rm=T))
    
  bench_AR1_eggslong <- pivot_longer(bench_AR1_eggs, cols = c(Smsy, Smsr), names_to = "par") |>
    filter(par=="Smsr") |>
    mutate(upper.bm = value)|>
    filter(value <= 35000)

  ggplot(bench_AR1_eggslong|> filter(period !="avg"), aes(value, fill = period, color = period)) +
    geom_density(alpha = 0.8, adjust = 4) +
    geom_vline(xintercept = 1.5) +
    facet_wrap(~CU, scales = "free_y") +
    theme_sleek()+
    xlab("S[MSR]") +
    ylab("") +
    theme(legend.position = "bottom",
          axis.ticks.y = element_blank(), 
          axis.text.y = element_blank()) +
    scale_color_viridis_d(aesthetics = c("fill", "color"))

my.ggsave(here("analysis/plots/demo_bench_compare.PNG"))
  
  
  ggplot(bench_AR1_eggslong, aes(y=value, x=  period, color = period)) +
    geom_violin() +
    facet_wrap(~CU, scales = "free_y") +
    theme(legend.position = "bottom") +
    theme_sleek()
    
  
  ggplot(bench_AR1_eggslong |> filter(period !="avg"), aes(y=value, x=  period, color = period)) +
    geom_boxplot(outlier.shape = NA) +
    facet_wrap(~CU, scales = "free_y") +
    theme(legend.position = "bottom")
 