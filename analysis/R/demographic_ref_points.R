# script to calculate reference points for demographic models
library(here)
library(tidyverse)
library(rstan)
library(ggsidekick)

source(here("analysis/R/data_functions.R"))

male_S_comps <- read.csv(here("analysis/data/raw/male_age_comps.csv"))[,3:6] # male spawner age composition
fem_S_comps <- read.csv(here("analysis/data/raw/female_age_comps.csv"))[,3:6] # female spawner age composition
fem_fec_comps <- read.csv(here("analysis/data/raw/female_length_comps_eggs.csv")) # female spawner eggs by age
fem_egg_mass_comps <- read.csv(here("analysis/data/raw/female_length_comps_egg_mass.csv")) # female spawner eggs by age

AR1.eggs.fits <- lapply(list.files(here("analysis/data/generated/model_fits/AR1_egg_mass"),
                                   full.names = T), 
                        readRDS)
names(AR1.eggs.fits) <- unique(sp_har$CU)[order(unique(sp_har$CU))]

TVA.eggs.fits <- lapply(list.files(here("analysis/data/generated/model_fits/TVA_egg_mass"), 
                              full.names = T), readRDS)
names(TVA.eggs.fits) <- unique(sp_har$CU)

# Demographic reference points ----
bench_AR1_eggs <- NULL
ER.preds <- NULL
brood.all <- NULL
a.yrs.all <- NULL

for(i in unique(sp_har$CU)){
  sub_dat <- filter(sp_har, CU==i)
  
  #load AR1 model with eggs parameters ---
  sub_AR1_eggs_pars <- rstan::extract(AR1.eggs.fits[[i]])

  #egg-recruit relationship
  egg.quant <- apply(exp(sub_AR1_eggs_pars$lnEM), 2, quantile, probs=c(0.1,0.5,0.9))[,1:(nyrs-a_min)]
  rec.quant <- apply(sub_AR1_eggs_pars$R, 2, quantile, probs=c(0.1,0.5,0.9))[,(A+a_min):nRyrs]
  
  brood_t <- as.data.frame(cbind(sub_dat$year[1:(nyrs-A)],t(egg.quant), t(rec.quant))) |>
    round(2)
  colnames(brood_t) <- c("BroodYear","EM_lwr","EM_med","EM_upr","R_lwr","R_med","R_upr")
  
  brood_t <- mutate(brood_t, CU = i)
  
  brood.all <- rbind(brood.all, brood_t)
  
  #SR relationship based on full posterior---
  em <- seq(0,max(brood_t$EM_upr),length.out=100)
  ER.pred <- matrix(NA,length(em), length(sub_AR1_eggs_pars$lnalpha))
  
  for(j in 1:length(sub_AR1_eggs_pars$lnalpha)){ 
    ln_a <- sub_AR1_eggs_pars$lnalpha[j]
    b <- sub_AR1_eggs_pars$beta[j]
    ER.pred[,j] <- (exp(ln_a)*em*exp(-b*em))
  }
  ER.pred <- as.data.frame(cbind(em,t(apply(ER.pred, 1, quantile,probs=c(0.1,0.5,0.9), na.rm=T))))|>
    round(2) |>
    mutate(CU = i)
  
  ER.preds <- rbind(ER.preds, ER.pred)
    
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
  
  #time varying alpha --------------------------------------------------------------------
  sub_pars_egg_TVA <- rstan::extract(TVA.eggs.fits[[i]])
  
  a.yrs <- apply(exp(sub_pars_egg_TVA$ln_alpha), 2, quantile, probs=c(0.1,0.5,0.9))
  a.yrs <- as.data.frame(cbind(sub_dat$year, t(a.yrs)))
  
  colnames(a.yrs) <- c("brood_year", "lwr", "mid", "upr")
  
  a.yrs.all <- rbind(a.yrs.all, data.frame(a.yrs, CU = i)) 
}

colnames(ER.preds) <- c("EM", "Rec_lwr","Rec_med","Rec_upr", "CU")

CU_order <- c("NorthernYukonR.andtribs.", "Whiteandtribs.", "Stewart",  
              "MiddleYukonR.andtribs.","Pelly", "Nordenskiold", "Big.Salmon", 
              "UpperYukonR.","YukonR.Teslinheadwaters")

ER.preds$CU_f <- factor(ER.preds$CU, levels = CU_order)
brood.all$CU_f <- factor(brood.all$CU, levels = CU_order)

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
  
  
# source(inference_figs.R)
bench.posts <- readRDS(here("analysis/data/generated/benchmark_posteriors.rds"))

bench_eggs <- bench_AR1_eggslong|> filter(period =="recent") |>
  mutate(unit = "egg mass") |>
  select(unit, CU, par, value)

bench_body <- pivot_longer(bench.posts, cols = c(Smsr), names_to = "par") |>
  mutate(unit = "spawners") |>
  select(unit, CU, par, value)

bench_post <- rbind(bench_eggs, bench_body) |>
  filter(value <= 35000)

ggplot(bench_post, aes(value, fill = unit, color = unit)) +
  geom_density(alpha = 0.8, adjust = 4) +
  geom_vline(xintercept = 1.5) +
  facet_wrap(~CU, scales = "free_y") +
  theme_sleek()+
  xlab(expression(S[MSR])) +
  ylab("") +
  theme(legend.position = "bottom",
        axis.ticks.y = element_blank(), 
        axis.text.y = element_blank()) +
  scale_color_viridis_d(aesthetics = c("fill", "color"))

my.ggsave(here("analysis/plots/demo_bench_compare_spawn-vs-recent.PNG"))

summary_bench_AR1_eggs <- bench_post |>
  group_by(CU,unit) |>
  summarize(upper.alt = median(value, na.rm=T))


# Demographic recruitment relationships ----

# SR fits ---
ggplot() +
  geom_abline(intercept = 0, slope = 1,col="dark grey") +
  geom_ribbon(data = ER.preds, aes(x = EM/1000, ymin = Rec_lwr/1000, ymax = Rec_upr/1000),
              fill = "grey80", alpha=0.5, linetype=2, colour="gray46") +
  geom_errorbar(data = brood.all, aes(x= EM_med/1000, y = R_med/1000, 
                                      ymin = R_lwr/1000, ymax = R_upr/1000),
                colour="grey", width=0, linewidth=0.3) +
  geom_errorbarh(data = brood.all, aes(y = R_med/1000, xmin = EM_lwr/1000, xmax = EM_upr/1000),
                 height=0, colour = "grey", linewidth = 0.3) +
  geom_point(data = brood.all,
             aes(x = EM_med/1000,
                 y = R_med/1000,
                 color=BroodYear),
             size = 1.5) +
  geom_line(data = ER.preds, aes(x = EM/1000, y = Rec_med/1000)) +
  facet_wrap(~CU_f, scales = "free") +
  scale_colour_viridis_c(name = "Brood Year")+
  labs(x = "Egg mass (kgs)",
       y = "Recruits (000s)") +
  theme_sleek()+
  theme(legend.position = "top",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.key.size = unit(0.4, "cm"),
        legend.title = element_text(size=7),
        legend.text = element_text(size=6))
my.ggsave(here("analysis/plots/EM-R_fits.PNG"))

# Demographic time-varying alpha (recruits/gram egg mass) relationships ----
ggplot(a.yrs.all
       |> filter(brood_year < 2018), aes(color = CU)) +
  geom_line(aes(x = brood_year , y = mid), lwd = 1.5) +
  scale_color_viridis_d() +
  theme_sleek() +
  labs(y = "Productivity (maximum R/EM)", x = "Brood year")

my.ggsave(here("analysis/plots/changing_demo_productivity.PNG"))
