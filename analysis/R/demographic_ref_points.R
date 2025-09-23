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
demo.par.summary.out <- NULL

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
  ER.pred <- as.data.frame(cbind(em,t(apply(ER.pred, 1, quantile, probs=c(0.1,0.5,0.9), na.rm=T))))|>
    round(2) |>
    mutate(CU = i)

  ER.preds <- rbind(ER.preds, ER.pred)

  #benchmarks based on early demographics
  bench_early <- matrix(NA,length(sub_AR1_eggs_pars$lnalpha),2,
                  dimnames = list(seq(1:length(sub_AR1_eggs_pars$lnalpha)), c("Smsy","Smsr")))

  for(j in 1:length(sub_AR1_eggs_pars$lnalpha)){
    ln_a <- sub_AR1_eggs_pars$lnalpha[j]
    beta <- sub_AR1_eggs_pars$beta[j]
    year_index = c(1:6)
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

  # parameter summary for research doc
  demo.par.summary <- rstan::summary(AR1.eggs.fits[[i]], pars=c("lnalpha", "beta", "phi"),
                                    probs=c(0.025, 0.5, 0.975))$summary
  demo.par.summary.out <- rbind(demo.par.summary.out, as.data.frame(demo.par.summary) |> mutate(CU = i))

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


bench_egg_mass_long <- pivot_longer(bench_AR1_eggs, cols = c(Smsy, Smsr), names_to = "par") |>
  group_by(CU,period, par) |>
  summarize(median = quantile(value,probs=c(0.5),na.rm=T),
            lower = quantile(value,probs=c(0.025),na.rm=T),
            upper = quantile(value,probs=c(0.975),na.rm=T))

write.csv(bench_egg_mass_long, here("analysis/data/generated/demographic_parameters.csv"),
          row.names = FALSE)

write.csv(demo.par.summary.out, here("analysis/data/generated/demographic_par_summary.csv"))

a.yrs.all$model<-"em"
write.csv(a.yrs.all, here("analysis/data/generated/demographic_TVA.csv"),
          row.names = FALSE)

  summary_bench_AR1_eggs <- bench_AR1_eggs |>
    group_by(CU,period) |>
    summarize(upper.BM = median(0.8*Smsy, na.rm=T),
              upper.alt = median(Smsr, na.rm=T))

  bench_AR1_eggslong <- pivot_longer(bench_AR1_eggs, cols = c(Smsy, Smsr), names_to = "par") |>
    filter(par=="Smsr") |>
    mutate(upper.bm = value)


  bs <- bench_AR1_eggslong |>
    filter(CU == "Big.Salmon",
           value <= 15000)

  n <- bench_AR1_eggslong |>
    filter(CU == "Nordenskiold",
           value <= 15000)

  s <- bench_AR1_eggslong |>
    filter(CU == "Stewart",
           value <= 15000)

  u <- bench_AR1_eggslong |>
    filter(CU == "UpperYukonR.",
           value <= 10000)


  w <- bench_AR1_eggslong |>
    filter(CU == "Whiteandtribs.",
           value <= 15000)

  t <- bench_AR1_eggslong |>
    filter(CU == "YukonR.Teslinheadwaters",
           value <= 20000)

  p <- bench_AR1_eggslong |>
    filter(CU == "Pelly",
           value <= 40000)

  m <- bench_AR1_eggslong |>
    filter(CU == "MiddleYukonR.andtribs.",
           value <= 40000)

  no <- bench_AR1_eggslong |>
    filter(CU == "NorthernYukonR.andtribs.",
           value <= 40000)

  custom.bench <- rbind(bs,n,s,u,w,t,p,m,no) |>
    mutate(CU_f = CU) |>
    left_join(CU_name_lookup, by="CU_f")
  custom.bench$CU_f <- factor(custom.bench$CU_f, levels = CU_order)

  custom.bench|> filter(period !="avg") |>
    mutate(period=str_to_sentence(period)) |>
    ggplot(aes(value, fill = period, color = period)) +
    geom_density(alpha = 0.8, adjust = 4) +
    facet_wrap(~CU_f, scales = "free", labeller=CU_labeller) +
    theme_sleek()+
    xlab(expression(S[MSR])) +
    labs(y="", fill="Period", color="Period") +
    theme(legend.position = "right",
          axis.ticks.y = element_blank(),
          axis.text.y = element_blank()) +
    scale_color_viridis_d(aesthetics = c("fill", "color")) +
    scale_x_continuous(limits = c(0, NA))

my.ggsave(here("analysis/plots/SR-models/demo_bench_compare_early-vs-recent.PNG"))


bench.posts <- readRDS(here("analysis/data/generated/benchmark_posteriors.rds"))

bench_eggs <- bench_AR1_eggslong|> filter(period =="recent") |>
  mutate(unit = "egg mass") |>
  select(unit, CU, par, value)

bench_body <- pivot_longer(bench.posts, cols = c(Smsr), names_to = "par") |>
  mutate(unit = "spawners") |>
  select(unit, CU, par, value)

bench_post <- rbind(bench_eggs, bench_body) |>
  filter(value <= 35000)|>
  mutate(unit=str_to_sentence(unit))


  ggplot(bench_post,aes(value, fill = unit, color = unit)) +
  geom_density(alpha = 0.8, adjust = 4) +
  facet_wrap(~CU, scales = "free_y", labeller=CU_labeller) +
  theme_sleek()+
  xlab(expression(S[MSR])) +
  labs(y="", fill="Period", color="Period") +
  theme(legend.position = "bottom",
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        strip.text = element_text(size=11),
        axis.title = element_text(size=13)) +
  scale_color_viridis_d(aesthetics = c("fill", "color"))+
  labs(color = "Reproductive measure",
       fill = "Reproductive measure")+
  theme(legend.position = c(0.88,0.925),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.key.size = unit(0.55, "cm"),
        legend.title = element_text(size=10, vjust=3),
        legend.text = element_text(size=8.5, angle=0, hjust=0),
        strip.text = element_text(size=10))

my.ggsave(here("analysis/plots/SR-models/demo_bench_compare.PNG"))
my.ggsave(here("csasdown/figure/demo_bench_compare.PNG"), height=7.20, width=8.10)


summary_bench_AR1_eggs <- bench_post |>
  group_by(CU,unit) |>
  summarize(upper.alt = median(value, na.rm=T))


# Demographic recruitment relationships ----

# SR fits ---
ggplot() +
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
  facet_wrap(~CU_f, scales = "free", labeller=CU_labeller) +
  scale_colour_viridis_c(name = "Brood Year")+
  labs(x = "Egg mass (kgs)",
       y = "Recruits (000s)") +
  theme_sleek()+
  theme(legend.position = c(0.94,0.925),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.key.size = unit(0.25, "cm"),
        legend.title = element_text(size=7, vjust=3),
        legend.text = element_text(size=6, angle=0, hjust=0),
        strip.text = element_text(size=10))

my.ggsave(here("analysis/plots/SR-models/EM-R_fits.PNG"))
ggsave(here("csasdown/figure/EM-R_fits.PNG"), height=800*2, width=900*2, units="px", dpi=240)

# Demographic time-varying alpha (recruits/gram egg mass) relationships ----
ggplot(a.yrs.all
       |> filter(brood_year < 2018), aes(color = CU)) +
  geom_line(aes(x = brood_year , y = mid*1000), lwd = 1.5) +
  scale_color_viridis_d() +
  theme_sleek() +
  labs(y = "Productivity (maximum R/EM[kg])", x = "Brood year")

my.ggsave(here("analysis/plots/SR-models/changing_demo_productivity.PNG"))

