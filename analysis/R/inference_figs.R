library(here)
library(tidyverse)
library(gsl)
#library(ggsidekick) #for theme_sleek() - doesn't work with some vs of R, hence the comment

source(here("analysis/R/data_functions.R"))

# read in model fits ---------------------------------------------------------------------
AR1.fits <- lapply(list.files(here("analysis/data/generated/model_fits/AR1"),
                              full.names = T), 
                   readRDS)
names(AR1.fits) <- unique(sp_har$CU)[order(unique(sp_har$CU))]

TVA.fits <- lapply(list.files(here("analysis/data/generated/model_fits/TVA"), 
                              full.names = T), 
                   readRDS)
names(TVA.fits) <- unique(sp_har$CU)[order(unique(sp_har$CU))]


# read in escapement estimates -----------------------------------------------------------

esc <- read.csv(here("analysis/data/raw/esc-data.csv")) 
esc$lwr <- as.numeric(esc$lower)
esc$upr <- as.numeric(esc$upper)
esc$year <- as.numeric(esc$year)

# process data and fits to make plots later ----------------------------------------------
bench.par.table <- NULL #empty objects to rbind CU's outputs to 
bench.posts <- NULL
SR.preds <- NULL
brood.all <- NULL
a.yrs.all <- NULL
AR1.resids <- NULL
TV.resids <- NULL
TV.SR.preds <- NULL
TV.spwn <- NULL
TV.harv <- NULL

for(i in unique(sp_har$CU)){
  sub_dat <- filter(sp_har, CU==i)
  sub_pars <- rstan::extract(AR1.fits[[i]])
  
  #latent states of spawners and recruits---
  spwn.quant <- apply(sub_pars$S, 2, quantile, probs=c(0.1,0.5,0.9))[,1:(nyrs-a_min)]
  rec.quant <- apply(sub_pars$R, 2, quantile, probs=c(0.1,0.5,0.9))[,(A+a_min):nRyrs]
  
  brood_t <- as.data.frame(cbind(sub_dat$year[1:(nyrs-A)],t(spwn.quant), t(rec.quant))) |>
    round(2)
  colnames(brood_t) <- c("BroodYear","S_lwr","S_med","S_upr","R_lwr","R_med","R_upr")
  
  brood_t <- mutate(brood_t, CU = i)
  
  brood.all <- rbind(brood.all, brood_t)
  
  #SR relationship based on full posterior---
  spw <- seq(0,max(brood_t$S_upr),length.out=100)
  SR.pred <- matrix(NA,length(spw), length(sub_pars$lnalpha))
  bench <- matrix(NA,length(sub_pars$lnalpha),6,
                  dimnames = list(seq(1:length(sub_pars$lnalpha)), c("Sgen", "Smsy.80", "Umsy", "Seq", "Smsr", "S.recent")))
  
  for(j in 1:length(sub_pars$lnalpha)){ 
    ln_a <- sub_pars$lnalpha[j]
    b <- sub_pars$beta[j]
    SR.pred[,j] <- (exp(ln_a)*spw*exp(-b*spw))
    
    bench[j,2] <- get_Smsy(ln_a, b) #S_MSY
    bench[j,1] <- get_Sgen(exp(ln_a),b,-1,1/b*2, bench[j,2]) #S_gen
    bench[j,3] <- (1 - lambert_W0(exp(1 - ln_a))) #U_MSY
    bench[j,4] <- ln_a/b #S_eq
    bench[j,5] <- 1/b #S_msr 
    bench[j,6] <- mean(sub_pars$S[j, (nyrs-4):nyrs]) #mean spawners in last generation 
  }
  
  SR.pred <- as.data.frame(cbind(spw,t(apply(SR.pred, 1, quantile,probs=c(0.1,0.5,0.9), na.rm=T))))|>
    round(2) |>
    mutate(CU = i)
  
  SR.preds <- rbind(SR.preds, SR.pred)
  
  # get benchmarks & pars ------------------------------------------------------------------
  bench[,2] <- bench[,2]*0.8 #make it 80% Smsy
  
  bench.posts <- rbind(bench.posts, as.data.frame(bench) |> mutate(CU = i))
  
  bench.quant <- apply(bench[,1:5], 2, quantile, probs=c(0.1,0.5,0.9), na.rm=T) |>
    t()
  
  mean <- apply(bench[,1:5],2,mean, na.rm=T) #get means of each
  
  sub_benchmarks <- cbind(bench.quant, mean) |>
    as.data.frame() |>
    mutate(CU = i) |>
    relocate('50%', 1)
  
  #other pars to report 
  alpha <- quantile(exp(sub_pars$lnalpha), probs = c(.1, .5, .9))
  beta <- quantile(sub_pars$beta, probs = c(.1, .5, .9))
  sigma <- quantile(sub_pars$sigma_R, probs = c(.1, .5, .9))
  phi <- quantile(sub_pars$phi, probs = c(.1, .5, .9))
  
  par.quants <- rbind(alpha, beta, sigma, phi)
  
  #make big table of bench and pars
  par.summary <- as.data.frame(rstan::summary(AR1.fits[[i]])$summary) |>
    select(mean, n_eff, Rhat)
  
  #summarise not other pars...
  par.summary <- filter(par.summary, row.names(par.summary) %in% c('lnalpha', 'beta',
                                                                   'sigma_R', 'phi')) |>
    mutate(CU = i)
  par.summary[1,1] <- exp(par.summary[1,1]) #exp ln_alpha
  
  pars <- cbind(par.quants, par.summary)
  
  sub.bench.par.table <- bind_rows(sub_benchmarks, pars) |>
    mutate(n_eff = round(n_eff, 0),
           Rhat = round(Rhat, 4))
  
  sub.bench.par.table <- mutate(sub.bench.par.table, bench.par = rownames(sub.bench.par.table)) 
  
  bench.par.table <- bind_rows(bench.par.table, sub.bench.par.table)
  
  # then residuals---
  resid.quant <- apply(sub_pars$lnresid, 2, quantile, probs=c(0.1,0.25,0.5,0.75,0.9))[,(A):nRyrs] 
  
  resids <- as.data.frame(cbind(sub_dat$year, t(resid.quant))) |>
    mutate(CU = i)
  colnames(resids) <- c("year","lwr","midlwr","mid","midupr","upr", "CU")
  
  AR1.resids <- rbind(AR1.resids, resids)
  
  #time varying alpha --------------------------------------------------------------------
  sub_pars_TVA <- rstan::extract(TVA.fits[[i]])
  
  a.yrs <- apply(exp(sub_pars_TVA$ln_alpha), 2, quantile, probs=c(0.1,0.5,0.9))
  a.yrs <- as.data.frame(cbind(sub_dat$year, t(a.yrs)))
  
  colnames(a.yrs) <- c("brood_year", "lwr", "mid", "upr")
  
  a.yrs.all <- rbind(a.yrs.all, data.frame(a.yrs, CU = i)) 
  
  #time varying alpha residuals  
  resid.quant <- apply(sub_pars_TVA$lnresid, 2, quantile, 
                       probs=c(0.1,0.25,0.5,0.75,0.9))[,(A):nRyrs]
  
  resids <- as.data.frame(cbind(sub_dat$year, t(resid.quant))) |>
    mutate(CU = i)
  TV.resids <- bind_rows(TV.resids, resids)
  
  # median SR pred fit by year
  for(j in 1:ncol(sub_pars_TVA$ln_alpha)){
    ln_alpha_yr <- sub_pars_TVA$ln_alpha[,j]
    pred <- data.frame(pred.R = (exp(median(ln_alpha_yr))*spw*exp(-median(sub_pars_TVA$beta)*spw))) |>
      mutate(CU = i, 
             year = unique(sp_har$year)[j]) |>
      cbind(spw)
    TV.SR.preds <- rbind(TV.SR.preds, pred)
  }
  
  #storing spawner quantiles for fwd sim plot
  TV.spwn.quant <- data.frame(t(apply(sub_pars_TVA$S, 2, quantile, probs=c(0.25,0.5,0.75)))) |>
    mutate(CU =i, 
           year = unique(sub_dat$year))
  
  TV.spwn <- bind_rows(TV.spwn, TV.spwn.quant)
  
  TV.harv.quant <- data.frame(t(apply(sub_pars_TVA$C, 2, quantile, probs=c(0.25,0.5,0.75)))) |>
    mutate(CU =i, 
           year = unique(sub_dat$year))
  
  TV.harv <- bind_rows(TV.harv, TV.harv.quant)
}
colnames(SR.preds) <- c("Spawn", "Rec_lwr","Rec_med","Rec_upr", "CU")
colnames(AR1.resids) <- c("year","lwr","midlwr","mid","midupr","upr", "CU")
colnames(TV.resids) <- c("year","lwr","midlwr","mid","midupr","upr", "CU")
colnames(TV.spwn) <- c("S.25", "S.50", "S.75", "CU", "year")
colnames(TV.harv) <- c("H.25", "H.50", "H.75", "CU", "year")

# write important tables to repo ---------------------------------------------------------
bench.par.table <- bench.par.table |>
  relocate(CU, 1) |>
  relocate(bench.par, .after = 1) |>
  relocate(mean, .after = 2) |>
  mutate_at(3:7, ~round(.,5)) |>
  arrange(bench.par, mean)

write.csv(bench.par.table, here("analysis/data/generated/bench_par_table.csv"), 
          row.names = FALSE) 

# make key plots for pub -----------------------------------------------------------------
# SR fits ---
ggplot() +
  geom_abline(intercept = 0, slope = 1,col="dark grey") +
  geom_ribbon(data = SR.preds, aes(x = Spawn/1000, ymin = Rec_lwr/1000, ymax = Rec_upr/1000),
              fill = "grey80", alpha=0.5, linetype=2, colour="gray46") +
  geom_errorbar(data = brood.all, aes(x= S_med/1000, y = R_med/1000, 
                                      ymin = R_lwr/1000, ymax = R_upr/1000),
                colour="grey", width=0, linewidth=0.3) +
  geom_errorbarh(data = brood.all, aes(y = R_med/1000, xmin = S_lwr/1000, xmax = S_upr/1000),
                 height=0, colour = "grey", linewidth = 0.3) +
  geom_point(data = brood.all,
             aes(x = S_med/1000,
                 y = R_med/1000,
                 color=BroodYear),
             size = 1.5) +
  geom_line(data = SR.preds, aes(x = Spawn/1000, y = Rec_med/1000)) +
  facet_wrap(~CU, scales = "free") +
  scale_colour_viridis_c(name = "Brood Year")+
  labs(x = "Spawners (000s)",
       y = "Recruits (000s)", 
       title = "Spawner-recruit fits") +
  theme(legend.position = "bottom",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.key.size = unit(0.4, "cm"),
        legend.title = element_text(size=9),
        legend.text = element_text(size=8))

my.ggsave(here("analysis/plots/SR_fits.PNG"))

# AR1 resids --- 
ggplot(AR1.resids, aes(x=year, y = mid)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr),  fill = "darkgrey", alpha = 0.5) +
  geom_ribbon(aes(ymin = midlwr, ymax = midupr),  fill = "black", alpha=0.2) +
  geom_line(lwd = 1.1) +
  coord_cartesian(ylim=c(-2,2)) +
  labs(x = "Return year",
       y = "Recruitment residuals", 
       title = "AR1 recruitment residuals") + 
  facet_wrap(~CU) +
  theme(legend.position = "none",
        panel.grid = element_blank()) +
  geom_abline(intercept = 0, slope = 0, col = "dark grey", lty = 2)

my.ggsave(here("analysis/plots/AR1_resids.PNG"))

# TV resids ---
ggplot(TV.resids, aes(x=year, y = mid)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr),  fill = "darkgrey", alpha = 0.5) +
  geom_ribbon(aes(ymin = midlwr, ymax = midupr),  fill = "black", alpha=0.2) +
  geom_line(lwd = 1.1) +
  labs(x = "Return year",
       y = "Recruitment residuals", 
       title = "Time-varying recruitment residuals") + 
  facet_wrap(~CU) +
  theme(legend.position = "none",
        panel.grid = element_blank()) +
  geom_abline(intercept = 0, slope = 0, col = "dark grey", lty = 2)

my.ggsave(here("analysis/plots/TV_resids.PNG"))

# TV alpha ---
ggplot(a.yrs.all
       |> filter(brood_year < 2018), aes(color = CU)) +
  geom_line(aes(x = brood_year , y = mid), lwd = 2) +
  scale_color_viridis_d() +
  labs(y = "Productivity (Ricker alpha)", x = "Brood year", 
       title = "Time-varying productivity")

my.ggsave(here("analysis/plots/changing_productivity.PNG"))

# TV SR fits --- 
ggplot() +
  geom_point(data = brood.all,
             aes(x = S_med/1000,
                 y = R_med/1000),
             size = 1.5) +
  geom_line(data = TV.SR.preds, aes(x = spw/1000, y = pred.R/1000, color = year, group = year)) +
  facet_wrap(~CU, scales = "free") +
  scale_colour_viridis_c(name = "Brood Year")+
  geom_abline(intercept = 0, slope = 1,col="dark grey") +
  labs(x = "Spawners (000s)",
       y = "Recruits (000s)", 
       title = "Time varying productivity spawner-recruit fits") +
  theme(legend.position = "bottom",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.key.size = unit(0.4, "cm"),
        legend.title = element_text(size=9),
        legend.text = element_text(size=8))

my.ggsave(here("analysis/plots/TV_SR_fits.PNG"))

# "status" plots ---
bench.long <- pivot_longer(bench.posts, cols = c(Sgen, Smsy.80, S.recent), names_to = "par") |>
  select(-Umsy, - Seq) |>
  arrange(CU, par, value) |>
  filter(value <= 10000) #hack to cut off fat tails to help with density visualization, also an IUCN cutoff... 

ggplot(bench.long, aes(value/1000, fill = par, color = par)) +
  geom_density(alpha = 0.3) +
  geom_vline(xintercept = 1.5) +
  facet_wrap(~CU, scales = "free_y") +
  theme(legend.position = "bottom") +
  scale_fill_manual(breaks = c("S.recent", "Sgen", "Smsy.80"),
                    values = c("black", "darkred", "forestgreen"), 
                    aesthetics = c("fill", "color"), 
                    labels = c(expression(italic(S[recent])), expression(italic(S[gen])), 
                               expression(italic(paste("80% ",S)[MSY])))) +
  theme(axis.ticks.y = element_blank(), 
        axis.text.y = element_blank(), 
        legend.title=element_blank()) +
  labs(x = "Spawners (000s)", y = "Posterior density", 
       title = "Recent spawners relative to benchmarks and 1500 cutoff")

# EXPERIMENTAL: benchmarks with Smsr ---
bench.long.Smsr <- pivot_longer(bench.posts, cols = c(Sgen, Smsy.80, Smsr, S.recent), names_to = "par") |>
  select(-Umsy, - Seq) |>
  arrange(CU, par, value) |>
  filter(value <= 15000)

ggplot(bench.long.Smsr, aes(value/1000, fill = par, color = par)) +
  geom_density(alpha = 0.3) +
  geom_vline(xintercept = 1.5) +
  facet_wrap(~CU, scales = "free_y") +
  theme(legend.position = "bottom") +
  theme(axis.ticks.y = element_blank(), 
        axis.text.y = element_blank(), 
        legend.title=element_blank()) +
  labs(x = "Spawners (000s)", y = "Posterior density", 
       title = "Recent spawners relative to benchmarks and 1500 cutoff")

# escapement plot ---

bench_plot <- bench.par.table |>
  filter(bench.par == "Smsr") |>
  mutate(upper = `50%`,
         lower = `50%`*0.2,
         stock = CU) |>
  select(stock, upper, lower)

ggplot(esc, aes(x = year, y = mean/1000)) + 
  geom_ribbon(aes(ymin = lwr/1000, ymax = upr/1000),  fill = "darkgrey", alpha = 0.5) +
  geom_line(lwd = 1.1) +
  xlab("Year") +
  ylab("Spawners (000s)") +
  facet_wrap(~stock, ncol=3, scales = "free_y") +
  geom_hline(data= bench_plot, aes(yintercept=upper/1000), lty=2, lwd=1,col="dark green") +
  geom_hline(data= bench_plot, aes(yintercept=lower/1000), lty=2, lwd=1,col="red") +
  theme_sleek()  
my.ggsave(here("analysis/plots/cu-escape.PNG"))

# forward simulations --------------------------------------------------------------------
S.fwd <- read.csv(here("analysis/data/generated/simulations/S_fwd.csv"))
H.fwd <- read.csv(here("analysis/data/generated/simulations/H_fwd.csv"))

ggplot(S.fwd) +
  geom_ribbon(aes(ymin = S.25/1000, ymax = S.75/1000, x = year, color = HCR, fill = HCR), 
              alpha = 0.2) +
  geom_ribbon(data = filter(TV.spwn, year >= max(TV.spwn$year)-7), 
              aes(ymin = S.25/1000, ymax = S.50/1000, 
                  x= year), #offset to return year 
              fill = "grey", color = "grey") +
  geom_line(data = filter(TV.spwn, year >= max(TV.spwn$year)-7), ##Should line up?
            aes(y=S.50/1000, x= year), color = "black") + 
  geom_line(aes(year, S.50/1000, color = HCR), lwd=1) +
  geom_hline(data = filter(bench.par.table, bench.par=="Smsr"), aes(yintercept = mean/1000), 
             color = "forestgreen", lty = 2) +
  geom_hline(data = filter(bench.par.table, bench.par=="Smsr"), aes(yintercept = (mean*0.2)/1000), 
             color = "darkred", lty = 2) +
  facet_wrap(~CU) +
  scale_x_continuous(expand = expansion(mult = c(0, .01))) +
  labs(title = "Forward simulation spawner trajectory with Smsr (green) and 20% Smsr (red)", 
       y = "Spawners (000s)") +
  theme_bw() +
  theme(legend.position = "bottom") +
  scale_color_viridis_d(aesthetics = c("fill", "color"))

my.ggsave(here("analysis/plots/S-fwd.PNG"))

ggplot(H.fwd) +
  geom_ribbon(aes(ymin = H.25, ymax = H.75, x = year, color = HCR, fill = HCR), 
              alpha = 0.2) +
  geom_ribbon(data = filter(TV.harv, year >= max(TV.harv$year)-7), 
              aes(ymin = H.25, ymax = H.75, x= year), #offset to return year 
              fill = "grey", color = "grey") +
  geom_line(data = filter(TV.harv, year >= max(TV.harv$year)-7), ##Should line up?
            aes(y=H.50, x= year), color = "black") + 
  geom_line(aes(year, H.50, color = HCR), lwd=1) +
  #geom_hline(data = filter(bench.par.table, bench.par=="Umsy"), aes(yintercept = mean), 
  #           color = "forestgreen", lty = 2) +
  facet_wrap(~CU) +
  scale_x_continuous(expand = expansion(mult = c(0, .01))) +
  labs(title = "Forward simulation harvest trajectory", 
       y = "Harvest (number caught)") +
  theme_bw() +
  theme(legend.position = "bottom") +
  scale_color_viridis_d(aesthetics = c("fill", "color"))

my.ggsave(here("analysis/plots/H-fwd.PNG"))
