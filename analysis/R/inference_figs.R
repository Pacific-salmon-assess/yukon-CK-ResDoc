library(here)
library(tidyverse)
library(gsl)
library(ggsidekick) #for theme_sleek() - doesn't work with some vs of R, hence the comment

source(here("analysis/R/data_functions.R"))

# read in data ---------------------------------------------------------------------------
# model fits --- 
AR1.fits <- lapply(list.files(here("analysis/data/generated/model_fits/AR1"),
                              full.names = T), readRDS)
names(AR1.fits) <- unique(sp_har$CU)

TVA.fits <- lapply(list.files(here("analysis/data/generated/model_fits/TVA"), 
                              full.names = T), readRDS)
names(TVA.fits) <- unique(sp_har$CU)

# escapement estimates ---
esc <- read.csv(here("analysis/data/raw/esc-data.csv")) |>
  mutate_at(2:6, as.numeric)

# process data and fits to make plots later ----------------------------------------------
bench.par.table <- NULL #empty objects to rbind CU's outputs to 
bench.posts <- NULL
par.posts <- NULL
SR.preds <- NULL
AR1.spwn <- NULL
brood.all <- NULL
a.yrs.all <- NULL
AR1.resids <- NULL
TV.resids <- NULL
TV.SR.preds <- NULL
TV.spwn <- NULL
TV.harv <- NULL


for(i in unique(sp_har$CU)){ # Loop over CUs to process model outputs 
  
  # AR1 (base S-R) models -----------------------------------------------------------------
  sub_dat <- filter(sp_har, CU==i)
  sub_pars <- rstan::extract(AR1.fits[[i]])
  
  AR1.spwn <- rbind(AR1.spwn, bind_cols(t(apply(sub_pars$S, 2, quantile, c(0.25, .5, .75))),
                                        unique(sub_dat$year),
                                        i))
  
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

  par <- matrix(NA,length(sub_pars$lnalpha),3,
                  dimnames = list(seq(1:length(sub_pars$lnalpha)), c("sample","ln_a","beta")))
  
  # get benchmarks & pars (AR1 model)------------------------------------------------------------------
  
  for(j in 1:length(sub_pars$lnalpha)){ 
    ln_a <- sub_pars$lnalpha[j]
    b <- sub_pars$beta[j]
    SR.pred[,j] <- (exp(ln_a)*spw*exp(-b*spw))
    
    bench[j,2] <- get_Smsy(ln_a, b) #S_MSY
    bench[j,1] <- get_Sgen(exp(ln_a),b,-1,1/b*2, bench[j,2]) #S_gen
    bench[j,3] <- (1 - lambert_W0(exp(1 - ln_a))) #U_MSY
    bench[j,4] <- ln_a/b #S_eq
    bench[j,5] <- 1/b #S_msr 
    bench[j,6] <- mean(sub_pars$S[j, (nyrs-4):nyrs]) #S recent - mean spawners in last generation 
    
    par[j,1] <- j
    par[j,2] <- ln_a
    par[j,3] <- b
  }
  
  SR.pred <- as.data.frame(cbind(spw,t(apply(SR.pred, 1, quantile,probs=c(0.1,0.5,0.9), na.rm=T))))|>
    round(2) |>
    mutate(CU = i)
  
  SR.preds <- rbind(SR.preds, SR.pred)
  
  bench[,2] <- bench[,2]*0.8 #make it 80% Smsy
  
  bench.posts <- rbind(bench.posts, as.data.frame(bench) |> mutate(CU = i))
  
  par.posts <- rbind(par.posts, as.data.frame(par) |> mutate(CU = i))
  
  bench.quant <- apply(bench[,1:6], 2, quantile, probs=c(0.1,0.5,0.9), na.rm=T) |>
    t()
  
  mean <- apply(bench[,1:6],2,mean, na.rm=T) #get means of each
  
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
}  # End data wrangling loop by CU 


colnames(SR.preds) <- c("Spawn", "Rec_lwr","Rec_med","Rec_upr", "CU")
colnames(AR1.spwn) <- c("S.25", "S.50", "S.75", "year", "CU")
colnames(AR1.resids) <- c("year","lwr","midlwr","mid","midupr","upr", "CU")
colnames(TV.resids) <- c("year","lwr","midlwr","mid","midupr","upr", "CU")
colnames(TV.spwn) <- c("S.25", "S.50", "S.75", "CU", "year")
colnames(TV.harv) <- c("H.25", "H.50", "H.75", "CU", "year")

CU_order <- c("NorthernYukonR.andtribs.", "Whiteandtribs.", "Stewart",  
              "MiddleYukonR.andtribs.","Pelly", "Nordenskiold", "Big.Salmon", 
              "UpperYukonR.","YukonR.Teslinheadwaters")

SR.preds$CU_f <- factor(SR.preds$CU, levels = CU_order)
AR1.resids$CU_f <- factor(AR1.resids$CU, levels = CU_order)
TV.resids$CU_f <- factor(TV.resids$CU, levels = CU_order)
TV.spwn$CU_f <- factor(TV.spwn$CU, levels = CU_order)
TV.harv$CU_f <- factor(TV.harv$CU, levels = CU_order)
brood.all$CU_f <- factor(brood.all$CU, levels = CU_order)
TV.SR.preds$CU_f <- factor(TV.SR.preds$CU, levels = CU_order)

# write important tables to repo ---------------------------------------------------------
bench.par.table.out <- bench.par.table |>
  relocate(CU, 1) |>
  relocate(bench.par, .after = 1) |>
  relocate(mean, .after = 2) |>
  mutate_at(3:7, ~round(.,5)) |> 
  arrange(bench.par, CU)

write.csv(bench.par.table.out, here("analysis/data/generated/bench_par_table.csv"), 
          row.names = FALSE)

write_rds(bench.posts, here("analysis/data/generated/benchmark_posteriors.rds"))

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
  facet_wrap(~CU_f, scales = "free") +
  scale_colour_viridis_c(name = "Brood Year")+
  labs(x = "Spawners (000s)",
       y = "Recruits (000s)", 
       title = "Spawner-recruit fits") +
  theme_sleek()+
  theme(legend.position = "top",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.key.size = unit(0.4, "cm"),
        legend.title = element_text(size=7),
        legend.text = element_text(size=6))

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
  facet_wrap(~CU_f) +
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
  facet_wrap(~CU_f) +
  theme(legend.position = "none",
        panel.grid = element_blank()) +
  geom_abline(intercept = 0, slope = 0, col = "dark grey", lty = 2)

my.ggsave(here("analysis/plots/TV_resids.PNG"))

# TV alpha ---
ggplot(a.yrs.all
       |> filter(brood_year < 2018), aes(color = CU)) +
  geom_line(aes(x = brood_year , y = mid), lwd = 1.5) +
  scale_color_viridis_d() +
  theme_sleek() +
  geom_hline(yintercept = 1, lty=2, col = "grey") +
  labs(y = "Productivity (maximum R/S)", x = "Brood year", 
       title = "Time-varying productivity")

my.ggsave(here("analysis/plots/changing_productivity.PNG"))

# TV SR fits --- 
ggplot() +
  geom_point(data = brood.all,
             aes(x = S_med/1000,
                 y = R_med/1000),
             size = 1.5) +
  geom_line(data = TV.SR.preds, aes(x = spw/1000, y = pred.R/1000, color = year, group = year)) +
  facet_wrap(~CU_f, scales = "free") +
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
  filter(value <= 15000) #hack to cut big tails of observations, worth looking without this line

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

# escapement plot ----
bench_plot <- bench.par.table |>
  filter(bench.par == "Smsr") |>
  mutate(upper = `50%`,
         lower = `50%`*0.2,
         stock = CU) |>
  select(stock, upper, lower)

esc$CU_f <- factor(esc$stock, levels = CU_order)
bench_plot$CU_f <- factor(bench_plot$stock, levels = CU_order)

ggplot(esc, aes(x = year, y = mean/1000)) + 
  geom_ribbon(aes(ymin = lower/1000, ymax = upper/1000),  fill = "darkgrey", alpha = 0.5) +
  geom_line(lwd = 1.1) +
  xlab("Year") +
  ylab("Spawners (000s)") +
  facet_wrap(~CU_f, ncol=3, scales = "free_y") +
  theme_sleek()  
my.ggsave(here("analysis/plots/cu-escape.PNG"))

# escapement plot all CUs plus aggregate ----
porcupine <- read.csv(here("analysis/data/raw/trib-spwn.csv")) |>
  filter(CU == "Porcupine") |>
  mutate(mean = estimate,
         lwr = NA,
         upr = NA) |>
  select(year, CU, mean,lwr, upr)
  
aggregrate <- read.csv(here("analysis/data/raw/cdn-esc.csv")) |>
  filter(year > 1984) |>
  mutate(stock = "Aggregate",
         mean = Escapement,
         CU = stock) |>
  select(year, CU, mean,lwr, upr)

esc_plus <- esc |>
  mutate(CU = stock,
         lwr = lower,
         upr = upper)|>
  select(year, CU, mean,lwr, upr)

esc_plus <- rbind(esc_plus, porcupine,aggregrate)
esc_plus$CU_f <- factor(esc_plus$CU, levels = c("Porcupine", CU_order, "Aggregate"))
  
ggplot(esc_plus, aes(x = year, y = mean/1000)) + 
  geom_ribbon(aes(ymin = lwr/1000, ymax = upr/1000),  fill = "darkgrey", alpha = 0.5) +
  geom_line(lwd = 1.1) +
  xlab("Year") +
  ylab("Spawners (000s)") +
  facet_wrap(~CU_f, ncol=4, scales = "free_y") +
  scale_y_continuous(limits = c(0, NA)) +
  theme_sleek()  

my.ggsave(here("analysis/plots/cu-agg-escape.PNG"))

# trib vs RR spawner relationships ----
esc_join <- esc |>
  mutate(CU = stock) |>
  select(CU, year, mean, CU_f)

tribs <- read.csv(here("analysis/data/raw/trib-spwn.csv")) |>
  filter(CU != "Porcupine") |>
  mutate(
    estimate = case_when(
      system == "whitehorse" ~ estimate*(1-hatch_contrib),
      .default = estimate)) |>
  unite(tributary, c("system", "type")) |>
  select(!hatch_contrib)

trib_rr <- left_join(tribs,esc_join,by = join_by("CU", "year")) |>
  drop_na()

ggplot(trib_rr, aes(x = mean, y = estimate)) +
  geom_smooth(method="lm", color="grey") +
  geom_point(size=2, color="dark grey")+ 
  xlab("CU spawners") +
  ylab("Tributary spawners") +
  theme_sleek() +
  facet_wrap(~tributary, scales = "free") 

my.ggsave(here("analysis/plots/RR-vs-trib-spawners.PNG"))


# forward simulations --------------------------------------------------------------------

# Use "standard" (non-TV) benchmarks
bench.par.table <- read.csv(here("analysis/data/generated/bench_par_table.csv"))

# Use alpha from which models? TVA = Reference set, recent 5 years; AR1 = robustness set, all yrs
#alpha_type = "TVA"
alpha_type = "AR1"

## Plots for TV alpha models - "Reference set" 
if(alpha_type == "TVA"){
    S.fwd <- read.csv(here("analysis/data/generated/simulations/S_fwd_tv.csv"))
    H.fwd <- read.csv(here("analysis/data/generated/simulations/H_fwd_tv.csv"))
    spwn.fwd <- TVA.spwn
    perf.metrics <- read.csv(here("analysis/data/generated/perf_metrics_tv.csv")) |>
      pivot_longer(-1, names_to = "metric") 
   } else if (alpha_type == "AR1"){
    S.fwd <- read.csv(here("analysis/data/generated/simulations/S_fwd_AR1.csv"))
    H.fwd <- read.csv(here("analysis/data/generated/simulations/H_fwd_AR1.csv"))
    spwn.fwd <- AR1.spwn 
    perf.metrics <- read.csv(here("analysis/data/generated/perf_metrics_AR1.csv")) |>
       pivot_longer(-1, names_to = "metric") 
}


S.fwd$CU_f <- factor(S.fwd$CU, levels = CU_order)
H.fwd$CU_f <- factor(H.fwd$CU, levels = CU_order)
bench.par.table$CU_f <- factor(bench.par.table$CU, levels = CU_order)
spwn.fwd$CU_f <- factor(spwn.fwd$CU, levels = CU_order)

# assign HCRs to groups
HCR_grps <- list(base = c("no.fishing", "fixed.ER.60", "status.quo"),
                 rebuilding = c("no.fishing",
                                unique(S.fwd$HCR[
                                  grepl("*rebuilding*", S.fwd$HCR)])),
                 status.quo = c("no.fishing", 
                                unique(S.fwd$HCR[grepl("*status.quo*", S.fwd$HCR)])),
                 simple = c("no.fishing", "status.quo", "rebuilding"),
                 fixed = unique(S.fwd$HCR[grepl("fixed.ER", S.fwd$HCR)]))

# next step is to make colours consistent
HCR_cols <- c("#B07300", "purple3", "grey25", "#CCA000", "#FEE106",  "#0F8A2E", "#3638A5")
names(HCR_cols) <- unique(S.fwd$HCR)[c(1, 14, 22:26)] # this will break easily

## Spawners projection
                  
for(i in 1:length(HCR_grps[1:4])) { # don't make this fig for all fixed exp rates
  S.fwd %>% filter(HCR %in% HCR_grps[[i]]) %>%
  ggplot() +
    geom_ribbon(aes(ymin = S.25/1000, ymax = S.75/1000, x = year, color=HCR, fill = HCR), 
                alpha = 0.2) +
    #geom_ribbon(data = filter(TV.spwn, year >= max(TV.spwn$year)-7), 
     #           aes(ymin = S.25/1000, ymax = S.75/1000, 
      #              x= year) #offset to return year 
       #         ) +
    #geom_line(data = filter(TV.spwn, year >= max(TV.spwn$year)-7), ##Should line up?
     #         aes(y=S.50/1000, x= year), color = "black") + 
    geom_line(aes(year, S.50/1000, color = HCR), lwd=1) +
    geom_hline(data = filter(bench.par.table, bench.par=="Smsr"), aes(yintercept = mean/1000), 
               color = "forestgreen", lty = 2) +
    geom_hline(data = filter(bench.par.table, bench.par=="Smsr"), aes(yintercept = (mean*0.2)/1000), 
               color = "darkred", lty = 2) +
    facet_wrap(~CU_f, scales = "free_y") +
    scale_x_continuous(expand = expansion(mult = c(0, .01))) +
    labs(title = "Forward simulation spawner trajectory with Smsr (green) and 20% Smsr (red)", 
         y = "Spawners (000s)") +
    theme_sleek() +
    theme(legend.position = "bottom") +
    scale_color_manual(values=HCR_cols, aesthetics = c("fill", "color"))
  
  my.ggsave(here(paste("analysis/plots/S-fwd", alpha_type, names(HCR_grps[i]), "grp", ".PNG", sep="_")))
}

## Harvest projection
for(i in 1:length(HCR_grps[1:4])) { # don't make this fig for all fixed exp rates
  H.fwd %>% filter(HCR %in% HCR_grps[[i]]) %>%
    ggplot() +
    geom_ribbon(aes(ymin = H.25/1000, ymax = H.75/1000, x = year, color=HCR, fill = HCR), 
                alpha = 0.2) +
    #geom_ribbon(data = filter(TV.spwn, year >= max(TV.spwn$year)-7), 
    #           aes(ymin = S.25/1000, ymax = S.75/1000, 
    #              x= year) #offset to return year 
    #         ) +
    #geom_line(data = filter(TV.spwn, year >= max(TV.spwn$year)-7), ##Should line up?
     #         aes(y=H.50/1000, x= year), color = "black") + 
    geom_line(aes(year, H.50/1000, color = HCR), lwd=1) +
    facet_wrap(~CU_f, scales = "free_y") +
    scale_x_continuous(expand = expansion(mult = c(0, .01))) +
    labs(title = "Forward simulation spawner trajectory with Smsr (green) and 20% Smsr (red)", 
         y = "Harvest (000s)") +
    theme_sleek() +
    theme(legend.position = "bottom") +
    scale_color_manual(values=HCR_cols, aesthetics = c("fill", "color"))
  
  my.ggsave(here(paste("analysis/plots/H-fwd", alpha_type, names(HCR_grps[i]), "grp", ".PNG", sep="_")))
}


# alternative forward projection of spawners (shorter time frame, only two scenarios)
ggplot(S.fwd |>
         filter(HCR != "status.quo",
                year < 2051)) +
  geom_ribbon(aes(ymin = S.25/1000, ymax = S.75/1000, x = year, color = HCR, fill = HCR), 
              alpha = 0.2) +
  geom_ribbon(data = filter(TV.spwn, year >= max(TV.spwn$year)-5), 
              aes(ymin = S.25/1000, ymax = S.75/1000, 
                  x= year), #offset to return year 
              fill = "grey", color = "grey") +
  geom_line(data = filter(TV.spwn, year >= max(TV.spwn$year)-5), ##Should line up?
            aes(y=S.50/1000, x= year), color = "black") + 
  geom_line(aes(year, S.50/1000, color = HCR), lwd=1) +
  facet_wrap(~CU_f, scales = "free_y") +
  scale_x_continuous(expand = expansion(mult = c(0, .01))) +
  labs(x = "Year", 
       y = "Spawners (000s)") +
  theme_sleek() +
  theme(legend.position = "bottom") +
  scale_color_viridis_d(aesthetics = c("fill", "color"))
my.ggsave(here(paste0("analysis/plots/S-fwd-bc-alternative_", alpha_type, ".PNG")))



## Spawners vs ER trade-off

spwn_v_ER <- S.fwd %>% filter(HCR %in% HCR_grps[["fixed"]]) %>%
  group_by(HCR, CU_f) %>% 
  summarize(mean_spwn = mean(S.50)) %>%
  mutate(ER = as.numeric(gsub("\\D", "", HCR))) %>% 
  ggplot() +
  geom_point(aes(y=ER, x=mean_spwn/1000, col=CU_f), shape='circle', size=2, alpha=0.7) +
  scale_colour_viridis_d() +
  theme_sleek() +
  theme(legend.position="none") +
  labs(x="Spawners (x1000)", y="Exploitation Rate", col="Conservation Unit")

harv_v_ER <- H.fwd %>% filter(HCR %in% HCR_grps[["fixed"]]) %>%
  group_by(HCR, CU_f) %>% 
  summarize(mean_harv = mean(H.50)) %>%
  mutate(ER = as.numeric(gsub("\\D", "", HCR))) %>% 
  ggplot() +
  geom_point(aes(y=ER, x=mean_harv/1000, col=CU_f), shape='circle', size=2, alpha=0.7) +
  scale_colour_viridis_d() +
  theme_sleek() +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank()) +
  labs(x="Harvest (x1000)", y="Exploitation Rate", col="Conservation Unit")

cowplot::plot_grid(spwn_v_ER, harv_v_ER, nrow=1, rel_widths=c(0.75,1))

my.ggsave(here(paste0("analysis/plots/fixed_ER_tradeoffs", alpha_type, ".PNG")))


# Performance metrics - ALL
perf.plot <- perf.metrics

perf.plot %>% filter(!(HCR %in% HCR_grps[["fixed"]])) %>%
  ggplot(aes(x=factor(HCR, levels=HCRs[1:6]), y = value, fill=factor(HCR, levels=HCRs[1:6]))) + 
  geom_col() +
  scale_fill_manual(values=HCR_cols) +
  facet_wrap(~metric, scales = "free_y") +
  theme_minimal() +
  theme(legend.position = "bottom", 
        axis.text.x = element_blank(), 
        legend.title = element_blank()) +
  guides(fill=guide_legend(nrow=2, byrow=T)) +
  labs(title = "Forward simulation performance metrics", x="", y="") 

my.ggsave(here(paste0("analysis/plots/perf_metrics_all_", alpha_type, ".PNG")))


# Performance metrics - basic
perf.plot <- filter(perf.metrics, metric %in% c("escapement", "ER", "harvest", "harv.stability"))

perf.plot %>% filter(!(HCR %in% HCR_grps[["fixed"]])) %>%
  ggplot(aes(x=HCR, y = value, fill=HCR)) + 
  geom_col() +
  scale_fill_manual(values=HCR_cols) +
  facet_wrap(~metric, scales = "free_y", nrow=1) +
  theme_bw() +
  theme(legend.position = "bottom", 
        axis.text.x = element_blank(), 
        legend.title = element_blank()) +
  labs(title = "Forward simulation performance metrics") 

my.ggsave(here(paste0("analysis/plots/perf_metrics_basic_", alpha_type, ".PNG")))


# Performance metrics - harvest
perf.plot <- filter(perf.metrics, metric %in% c("harvest", "harv.stability", "pr.basic.needs", "pr.no.harvest"))

perf.plot %>% filter(!(HCR %in% HCR_grps[["fixed"]])) %>%
  ggplot(aes(x=HCR, y = value, fill=HCR)) + 
  geom_col() +
  scale_fill_manual(values=HCR_cols) +
  facet_wrap(~metric, scales = "free_y", nrow=1) +
  theme_bw() +
  theme(legend.position = "bottom", 
        axis.text.x = element_blank(), 
        legend.title = element_blank()) +
  labs(title = "Forward simulation harvest performance metrics") 

my.ggsave(here(paste0("analysis/plots/perf_metrics_harv_", alpha_type, ".PNG")))


# Performance status
perf.status <- perf.metrics |>
  filter(metric %in% c("n.above.USR", "n.between.ref", "n.below.LSR", "n.extinct")) |>
  mutate(status = factor(gsub("^n ", "", gsub("\\.", " ", metric)), 
    levels=c("above USR", "between ref", "below LSR", "extinct")))

perf.status %>% filter(!(HCR %in% HCR_grps[["fixed"]])) %>%
  ggplot(aes(x = HCR, y = value, fill = status)) + 
  geom_col() +
  scale_fill_discrete(type = c("forestgreen", "darkorange", "darkred", "black")) +
  scale_y_continuous(breaks = c(2,4,6,8)) +
  labs(y = "Number of CUs", title = "CU status at the end of forward simulation") +
  theme_bw()

my.ggsave(here(paste0("analysis/plots/perf_status_", alpha_type, ".PNG")))


# Performance status - Fixed ER

perf.status %>% filter(HCR %in% HCR_grps[["fixed"]]) %>%
  mutate(ER = as.numeric(gsub("\\D", "", HCR))) %>%
  ggplot(aes(x=value, y=factor(ER), fill=status)) +
  geom_col() +
  scale_fill_discrete(type = c("forestgreen", "darkorange", "darkred", "black")) +
  scale_x_continuous(breaks = c(2,4,6,8)) +
  scale_y_discrete(breaks = seq(0,100,20)) +
  labs(y = "Exploitation Rate", title = "CU status at the end of forward simulation", fill="") +
  theme_bw()

my.ggsave(here(paste0("analysis/plots/fixed_ER_status_", alpha_type, ".PNG")))

# covariance plots ----

library(ggcorrplot)

colnames(sig.R.samps) <- names(AR1.fits)
sig.R.samps.order <- sig.R.samps[,c(4,8,6,2,5,3,1,7,9)]
Sig.R <- cov(sig.R.samps.order) 

colnames(Sig.R) <- colnames(sig.R.samps.order)
rownames(Sig.R) <- colnames(sig.R.samps.order)

ggcorrplot(Sig.R, hc.order = TRUE, type = "lower",
           outline.col = "white",
           lab=TRUE)
my.ggsave(here("analysis/plots/recruit-corr-matrix.PNG"))
