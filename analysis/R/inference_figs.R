library(here)
library(tidyverse)
library(gsl)
library(ggsidekick) #for theme_sleek() - doesn't work with some vs of R, hence the comment
library(ggcorrplot)
library(ggpubr)
library(viridis)
source(here("analysis/R/data_functions.R"))

# read in data ---------------------------------------------------------------------------
# model fits --- 
AR1.fits <- lapply(list.files(here("analysis/data/generated/model_fits/AR1"),
                              full.names = T), readRDS)
names(AR1.fits) <- unique(sp_har$CU)

TVA.fits <- lapply(list.files(here("analysis/data/generated/model_fits/TVA"), 
                              full.names = T), readRDS)
names(TVA.fits) <- unique(sp_har$CU)

load(here("analysis/R/run-reconstructions/fittedMod/rpt.Rdata")) # RR model fit

# escapement estimates ---
esc <- read.csv(here("analysis/data/generated/esc-data.csv")) |>
  mutate_at(2:6, as.numeric)

# process data and fits to make plots later ----------------------------------------------
bench.par.table <- NULL #empty objects to rbind CU's outputs to 
bench.posts <- NULL
par.posts <- NULL
par.posts.tva <- NULL
SR.preds <- NULL
AR1.spwn <- NULL
AR1.harv <- NULL
AR1.resids <- NULL
brood.all <- NULL
brood.all.long <- NULL
a.yrs.all <- NULL
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

  AR1.harv <- rbind(AR1.harv, bind_cols(t(apply(sub_pars$H, 2, quantile, c(0.25, .5, .75))),
                                        unique(sub_dat$year),
                                        i))
  
  #latent states of spawners and recruits---
  spwn.quant <- apply(sub_pars$S, 2, quantile, probs=c(0.1,0.5,0.9))[,1:(nyrs-a_min)]
  rec.quant <- apply(sub_pars$R, 2, quantile, probs=c(0.1,0.5,0.9))[,(A+a_min):nRyrs]

  brood_t <- as.data.frame(cbind(sub_dat$year[1:(nyrs-A)], t(spwn.quant), t(rec.quant))) |>
    round(2)
  colnames(brood_t) <- c("BroodYear","S_lwr","S_med","S_upr","R_lwr","R_med","R_upr")
  
  brood_t <- mutate(brood_t, CU = i)
  
  brood.all <- rbind(brood.all, brood_t)
  
  spwn.quant.long <- apply(sub_pars$S, 2, quantile, probs=c(0.1,0.5,0.9))[,1:(nyrs)]
  rec.quant.long <- apply(sub_pars$R, 2, quantile, probs=c(0.1,0.5,0.9))[,(A+a_min):nRyrs]
  
  brood_t.long <- as.data.frame(cbind(sub_dat$year[1:nyrs],t(spwn.quant.long), rbind(t(rec.quant.long),matrix(NA,4,3)))) |>
    round(2)
  colnames(brood_t.long) <- c("BroodYear","S_lwr","S_med","S_upr","R_lwr","R_med","R_upr")
  
  brood_t.long <- mutate(brood_t.long, CU = i)
  
  brood.all.long <- rbind(brood.all.long, brood_t.long)
  
  #SR relationship based on full posterior---
  spw <- seq(0,max(brood_t$S_upr),length.out=100)
  SR.pred <- matrix(NA,length(spw), length(sub_pars$lnalpha))
  bench <- matrix(NA,length(sub_pars$lnalpha),8,
                  dimnames = list(seq(1:length(sub_pars$lnalpha)), c("Sgen", "Smsy", "Umsy", "Seq", "Smsr", "S.recent","Smsr.20","Smsr.40")))

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
    bench[j,6] <- exp(mean(log(sub_pars$S[j, (nyrs-5):nyrs]))) #S recent - mean spawners in last generation 
    bench[j,7] <- (1/b)*0.2 
    bench[j,8] <- (1/b)*0.4 
    
    par[j,1] <- j
    par[j,2] <- ln_a
    par[j,3] <- b
  }
  
  SR.pred <- as.data.frame(cbind(spw,t(apply(SR.pred, 1, quantile,probs=c(0.1,0.5,0.9), na.rm=T))))|>
    round(2) |>
    mutate(CU = i)
  
  SR.preds <- rbind(SR.preds, SR.pred)
  
  bench.posts <- rbind(bench.posts, as.data.frame(bench) |> mutate(CU = i))
  
  par.posts <- rbind(par.posts, as.data.frame(par) |> mutate(CU = i))
  
  bench.quant <- apply(bench[,1:8], 2, quantile, probs=c(0.1,0.5,0.9), na.rm=T) |>
    t()
  
  mean <- apply(bench[,1:8],2,mean, na.rm=T) #get means of each
  
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
  
  
  
  # TVA (time varying alpha) models ---------------------------------------------
  sub_pars_TVA <- rstan::extract(TVA.fits[[i]])
  par_TVA <- matrix(NA,length(sub_pars_TVA$beta),42,
                    dimnames = list(seq(1:length(sub_pars_TVA$beta)),c("sample", seq(1:40),"beta")))
  par_TVA[,1] <- seq(1:length(sub_pars_TVA$beta))
  par_TVA[,c(2:41)] <- sub_pars_TVA$ln_alpha
  par_TVA[,c(42)] <- sub_pars_TVA$beta
  par.posts.tva <- rbind(par.posts.tva, as.data.frame(par_TVA) |> mutate(CU = i))

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
colnames(AR1.resids) <- c("year","lwr","midlwr","mid","midupr","upr", "CU")
colnames(AR1.spwn) <- c("S.25", "S.50", "S.75", "year", "CU")
colnames(AR1.harv) <- c("H.25", "H.50", "H.75", "year", "CU")
colnames(TV.resids) <- c("year","lwr","midlwr","mid","midupr","upr", "CU")
colnames(TV.spwn) <- c("S.25", "S.50", "S.75", "CU", "year")
colnames(TV.harv) <- c("H.25", "H.50", "H.75", "CU", "year")

CU_order <- c("NorthernYukonR.andtribs.", 
              "Whiteandtribs.", "Stewart",  
              "MiddleYukonR.andtribs.","Pelly",
              "Nordenskiold", "Big.Salmon", 
              "UpperYukonR.","YukonR.Teslinheadwaters")
CU_prettynames <- c("Northern Yukon R. and tribs.", 
                    "White and tribs.", "Stewart",  
                    "Middle Yukon R. and tribs.","Pelly",
                    "Nordenskiold", "Big Salmon", 
                    "Upper Yukon R.","Yukon R. Teslin Headwaters")
CU_name_lookup <- data.frame(CU_f = factor(CU_order, levels=CU_order), 
                             CU_pretty = factor(CU_prettynames, levels=CU_prettynames))

  
SR.preds$CU_f <- factor(SR.preds$CU, levels = CU_order)
AR1.spwn$CU_f <- factor(AR1.spwn$CU, levels = CU_order)
AR1.resids$CU_f <- factor(AR1.resids$CU, levels = CU_order)
AR1.harv$CU_f <- factor(AR1.harv$CU, levels = CU_order)
TV.resids$CU_f <- factor(TV.resids$CU, levels = CU_order)
TV.spwn$CU_f <- factor(TV.spwn$CU, levels = CU_order)
TV.harv$CU_f <- factor(TV.harv$CU, levels = CU_order)
brood.all$CU_f <- factor(brood.all$CU, levels = CU_order)
brood.all.long$CU_f <- factor(brood.all.long$CU, levels = CU_order)
TV.SR.preds$CU_f <- factor(TV.SR.preds$CU, levels = CU_order)
esc$CU_f <- factor(esc$stock, levels = CU_order)

# write important tables to repo ---------------------------------------------------------
bench.par.table.out <- bench.par.table |>
  relocate(CU, 1) |>
  relocate(bench.par, .after = 1) |>
  relocate(mean, .after = 2) |>
  mutate_at(3:7, ~round(.,5)) |> 
  arrange(bench.par, CU)

bench.par.table.out[82:90,3:6] <- bench.par.table.out[82:90,3:6]*10000 # makes beta on scale of 10k so not crazu sig. digits in tables

write.csv(bench.par.table.out, here("analysis/data/generated/bench_par_table.csv"), 
          row.names = FALSE)

write_rds(bench.posts, here("analysis/data/generated/benchmark_posteriors.rds"))

write.csv(par.posts, here("analysis/data/generated/AR1_posteriors.csv"))

write.csv(par.posts.tva, here("analysis/data/generated/TVA_posteriors.csv"))

write.csv(brood.all.long, here("analysis/data/generated/brood_table_long.csv"), 
          row.names = FALSE)

a.yrs.all$model<-"spw"
write.csv(a.yrs.all, here("analysis/data/generated/spw_TVA.csv"), 
          row.names = FALSE)
# make key plots for pub -----------------------------------------------------------------

# SR fits ----
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
  facet_wrap(~CU_f, scales = "free", labeller=CU_labeller) +
  scale_colour_viridis_c(name = "Brood Year")+
  labs(x = "Spawners (000s)",
       y = "Recruits (000s)") +
  theme_sleek()+
  theme(legend.position = c(0.94,0.925),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.key.size = unit(0.25, "cm"),
        legend.title = element_text(size=7, vjust=3),
        legend.text = element_text(size=6, angle=0, hjust=0),
        strip.text = element_text(size=10))

my.ggsave(here("analysis/plots/SR-models/SR_fits_AR1.PNG"))
ggsave(here("csasdown/figure/SR_fits_AR1.PNG"), height=800*2, width=900*2, units="px", dpi=240)


# AR1 resids ---- 
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

my.ggsave(here("analysis/plots/SR-models/AR1_resids.PNG"))

# TV resids ----
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

my.ggsave(here("analysis/plots/SR-models/TV_resids.PNG"))

# TV alpha ----
a.yrs.all |> 
  filter(brood_year < 2018) |>
  left_join(CU_name_lookup, by=c("CU"= "CU_f")) |>
  ggplot(aes(color = CU_pretty)) +
  geom_line(aes(x = brood_year , y = mid), lwd = 1.5) +
  scale_color_viridis_d() +
  theme_sleek() +
  geom_hline(yintercept = 1, lty=2, col = "grey") +
  labs(y ="Productivity (\U03B1)", x = "Brood year")+ 
  guides(color=guide_legend(title="Conservation Unit")) +
  theme(legend.position = c(0.75,0.75),
        plot.margin = margin(60,60,10,60),
        legend.text = element_text(size=7),
        axis.title = element_text(size=8))

my.ggsave(here("analysis/plots/SR-models/changing_productivity.PNG"))
ggsave(here("csasdown/figure/changing_productivity.PNG"), height = 550*2, 
       width = 700*2, units="px", dpi=200)


# TV SR fits ---- 
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

my.ggsave(here("analysis/plots/SR-models/TV_SR_fits.PNG"))

# "status" plots ----
bench.long <- pivot_longer(bench.posts, cols = c(Smsr.20, Smsr.40, S.recent), names_to = "par") |>
  arrange(CU, par, value) |>
  filter(value <= 10000) #hack to cut off fat tails to help with density visualization, also an IUCN cutoff... 

bs <- bench.long |>
  filter(CU == "Big.Salmon",
         value <= 6000)

n <- bench.long |>
  filter(CU == "Nordenskiold",
         value <= 2500)

s <- bench.long |>
  filter(CU == "Stewart",
         value <= 4000)

u <- bench.long |>
  filter(CU == "UpperYukonR.",
         value <= 4000)


w <- bench.long |>
  filter(CU == "Whiteandtribs.",
         value <= 5000)

t <- bench.long |>
  filter(CU == "YukonR.Teslinheadwaters",
         value <= 5000)

p <- bench.long |>
  filter(CU == "Pelly",
         value <= 9000)

m <- bench.long |>
  filter(CU == "MiddleYukonR.andtribs.",
         value <= 10000)

no <- bench.long |>
  filter(CU == "NorthernYukonR.andtribs.",
         value <= 10000)

custom.bench <- rbind(bs,n,s,u,w,t,p,m,no) |>
  mutate(CU_f = CU) |>
  left_join(CU_name_lookup, by="CU_f")

ggplot(custom.bench, aes(value/1000, fill = par, color = par)) +
  geom_density(alpha = 0.3) +
  geom_vline(xintercept = 1.5, lty=2, col="grey") +
  facet_wrap(~CU_pretty, scales = "free") +
  theme(legend.position = "bottom") +
  scale_fill_manual(breaks = c("S.recent", "Smsr.20", "Smsr.40"),
                    values = c("black", "darkred", "forestgreen"), 
                    aesthetics = c("fill", "color"), 
                    labels = c(expression(italic(S[recent])), expression(italic(paste("20% ",S)[MSR])), 
                               expression(italic(paste("40% ",S)[MSR])))) +
  labs(x = "Spawners (000s)", y = "Posterior density") +
  theme_sleek()   +
  theme(axis.ticks.y = element_blank(), 
        axis.text.y = element_blank(), 
        legend.title=element_blank()) +
  scale_x_continuous(limits = c(0, NA))
my.ggsave(here("analysis/plots/SR-models/status.PNG"))

# escapement plot ----
bench_plot <- bench.par.table |>
  filter(bench.par == "Smsr") |>
  mutate(rebuilding = `50%`,
         upper = `50%`*0.4,
         lower = `50%`*0.2,
         stock = CU) |>
  select(stock, rebuilding, upper, lower) |> 
  mutate(CU_f = factor(stock, levels=CU_order)) |>
  left_join(CU_name_lookup, by="CU_f") |>
  left_join(data.frame(summarize(esc, plot_lim=round(max(upper)/1000), .by="stock")),
            by="stock") |> 
  pivot_longer(cols=c("rebuilding", "upper", "lower"), values_to="value", names_to="bench")
esc$CU_f <- factor(esc$stock, levels = CU_order)

esc |>
  left_join(CU_name_lookup, by="CU_f") |>
  ggplot(aes(x = year, y = mean/1000)) + 
  geom_ribbon(aes(ymin = lower/1000, ymax = upper/1000),  fill = "darkgrey", alpha = 0.5) +
  geom_line(lwd = 1.1, col="grey30") +
  xlab("Year") +
  ylab("Spawners (000s)") +
  facet_wrap(~CU_pretty, ncol=3, scales = "free_y") +
  theme_sleek() +
  theme(strip.text = element_text(size=10))

my.ggsave(here("analysis/plots/trib-rr/cu-escape.PNG"), width = 11)
ggsave(here("csasdown/figure/cu-escape.PNG"), width=900*2, height=800*2, units="px",
       dpi=240)

# Add benchmarks
esc |>
  left_join(CU_name_lookup, by="CU_f") |>
  ggplot() + 
  geom_hline(data=bench_plot, aes(yintercept=value/1000, col=label), lty="dashed") +
  scale_colour_manual(values=c("darkred", "pink", "forestgreen"), 
                      labels=c(expression(paste("20%", S[MSR])), expression(paste("40%", S[MSR])), 
                               expression(S[MSR]))) +
   #          col="pink3") +
  #geom_rect(data=bench_plot, aes(xmin=1985, xmax=2025, ymin=0, ymax=lower/1000), fill="darkred", alpha=0.2) +
  #geom_rect(data=bench_plot, aes(xmin=1985, xmax=2025, ymin=lower/1000, ymax=upper/1000),  fill="yellow", alpha=0.2) +
  #geom_rect(data=bench_plot, aes(xmin=1985, xmax=2025, ymin=upper/1000, ymax=plot_lim),  fill="forestgreen", alpha=0.2) +
  geom_ribbon(aes(x=year, ymin = lower/1000, ymax = upper/1000),  fill = "darkgrey", alpha = 0.5) +
  geom_line(aes(x = year, y = mean/1000), lwd = 1.1, col="grey30") +
  labs(x="Year", y="Spawners (000s)", col="") +
  facet_wrap(~CU_pretty, ncol=3, scales = "free_y") +
  theme_sleek() +
  theme(strip.text = element_text(size=10),
        legend.position="right",
        legend.text = element_text(size=9))



ggsave(here("csasdown/figure/cu-escape-bench.PNG"), width=900*2, height=800*2, units="px",
       dpi=240)



# trib vs RR spawner relationships ----
esc_join <- esc |>
  mutate(CU = stock,
         mean = mean/1000) |>
  select(CU, year, mean, CU_f)

tribs <- read.csv(here("analysis/data/raw/trib-spwn.csv")) |>
  filter(CU != "Porcupine")|>
  filter(system != "teslin") |>
  mutate(
    estimate = case_when(
      system == "whitehorse" ~ estimate*(1-hatch_contrib),
      .default = estimate)) |>
  unite(tributary, c("system", "type")) |>
  select(!hatch_contrib)

trib_rr <- left_join(tribs,esc_join,by = join_by("CU", "year")) |>
  drop_na(c("mean", "estimate")) |>
  filter(! tributary %in% c("morley_aerial", "chandindu_weir","nisutlin_sonar", "pelly_aerial", "ross_aerial"))

trib_rr |> 
  mutate(tribs_name = gsub("creek", " Creek", 
                           gsub("salmon", " Salmon", 
                                gsub("_", " ", str_to_sentence(tributary))))) |>
  ggplot(aes(x = mean, y = estimate)) +
  geom_smooth(method="lm", color="grey") +
  geom_point(size=2, color="dark grey")+ 
  xlab("CU spawners (000s)") +
  ylab("Tributary spawners") +
  theme_sleek() +
  facet_wrap(~tribs_name, scales = "free", ncol = 4) +
  theme(axis.title = element_text(size=12))

my.ggsave(here("analysis/plots/trib-rr/RR-vs-trib-spawners.PNG"), width = 11)
ggsave(here("csasdown/figure/RR-vs-trib-spawners.PNG"), width=777*2, height=800*2, units="px",
       dpi=240)

# trib time series ----
tribs.all <- read.csv(here("analysis/data/raw/trib-spwn.csv")) |>
  mutate(
    estimate = case_when(
      system == "whitehorse" ~ estimate*(1-hatch_contrib),
      .default = estimate)) |>
  unite(tributary, c("system_alt", "type"), sep = " ") |>
  select(!hatch_contrib)


dat_text <- tribs.all |>
  group_by(tributary) |>
  slice_head() |>
  select(tributary, CU)

tribs.all |>
  mutate(tribs_name = gsub("salmon", "Salmon", gsub("_", "-", str_to_sentence(tributary)))) |>
  ggplot(aes(x = year, y = estimate/1000)) + 
  geom_line(lwd = 0.8, col="grey") +
  xlab("Year") +
  ylab("Spawners (000s)") +
  facet_wrap(~tribs_name, ncol=4, scales = "free_y") +
  scale_y_continuous(limits = c(0, NA)) +
  theme_sleek() + 
  theme(axis.title = element_text(size=12),
        strip.text = element_text(size=8))

my.ggsave(here("analysis/plots/trib-rr/trib-escape.PNG"), width = 11, dpi= 180)
ggsave(here("csasdown/figure/trib-escape.PNG"), height = 900*2, 
       width=720*2, units="px", dpi= 240)

  
# forward simulations ----

## reference vs robustness productivity ----  
AR1.par.posts <- read.csv(here("analysis/data/generated/AR1_posteriors.csv"))
TVA.par.posts <- read.csv(here("analysis/data/generated/TVA_posteriors.csv"))
  
TV.pp.ref.long <- pivot_longer(TVA.par.posts, cols = c(24:34), names_to = "par") |>
  select(CU, par, value) |>
  mutate(scenario = "reference")       

TV.pp.rob.long <- pivot_longer(TVA.par.posts, cols = c(35), names_to = "par") |>
  select(CU, par, value) |>
  mutate(scenario = "robustness")   

AR.pp.rob <- AR1.par.posts |>
  mutate(scenario = "stationary",
         par= 1,
         value = ln_a) |>
  select(CU, par, value, scenario) 

alpha.posts <- rbind(TV.pp.ref.long, TV.pp.rob.long, AR.pp.rob)

alpha.stats <- alpha.posts |>
  group_by(CU, scenario) |>
  summarise(mean.prod = mean(value)) |>
  filter(scenario != "stationary")

alpha.posts |> filter(scenario != "stationary") |>
  mutate(CU_f = factor(CU, levels=CU_order)) |>
  mutate(scenario = str_to_sentence(scenario)) |>
  ggplot(aes(value, fill = scenario, color = scenario)) +
  geom_density(alpha = 0.3) +
  facet_wrap(~CU_f, scales = "free_y", labeller=CU_labeller) +
  theme_sleek() +
  theme(legend.position = "bottom") +
  geom_vline(xintercept = 0, lty=2, col="grey") +
  scale_colour_grey(aesthetics = c("colour", "fill"),start = 0.3, end = 0.6) +
  labs(y = "", x = expression(Log(alpha)), fill="Productivity Scenario", color="Productivity Scenario") 
  
my.ggsave(here("analysis/plots/fwd-sim/OM-productivity-scenarios.PNG"))
ggsave(here("csasdown/figure/OM-productivity-scenarios.PNG"), height=600*2, width=777*2, 
       units = "px", dpi=240)


# Use demographic benchmarks
#bench.par.table <- read.csv(here("analysis/data/generated/bench_par_table.csv"))
demo.bench <- read.csv(here("analysis/data/generated/demographic_parameters.csv"))

# Generate plots for which set of fwd simulations?
fit_type <- c("TVA", "TVA2")#, "AR1") # Can omit one to avoid re-generating figures

for(k in fit_type) { # generate Fwd-sim figures for reference set (TVA) & robustness set (AR1)
  
if(k == "TVA"){ 
    S.fwd <- read.csv(here("analysis/data/generated/simulations/S_fwd_TVA.csv"))
    H.fwd <- read.csv(here("analysis/data/generated/simulations/H_fwd_TVA.csv"))
    perf.metrics <- read.csv(here("analysis/data/generated/perf_metrics_TVA.csv")) |>
      pivot_longer(-c(1,11), names_to = "metric") |>
      pivot_wider(names_from=prob, values_from=value)
    Sig.R <- read.csv(here("analysis/data/generated/simulations/var_covar_TVA.csv"), 
                      row.names = 1)
    spwn.obs = TV.spwn; harv.obs = TV.harv
  } else if (k == "TVA2"){
    S.fwd <- read.csv(here("analysis/data/generated/simulations/S_fwd_TVA2.csv"))
    H.fwd <- read.csv(here("analysis/data/generated/simulations/H_fwd_TVA2.csv"))
    perf.metrics <- read.csv(here("analysis/data/generated/perf_metrics_TVA2.csv")) |>
      pivot_longer(-c(1,11), names_to = "metric") |>
      pivot_wider(names_from=prob, values_from=value)
    Sig.R <- read.csv(here("analysis/data/generated/simulations/var_covar_TVA2.csv"),
                    row.names = 1)
    spwn.obs = TV.spwn; harv.obs = TV.harv
   } else if (k == "AR1"){
    S.fwd <- read.csv(here("analysis/data/generated/simulations/S_fwd_AR1.csv"))
    H.fwd <- read.csv(here("analysis/data/generated/simulations/H_fwd_AR1.csv"))
    perf.metrics <- read.csv(here("analysis/data/generated/perf_metrics_AR1.csv")) |>
      pivot_longer(-c(1,11), names_to = "metric") |>
      pivot_wider(names_from=prob, values_from=value)
    Sig.R <- read.csv(here("analysis/data/generated/simulations/var_covar_AR1.csv"),
                      row.names = 1)
    spwn.obs = AR1.spwn; harv.obs = AR1.harv
    
}

# Variables as factors for plotting
S.fwd$CU_f <- factor(S.fwd$CU, levels = CU_order)
H.fwd$CU_f <- factor(H.fwd$CU, levels = CU_order)
demo.bench$CU_f <- factor(demo.bench$CU, levels = CU_order)

# assign HCRs to groups for subsetting harvest scenarios / HCRs
ER_seq <- seq(5, 100, 5) # Must match ER_seq in "fwd_sim.R"
HCRs <- c("no.fishing", "moratorium", "IMEG", "moratorium.cap", "IMEG.cap", "PA.alternative", paste0("fixed.ER.", ER_seq)) 
HCR_grps <- list(base = c("no.fishing", "fixed.ER.60", "IMEG"),
                 moratorium = c("no.fishing", "moratorium", "moratorium.cap", "PA.alternative"),
                 IMEG = c("no.fishing", 
                                unique(S.fwd$HCR[grepl("*IMEG*", S.fwd$HCR)])),
                 simple = c("no.fishing", "IMEG", "moratorium"),
                 fixed = unique(S.fwd$HCR[grepl("fixed.ER", S.fwd$HCR)]))
# colours 
HCR_cols <- c("#B07300", "purple3", "grey25", "#CCA000", "#FEE106",  "#0F8A2E", "#3638A5")
names(HCR_cols) <- c("PA Alternative", "Fixed ER 60", "No fishing", "Moratorium", "Moratorium cap", "IMEG", "IMEG cap")
HCR_lookup <- data.frame(HCR=c(HCRs[1:6], "fixed.ER.60"), HCR_name = names(HCR_cols)[c(3:4,6,5,7,1:2)])
# CHECK HCR_lookup to make sure names match.

## Spawners projection ----
                  
for(i in 1:length(HCR_grps[1:4])) { # don't make this fig for all fixed exp rates
  S.fwd |> filter(HCR %in% HCR_grps[[i]]) |>
    left_join(HCR_lookup, by="HCR") |>
    ggplot() +
    # Observations:
    geom_ribbon(data = filter(spwn.obs, year >= max(spwn.obs$year)-7), 
                aes(ymin = S.25/1000, ymax = S.75/1000, 
                    x= year), #offset to return year 
                fill = "grey", color = "grey") +
    geom_line(data = filter(spwn.obs, year >= max(spwn.obs$year)-7), 
              aes(y=S.50/1000, x= year), color = "black") + 
    # Projections: 
    geom_ribbon(aes(ymin = S.25/1000, ymax = S.75/1000, x = year, color=HCR_name, fill = HCR_name), 
                alpha = 0.2) +
    geom_line(aes(year, S.50/1000, color = HCR_name), lwd=1) +
    # Benchmarks:
    geom_hline(data = filter(demo.bench, par=="Smsr", period=="recent"), aes(yintercept = median/1000), 
               color = "forestgreen", lty = 2) +
    geom_hline(data = filter(demo.bench, par=="Smsr", period=="recent"), aes(yintercept = 
                             (median*0.2)/1000), color = "darkred", lty = 2) +
    facet_wrap(~CU_f, scales = "free_y", labeller=CU_labeller) +
    scale_x_continuous(expand = expansion(mult = c(0, .01))) +
    labs(y = "Spawners (000s)", x="Year", col="", fill="") +
    theme_sleek() +
    theme(legend.position = "bottom",
          strip.text = element_text(size=10),
          legend.text = element_text(size=10)) +
    scale_color_manual(values=HCR_cols, aesthetics = c("fill", "color"))
  
  my.ggsave(here(paste("analysis/plots/fwd-sim/S-fwd", names(HCR_grps[i]), "grp", paste0(k, ".PNG"), sep="_")))
  if(names(HCR_grps[i])=="simple" && k=="TVA"){
    ggsave(here(paste("csasdown/figure/S-fwd", names(HCR_grps[i]), "grp", paste0(k, ".PNG"), sep="_")),
           height=650*2, width=810*2, units="px", dpi=240)
  }
}


## Harvest projection ----
for(i in 1:length(HCR_grps[1:4])) { # don't make this fig for all fixed exp rates
  H.fwd |> filter(HCR %in% HCR_grps[[i]], HCR != "no.fishing") |>
    left_join(HCR_lookup, by="HCR") |>
    ggplot() +
    # Observations:
    geom_ribbon(data = filter(harv.obs, year >= max(harv.obs$year)-7), 
                aes(ymin = H.25/1000, ymax = H.75/1000, 
                    x= year), #offset to return year 
                fill = "grey", color = "grey") +
    geom_line(data = filter(harv.obs, year >= max(harv.obs$year)-7), 
              aes(y=H.50/1000, x= year), color = "black") + 
    # Projections:
    geom_ribbon(aes(ymin = H.25/1000, ymax = H.75/1000, x = year, color=HCR_name, fill = HCR_name), 
                alpha = 0.2) +
    geom_line(aes(year, H.50/1000, color = HCR_name), lwd=1) +
    facet_wrap(~CU_f, scales = "free_y", labeller=CU_labeller) +
    scale_x_continuous(expand = expansion(mult = c(0, .01))) +
    labs(y = "Harvest (000s)", col="", fill="") +
    theme_sleek() +
    theme(legend.position = "bottom") +
    scale_color_manual(values=HCR_cols, aesthetics = c("fill", "color"))
  
  my.ggsave(here(paste("analysis/plots/fwd-sim/H-fwd", names(HCR_grps[i]), "grp", paste0(k, ".PNG"), sep="_")))
}


## -- performance metrics and status multipanel (all HCR excl. fixed ER, all PMs) ---- 
#metrics
perf.metrics <- perf.metrics |> left_join(HCR_lookup, by="HCR") |>
  mutate(HCR_name = factor(HCR_name, levels=HCR_lookup$HCR_name))

pm_plot <- perf.metrics |> 
  filter(!(HCR %in% HCR_grps[["fixed"]])) |> 
  filter(!(metric %in% c("n.below.LSR", "n.between.ref", "n.above.USR", "n.extinct"))) |>
  mutate(metric_name = case_when(metric == "ER" ~ "Exploitation rate",
                                 metric == "pr.basic.needs" ~ "Prop. years basic needs met",
                                 metric == "pr.no.harv" ~ "Prop. years with no harvest",
                                 .default = str_to_sentence(metric)))|>
  mutate(metric_name = factor(metric_name, levels=c("Escapement", "Harvest", "Prop. years with no harvest", "Prop. years basic needs met", "Exploitation rate"))) |>
  ggplot() + 
  #geom_point(aes(x=factor(HCR, levels=HCRs[1:6]), y=mean), col="grey75", alpha=0.4) + #see means?
  geom_col(aes(x=HCR_name, y = median, fill=HCR_name)) +
  geom_segment(aes(x=HCR_name,
                   xend=HCR_name,
                   y=q_25, yend=q_75), col="grey30") +
  scale_fill_manual(values=HCR_cols) +
  facet_wrap(~metric_name, scales = "free_y", nrow=4) +
  theme_sleek() +
  scale_y_continuous() +
  theme(legend.position = c(0.8,0.15), 
        axis.text.x = element_blank(), 
        legend.title = element_blank()) +
  guides(fill=guide_legend(ncol=2, byrow=T)) +
  labs(x="", y="") 

# status
perf.status <- perf.metrics |>
  filter(metric %in% c("n.above.USR", "n.between.ref", "n.below.LSR", "n.extinct")) |>
  mutate(status = factor(gsub("^n ", "", gsub("\\.", " ", metric)), 
    levels=c("above USR", "between ref", "below LSR", "extinct")))

status_plot <- perf.status |>
  filter(!(HCR %in% HCR_grps[["fixed"]])) |>
  ggplot(aes(x = HCR_name, y = mean, fill = status)) + 
  geom_col() +
  scale_fill_discrete(type = c("forestgreen", "darkorange", "darkred", "black")) +
  scale_y_continuous(breaks = c(2,4,6,8)) +
  labs(x="Harvest control rule", y = "Average # of CUs (over 1000 simulations)", fill="Status") +
  theme_sleek()

cowplot::plot_grid(pm_plot, status_plot, nrow=2, labels="auto", rel_heights = c(1.5,1))

my.ggsave(here(paste0("analysis/plots/fwd-sim/perf_metrics_status_", k, ".PNG")))
ggsave(here(paste0("csasdown/figure/perf_metrics_status_", k, ".PNG")), height=900*2, 
       width=800*2, units="px", dpi=240,bg = "white")



## fixed ER trade-off multipanel ----

spwn_v_ER <- S.fwd %>% filter(HCR %in% HCR_grps[["fixed"]]) %>%
  group_by(HCR, CU_f) %>% 
  summarize(mean_spwn = mean(S.50)) %>%
  mutate(ER = as.numeric(gsub("\\D", "", HCR))) %>% 
  ggplot() +
  geom_point(aes(y=ER, x=mean_spwn/1000, col=CU_f), shape='circle', size=2, alpha=0.7) +
  scale_colour_viridis_d() + scale_y_continuous(breaks = seq(0,100,20)) +
  theme_sleek() +
  theme(legend.position="none",
        axis.title = element_text(size=12)) +
  labs(x="Spawners (000s)", y="Exploitation rate", col="Conservation Unit")

harv_v_ER <- H.fwd |> filter(HCR %in% HCR_grps[["fixed"]]) |>
  group_by(HCR, CU_f) |> 
  summarize(mean_harv = mean(H.50)) |>
  mutate(ER = as.numeric(gsub("\\D", "", HCR))) |> 
  left_join(CU_name_lookup) |>
  ggplot() +
  geom_point(aes(y=ER, x=mean_harv/1000, col=CU_pretty), shape='circle', size=2, alpha=0.7) +
  scale_colour_viridis_d() + scale_y_continuous(breaks = seq(0,100,20)) +
  theme_sleek() +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        legend.text = element_text(size=10),
        axis.title = element_text(size=12)) +
  labs(x="Harvest (000s)", y="Exploitation rate", col="Conservation Unit")

status_ER <- perf.status %>% filter(HCR %in% HCR_grps[["fixed"]]) %>%
  mutate(ER = as.numeric(gsub("\\D", "", HCR))) %>%
  filter(ER != 100) %>% 
  ggplot(aes(x=mean, y=factor(ER), fill=status)) +
  geom_col() +
  geom_vline(data=data.frame(x=seq(1:9)), aes(xintercept=x), col="white", linewidth=0.05) +
  scale_fill_discrete(type = c("forestgreen", "darkorange", "darkred", "black")) +
  scale_x_continuous(breaks = c(2,4,6,8)) +
  scale_y_discrete(breaks = seq(0,90,20)) +
  labs(y = "Exploitation Rate", fill="", x="Average # of CUs (over 1000 simulations)") +
  theme_sleek() + 
  theme(legend.margin = margin(l=20, r=20),
        legend.text = element_text(size=10),
        axis.title = element_text(size=12))

b <- cowplot::plot_grid(spwn_v_ER, harv_v_ER, nrow=1, rel_widths=c(0.55,1), labels=c("b","c"), label_x = c(0,-0.015))
cowplot::plot_grid(status_ER, b, nrow=2, rel_heights=c(1,1), labels=c(NULL,"a"))

my.ggsave(here(paste0("analysis/plots/fwd-sim/fixed_ER_tradeoffs_", k, ".PNG")))
ggsave(here(paste0("csasdown/figure/fixed_ER_tradeoffs_", k, ".PNG")), height=600*2, 
       width=800*2, units="px", dpi=240)


## visualize HCRs ----
HCR_order <- c("Moratorium", "IMEG", "Moratorium cap", "IMEG cap", "PA Alternative")
out <- visualize_HCR(HCRs=HCRs[2:6]) # get simulated HRs
out <- left_join(out, HCR_lookup, by="HCR")
out$HCR_name <- factor(out$HCR_name, levels=HCR_order)

ggplot(out) + geom_line(aes(x=run_size/1000, y=HR*100, col=HCR_name), linewidth=0.75) +
  scale_colour_manual(values=HCR_cols, guide="none") +
  scale_fill_manual(values="grey70", guide="legend") +
  facet_wrap(~HCR_name) + 
  labs(x="Run Size (000s)", y="Harvest Rate (%)") +
  theme_minimal() + theme(strip.text = element_text(size=10),
                          axis.title = element_text(size=10)) +
  lims(x=c(0,200)) +
  scale_y_continuous(breaks=seq(0,100,20), limits=c(0,100)) 

my.ggsave(here("analysis/plots/fwd_sim/HCR_visualize.PNG"))
my.ggsave(here("csasdown/figure/HCR_visualize.PNG"), height=4, 
          width=6)


## covariance plots ----

Sig.R.order <- Sig.R[c(4,8,6,2,5,3,1,7,9),c(4,8,6,2,5,3,1,7,9)]

ggcorrplot(Sig.R.order, hc.order = TRUE, type = "lower",
           outline.col = "white",
           lab=TRUE) 

my.ggsave(here(paste0("analysis/plots/fwd-sim/recruit-corr-matrix_", k, ".PNG")))


} # end k loop


# SMU run and escapement ----

SMU_RR <- read.csv(here("analysis/data/raw/rr_95_table.csv")) |>
  filter(Stock == "Canada", 
         Counts != "Harvest") 

a<- ggplot(SMU_RR |>
         filter(Counts != "Exploitation")) + 
  geom_hline(yintercept = 19, col = "darkred", lty=2) +
  geom_hline(yintercept = 95, col = "dark green", lty=2) +
  geom_ribbon(aes(x = Year, ymin = Lower95./1000, ymax = Upper95./1000, col = Counts, fill = Counts), alpha=0.5) +
  geom_line(aes(x = Year, y = Median50./1000, col = Counts), size = 1) + 
  ylab("Fish (000s)") +
  xlab("Year") +
  scale_color_manual(values=c('#999999','#E69F00')) +
  scale_fill_manual(values=c('#999999', '#E69F00')) +
  theme_sleek() +
  theme(legend.position = c(0.75,0.85),
        legend.title = element_blank(),
        legend.text = element_text(size=11),
        axis.title = element_text(size=12),
        plot.margin = margin(0.5,20,0.5,0.5))

b<- ggplot(SMU_RR |>
             filter(Counts == "Exploitation")) + 
  geom_hline(yintercept = 40, col = "darkred", lty=2) +
  geom_ribbon(aes(x = Year, ymin = Lower95., ymax = Upper95.), fill="darkblue", col="darkblue", alpha=0.4) +
  geom_line(aes(x = Year, y = Median50.), size = 1, col = "darkblue") + 
  ylab("Harvest rate (%)") +
  xlab("Year") +
  theme_sleek() + theme(plot.margin = margin(0.5,20,0.5,0.5),
                        axis.title = element_text(size=12))
 
cowplot::plot_grid(a, b, labels="auto", ncol=2)


my.ggsave(here("analysis/plots/trib-rr/SMU-run-esc.PNG"), width = 13, height = 6)
ggsave(here("csasdown/figure/SMU-run-esc.PNG"), width = 975*2, height = 350*2, 
       units="px", dpi=240)


# CU run-reconstructions plot ----

border <- read.csv(here("analysis/data/raw/border-passage.csv")) |>
  filter(year>1984)

par( mar=c(5,5,0,0) )

#t <- !is.na(rpt$I_t)
t <- rep(TRUE,rpt$nT)
yr  <- rpt$years[t]
I_t <- rpt$I_t[t]*1e-3/exp(rpt$lnqI_s[1])
E_t <- colSums(exp(rpt$lnRunSize_st))[t]*1e-3
sonarN_t <- colSums(rpt$E_dtg[ ,t,1])*1e-3
ymax <- max(I_t,E_t,sonarN_t,na.rm=TRUE)
fw_t <- 1e-3*colSums(rpt$E_dtg[,,2])/exp(rpt$lnqE_sg[1,2])

png(file=here("csasdown/figure/CU-RR-fits.PNG"), width= 8, height = 6,units="in", res =700 )

plot( x=yr, y=I_t, type="n", las=1, yaxs="i", xlab="Year",
      ylab="Total border passage (1000s)", ylim=c(0,1.3*ymax) )
grid()
box()

if( is.finite(rpt$sdrpt[1,5]) )
{
  Ese <- filter(rpt$sdrpt,par=="runSize_t")[t, ]
  segments( x0=yr+0.2, y0=Ese$lCI*1e-3, y1=Ese$uCI*1e-3, col="grey70", lwd=2 )
  segments( x0=yr, y0=border$lwr*1e-3, y1=border$upr*1e-3, col="black", lwd=1.5 )
  
}

points( x=yr+0.2, y=E_t, pch=16, col="grey40" )
points( x=yr, y=I_t, pch=0, lwd=1.5 )
points( x=yr, y=sonarN_t, pch=16, lwd=1.5, col="red" )
points( x=yr, y=fw_t, pch=2, lwd=1.5, col="green" )

legend( x="bottomleft", bty="n",
        legend=c("CU run reconstruction estimates","Aggregrate run reconstruction estimates","Sonar counts", "Fish wheel counts"),
        pch=c(NA,NA,NA,NA), lwd=c(1,1,NA,NA), col=c("grey70","black",NA,NA), lty=c(1,1,NA,NA), cex=0.75 )
legend( x="bottomleft", bty="n",
        legend=c("CU run reconstruction estimates","Aggregrate run reconstruction estimates","Sonar counts", "Fish wheel counts"),
        pch=c(16,0,16,2), lwd=c(1.5), col=c("grey40","black","red","green"), lty=c(0), cex=0.75 )

dev.off()

# Compare CU and aggregrate border passage ----

Ese <- filter(rpt$sdrpt,par=="runSize_t")
Ese$year <- seq(1985,2024)

border <- read.csv(here("analysis/data/raw/border-passage.csv")) |>
  filter(year>1984)

cu_rr <- Ese[,c(6,2,4:5)]; colnames(cu_rr) <- c("year", "est", "lwr", "upr"); cu_rr$model <- "CU_rr"
agg_rr <- border[,c(1:2,4:5)]; colnames(agg_rr) <- c("year", "est", "lwr", "upr"); agg_rr$model <- "Agg_rr"


bp_models <- rbind(cu_rr,agg_rr)

ggplot(bp_models,aes(x = year, y = est/1000, fill=model)) + 
  geom_ribbon(aes(ymin = lwr/1000, ymax = upr/1000),   alpha = 0.5) +
  geom_line(lwd = 1.1, ) +
  xlab("Year") +
  ylab("Border passage (000s)")+
  theme_sleek() +
  theme(strip.text = element_text(size=10))

my.ggsave(here("analysis/R/run-reconstructions/fittedMod/cu-vs-agg-rr-border.PNG"))


# CU run-timing plot ----
CU_order <- c(1,2,5,3,6,9,4,8)
cols <- viridis(rpt$nS+1)
cols_ordered <- cols[CU_order]
x <- rpt$day_d 
x_trunc <- x[x>=170 & x <=270]

y_dpt <- rpt$rho_dst
y_dp  <- apply( y_dpt, 1:2, mean )

png( file=here("csasdown/figure/CU-run-timing.PNG"), width= 9, height = 5.562,units="in", res =700 )

par(mar=c(5,15,1,1),oma=c(0,0,0,0),
    col.lab = "grey40")

plot( x=range(x_trunc), c(1,rpt$nS+1.5), type="n", axes=FALSE,
      xlab="Day of year", ylab="" )
axis( side=1 , col="grey40", col.axis="grey40") 
axis( side=2, at=rpt$nS:1, labels=CU_name_lookup$CU_pretty[match(rpt$stocks, CU_name_lookup$CU_f)], 
      las=1, cex.axis=1, col="grey40", col.axis="grey40" ) 
for( p in 1:rpt$nS )
{
  y <- 1.5*(y_dp[ ,p]/max(y_dp[ ,p]))
  polygon( x=c(x,rev(x)), y=rpt$nS-p+1+c(y,rep(0,length(x))),
           border=NA, col=cols_ordered[p] )
}

abline(v=205, col="grey",lty=2)
dev.off()



# escapement stats ----

esc_summary <- esc |>
  group_by(stock) |>
  summarize(avg_spwn = mean(mean[1:33]),
            recent_spwn = mean(mean[34:40]),
            change_spwn = (1-(recent_spwn/avg_spwn))*-100)

mean(esc_summary$change_spwn)


# enhancement Figures ------------------------------------------------------

## hatchery releases ---- 
rel_rep <- read.csv("./analysis/data/raw/yukon_releases.csv")

yukn_rel <- rel_rep %>%
  filter(Status != "Exclude") %>%
  filter(STOCK_PROD_AREA_CODE == "YUKN") %>% #Exclude - not part of the Yukon watershed
  filter(!is.na(REL_CU_INDEX)) %>%
  mutate(REL_CU = paste(REL_CU_INDEX, REL_CU_NAME)) %>%
  mutate(FACILITY_NAME = case_when(FACILITY_NAME == "Yukon River H" ~ "Whitehorse Rapids Fish H",
                                   FACILITY_NAME == "McIntyre Creek H" ~ "McIntyre Creek Fish Incubation Facility",
                                   FACILITY_NAME == "Klondike River, North H" ~ "North Klondike River H",
                                   grepl("Schools", FACILITY_NAME) ~ "School Programs",
                                   .default = FACILITY_NAME)) %>%
  mutate(FACILITY_NAME = gsub(" H", " Hatchery", FACILITY_NAME)) %>%
  group_by(RELEASE_YEAR,FACILITY_NAME, REL_CU) %>%
  dplyr::summarise(TotalRelease = sum(TotalRelease, na.rm = TRUE),
                   .groups = 'drop') %>%
  arrange(REL_CU)

breakV <- seq(min(yukn_rel$RELEASE_YEAR), max(yukn_rel$RELEASE_YEAR), by = 5) # breaks for fig

ggplot(yukn_rel, aes(x=RELEASE_YEAR, y=TotalRelease/1000, fill = REL_CU))+
  geom_bar(stat = "identity", width = 1, colour="white", linewidth = 0.1) +
  facet_wrap(~FACILITY_NAME, ncol=2,
             scales = "free_y") +
  scale_x_continuous(name = "Release Year", breaks = breakV) +
  scale_y_continuous(name = "Total Releases (000s)") +
  scale_fill_viridis_d() + labs(fill="Conservation Unit of release") +
  theme_sleek() +
  theme(legend.position = c(0.75, 0.15),
        strip.text = element_text(size=11),
        axis.title = element_text(size=13))

my.ggsave(here("analysis/plots/trib-rr/hatch_bar.PNG"))
my.ggsave(here("csasdown/figure/hatch_bar.PNG"))



## proportion hatchery fish at whitehorse fishway ---- 

trib.spwn <- read.csv(here('analysis/data/raw/trib-spwn.csv'))
wh.hatch <- trib.spwn[!is.na(trib.spwn$hatch_contrib),]

wh.hatch <- wh.hatch %>% mutate(Hatchery = estimate*hatch_contrib,
                    Wild = estimate*(1-hatch_contrib)) %>% 
  pivot_longer(cols=c("Hatchery", "Wild"), names_to="Origin", values_to="Returns")

ggplot(wh.hatch, aes(x=year)) + 
  geom_bar(aes(y=Returns, fill=Origin), stat="identity") +
  scale_y_continuous(sec.axis = sec_axis(~./4000, name="Proportion")) +
  geom_line(aes(y=hatch_contrib*4000, col="Proportion hatchery-origin \nspawners"), linewidth=0.6, alpha=0.6) +
  geom_line(aes(y=pni*4000, col="Proportionate natural \ninfluence (PNI)"), linewidth=0.6, alpha=0.6) +
  scale_color_manual(name = "", values=c("Proportion hatchery-origin \nspawners" = "darkblue", "Proportionate natural \ninfluence (PNI)" = "darkred"), guide ="legend") +
  scale_fill_manual(name = "Returns", values=c("Hatchery" = "pink2", "Wild" = "green4")) +
  theme_sleek() + labs(x="Year" ) +
  theme(axis.text.y.right = element_text(margin = margin(r=9)),
        legend.position = c(0.27,0.7),
        legend.text = element_text(size=8),
        legend.title = element_text(size=9),
        axis.title = element_text(size=10),
        legend.spacing.y = unit(0, "pt"))

my.ggsave(here("analysis/plots/trib-rr/hatch_prop.PNG"))
ggsave(here("csasdown/figure/hatch_prop.PNG"), width=700*2, height=400*2, dpi=240,
       units="px")


# time-varying productivity SR-EMR ----
sp.TVA <- read.csv(here('analysis/data/generated/spw_TVA.csv'))
em.TVA <- read.csv(here('analysis/data/generated/demographic_TVA.csv'))

TVA <- rbind(sp.TVA,em.TVA) |>
  group_by(model, CU) |>
  mutate(scale_prod = scale(mid)[,1],
         model.type = case_when(
                        model == "spw" ~ "Spawners",
                        model == "em" ~ "Egg mass"
         )) 

TVA$CU_f <- factor(TVA$CU, levels = CU_order)


ggplot(TVA
       |> filter(brood_year < 2018), aes(col = model.type)) +
  geom_line(aes(x = brood_year , y = scale_prod), lwd = 1) +
  facet_wrap(~CU_f, scales = "free_y",nrow = 3, 
             labeller=CU_labeller) +
  theme_sleek() +
  labs(y = "Scaled productivity index", x = "Brood year") +
  scale_color_viridis_d(end=0.9) +
  labs(color = "SR model ") +
  theme(legend.position = c(0.92,0.925),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.key.size = unit(0.35, "cm"),
        legend.title = element_text(size=12),
        legend.text = element_text(size=11, angle=0, hjust=0),
        strip.text = element_text(size=12),
        axis.title = element_text(size=12))

my.ggsave(here("analysis/plots/SR-models/spw-vs-em-SR-TVA.PNG"))
my.ggsave(here("csasdown/figure/spw-vs-em-SR-TVA.PNG"), height=6.60, width=7.55)

# CU contribution timeseries ----

esc |>
  left_join(CU_name_lookup, by="CU_f") 

total_spawn <- esc |>
  left_join(CU_name_lookup, by="CU_f") |>
  group_by(year) |>
  summarise(total_spwn = sum(mean)) |>
  ungroup()

prop_CU <- esc |>
  left_join(total_spawn, by="year") |>
  left_join(CU_name_lookup, by="CU_f")  |>
  mutate(contr = (mean/total_spwn)*100) |>
  select(CU_f,CU_pretty,year,contr)

ggplot(prop_CU, aes(x=year, y=contr, fill=CU_pretty)) + 
  geom_area() +
  scale_fill_viridis_d() + 
  labs(fill="Conservation Unit", x = "Year", y = "Contribution to total spawners (%)") +
  theme_sleek() 

my.ggsave(here("csasdown/figure/percent-cc-contribution.PNG"), height=4.25, width=8)

# CU contribution from raw GSI

cugsitable <- read.csv(here::here("csasdown/data/CU-gsi-annual-summary.csv"))
names(cugsitable)[2:ncol(cugsitable)] <- c("Northern Yukon R. and tribs.", 
                                           "White and tribs.", "Stewart",  
                                           "Middle Yukon R. and tribs.","Pelly",
                                           "Nordenskiold",
                                           "Upper Yukon R.","Yukon R. Teslin Headwaters")

cugsitable.long <- cugsitable |>
  pivot_longer(!Year, names_to = "Conservation Unit", values_to = "contr" )

cugsitable.long$CU_f <- factor(cugsitable.long$'Conservation Unit', levels = CU_prettynames)

a <- ggplot(prop_CU, aes(x=year, y=contr, fill=CU_f)) + 
  geom_area() +
  scale_fill_viridis_d() + 
  labs(fill="Conservation Unit", x = "Year", y = "Contribution to total spawners (%)") +
  theme_sleek() +
  theme(legend.position="none")

b<- ggplot(cugsitable.long, aes(x=Year, y=contr, fill=CU_f )) + 
  geom_area() +
  scale_fill_viridis_d() + 
  labs(fill="Conservation Unit", x = "Year", y = "Contribution to total spawners (%)") +
  theme_sleek() +
  theme(legend.position="none")

cowplot::plot_grid(a, b, nrow=2, labels="auto", rel_heights = c(1,1))

my.ggsave(here("analysis/plots/trib-rr/compare-percent-cc-contribution.PNG"), height=11, width=9)

# GSI summary stats ----

gsi <- read.csv(here("analysis/data/raw/border-gsi-table-2024-update-full.csv")) %>%
  rename( sample_num=fish )

gsi_summary <- gsi |>
  group_by(year, CU) |>
  summarize(sum_prob = sum(prob),
            samples = n_distinct(sample_num)) |>
  mutate(cu_percent = round(sum_prob/samples, 5)) |>
  select(year, CU, cu_percent) |>
  pivot_wider(names_from=CU, values_from=cu_percent)

write.csv(gsi_summary, here("analysis/data/generated/CU-gsi-annual-summary.csv"), row.names = FALSE)

# FSAR 4 panel plot ----


SMU_RR <- read.csv(here("analysis/data/raw/rr_95_table.csv")) |>
  filter(Stock == "Canada")


a <- ggplot(SMU_RR |>
             filter(Counts == "Harvest")) + 
  geom_ribbon(aes(x = Year, ymin = Lower95./1000, ymax = Upper95./1000), col = c('#999999'),alpha=0.25) +
  geom_line(aes(x = Year, y = Median50./1000, col = Counts), size = 1, col="black") + 
  ylab("Catch (000s)") +
  xlab("Year") +
  theme_sleek() +
  theme(legend.position = c(0.75,0.85),
        legend.title = element_blank(),
        legend.text = element_text(size=11),
        axis.title = element_text(size=12))

b <- ggplot(SMU_RR |>
             filter(Counts == "Escapement")) + 
  geom_ribbon(aes(x = Year, ymin = Lower95./1000, ymax = Upper95./1000), col = c('#999999'),alpha=0.25) +
  geom_line(aes(x = Year, y = Median50./1000, col = Counts), size = 1, col="black") + 
  geom_hline(yintercept = 19, col = "darkred", lty=2, lwd = 1) +
  geom_hline(yintercept = 95, col = "darkgreen", lty=2, lwd = 1) +
  ylab("Spawners (000s)") +
  xlab("Year") +
  theme_sleek() +
  theme(legend.position = c(0.75,0.85),
        legend.title = element_blank(),
        legend.text = element_text(size=11),
        axis.title = element_text(size=12))

c <- ggplot(SMU_RR |>
             filter(Counts == "Exploitation")) + 
  geom_hline(yintercept = 35, col = "darkred", lty=2, lwd=1) +
  geom_ribbon(aes(x = Year, ymin = Lower95., ymax = Upper95.), col = c('#999999'),alpha=0.25) +
  geom_line(aes(x = Year, y = Median50.), size = 1, col = "black") + 
  ylab("Exploitation rate (%)") +
  xlab("Year") +
  theme_sleek() + theme(axis.title = element_text(size=12))

d <- ggplot(SMU_RR |>
              filter(Counts == "Run")) + 
  geom_ribbon(aes(x = Year, ymin = Lower95./1000, ymax = Upper95./1000), col = c('#999999'),alpha=0.25) +
  geom_line(aes(x = Year, y = Median50./1000, col = Counts), size = 1, col="black") + 
  ylab("Returns (000s)") +
  xlab("Year") +
  theme_sleek() +
  theme(legend.position = c(0.75,0.85),
        legend.title = element_blank(),
        legend.text = element_text(size=11),
        axis.title = element_text(size=12))

cowplot::plot_grid(a, b, c, d, labels="auto", ncol=2)


ggsave(here("csasdown/figure/fsar-4-panel.PNG"), width = 975*2, height = 675*2, 
       units="px", dpi=240)
