# fit SR models for all CUs, run diagnostics, make initial plots
library(here)
library(tidyverse)
library(ggplotify) #for as.ggplot() to help mcmc_combo() plotting
library(rstan)
library(bayesplot)
library(shinystan)
library(gsl) #for lambertw0() to calc U_MSY
source(here("analysis/R/data_functions.R"))

refit <- FALSE #toggle T/F if you want to refit models
#refit <- TRUE

# fit AR1 and time varying productivity (TVA) models--------------------------------------
if(refit == TRUE){
  for(i in unique(sp_har$cu)){
    
    sp_har1 <- filter(sp_har, cu == i) 
    
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
    
    #AR1.fit <- stan(file = here("analysis/Stan/SS-SR_AR1.stan"), 
    #                data = stan.data,
    #                iter = 8000)
    
    #saveRDS(AR1.fit, here("analysis/data/generated/model_fits/AR1/", 
    #                      paste0(i, "_AR1.rds")))
    
    TV.fit <- stan(file = here("analysis/Stan/SS-SR_TVA2.stan"), 
                    data = stan.data,
                    iter = 8000)
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
  names(AR1.fits) <- unique(sp_har$cu)[order(unique(sp_har$cu))]
  
  TVA.fits <- lapply(list.files(here("analysis/data/generated/model_fits/TVA"), 
                                full.names = T), 
                     readRDS)
  names(TVA.fits) <- unique(sp_har$cu)[order(unique(sp_har$cu))]
}

# describe diagnostics in loop -----------------------------------------------------------
min_ESS_AR1 <- NULL
max_Rhat_AR1 <- NULL
min_ESS_TVA <- NULL
max_Rhat_TVA <- NULL
for(i in unique(sp_har$cu)){
  sub_dat <- filter(sp_har, cu==i)
  
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
  
#modelling results -----------------------------------------------------------------------
bench.par.table <- NULL 
bench.posts <- NULL
a.yrs.all <- NULL

for(i in unique(sp_har$cu)){
  sub_dat <- filter(sp_har, cu==i)
  sub_pars <- rstan::extract(AR1.fits[[i]])
  
  #latent states of spawners and recruits---
  spwn.quant <- apply(sub_pars$S, 2, quantile, probs=c(0.1,0.5,0.9))[,1:(nyrs-a_min)]
  rec.quant <- apply(sub_pars$R, 2, quantile, probs=c(0.1,0.5,0.9))[,(A+a_min):nRyrs]
  
  brood_t <- as.data.frame(cbind(sub_dat$year[1:(nyrs-A)],t(spwn.quant), t(rec.quant))) |>
    round(2)
  colnames(brood_t) <- c("BroodYear","S_lwr","S_med","S_upr","R_lwr","R_med","R_upr")
  
  #SR relationship based on full posterior---
  spw <- seq(0,max(brood_t$S_upr),length.out=100)
  iter <- length(sub_pars$lnalpha)
  SR_pred <- matrix(NA,length(spw), iter)
  
  bench <- matrix(NA,iter,4,
                  dimnames = list(seq(1:iter), c(paste0("Sgen_", i), paste0("Smsy_", i), 
                                                 paste0("Umsy_", i), paste0("Seq_", i))))
  
  for(j in 1:iter){ 
    ln_a <- sub_pars$lnalpha[j]
    b <- sub_pars$beta[j]
    SR_pred[,j] <- (exp(ln_a)*spw*exp(-b*spw))
    
    bench[j,2] <- get_Smsy(ln_a, b) #S_MSY
    bench[j,1] <- get_Sgen(exp(ln_a),b,-1,1/b*2, bench[j,2]) #S_gen
    bench[j,3] <- (1 - lambert_W0(exp(1 - ln_a))) #U_MSY
    bench[j,4] <- ln_a/b #S_eq
  }
  
  SR_pred <- as.data.frame(cbind(spw,t(apply(SR_pred, 1, quantile,probs=c(0.1,0.5,0.9), na.rm=T))))|>
    round(2)
  colnames(SR_pred) <- c("Spawn", "Rec_lwr","Rec_med","Rec_upr")
  
  # get benchmarks & pars ------------------------------------------------------------------
  bench[,2] <- bench[,2]*0.8 #make it 80% Smsy
  
  bench.posts <- cbind(bench.posts, bench)
  
  bench.quant <- apply(bench, 2, quantile, probs=c(0.1,0.5,0.9), na.rm=T) |>
    t()
  
  mean <- apply(bench,2,mean, na.rm=T) #get means of each
  
  sub_benchmarks <- cbind(bench.quant, mean) |>
    as.data.frame() |>
    mutate(cu = i) |>
    relocate('50%', 1)
  
  #other pars to report 
  alpha <- exp(quantile(sub_pars$lnalpha, probs = c(.1, .5, .9)))
  beta <- quantile(sub_pars$beta, probs = c(.1, .5, .9))
  sigma <- quantile(sub_pars$sigma_R, probs = c(.1, .5, .9))
  phi <- quantile(sub_pars$phi, probs = c(.1, .5, .9))
  
  par.quants <- rbind(alpha, sigma, beta, phi)
  
  #make big table of bench and pars
  par.summary <- as.data.frame(rstan::summary(AR1.fits[[i]])$summary) |>
    select(mean, n_eff, Rhat)
  
  #summarise not other pars...
  par.summary <- filter(par.summary, row.names(par.summary) %in% c('lnalpha', 'beta',
                                                                   'phi', 'sigma_R')) |>
    slice(1,4,3,2) |>
    mutate(cu = i)
  
  pars <- cbind(par.quants, par.summary)
  
  sub.bench.par.table <- bind_rows(sub_benchmarks, pars) |>
    mutate(n_eff = round(n_eff, 0),
           Rhat = round(Rhat, 4))
  
  sub.bench.par.table <- mutate(sub.bench.par.table, bench.par = rownames(sub.bench.par.table)) 
  
  bench.par.table <- bind_rows(bench.par.table, sub.bench.par.table)
  
  #Plotting ------------------------------------------------------------------------------
  # plot SR relationship ---
  ggplot() +
    geom_abline(intercept = 0, slope = 1,col="dark grey") +
    geom_ribbon(data = SR_pred, aes(x = Spawn, ymin = Rec_lwr, ymax = Rec_upr),
                fill = "grey80", alpha=0.5, linetype=2, colour="gray46") +
    geom_line(data = SR_pred, aes(x = Spawn, y = Rec_med)) +
    geom_errorbar(data = brood_t, aes(x= S_med, y = R_med, ymin = R_lwr, ymax = R_upr),
                  colour="grey", width=0, linewidth=0.3) +
    geom_errorbarh(data = brood_t, aes(y = R_med, xmin = S_lwr, xmax = S_upr),
                   height=0, colour = "grey", linewidth = 0.3) +
    geom_point(data = brood_t,
               aes(x = S_med,
                   y = R_med,
                   color=BroodYear),
               size = 3) +
    scale_colour_viridis_c(name = "Brood Year")+
    labs(x = "Spawners",
         y = "Recruits", 
         title = paste(i, "S-R fit")) +
    theme(legend.position = "bottom",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.key.size = unit(0.4, "cm"),
          legend.title = element_text(size=9),
          legend.text = element_text(size=8))
  my.ggsave(here("analysis/plots/", paste0("SRR_", i, ".PNG")))
  
  # then residuals---
  ## NEED TO CHECK HOW THESE YEARS LINE UP
  resid.quant <- apply(sub_pars$lnresid, 2, quantile, probs=c(0.1,0.25,0.5,0.75,0.9))[,(A):nRyrs] ##CHECK INDEX
  
  resids <- as.data.frame(cbind(sub_dat$year, t(resid.quant))) ##CHECK INDEX
  colnames(resids) <- c("year","lwr","midlwr","mid","midupr","upr")
  
  ggplot(resids, aes(x=year, y = mid)) +
    geom_ribbon(aes(ymin = lwr, ymax = upr),  fill = "darkgrey", alpha = 0.5) +
    geom_ribbon(aes(ymin = midlwr, ymax = midupr),  fill = "black", alpha=0.2) +
    geom_line(lwd = 1.1) +
    coord_cartesian(ylim=c(-2,2)) +
    labs(x = "Return year",
         y = "Recruitment residuals", 
         title = paste(i, "recruitment residuals")) +
    theme(legend.position = "none",
          panel.grid = element_blank()) +
    geom_abline(intercept = 0, slope = 0, col = "dark grey", lty = 2)
  my.ggsave(here("analysis/plots/", paste0("rec_resids_", i, ".PNG")))
  
  #KOBE PLOT - sort out the SR links first
  
  #time varying alpha plot 
  sub_pars_TVA <- rstan::extract(TVA.fits[[i]])
  
  a.yrs <- apply(sub_pars_TVA$ln_alpha, 2, quantile, probs=c(0.1,0.5,0.9))
  a.yrs <- as.data.frame(cbind(sub_dat$year, t(a.yrs)))
  
  colnames(a.yrs) <- c("brood_year", "lwr", "mid", "upr")
  
  ggplot(a.yrs) +
    geom_ribbon(aes(x = brood_year, ymin = lwr, ymax = upr), fill = "darkgrey", alpha = 0.5) +
    geom_line(aes(x = brood_year, y = mid), lwd = 2,  color = "black") +
    labs(y = "Productivity (Ricker alpha 80th percentiles)", x = "Brood year", 
         title = paste(i, "time-varying productivity"))
  my.ggsave(here("analysis/plots/", paste0("TV_alpha_", i, ".PNG")))
  
  a.yrs.all <- rbind(a.yrs.all, data.frame(a.yrs, cu = i)) #store all alpha trends for plotting outside loop 
  
  #time varying alpha residuals  
  resid.quant <- apply(sub_pars_TVA$lnresid, 2, quantile, 
                       probs=c(0.1,0.25,0.5,0.75,0.9))[,(A):nRyrs] ##CHECK INDEX
  
  resids <- as.data.frame(cbind(sub_dat$year, t(resid.quant))) ##CHECK INDEX
  colnames(resids) <- c("year","lwr","midlwr","mid","midupr","upr")
  
  ggplot(resids, aes(x=year, y = mid)) +
    geom_ribbon(aes(ymin = lwr, ymax = upr),  fill = "darkgrey", alpha = 0.5) +
    geom_ribbon(aes(ymin = midlwr, ymax = midupr),  fill = "black", alpha=0.2) +
    geom_line(lwd = 1.1) +
    coord_cartesian(ylim=c(-2,2)) +
    labs(x = "Return year",
         y = "Recruitment residuals", 
         title = paste(i, "recruitment residuals (time-varying productivity)")) +
    theme(legend.position = "none",
          panel.grid = element_blank()) +
    geom_abline(intercept = 0, slope = 0, col = "dark grey", lty = 2)
  my.ggsave(here("analysis/plots/", paste0("TV_rec_resids_", i, ".PNG"))) 
}

ggplot(a.yrs.all, aes(color = cu)) +
  geom_line(aes(x = brood_year, y = mid), lwd = 2) +
  scale_color_viridis_d() +
  labs(y = "Productivity (Ricker alpha 80th percentiles)", x = "Brood year", 
       title = "Time-varying productivity across CUs")
my.ggsave(here("analysis/plots/TVA_all.PNG"))

write_rds(bench.posts, here("analysis/data/generated/benchmark_posteriors.rds")) 

bench.par.table <- bench.par.table |>
  relocate(cu, 1) |>
  relocate(bench.par, .after = 1) |>
  relocate(mean, .after = 2)

write.csv(bench.par.table, here("analysis/data/generated/bench_par_table.csv")) 
