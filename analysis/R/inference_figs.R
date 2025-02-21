#old code for looping through and making individual plots 1 by 1
library(here)
library(tidyverse)
library(ggplotify) #for as.ggplot() to help mcmc_combo() plotting
library(rstan)
library(bayesplot)
library(shinystan)
library(gsl) #for lambertw0() to calc U_MSY
source(here("analysis/R/data_functions.R"))

# read in model fits ---------------------------------------------------------------------
AR1.fits <- lapply(list.files(here("analysis/data/generated/model_fits/AR1"),
                              full.names = T), 
                   readRDS)
names(AR1.fits) <- unique(sp_har$cu)[order(unique(sp_har$cu))]

TVA.fits <- lapply(list.files(here("analysis/data/generated/model_fits/TVA"), 
                              full.names = T), 
                   readRDS)
names(TVA.fits) <- unique(sp_har$cu)[order(unique(sp_har$cu))]

# processing data and fits to make plots later -------------------------------------------
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
  
  bench <- matrix(NA,iter,5,
                  dimnames = list(seq(1:iter), c("Sgen", "80.Smsy", "Umsy", "Seq", "S.recent")))
  
  for(j in 1:iter){ 
    ln_a <- sub_pars$lnalpha[j]
    b <- sub_pars$beta[j]
    SR_pred[,j] <- (exp(ln_a)*spw*exp(-b*spw))
    
    bench[j,2] <- get_Smsy(ln_a, b) #S_MSY
    bench[j,1] <- get_Sgen(exp(ln_a),b,-1,1/b*2, bench[j,2]) #S_gen
    bench[j,3] <- (1 - lambert_W0(exp(1 - ln_a))) #U_MSY
    bench[j,4] <- ln_a/b #S_eq
    bench[j,5] <- mean(sub_pars$S[j, (nyrs-4):nyrs]) #mean spawners in last generation 
  }
  
  SR_pred <- as.data.frame(cbind(spw,t(apply(SR_pred, 1, quantile,probs=c(0.1,0.5,0.9), na.rm=T))))|>
    round(2)
  colnames(SR_pred) <- c("Spawn", "Rec_lwr","Rec_med","Rec_upr")
  
  # get benchmarks & pars ------------------------------------------------------------------
  bench[,2] <- bench[,2]*0.8 #make it 80% Smsy
  
  bench.posts <- rbind(bench.posts, as.data.frame(bench) |> mutate(CU = i))
  
  bench.quant <- apply(bench[,1:4], 2, quantile, probs=c(0.1,0.5,0.9), na.rm=T) |>
    t()
  
  mean <- apply(bench[,1:4],2,mean, na.rm=T) #get means of each
  
  sub_benchmarks <- cbind(bench.quant, mean) |>
    as.data.frame() |>
    mutate(cu = i) |>
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
    mutate(cu = i)
  par.summary[1,1] <- exp(par.summary[1,1]) #exp ln_alpha
  
  pars <- cbind(par.quants, par.summary)
  
  sub.bench.par.table <- bind_rows(sub_benchmarks, pars) |>
    mutate(n_eff = round(n_eff, 0),
           Rhat = round(Rhat, 4))
  
  sub.bench.par.table <- mutate(sub.bench.par.table, bench.par = rownames(sub.bench.par.table)) 
  
  bench.par.table <- bind_rows(bench.par.table, sub.bench.par.table)
  
  #Plotting ------------------------------------------------------------------------------

  # then residuals---
  resid.quant <- apply(sub_pars$lnresid, 2, quantile, probs=c(0.1,0.25,0.5,0.75,0.9))[,(A):nRyrs] 
  
  resids <- as.data.frame(cbind(sub_dat$year, t(resid.quant))) ##CHECK INDEX
  colnames(resids) <- c("year","lwr","midlwr","mid","midupr","upr")
  
  #time varying alpha plot 
  sub_pars_TVA <- rstan::extract(TVA.fits[[i]])
  
  a.yrs <- apply(exp(sub_pars_TVA$ln_alpha), 2, quantile, probs=c(0.1,0.5,0.9))
  a.yrs <- as.data.frame(cbind(sub_dat$year, t(a.yrs)))
  
  colnames(a.yrs) <- c("brood_year", "lwr", "mid", "upr")
  
  
  a.yrs.all <- rbind(a.yrs.all, data.frame(a.yrs, cu = i)) #store all alpha trends for plotting outside loop 
  
  #time varying alpha residuals  
  resid.quant <- apply(sub_pars_TVA$lnresid, 2, quantile, 
                       probs=c(0.1,0.25,0.5,0.75,0.9))[,(A):nRyrs]
  
  resids <- as.data.frame(cbind(sub_dat$year, t(resid.quant))) 
  colnames(resids) <- c("year","lwr","midlwr","mid","midupr","upr")
  

}

ggplot(a.yrs.all
       |> filter(brood_year < 2018), aes(color = cu)) +
  geom_line(aes(x = brood_year , y = mid), lwd = 2) +
  scale_color_viridis_d() +
  labs(y = "Productivity (Ricker alpha 80th percentiles)", x = "Brood year", 
       title = "Time-varying productivity across CUs") #+
# geom_hline(yintercept = c(2.51, 3.88), lwd = 1, lty = 2) #plugging in range of AR1 means

my.ggsave(here("analysis/plots/TVA_all.PNG"))

write_rds(bench.posts, here("analysis/data/generated/benchmark_posteriors.rds")) 

bench.par.table <- bench.par.table |>
  relocate(cu, 1) |>
  relocate(bench.par, .after = 1) |>
  relocate(mean, .after = 2) |>
  mutate_at(3:6, ~round(.,5)) |>
  arrange(bench.par, mean)

write.csv(bench.par.table, here("analysis/data/generated/bench_par_table.csv"), 
          row.names = FALSE) 
