# using different periods of time to estimate time-varying (TV) benchmarks
library(here)
library(tidyverse)
library(gsl)
source(here("analysis/R/data_functions.R"))

TVA.fits <- lapply(list.files(here("analysis/data/generated/model_fits/TVA"), 
                              full.names = T), 
                   readRDS)
names(TVA.fits) <- unique(sp_har$CU)[order(unique(sp_har$CU))]

candidate.bench.par.table <- NULL #empty objects to rbind CU's outputs to 
candidate.bench.posts <- NULL
period <- (nyrs-A):nyrs #last generation - can change if desired

for(i in unique(sp_har$CU)){
  sub_dat <- filter(sp_har, CU==i)
  sub_pars <- rstan::extract(TVA.fits[[i]])
  
  bench <- matrix(NA,length(sub_pars$beta),5,
                  dimnames = list(seq(1:length(sub_pars$beta)), c("Sgen", "80.Smsy", "Umsy", "Seq", "S.recent")))
  
  for(j in 1:length(sub_pars$beta)){ 
    ln_a <- median(sub_pars$ln_alpha[j, period]) #TOGGLE TO PLAY WITH DIFFERENT PERIODS of alpha
    b <- sub_pars$beta[j]

    bench[j,2] <- get_Smsy(ln_a, b) #S_MSY
    bench[j,1] <- get_Sgen(exp(ln_a),b,-1,1/b*2, bench[j,2]) #S_gen
    bench[j,3] <- (1 - lambert_W0(exp(1 - ln_a))) #U_MSY
    bench[j,4] <- ln_a/b #S_eq
    bench[j,5] <- mean(sub_pars$S[j, (nyrs-4):nyrs]) #mean spawners in last generation 
  }
  
  # get benchmarks & pars ------------------------------------------------------------------
  bench[,2] <- bench[,2]*0.8 #make it 80% Smsy
  
  candidate.bench.posts <- rbind(candidate.bench.posts, as.data.frame(bench) |> mutate(CU = i))
  
  bench.quant <- apply(bench[,1:4], 2, quantile, probs=c(0.1,0.5,0.9), na.rm=T) |>
    t()
  
  mean <- apply(bench[,1:4],2,mean, na.rm=T) #get means of each
  
  sub_benchmarks <- cbind(bench.quant, mean) |>
    as.data.frame() |>
    mutate(CU = i) |>
    relocate('50%', 1)
  
  #other pars to report 
  alpha <- quantile(exp(sub_pars$ln_alpha[, period]), probs = c(.1, .5, .9))
  beta <- quantile(sub_pars$beta, probs = c(.1, .5, .9))
  sigma <- quantile(sub_pars$sigma_R, probs = c(.1, .5, .9))

  par.quants <- as.data.frame(rbind(alpha, beta, sigma)) |>
    mutate(CU = i)
  
  #make big table of bench and pars
  #par.summary <- as.data.frame(rstan::summary(TVA.fits[[i]])$summary) |>
  #  select(mean, n_eff, Rhat)
  
  #summarise not other pars...
  #par.summary <- filter(par.summary, row.names(par.summary) %in% c('ln_alpha', 'beta',
                                                                   #'sigma_R')) |>
  #  mutate(CU = i)
  #par.summary[1,1] <- exp(par.summary[1,1]) #exp ln_alpha
  #^ doesn't work nicely with the multiple alphas; turned off for now
  
  #pars <- cbind(par.quants, par.summary)
  
  sub.bench.par.table <- bind_rows(sub_benchmarks, par.quants)# |>
   # mutate(n_eff = round(n_eff, 0),
    #       Rhat = round(Rhat, 4))
  
  sub.bench.par.table <- mutate(sub.bench.par.table, bench.par = rownames(sub.bench.par.table)) 
  
  candidate.bench.par.table <- bind_rows(candidate.bench.par.table, sub.bench.par.table)
}

bench.long <- pivot_longer(candidate.bench.posts, cols = c(Sgen, '80.Smsy', S.recent), names_to = "par") |>
  select(-Umsy, - Seq) |>
  arrange(CU, par, value) |>
  filter(value <= 10000)

ggplot(bench.long, aes(value/1000, fill = par, color = par)) +
  geom_density(alpha = 0.3) +
  geom_vline(xintercept = 1.5) +
  facet_wrap(~CU, scales = "free_y") +
  theme(legend.position = "bottom") +
  scale_fill_manual(values = c("black", "darkred", "forestgreen"), 
                    aesthetics = c("fill", "color"), 
                    labels = c(expression(italic(S[recent])), expression(italic(S[gen])), 
                               expression(italic(paste("80% ",S)[MSY])))) +
  theme(axis.ticks.y = element_blank(), 
        axis.text.y = element_blank(), 
        legend.title=element_blank()) +
  labs(x = "Spawners (thousands)", y = "Posterior density", 
       title = "Recent spawners relative to TIME- VARYING benchmarks and 1500 cutoff")
my.ggsave(here("analysis/plots/TV_status.PNG"))

candidate.bench.par.table <- candidate.bench.par.table |>
  relocate(CU, 1) |>
  relocate(bench.par, .after = 1) |>
  relocate(mean, .after = 2) |>
  mutate_at(3:6, ~round(.,5)) |>
  arrange(bench.par, mean)