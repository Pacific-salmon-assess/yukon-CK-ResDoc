# fwd sim 
library(here)
library(tidyverse)
library(mvtnorm) #for rmvnorm() nested in the functions
source(here("analysis/R/data_functions.R"))

set.seed(2)

# wrangle fits into same structure as the "samps" matrix -----------------------------
#here(https://github.com/DylanMG/Kusko-harvest-diversity-tradeoffs/blob/master/load.R)

#fit_type <- "AR1" # delete?

AR1.fits <- lapply(list.files(here("analysis/data/generated/model_fits/AR1"), 
                              full.names = T), readRDS)
names(AR1.fits) <- unique(sp_har$CU)

AR1.fits <- lapply(AR1.fits, rstan::extract)

bench.posts <- read_rds(here("analysis/data/generated/benchmark_posteriors.rds"))

#populating the "samps" object  -------------------------------------
samps <- NULL
pi.samps <- array(NA, dim = c(nrow(AR1.fits[[1]]$beta), A, length(AR1.fits)))
p.samps <- array(NA, dim = c(nrow(AR1.fits[[1]]$beta), A, length(AR1.fits), 3))
sig.R.samps <- NULL
for(i in 1:length(names(AR1.fits))){
  sub_samps <- cbind(exp(AR1.fits[[i]]$lnalpha),
                     AR1.fits[[i]]$beta,
                     filter(bench.posts, CU == unique(bench.posts$CU)[i])$Umsy,
                     filter(bench.posts, CU == unique(bench.posts$CU)[i])$Smsy.80,
                     AR1.fits[[i]]$S[,(nyrs-A+1):nyrs], #last 4 spawner states 
                     AR1.fits[[i]]$R[,(nRyrs-A+2):nRyrs], #last 3 rec states 
                     AR1.fits[[i]]$lnresid[,nRyrs]) #last resid
  colnames(sub_samps) <- c(paste0("alpha_", i), 
                           paste0("beta_", i),
                           paste0("Umsy_", i),
                           paste0("Smsy_", i), 
                           paste0("S_", (nyrs-A+1):nyrs, "_", i), 
                           paste0("R_", (nRyrs-A+2):nRyrs, "_", i), 
                           paste0("last_resid_", i))
  samps <- cbind(samps, sub_samps) #cbind posteriors of parms we care about
  
  pi.samps[,,i] <- AR1.fits[[i]]$pi #store pis to summarise later
  
  for(j in 1:3){
    p.samps[,,i,j] <- AR1.fits[[i]]$p[,nyrs+j, ] #store ps for last 3 nRyrs to summarise later
  }
  AR1.fits[[i]]$lnresid
  sig.R.samps <- cbind(sig.R.samps, apply(AR1.fits[[i]]$lnresid, 2, median)[8:nRyrs])
}

#get median of pis and ps for all pops 
#(i.e. take the median across posterior slices by "page" of array)
median.pi.samps <- apply(pi.samps, c(1,2), median) # Is this supposed to avg across CUs?
colnames(median.pi.samps) <- paste0("pi_", 1:4) 

median.p.samps <- NULL
for(i in 1:3){
  median.p.samps <- cbind(median.p.samps, apply(p.samps[,,,i], c(1,2), median)) 
}
median.p.samps <- median.p.samps[,  c(1,5,9, 2,6,10, 3,7,11, 4,8,12)] #rearrange order of p match kusko
colnames(median.p.samps) <- paste("p", (nyrs+1):nRyrs, rep(1:4, each=3), sep = "_")

Sig.R <- cov(sig.R.samps) ##look at this object and double check all good; this should calculate var-cov matrix based on median estimates of recruitment residuals by CU and year

#bind all samps into one object 
samps <- cbind(samps, median.p.samps, median.pi.samps)

#Set common conditions for simulations----------------------------------------------------
num.sims = 500 # number of Monte Carlo trials
ny = 26 # number of years in forward simulation (through 2050)
pm.yr <- ny-20 # nyrs that we evaluate pms across
for.error <- 0.79 # empirical estimated based on forecast vs true run 2000-present  
OU <- 0.1  ## could also base this off something else from fisheries management 

# --- Create array to store outcomes -----------------------------------------------------
ER_seq <- seq(5, 100, 5) # how many fixed ERs to test?
HCRs <- c("no.fishing", "status.quo", "status.quo.cap", "rebuilding", "rebuilding.cap", "alt.rebuilding", paste0("fixed.ER.", ER_seq))
sim.outcomes <- NULL
S.time <- NULL #null objects to bind to - because need dataframes for ggplot
H.time <- NULL
# --- Stationary Ricker SR dynamics ----------------------------------------------------

# run simulation
for(i in 1:length(HCRs)){
  for(j in 1:num.sims){
    HCR <- HCRs[i]
    draw <- sample(nrow(samps),1)
    alpha <- process.iteration(samps[draw,])$alpha
    beta <- process.iteration(samps[draw,])$beta
    vcov.matrix <- Sig.R
    mat <- process.iteration(samps[draw,])$pis
    Rec <- process.iteration(samps[draw,])$R
    Spw <- process.iteration(samps[draw,])$S
    lst.resid <- process.iteration(samps[draw,])$last_resid
    phi <- 0.75 # mean across CUs
    if(grepl('fixed.ER', HCR)) {
      ER <- as.numeric(gsub("\\D", "", HCR))/100
      } else {ER <- NULL} # Set ER for fixed ER sims
      
    out <- process(HCR,ny,vcov.matrix,mat,alpha,beta,pm.yr,for.error,OU,Rec,Spw,lst.resid,phi, ER)
  
    

    sim.outcomes <- rbind(sim.outcomes, cbind(rep(HCR, nrow(out$PMs)), rep(j, nrow(out$PMs)), out$PMs))
    S.time <- rbind(S.time, cbind(out$S[7:ny,], rep(HCR, ny-a_max+1), 
                                  (max(sp_har$year)):(max(sp_har$year)+ny-a_max), #store trajectory while clipping out observed states
                                  rep(j, ny-a_max+1))) 
    H.time <- rbind(H.time, cbind(out$H[7:ny,], rep(HCR, ny-a_max+1), 
                                  (max(sp_har$year)):(max(sp_har$year)+ny-a_max),
                                  rep(j, ny-a_max+1)))
    
  }
} # end of simulation loop

pms <- c("escapement", "harvest", "ER", "pr.no.harv", "pr.basic.needs", "harv.stability", "n.below.LSR", "n.between.ref", "n.above.USR", "n.extinct")
colnames(sim.outcomes) <- c("HCR", "sim", pms)

sim.outcome.summary <- as.data.frame(sim.outcomes) |>
  mutate_at(2:ncol(sim.outcomes), as.numeric) |>
  group_by(HCR) |>
  summarise(escapement = mean(escapement), 
            harvest = mean(harvest), 
            ER = mean(ER), 
            pr.no.harv = mean(pr.no.harv),
            pr.basic.needs = mean(pr.basic.needs),
            harv.stability = mean(harv.stability, na.rm=TRUE), # how should NAs be treated? there should be NAs where harv=0 
            n.below.LSR = mean(n.below.LSR), 
            n.between.ref = mean(n.between.ref), 
            n.above.USR = mean(n.above.USR), 
            n.extinct = mean(n.extinct))

write.csv(sim.outcome.summary, here("analysis/data/generated/perf_metrics_AR1.csv"), 
          row.names = FALSE)

colnames(S.time) <- c(names(AR1.fits), "HCR", "year", "sim")
S.time <- as.data.frame(S.time) |>
  pivot_longer(1:9, names_to = "CU") |>
  rename(Spawners = value) |>
  mutate(Spawners = round(as.numeric(Spawners), 0), 
         year = as.numeric(year))

colnames(H.time) <- c(names(AR1.fits), "HCR", "year", "sim")
H.time <- as.data.frame(H.time) |>
  pivot_longer(1:9, names_to = "CU") |>
  rename(Harvest = value) |>
  mutate(Harvest = round(as.numeric(Harvest), 0), 
         year = as.numeric(year))

S.fwd.summmary <- S.time |>
  group_by(HCR, CU, year) |>
  summarise(S.50 = median(Spawners), 
            S.25 = quantile(Spawners, 0.25), 
            S.75 = quantile(Spawners, 0.75))
write.csv(S.fwd.summmary, here("analysis/data/generated/simulations/S_fwd_AR1.csv"), 
          row.names = FALSE)

H.fwd.summmary <- H.time |>
  group_by(HCR, CU, year) |>
  summarise(H.50 = median(Harvest), 
            H.25 = quantile(Harvest, 0.25, na.rm = T), 
            H.75 = quantile(Harvest, 0.75, na.rm = T))
write.csv(H.fwd.summmary, here("analysis/data/generated/simulations/H_fwd_AR1.csv"), 
          row.names = FALSE)

