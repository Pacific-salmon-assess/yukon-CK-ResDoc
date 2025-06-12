# fwd sim 
library(here)
library(tidyverse)
library(mvtnorm) #for rmvnorm() nested in the functions
source(here("analysis/R/data_functions.R"))


fit_type <- c("TVA") #,"TVA2")#, "AR1") # Can omit one to avoid re-simulating it
set.seed(2)


# Load SR fits ---------------------------------------------------------------

# Load TVA fits
TVA.fits <- lapply(list.files(here("analysis/data/generated/model_fits/TVA"), 
                              full.names = T), readRDS)
names(TVA.fits) <- unique(sp_har$CU)
TVA.fits <- lapply(TVA.fits, rstan::extract)

# Load AR1 fits
AR1.fits <- lapply(list.files(here("analysis/data/generated/model_fits/AR1"), 
                              full.names = T), readRDS)
names(AR1.fits) <- unique(sp_har$CU)
AR1.fits <- lapply(AR1.fits, rstan::extract)

# Load benchmarks
#bench.posts <- read_rds(here("analysis/data/generated/benchmark_posteriors.rds"))
demo.bench.posts <- read.csv(here("analysis/data/generated/demographic_parameters.csv")) |> select(1:4) |> pivot_wider(values_from= "median", names_from = "par") |> filter(period == "recent")


## Loop over TVA and AR1 fit types to create a 'samps' object for each
for(k in fit_type){

# Populate "samps" object ----------------------------------------------

  if(k=="AR1") {
    fits = AR1.fits} else if(k == "TVA") {
      fits = TVA.fits} else if(k == "TVA2") {
        fits = TVA.fits}

samps <- NULL
pi.samps <- array(NA, dim = c(nrow(fits[[1]]$beta), A, length(fits)))
p.samps <- array(NA, dim = c(nrow(fits[[1]]$beta), A, length(fits), 3))
sig.R.samps <- NULL
comp.brood <- nyrs-6 # remove incompletely observed brood years
a_yrs = 10 # number of years to average productivity over for reference test
rob.yr <- 35 # year(s) index for robustness test

for(i in 1:length(names(fits))){ # loop over CUs
  if(k=="TVA") {
    alpha <- exp(apply(fits[[i]]$ln_alpha[, (comp.brood-a_yrs):comp.brood], 1, median))
  } else if(k=="TVA2") {
    alpha <- exp(fits[[i]]$ln_alpha[, rob.yr])
  }else if(k=="AR1") {
    alpha <- exp(fits[[i]]$lnalpha)
  }
  sub_samps <- cbind(alpha,
                     fits[[i]]$beta,
                     #filter(bench.posts, CU == names(fits)[i])$Umsy,
                     #filter(bench.posts, CU == names(fits[i]))$Smsy.80,
                     fits[[i]]$S[,(nyrs-A+1):nyrs], #last 4 spawner states 
                     fits[[i]]$R[,(nRyrs-A+2):nRyrs], #last 3 rec states 
                     fits[[i]]$lnresid[,nRyrs]) #last resid
  colnames(sub_samps) <- c(paste0("alpha_", i), 
                           paste0("beta_", i),
                           #paste0("Umsy_", i),
                           #paste0("Smsy_", i), 
                           paste0("S_", (nyrs-A+1):nyrs, "_", i), 
                           paste0("R_", (nRyrs-A+2):nRyrs, "_", i), 
                           paste0("last_resid_", i))
  samps <- cbind(samps, sub_samps) #cbind posteriors of parms we care about
  
  pi.samps[,,i] <- fits[[i]]$pi #store pis to summarise later
  
  for(j in 1:3){
    p.samps[,,i,j] <- fits[[i]]$p[,nyrs+j, ] #store ps for last 3 nRyrs to summarise later
  }
  fits[[i]]$lnresid
  sig.R.samps <- cbind(sig.R.samps, apply(fits[[i]]$lnresid, 2, median)[8:nRyrs])
} # end CU loop

#get median of pis and ps for all pops 
#(i.e. take the median across posterior slices by "page" of array)
median.pi.samps <- apply(pi.samps, c(1,2), median)
colnames(median.pi.samps) <- paste0("pi_", 1:4) 

median.p.samps <- NULL
for(i in 1:3){
  median.p.samps <- cbind(median.p.samps, apply(p.samps[,,,i], c(1,2), median)) 
}
median.p.samps <- median.p.samps[,  c(1,5,9, 2,6,10, 3,7,11, 4,8,12)] #rearrange order of p
colnames(median.p.samps) <- paste("p", (nyrs+1):nRyrs, rep(1:4, each=3), sep = "_")


Sig.R <- cov(sig.R.samps) ##look at this object and double check all good; this should calculate var-cov matrix based on median estimates of recruitment residuals by CU and year
colnames(Sig.R) = rownames(Sig.R) = names(fits) 
write.csv(Sig.R, here('analysis', 'data', 'generated', 'simulations', paste0('var_covar_', k, '.csv')))

#bind all samps into one object 
samps <- cbind(samps, median.p.samps, median.pi.samps)


#Set common conditions for simulations----------------------------------

num.sims = 1000  # number of Monte Carlo trials
ny = 34 # number of years in forward simulation (complete years through 2050; 26+8)
pm.yr <- ny-20 # nyrs that we evaluate pms across
for.error <- 0.79 # empirical estimated based on forecast vs true run 2000-present  
OU <- 0.1  ## could also base this off something else from fisheries management 

# --- Create array to store outcomes --------------------------------------
ER_seq <- seq(5, 100, 5) # how many fixed ERs to test?
HCRs <- c("no.fishing", "IMEG", "IMEG.cap", "moratorium", "moratorium.cap", "PA.alternative", "realistic", paste0("fixed.ER.", ER_seq))
sim.outcomes <- NULL
S.time <- NULL #null objects to bind to - because need dataframes for ggplot
H.time <- NULL

# -- Simulation loop -------------------------------------------------------

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
    Smax <- demo.bench.posts$Smsr
    if(grepl('fixed.ER', HCR)) {
      ER <- as.numeric(gsub("\\D", "", HCR))/100
      } else {ER <- NULL} # Set ER for fixed ER sims
      
    out <- process(HCR,ny,vcov.matrix,mat,alpha,beta,Smax, pm.yr,for.error,OU,Rec,Spw,lst.resid,phi, ER)
  
    
    if(anyNA(out$S) | any(is.nan(out$S[,]))) { print("Stop: NA/NAN spawners produced")}
    if(anyNA(out$H) | any(is.nan(out$H[,]))) { print("Stop: NA/NAN harvest produced")}

     sim.outcomes <- rbind(sim.outcomes, cbind(rep(HCR, nrow(out$PMs)), rep(j, nrow(out$PMs)), out$PMs))
    S.time <- rbind(S.time, cbind(out$S[7:ny,], rep(HCR, ny-a_max+1), 
                                  (max(sp_har$year)):(max(sp_har$year)+ny-a_max), #store trajectory while clipping out observed states
                                  rep(j, ny-a_max+1))) 
    H.time <- rbind(H.time, cbind(out$H[7:ny,], rep(HCR, ny-a_max+1), 
                                  (max(sp_har$year)):(max(sp_har$year)+ny-a_max),
                                  rep(j, ny-a_max+1)))
    
  }
} # end of simulation loop

pms <- c("escapement", "harvest", "ER", "pr.no.harv", "pr.basic.needs",  "n.below.LSR", "n.between.ref", "n.above.USR", "n.above.reb", "n.extinct")
colnames(sim.outcomes) <- c("HCR", "sim", pms)
qmean <- function(x){
  c(quantile(x, c(0.25, 0.5, 0.75)), mean(x))
}
sim.outcome.summary <- as.data.frame(sim.outcomes) |>
  mutate_at(2:ncol(sim.outcomes), as.numeric) |>
  group_by(HCR) |>
  reframe(escapement = qmean(escapement),
          harvest = qmean(harvest),
          ER = qmean(ER), 
          pr.no.harv = qmean(pr.no.harv),
          pr.basic.needs = qmean(pr.basic.needs),
          n.below.LSR = qmean(n.below.LSR), 
          n.between.ref = qmean(n.between.ref), 
          n.above.USR = qmean(n.above.USR), 
          n.above.reb = qmean(n.above.reb),
          n.extinct = qmean(n.extinct),
          prob=c("q_25","median","q_75","mean")
          )

write.csv(sim.outcome.summary, here('analysis', 'data', 'generated', paste0('perf_metrics_', k, '.csv')), 
          row.names = FALSE)

colnames(S.time) <- c(names(fits), "HCR", "year", "sim")
S.time <- as.data.frame(S.time) |>
  pivot_longer(1:9, names_to = "CU") |>
  rename(Spawners = value) |>
  mutate(Spawners = round(as.numeric(Spawners), 0), 
         year = as.numeric(year))

colnames(H.time) <- c(names(fits), "HCR", "year", "sim")
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
write.csv(S.fwd.summmary, here('analysis', 'data', 'generated', 'simulations', 
                               paste0('S_fwd_', k, '.csv')),  row.names = FALSE)

H.fwd.summmary <- H.time |>
  group_by(HCR, CU, year) |>
  summarise(H.50 = median(Harvest), 
            H.25 = quantile(Harvest, 0.25), 
            H.75 = quantile(Harvest, 0.75))
write.csv(H.fwd.summmary, here('analysis', 'data', 'generated', 'simulations', 
                               paste0('H_fwd_', k, '.csv')),  row.names = FALSE)

} # end k loop
