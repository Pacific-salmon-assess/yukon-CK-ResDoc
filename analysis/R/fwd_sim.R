# fwd sim "##" comments are Q's for BC to address
library(here)
library(tidyverse)
library(mvtnorm) #for rmvnorm()
source(here("analysis/R/data_functions.R"))

# wrangle TVA fits into same structure as the "samps" matrix -----------------------------
#here(https://github.com/DylanMG/Kusko-harvest-diversity-tradeoffs/blob/master/load.R)
TVA.fits <- lapply(list.files(here("analysis/data/generated/model_fits/TVA"), 
                              full.names = T), 
                   readRDS)
names(TVA.fits) <- unique(sp_har$CU)[order(unique(sp_har$CU))]

TVA.fits <- lapply(TVA.fits, rstan::extract)

bench.posts <- read_rds(here("analysis/data/generated/benchmark_posteriors.rds"))

#infilling the "samps" object from the kusko example -------------------------------------
samps <- NULL
pi.samps <- array(NA, dim = c(nrow(TVA.fits$MiddleYukonR.andtribs.$beta), A, length(TVA.fits)))
p.samps <- array(NA, dim = c(nrow(TVA.fits$MiddleYukonR.andtribs.$beta), A, length(TVA.fits), 3))
sig.R.samps <- NULL
for(i in 1:length(names(TVA.fits))){
  sub_samps <- cbind(exp(apply(TVA.fits[[i]]$ln_alpha[, (nyrs-a_max+1):nyrs], 1, median)),#one previous gen median alpha
                     TVA.fits[[i]]$beta,
                     filter(bench.posts, CU == unique(bench.posts$CU)[i])$Umsy,
                     filter(bench.posts, CU == unique(bench.posts$CU)[i])$`80.Smsy`,
                     #latent S&R states. Kusko has same ages, but nyrs=42 
                     TVA.fits[[i]]$S[,(nyrs-A+1):nyrs], #last 4 spawner states - 39:42 in kusko 
                     TVA.fits[[i]]$R[,(nRyrs-A+2):nRyrs], #last 3 rec states - 43:45 in kusko
                     TVA.fits[[i]]$lnresid[,nRyrs]) #last resid
  colnames(sub_samps) <- c(paste0("alpha_", i), 
                           paste0("beta_", i),
                           paste0("Umsy_", i),paste0("Smsy_", i), 
                           paste0("S_", (nyrs-A+1):nyrs, "_", i), 
                           paste0("R_", (nRyrs-A+2):nRyrs, "_", i), 
                           paste0("last_resid_", i))
  samps <- cbind(samps, sub_samps) #cbind posteriors of parms we care about
  
  pi.samps[,,i] <- TVA.fits[[i]]$pi #store pis to summarise later
  
  for(j in 1:3){
    p.samps[,,i,j] <- TVA.fits[[i]]$p[,nyrs+j, ] #store ps for last 3 nRyrs to summarise later
  }
  sig.R.samps <- cbind(sig.R.samps, TVA.fits[[i]]$sigma_R)
}

#get median of pis and ps for all pops 
#(i.e. take the median across posterior slices by "page" of array)
median.pi.samps <- apply(pi.samps, c(1,2), median)
colnames(median.pi.samps) <- paste0("pi_", 1:4)

median.p.samps <- NULL
for(i in 1:3){
  median.p.samps <- cbind(median.p.samps, apply(p.samps[,,,i], c(1,2), median)) 
}
median.p.samps <- median.p.samps[,  c(1,5,9, 2,6,10, 3,7,11, 4,8,12)] #rearrange order of p match kusko
colnames(median.p.samps) <- paste("p", (nyrs+1):nRyrs, rep(1:4, each=3), sep = "_")

Sig.R <- cov(sig.R.samps)

#bind all samps into one object 
samps <- cbind(samps, median.p.samps, median.pi.samps)

#Set common conditions for simulations----------------------------------------------------
num.sims = 50 # number of Monte Carlo trials #originally 500
ny = 50 # number of years in forward simulation #originally 50
pm.yr <- ny-20
for.error <- 0.27 ## base this off something observed  
OU <- 0.1 

# --- Create array to store outcomes -----------------------------------------------------
HCRs <- c("no.fishing", "status.quo")
sim.outcomes <- array(NA, dim=c(length(HCRs), 5, num.sims))  #5 perf metrics
sim.outcomes.time <- array(NA, dim=c(length(HCRs), 3, num.sims)) #3 states: catch, esc, state 
# --- Time-varying Ricker SR dynamics ----------------------------------------------------

# run simulations
ptm <- proc.time() #w & k deleted
for(i in 1:length(HCRs)){
  for(j in 1:num.sims){ #l been changed to j
    draw <- sample(nrow(samps),1)
    alpha <- process.iteration(samps[draw,])$alpha
    beta <- process.iteration(samps[draw,])$beta
    vcov.matrix <- Sig.R # had to change this, ##check w/ BC if appropriate to leave it fixed.
    mat <- process.iteration(samps[draw,])$pis
    Rec <- process.iteration(samps[draw,])$R
    Spw <- process.iteration(samps[draw,])$S
    lst.resid <- process.iteration(samps[draw,])$last_resid
    expan <- 1/(rnorm(1,0.56,0.05))
    
    out <- process(ny,vcov.matrix,phi=NULL,mat,alpha,beta,sub=NULL,com=NULL,egfloor[k]=NULL,pm.yr,
                   for.error,OU,Rec,Spw,lst.resid,SR_rel,BH.alpha.CV=NULL,
                   period=NULL,dir.SR, SR_devs=NULL, expan)
    sim.outcomes[] <- out$PMs
    #sim.outcomes.spw.time[,,i] <- out$S
  }
}

(proc.time() - ptm)/60

saveRDS(sim.outcomes, here("analysis/data/generated/base_sims.rickerTV"))  
saveRDS(sim.outcomes.spw.time, here("analysis/data/generated/base_sims_projections.rickerTV"))  

#sim.outcomes <- read_rds(here("analysis/data/generated/base_sims.rickerTV"))
#sim.outcomes.spw.time <-  read_rds(here("analysis/data/generated/base_sims_projections.rickerTV"))
#}


#sim inference?