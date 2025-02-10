library(here)
library(tidyverse)
source(here("analysis/R/data_functions.R"))

#wrangle TVA fits into same structure as the "samps" matrix here(https://github.com/DylanMG/Kusko-harvest-diversity-tradeoffs/blob/master/load.R)
TVA.fits <- lapply(list.files(here("analysis/data/generated/model_fits/TVA"), 
                              full.names = T), 
                   readRDS)
names(TVA.fits) <- unique(sp_har$cu)[order(unique(sp_har$cu))]

TVA.fits <- lapply(TVA.fits, rstan::extract)

bench.posts <- read_rds(here("analysis/data/generated/benchmark_posteriors.rds"))

#infilling the "samps" object from the kusko example, not I am NOT SURE if these need to be in order!
  #- especially ordering on the things that will be matricies (i.e., Sigma_R, R, S)
samps <- NULL
pi.samps <- array(NA, dim = c(nrow(TVA.fits$MiddleYukonR.andtribs.$beta), A, length(TVA.fits)))
p.samps <- array(NA, dim = c(nrow(TVA.fits$MiddleYukonR.andtribs.$beta), A, length(TVA.fits), 3))

for(i in 1:length(names(TVA.fits))){
  sub_samps <- cbind(exp(apply(TVA.fits[[i]]$ln_alpha[, (nyrs-a_max+1):nyrs], 1, median)),#one previous gen median alpha
                     TVA.fits[[i]]$beta,
                     bench.posts[, grep("Umsy", colnames(bench.posts))[i]],
                     bench.posts[, grep("Smsy", colnames(bench.posts))[i]], 
                     #latent S&R states. Kusko has same ages, but nyrs=42 
                     TVA.fits[[i]]$S[,(nyrs-A+1):nyrs], #last 4 spawner states - 39:42 in kusko 
                     TVA.fits[[i]]$R[,(nRyrs-A+2):nRyrs], #last 3 rec states - 43:45 in kusko
                     TVA.fits[[i]]$lnresid[,nRyrs]) #last resid
                     #TVA.fits[[i]]$pi) #prob at age posterior matrix [iter, ages] ## difficulty because fit 1 at a time
  colnames(sub_samps) <- c(paste0("alpha_", i), 
                           paste0("beta_", i),
                           paste0("Umsy_", i),paste0("Smsy_", i), 
                        #paste0("Sigma_R_", 1:length(names(TVA.fits)), "_", i),                        
                        paste0("S_", (nyrs-A+1):nyrs, "_", i), 
                        paste0("R_", (nRyrs-A+2):nRyrs, "_", i), 
                        paste0("last_resid_", i))#, 
                        #stuff for shared p
                        #paste0("pi_", 1:4))
                  
  pi.samps[,,i] <- TVA.fits[[i]]$pi #store pis to summarise later
  
  for(j in 1:3){
    p.samps[,,i,j] <- TVA.fits[[i]]$p[,nyrs+j, ] #store ps for last 3 nRyrs to summarise later
  }
  samps <- cbind(samps, sub_samps) #cbind posteriors of parms we care about
}

#get median of pis and ps for all pops 
  #(i.e. take the median across posterior slices by "page" of array)
median.pi.samps <- apply(pi.samps, c(1,2), median)
colnames(median.pi.samps) <- paste0("pi_", 1:4)

median.p.samps <- NULL
for(i in 1:3){
  median.p.samps <- cbind(median.p.samps, apply(p.samps[,,,i], c(1,2), median))
}
median.p.samps <- median.p.samps[,  c(1,5,9, 2,6,10, 3,7,11, 4,8,12)] #rearrange to match kusko
colnames(median.p.samps) <- paste("p", (nyrs+1):nRyrs, rep(1:4, each=3), sep = "_")




#get sigma_R varcov matrix 
for(i in 1:length(names(TVA.fits))){
  samps <- cbind(samps, TVA.fits[[i]]$sigma_R) #bind un-correlated sigma_R to samps
}
colnames(samps)[85:92] <- paste0("Sigma_R_", 1:8, "_", 1:8) #just the diagonal...

#populate the off diagonal of the sigma_R "matrix"
samps <- cbind(samps, matrix(0, nrow(samps), 8*8-8)) 

colnames(samps)[93:ncol(samps)] #need to figure out how to name these MFers and if order matters

paste0("Sigma_R_", 1:8, "_", rep(1:8, 8)) #hack to get all the sigma R names...
  #^ could create full empty matrix of sigma_R's as 0's then populate them in a loop with i = i+8 or something

#bind all samps into one object 

#------------------------------------------------------------------------------#
# Multi-stock simulation function with alternative structural forms
#------------------------------------------------------------------------------#
# ny <- the number of years
# vcov.matrix <- process error variance-covariance matrix
# phi <- the expected correlation through time
# mat <- maturation schedule
# alpha <- sub-stock productivity (NOT in log space)
# beta <- sub-stock density dependence 
# sub <- subsistence requirement
# com <- maximum commercial harvest
# egfloor <- escapement goal
# pm.yr <- year of simulation that pms start to be calculated over
# for.error <- forecast error (CV)
# OU <- outcome uncertainty (CV)
# Rec <- estimated recruits from last years of empirical data 
# Spw <- estimated spawners from last years of empirical data
# lst.resid <- estimated recruitment deviation from last year of empirical data
# SR_rel <- structural form of the SR relationship ("Ricker" or "Beverton-Holt")
# BH.alpha.CV <- magnitude (amplitude) of environmental forcing on alpha if SR_rel = "Beverton-Holt"
# period <- period of enviro forcing cycle if SR_rel = "Beverton-Holt"
# dir.SR <- flag for directional change in SR parameters ("Y" or "N")
# SR.devs <- deviations in SR parameters by time step if dir.SR == "Y"
# expan <- expansion of system to account for unmodelled spawning populations

process = function(ny,vcov.matrix,phi,mat,alpha,beta,sub,com,egfloor,pm.yr,for.error,
                   OU,Rec,Spw,lst.resid,SR_rel,BH.alpha.CV,period,dir.SR,SR_devs,expan){
  ns = length(alpha) #number of sub-stocks
  for.error = for.error
  OU = OU
  m.alpha <- alpha
  m.beta <- beta
  
  # create vectors of time varying alpha
  if (SR_rel == "Beverton-Holt"){ 
    beta.tim <- (alpha/beta)*exp(-1)
    alpha.time <- matrix(NA,ny,length(alpha))
    for (t in 1:ny){
      alpha.time[t,] <- sin(2*pi*(t/period))*((alpha + (alpha * BH.alpha.CV)) - alpha) + alpha
    }
  }
  
  #Create recruitment deviations that are correlated among stocks 
  epi = rmvnorm(ny, sigma= vcov.matrix)
  
  #Build time series of Spawners (S), abundance of returning spawners pre-harvest
  # (N), and the component of the residual that is correlated throught time (v)
  R = t(matrix(0,ns,ny))
  S = R * (1-0)
  v = R; v[,]=0
  R[1:3,]=Rec
  N = array(0,dim=c(ny,4,ns))
  Ntot = R; Ntot[,]=0
  H = Ntot; S = Ntot
  S[4:7,] = Spw
  predR = Ntot
  
  # populate first few years with realized states
  R[4,] = alpha[]*S[4,]*exp(-beta[]*S[4,]+(phi*lst.resid)+epi[4,])
  predR[4,] = alpha[]*S[4,]*exp(-beta[]*S[4,])
  v[4,] = log(R[4,])-log(predR[4,])
  v[v[,]=='NaN'] <- 0
  
  for(i in 5:7){
    R[i,] = alpha[]*S[i,]*exp(-beta[]*S[i,]+phi*v[i-1,]+epi[i,])
    predR[i,] = alpha[]*S[i,]*exp(-beta[]*S[i,])
    v[i,] = log(R[i,])-log(predR[i,])
    v[v[,]=='NaN'] <- 0	
  }
  
  N[4:7,1,]=R[4:7-(3),] * mat[1]
  N[5:7,2,]=R[5:7-(4),] * mat[2]
  N[6:7,3,]=R[6:7-(5),] * mat[3]
  N[7,4,]=R[7-(6),] * mat[4]
  
  # Loop through years of simulation	
  for(i in (7+1):ny){
    N[i,1,] = R[i-(4),] * mat[1]
    N[i,2,] = R[i-(5),] * mat[2]
    N[i,3,] = R[i-(6),] * mat[3]
    N[i,4,] = R[i-(7),] * mat[4]
    Ntot[i,] = colSums(N[i,,])
    
    # apply harvest control rule
    run.size <- sum(Ntot[i,])
    if(is.na(run.size)==TRUE){run.size <- 0}
    if(run.size > 999000) {run.size <- 1000000}
    HR.all = sub_hcr(sub,com,egfloor, run.size,for.error)
    HR_adj = 1
    realized.HR <- (HR.all*HR_adj); realized.HR[realized.HR < 0] <- 0; realized.HR[realized.HR > 1] <-1
    outcome_error <- (1+rnorm(1,0,OU))
    H[i,] =  realized.HR*Ntot[i,]*ifelse(outcome_error<0, 0, outcome_error) 
    S_exp = Ntot[i,]-H[i,]
    S_exp[S_exp<0] = 0
    S_exp[S_exp<50] = 0
    S[i,] = S_exp
    if (dir.SR == "T") {
      alpha <- m.alpha* SR_devs[i,1,]
      beta <- m.beta* SR_devs[i,2,]
    }
    
    # predict recruitment
    if (SR_rel == "Ricker"){
      R[i,] = alpha[]*S[i,]*exp(-beta[]*S[i,]+phi*v[i-1,]+epi[i,])
      predR[i,] = alpha[]*S[i,]*exp(-beta[]*S[i,])
      v[i,] = log(R[i,])-log(predR[i,])
      v[v[,]=='NaN'] <- 0
    }
    
    if (SR_rel == "Beverton-Holt"){
      R[i,] = alpha.time[i,]*S[i,]/(1+(alpha.time[i,]/beta.tim[])*S[i,])*exp(phi*v[i-1,]+epi[i,])
      predR[i,] = alpha.time[i,]*S[i,]/(1+(alpha.time[i,]/beta.tim[])*S[i,])
      v[i,] = log(R[i,])-log(predR[i,])
      v[v[,]=='NaN'] <- 0
    }
  }
  
  # Output
  # Performance measures:
  #	1: escapement
  #	2: harvest
  #	3: harvest rate
  #	4: predicted overfished
  #	5: predicted trending towards extinction
  #	6: empirical extinction
  #	7: proportion of years failed to meet subsistence goal
  #	8: CV in harvest
  #	9: proportion of tributary goals met
  
  pms <- matrix(NA,1,9) 
  
  S[S[,]=='NaN'] <- 0
  Ntot[Ntot[,]=='NaN'] <- 0
  over <- matrix(NA,length(alpha))
  ext <- matrix(NA,length(alpha))
  ext.emp <-ext
  trib.gl <-ext
  harvest_rate <- (H[pm.yr:ny,]/Ntot[pm.yr:ny,])[,1]
  harvest_rate[harvest_rate>1] <- 1
  harvest_rates <- (H[pm.yr:ny,]/Ntot[pm.yr:ny,])
  harvest_rates[harvest_rates>1] <- 1
  harvest_rate[harvest_rate=='NaN']<-1
  harvest_rates[harvest_rates=='NaN']<-1
  Smax <- round((m.alpha/m.beta)/m.alpha,digits=0)  
  ln.alpha <- log(m.alpha)
  Smsy <- round((ln.alpha*(0.5-0.07* ln.alpha))/m.beta)
  for(j in 1:length(alpha)){
    over[j] <- SC.eq(median(harvest_rates[,j]),alpha[j],beta[j])[3]
    ext[j] <- SC.eq(median(harvest_rates[,j]),alpha[j],beta[j])[4]
    ext.emp[j] <- ifelse(median(S[(ny-pm.yr):ny,j]) < ((log(alpha)/beta)*0.05)[j],1,0) # less than 5% of unfished biomass/abundance
    trib.gl[j] <- ifelse(median(S[(ny-pm.yr):ny,j]) >= (Smsy[j]),1,0)
  }
  
  pms[,1] <- (sum(S[pm.yr:ny,])/(ny - pm.yr +1)) * expan
  pms[,2] <- (sum(H[pm.yr:ny,])/(ny - pm.yr +1)) * expan
  pms[,3] <- median(harvest_rate)
  pms[,4] <- sum(over)/length(alpha)
  pms[,5] <- sum(ext)/length(alpha)
  pms[,6] <- sum(ext.emp)/length(alpha)
  pms[,7] <- sum(rowSums(H[pm.yr:ny,]) < (sub*0.90))/(ny - pm.yr +1)
  pms[,8] <- sd(H[pm.yr:ny,])/mean(H[pm.yr:ny,]) 
  pms[,9] <- sum(trib.gl)/length(alpha) 
  
  if (SR_rel == "Beverton-Holt"){
    list(S=S[,],R=R[,], N=Ntot[,],H=H[,],BH_alpha = alpha.time, BH_beta = beta.tim, PMs=pms)}
  else{
    list(S=S[,],R=R[,], N=Ntot[,],H=H[,],PMs=pms)}
  
}
