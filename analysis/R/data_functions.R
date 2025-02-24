#common data to be read in ---------------------------------------------------------------
harvest <- read.csv(here("analysis/data/raw/harvest-data.csv")) |>
  dplyr::rename(stock = population, 
                harv_cv = cv)

sp_har <- read.csv(here("analysis/data/raw/esc-data.csv")) |>
  dplyr::rename(spwn = mean, 
                spwn_cv = cv) |>
  select(-obs, - se) |>
  left_join(harvest, by = c("stock", "year")) |>
  dplyr::rename(cu = stock) |>
  mutate(N = spwn+harv)

ages <- read.csv(here("analysis/data/raw/run-age-comp.csv")) |>
  filter(Year > 1984,
         Year < 2025)

A_obs <- ages |>
  select(a4:a7) |>
  as.matrix()

a_min <- 4
a_max <- 7 
nyrs <- max(sp_har$year)-min(sp_har$year)+1 #number of years of observations
A <- a_max - a_min + 1 #total age classes
nRyrs <- nyrs + A - 1 #number of recruitment years: unobserved age classes at ? to predict last year of spawners

rm(harvest, ages)

#functions -------------------------------------------------------------------------------
my.ggsave <- function(filename = default_name(plot), plot = last_plot(), 
                      width= 9, height = 5.562, dpi= 180){
  ggsave(filename=filename, plot = last_plot(), width=width, height=height, dpi=dpi, bg="white")
}

# benchmark functions ---
get_Smsy <- function(a, b){
  Smsy <- (1 - lambert_W0(exp(1 - a))) / b
  if(Smsy <0){Smsy <- 0.001} #dumb hack for low draws so Smsy doesnt go negative
  return(Smsy)
}

get_Sgen <- function(a, b, int_lower, int_upper, Smsy){
  fun_Sgen <- function(Sgen, a, b, Smsy) {Sgen * a * exp(-b * Sgen) - Smsy}
  Sgen <- uniroot(fun_Sgen, interval=c(int_lower, int_upper), a=a, b=b, Smsy=Smsy)$root
  return(Sgen)
}

#from BC's kusko code (https://github.com/brendanmichaelconnors/Kusko-harvest-diversity-tradeoffs/blob/master/functions.R#L237)

#------------------------------------------------------------------------------#
# Subsistence harvest control rule function
#------------------------------------------------------------------------------#
# sub < - subsistence requirement
# egfloor <- escapement goal
# run <- run size (i.e., pre-harvest recruitment) 
# com <- maximum commercial harvest
# for.error <- forecast error 
# OU <- outcome uncertainty
sub_hcr = function(sub, com, egfloor, run,for.error){
  run.est <- run * rlnorm(1, 0, for.error); if(is.na(run.est)==TRUE){run.est <- run};if(is.na(run)==TRUE){run <- 0}
  if(run.est - egfloor <= 0){hr = 0}
  if(run.est > egfloor){
    if((run.est - egfloor) <= (sub)){hr = (run - egfloor)/run}
    if((run.est - egfloor) > (sub)){
      if((run.est - egfloor) > (sub + com)){hr = (sub + com)/run}
      if((run.est - egfloor) <= (sub + com)){hr = (sub + (run - egfloor - sub))/run}
    }
  }
  if(hr < 0){hr=0}
  if(hr >1 ){hr=1}
  return(hr)
}

#------------------------------------------------------------------------------#
# Status function (to estimate whether stock is overfished or predicted to go 
#  extinct at a given harvest rate, over the long-term)
#------------------------------------------------------------------------------#
# U <- harvest rate
# a <- productivity (Ricker a parameter)
# b <- density dependence (Ricker beta parameter)
SC.eq <- function(U,a,b){
  a <- log(a)
  S.eq <- max(0,(a-(-log(1-U)))/b)
  C.eq <- max(0,((a-(-log(1-U)))/b)*exp(a-b*((a-(-log(1-U)))/b))-((a-(-log(1-U)))/b))
  OF <- ifelse(U>0.5*a-0.07*a^2,1,0)
  EX <- ifelse(S.eq==0,1,0)
  return(c(S.eq,C.eq,OF,EX))
}

#take a slice of the posterior ---
process.iteration = function(samp) {
  # 1.) extract names
  nms = names(samp)
  A = 4
  ns = 8
  # 2.) extract elements according to the names and put them into the appropriate data structure
  
  # parameters
  alpha = unname(samp[substr(nms, 1, 5) == "alpha"])
  beta = unname(samp[substr(nms, 1, 4) == "beta"])
  last_resid = unname(samp[substr(nms, 1, 10) == "last_resid"])
  #phi = unname(samp["phi"])
  #Sigma_R = matrix(samp[substr(nms, 1, 7) == "Sigma_R"], ns, ns) 
  pis = c(as.numeric(samp["pi_1"]), as.numeric(samp["pi_2"]), as.numeric(samp["pi_3"]), as.numeric(samp["pi_4"]))
  
  # states
  S = matrix(samp[substr(nms, 1, 2) == "S_"], A, ns)
  R = matrix(samp[substr(nms, 1, 2) == "R_"], A - 1, ns)
  
  # 3.) create output list
  output = list(
    alpha = alpha,
    beta = beta,
    last_resid = last_resid,
    S = S,
    R = R,
    pis = pis
  )
  
  # 4.) return output
  return(output)
  
}

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

process = function(ny,vcov.matrix,phi=NULL,mat,alpha,beta,sub,com,egfloor,pm.yr,
                   for.error,OU,Rec,Spw,lst.resid,SR_rel,BH.alpha.CV=NULL,
                   period=NULL,dir.SR, SR_devs=NULL, expan){
  ns <- length(alpha) #number of sub-stocks
  for.error <- for.error
  OU <- OU
  m.alpha <- alpha
  m.beta <- beta

  
  #Create recruitment deviations that are correlated among stocks 
  epi <- rmvnorm(ny, sigma= vcov.matrix)
  #Build time series of Spawners (S), abundance of returning spawners pre-harvest
  # (N), and the component of the residual that is correlated throught time (v)
  R <- t(matrix(0,ns,ny))
  S <- R * (1-0)
  v <- R; v[,] <- 0
  R[1:3,] <- Rec
  N <- array(0,dim=c(ny,4,ns))
  Ntot <- R; Ntot[,]<-0
  H <- Ntot; S <- Ntot
  S[4:7,] <- Spw
  predR <- Ntot
  
  # populate first few years with realized states
  R[4,] <- alpha[]*S[4,]*exp(-beta[]*S[4,]+(lst.resid)+epi[4,])
  predR[4,] <- alpha[]*S[4,]*exp(-beta[]*S[4,])
  v[4,] <- log(R[4,])-log(predR[4,])
  v[v[,]=='NaN'] <- 0
  
  for(i in 5:7){
    R[i,] <- alpha[]*S[i,]*exp(-beta[]*S[i,]+v[i-1,]+epi[i,])
    predR[i,] <- alpha[]*S[i,]*exp(-beta[]*S[i,])
    v[i,] <- log(R[i,])-log(predR[i,])
    v[v[,]=='NaN'] <- 0	
  }
  
  N[4:7,1,] <- R[4:7-(3),] * mat[1]
  N[5:7,2,] <- R[5:7-(4),] * mat[2]
  N[6:7,3,] <- R[6:7-(5),] * mat[3]
  N[7,4,] <- R[7-(6),] * mat[4]
  
  # Loop through years of simulation	
  for(i in (7+1):ny){
    N[i,1,] <- R[i-(4),] * mat[1]
    N[i,2,] <- R[i-(5),] * mat[2]
    N[i,3,] <- R[i-(6),] * mat[3]
    N[i,4,] <- R[i-(7),] * mat[4]
    Ntot[i,] <- colSums(N[i,,])
    
    # apply harvest control rule ## needs changing based on Yukon HCR (no fishing, 42k esc goal, havest above)
    run.size <- sum(Ntot[i,])
    if(is.na(run.size)==TRUE){run.size <- 0}
    if(run.size > 999000) {run.size <- 1000000}
    HR.all <- sub_hcr(sub,com,egfloor, run.size,for.error)
    HR_adj <- 1
    realized.HR <- (HR.all*HR_adj); realized.HR[realized.HR < 0] <- 0; realized.HR[realized.HR > 1] <-1
    outcome_error <- (1+rnorm(1,0,OU))
    H[i,] <- realized.HR*Ntot[i,]*ifelse(outcome_error<0, 0, outcome_error) 
    S_exp <- Ntot[i,]-H[i,]
    S_exp[S_exp<0] <- 0
    S_exp[S_exp<50] <- 0
    S[i,] <- S_exp
    if (dir.SR == "T") {
      alpha <- m.alpha* SR_devs[i,1,]
      beta <- m.beta* SR_devs[i,2,]
    }
    
    # predict recruitment
    if (SR_rel == "Ricker"){ ##dump if statement - this is the only one
      R[i,] <- alpha[]*S[i,]*exp(-beta[]*S[i,]+v[i-1,]+epi[i,]) 
      predR[i,] <- alpha[]*S[i,]*exp(-beta[]*S[i,])
      v[i,] <- log(R[i,])-log(predR[i,])
      v[v[,]=='NaN'] <- 0
    }
  }
  
  # Output
  # Performance measures: ##set up to describe FULL sim, not like each time step for pinks 
  #	1: escapement
  #	2: harvest
  #	3: harvest rate (associated with REALIZED harvest, i.e. including outcome uncertainty)
  ## new 4: which zone is the CU in? vector (length(CU)) of what "zone" the CU is in (below LRP, between, or above USR)
    #^at the END of the simulation, what is your status (sim.yrs-5):sim.yrs
  #	8: CV in harvest
  
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
  
 ## if (SR_rel == "Beverton-Holt"){
##    list(S=S[,],R=R[,], N=Ntot[,],H=H[,],BH_alpha = alpha.time, BH_beta = beta.tim, PMs=pms)}
  #else{
    list(S=S[,],R=R[,], N=Ntot[,],H=H[,],PMs=pms)#}
  
}