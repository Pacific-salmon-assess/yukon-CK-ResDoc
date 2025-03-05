#common data to be read in ---------------------------------------------------------------
harvest <- read.csv(here("analysis/data/raw/harvest-data.csv")) |>
  dplyr::rename(stock = population, 
                harv_cv = cv)

sp_har <- read.csv(here("analysis/data/raw/esc-data.csv")) |>
  dplyr::rename(spwn = mean, 
                spwn_cv = cv) |>
  left_join(harvest, by = c("stock", "year")) |>
  dplyr::rename(CU = stock) |>
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

fem_S_comps <- read.csv(here("analysis/data/raw/female_age_comps.csv"))[,3:6] # female spawner age composition
fem_fec_comps <- read.csv(here("analysis/data/raw/female_length_comps_eggs.csv")) # female spawner eggs by age
fem_egg_mass_comps <- read.csv(here("analysis/data/raw/female_length_comps_egg_mass.csv")) # female spawner eggs by age

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

# adapted from BC's kusko code (https://github.com/brendanmichaelconnors/Kusko-harvest-diversity-tradeoffs/blob/master/functions.R#L237)


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
  ns = length(unique(sp_har$CU))
  
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
# HCR <- which pre-determined HCR are we using
# ny <- the number of years
# vcov.matrix <- process error variance-covariance matrix
# phi <- the expected correlation through time
# mat <- maturation schedule
# alpha <- sub-stock productivity (NOT in log space)
# beta <- sub-stock density dependence 
# pm.yr <- year of simulation that pms start to be calculated over
# for.error <- forecast error (CV)
# OU <- outcome uncertainty (CV)
# Rec <- estimated recruits from last years of empirical data 
# Spw <- estimated spawners from last years of empirical data
# lst.resid <- estimated recruitment deviation from last year of empirical data
# ER <- fixed exploitation rate (illustrative)

process = function(HCR,ny,vcov.matrix,mat,alpha,beta,pm.yr,for.error,OU,Rec,Spw,
                   lst.resid, ER){
  ns <- length(alpha) #number of sub-stocks
  for.error <- for.error
  OU <- OU
  m.alpha <- alpha
  m.beta <- beta
  
  #Create recruitment deviations that are correlated among stocks 
  epi <- rmvnorm(ny, sigma= vcov.matrix)
  #Build time series of Spawners (S), abundance of returning spawners pre-harvest
  # (N), and the component of the residual that is correlated through time (v)
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
    if(HCR == "no.fishing"){HR.all <- 0} ##BC double check this should be HR.all like the sub_hcr() in kusko
    if(HCR == "status.quo"){
      catch <- ifelse(run.size<=42000, 0, run.size-42000)
      if(run.size==0){HR.all <- 0}else{
        HR.all <- catch/run.size}
      }
    if(HCR == "fixed.ER"){
      if(run.size==0){ER <- 0}
      catch <- run.size*ER
      HR.all <- ER
    }
    
    HR_adj <- 1 ##what is this? harvest adjuster? omit if not necessary?
    realized.HR <- (HR.all*HR_adj); realized.HR[realized.HR < 0] <- 0; realized.HR[realized.HR > 1] <-1
    outcome_error <- (1+rnorm(1,0,OU))
    H[i,] <- realized.HR*Ntot[i,]*ifelse(outcome_error<0, 0, outcome_error) 
    S_exp <- Ntot[i,]-H[i,]
    S_exp[S_exp<0] <- 0  ##cutting out small and negative spawner obs? add comment
    S_exp[S_exp<50] <- 0
    S[i,] <- S_exp
    
    # predict recruitment
    R[i,] <- alpha[]*S[i,]*exp(-beta[]*S[i,]+v[i-1,]+epi[i,]) 
    predR[i,] <- alpha[]*S[i,]*exp(-beta[]*S[i,])
    v[i,] <- log(R[i,])-log(predR[i,])
    v[v[,]=='NaN'] <- 0
  }
  
  # Output
  # Performance measures:
  #	1: escapement
  #	2: harvest
  #	3: harvest rate (associated with REALIZED harvest, i.e. including outcome uncertainty)
  #	4: CV in harvest
  # 5: number of CUs below LRP at end of sim
  # 6: number of CUs between RPs " "
  # 7: number of CUs above USR at ""

  pms <- matrix(NA,1,7) 
  
  S[S[,]=='NaN'] <- 0
  Ntot[Ntot[,]=='NaN'] <- 0
#  over <- matrix(NA,length(alpha)) ##don't think I need any of this stuff I just commented out?
#  ext <- matrix(NA,length(alpha))
#  ext.emp <-ext
#  trib.gl <-ext
  harvest_rate <- (H[pm.yr:ny,]/Ntot[pm.yr:ny,])[,1] ##not really sure what's going on here, why take only the first CU (i.e. column)?
  harvest_rate[harvest_rate>1] <- 1
  harvest_rates <- (H[pm.yr:ny,]/Ntot[pm.yr:ny,])
  harvest_rates[harvest_rates>1] <- 1
  harvest_rate[harvest_rate=='NaN']<-1 ## i.e. fully harvested because Nan means run.size = 0?
  harvest_rates[harvest_rates=='NaN']<-1
    ##^why make these = 1 and not 0? if extinct you arent harvesting at all... 
  Smax <- round((m.alpha/m.beta)/m.alpha,digits=0)  
  ln.alpha <- log(m.alpha)
  Smsy <- round((ln.alpha*(0.5-0.07* ln.alpha))/m.beta)
 # for(j in 1:length(alpha)){
#    over[j] <- SC.eq(median(harvest_rates[,j]),alpha[j],beta[j])[3] #if overfished == 1 
#    ext[j] <- SC.eq(median(harvest_rates[,j]),alpha[j],beta[j])[4] #if extinct == 1
#    ext.emp[j] <- ifelse(median(S[(ny-pm.yr):ny,j]) < ((log(alpha)/beta)*0.05)[j],1,0) # less than 5% of unfished biomass/abundance
#    trib.gl[j] <- ifelse(median(S[(ny-pm.yr):ny,j]) >= (Smsy[j]),1,0)
#  }
  
  pms[,1] <- (sum(S[pm.yr:ny,])/(ny - pm.yr +1))# * expan ## got rid of this. unnecessary?
  pms[,2] <- (sum(H[pm.yr:ny,])/(ny - pm.yr +1))# * expan
  pms[,3] <- median(harvest_rate)
  pms[,4] <- sd(H[pm.yr:ny,])/mean(H[pm.yr:ny,])
  #"status" - how many CUs are in each zone IN THE FINAL YEAR? ##OPEN FOR ADJUSTMENTS! 
  pms[,5] <- sum(S[ny,] < 0.2*Smax) #1 #below 20% Smax (i.e. Smsr) 
  pms[,6] <- sum(S[ny,] > 0.2*Smax & S[ny,] < Smax)#between 20% and Smax
  pms[,7] <- sum(S[ny,] > Smax) #above Smax
  
  list(S=S[,],R=R[,], N=Ntot[,],H=H[,],PMs=pms) ##I can just pull H right here for the harvest timeseries
}

#------------------------------------------------------------------------------#
# functions to estimate spawners associated with maximum yield or recruits, 
#   based on yield per recruit under fishing mortality Fmax under a give age/sex 
#   class-specific reproductive output, and  probability of returning by age/sex
#------------------------------------------------------------------------------#
calcCeq <- function(Fmax)
{
  U_as <- 1-exp(-Fmax) # assume all ages/sexes equally vulnerable sp set v_as = 1
  reproOutputPerSpawner <- sum((1-U_as)*z_as*eta_as)
  alpha <- exp(ln_a)
  Req <- log(alpha*reproOutputPerSpawner) / (beta*reproOutputPerSpawner)# equilibrium population fished at Fmax
  Neq_as <- Req * eta_as # age/sex structured abundance
  Ceq_as <- Neq_as * U_as # age/sex structured harvest
  Seq_as <<- Neq_as * (1-U_as) # age/sex structured escapement
  sum(Ceq_as)
}

calcNeq <- function(Fmax)
{
  U_as <- 1-exp(-Fmax) # assume all ages/sexes equally vulnerable sp set v_as = 1
  reproOutputPerSpawner <- sum((1-U_as)*z_as*eta_as)
  alpha <- exp(ln_a)
  Req <- log(alpha*reproOutputPerSpawner) / (beta*reproOutputPerSpawner)# equilibrium population fished at Fmax
  Neq_as <- Req * eta_as # age/sex structured abundance
  Ceq_as <- Neq_as * U_as # age/sex structured harvest
  Seq_as <<- Neq_as * (1-U_as) # age/sex structured escapement
  sum(Neq_as)
}