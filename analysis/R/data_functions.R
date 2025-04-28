#common data to be read in ---------------------------------------------------------------
harvest <- read.csv(here("analysis/data/raw/harvest-data.csv")) |>
  dplyr::rename(stock = population, 
                harv_cv = cv)

sp_har <- read.csv(here("analysis/data/raw/esc-data.csv")) |>
  dplyr::rename(spwn = mean, 
                spwn_cv = cv) |>
  left_join(harvest, by = c("stock", "year")) |>
  dplyr::rename(CU = stock) |>
  mutate(N = spwn+harv) |>
  arrange(CU)

ages <- read.csv(here("analysis/data/raw/run-age-comp.csv")) |>
  filter(Year >= min(sp_har$year),
         Year <= max(sp_har$year))

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
  Smsy <- (1-lambert_W0(exp(1-a)))/b
  if(Smsy <0){Smsy <- 0.001} #dumb hack for low draws so Smsy doesnt go negative
  return(Smsy)
}

get_Sgen <- function(a, b, int_lower, int_upper, Smsy){
  fun_Sgen <- function(Sgen, a, b, Smsy){Sgen*a*exp(-b*Sgen)-Smsy}
  Sgen <- uniroot(fun_Sgen, interval=c(int_lower, int_upper), a=a, b=b, Smsy=Smsy)$root
  return(Sgen)
}

# below funs adapted from BC's kusko code (https://github.com/brendanmichaelconnors/Kusko-harvest-diversity-tradeoffs/blob/master/functions.R#L237)

#------------------------------------------------------------------------------#
# Status function (to estimate whether stock is overfished or predicted to go 
#  extinct at a given harvest rate, over the long-term)
#------------------------------------------------------------------------------#
# U <- harvest rate
# a <- productivity (Ricker a parameter)
# b <- density dependence (Ricker beta parameter)
SC.eq <- function(U,a,b){ ## Think we can delete this?
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
# Smax <- Smax calculated externally (or other benchmark)
# pm.yr <- year of simulation that pms start to be calculated over
# for.error <- forecast error (lognormal SD)
# OU <- outcome uncertainty (CV)
# Rec <- estimated recruits from last years of empirical data 
# Spw <- estimated spawners from last years of empirical data
# lst.resid <- estimated recruitment deviation from last year of empirical data
# phi <- expected correlation in recruitment deviation from one year to next
# ER <- fixed exploitation rate (illustrative)

process = function(HCR,ny,vcov.matrix,mat,alpha,beta,Smax,pm.yr,for.error,OU,Rec,Spw,
                   lst.resid, phi, ER){

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
  R[1:3,] <- Rec # observed recruitment indexed by brood/spawning year 
  N <- array(0,dim=c(ny,4,ns))
  Ntot <- R; Ntot[,]<-0
  H <- Ntot; S <- Ntot
  S[4:7,] <- Spw
  predR <- Ntot
  error <- matrix(NA, nrow=ny, ncol=4)
  
  # populate first few years with realized states
  R[4,] <- alpha[]*S[4,]*exp(-beta[]*S[4,]+(phi*lst.resid)+epi[4,])
  predR[4,] <- alpha[]*S[4,]*exp(-beta[]*S[4,])
  v[4,] <- log(R[4,])-log(predR[4,])
  v[v[,]=='NaN'] <- 0
  
  for(i in 5:7){
    R[i,] <- alpha[]*S[i,]*exp(-beta[]*S[i,]+phi*v[i-1,]+epi[i,])
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
    
    # apply harvest control rules
    run.size.true <- sum(Ntot[i,])
    if(is.na(run.size.true)){run.size.true <- 0}
    run.size <- rlnorm(1,log(run.size.true),for.error) # forecasted run-size
    if(is.na(run.size)==TRUE){run.size <- 0}
    if(HCR == "no.fishing"){HR.all <- 0}
    if(HCR == "status.quo"){ 
      catch <- ifelse(run.size<=42500, 0, run.size-42500)
      HR.all <- ifelse(run.size==0, 0, catch/run.size)
      if(HR.all > 0.8){       ## Add ER cap (80%) 
          catch <- run.size*0.8
          HR.all <- catch/run.size }}
    if(HCR == "status.quo.cap"){
      catch <- ifelse(run.size<=42500, 0, run.size-42500)
      HR.all <- ifelse(run.size==0, 0, catch/run.size)
      if(HR.all > 0.4){       ## Lower ER cap (40%) 
        catch <- run.size*0.4
        HR.all <- catch/run.size }}
    if(HCR == "moratorium"){
      catch <- ifelse(run.size<=71000, 0, run.size-71000)
      HR.all <- ifelse(run.size==0, 0, catch/run.size)
      if(HR.all > 0.8){       ## ER cap (80%) 
        catch <- run.size*0.8
        HR.all <- catch/run.size }}
    if(HCR == "moratorium.cap"){
      catch <- ifelse(run.size<=71000, 0, run.size-71000)
      HR.all <- ifelse(run.size==0, 0, catch/run.size)
      if(HR.all > 0.4){       ## lower ER cap (40%) 
        catch <- run.size*0.4
        HR.all <- catch/run.size }}
    if(grepl("fixed.ER", HCR)){
      if(run.size==0){ER <- 0}
      catch <- run.size*ER
      HR.all <- catch/run.size}
    if(HCR == "alternative"){
      if(run.size <= 19000) catch <- 0
      if(run.size >= 158333) catch <- run.size*0.4
      if(run.size > 19000 & run.size < 158333){
        dat <- data.frame(R=c(19000,158333), ER=c(0,0.4))
        lin <- lm(ER ~ R, data=dat)
        ER <- coef(lin)[1] + coef(lin)[2]*run.size
        catch <- run.size*ER
      }
      HR.all <- catch/run.size
    }


    
    HR_adj <- 1 ##what is this? harvest adjuster? omit if not necessary? BC: can remove
    realized.HR <- (HR.all*HR_adj); realized.HR[realized.HR < 0] <- 0; realized.HR[realized.HR > 1] <-1
    outcome_error <- (1+rnorm(1,0,OU))
    H[i,] <- realized.HR*Ntot[i,]*ifelse(outcome_error<0, 0, outcome_error) 
    S_exp <- Ntot[i,]-H[i,]
    S_exp[S_exp<0] <- 0  ##cutting out small and negative spawner obs
    S_exp[S_exp<50] <- 0
    S[i,] <- S_exp
    error[i,1] <- HR.all*ifelse(outcome_error<0, 0, outcome_error); error[i,2] <- realized.HR; error[i,3] <- run.size; error[i,4] <- run.size.true
    
    # predict recruitment
    R[i,] <- alpha[]*S[i,]*exp(-beta[]*S[i,]+phi*v[i-1,]+epi[i,]) 
    predR[i,] <- alpha[]*S[i,]*exp(-beta[]*S[i,])
    v[i,] <- log(R[i,])-log(predR[i,])
    v[v[,]=='NaN'] <- 0
    
  } # end years loop
  
  # Output
  # SMU-level Performance measures:
  #	1: avg annual escapement
  #	2: avg annual harvest
  #	3: harvest rate (associated with REALIZED harvest, i.e. including outcome uncertainty) (ER)
  # 4: % of years with no harvest
  # 5: % of yrs harv > 10k (basic needs)
  # 6: number of CUs below LRP at end of sim
  # 7: number of CUs between RPs " "
  # 8: number of CUs above USR at ""
  # 9: number of extinct pops 

  #browser()
    pms <- matrix(NA,1,9) 
  
  # SMU level PMs 
  S[S[,]=='NaN'] <- 0
  H[H[,] == 'NAN'] <- 0 # not sure if this is an appropriate fix
  Ntot[Ntot[,]=='NaN'] <- 0
  harvest_rates <- (H[pm.yr:ny,]/Ntot[pm.yr:ny,])
  harvest_rates[harvest_rates[,]=='NaN'] <- 0
  #Smax <- round((m.alpha/m.beta)/m.alpha,digits=0)
  Smax <- Smax
  ln.alpha <- log(m.alpha)
  Smsy <- round((ln.alpha*(0.5-0.07* ln.alpha))/m.beta)
  pms[,1] <- (sum(S[pm.yr:ny,])/(ny - pm.yr +1)) 
  pms[,2] <- (sum(H[pm.yr:ny,])/(ny - pm.yr +1))
  pms[,3] <- mean(harvest_rates) #was median before
  pms[,4] <- sum(rowSums(H[pm.yr:ny,])==0)/(ny - pm.yr +1) 
  pms[,5] <- sum(rowSums(H[pm.yr:ny,])> 10000)/(ny - pm.yr +1)
  #"status" - how many CUs are in each zone IN THE FINAL YEAR? 
  pms[,6] <- sum(S[ny,] < 0.2*Smax & S[ny,] !=0) 
  pms[,7] <- sum(S[ny,] > 0.2*Smax & S[ny,] < Smax)
  pms[,8] <- sum(S[ny,] > Smax)
  pms[,9] <- sum(S[ny,] ==0)
  
  
  # CU-level performance measures
  # 1: % yrs above LRP (20% Smsr/Smax)
  # 2: % yrs above Smsr/Smax
  
  pms_cu <- matrix(NA,length(beta),2)
  pms_cu[,1] <- colSums(S[pm.yr:ny,] > matrix(0.2*Smax, 21,9, byrow = T))/(ny - pm.yr +1)
  pms_cu[,2] <- colSums(S[pm.yr:ny,] > matrix(Smax, 21,9, byrow = T))/(ny - pm.yr +1)
  
  
  
  list(S=S[,],R=R[,], N=Ntot[,],H=H[,],PMs=pms, PMs_cu=pms_cu, error)
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



# ------------------------------------------------------------------------------#
# Simple simulation to visualize HCRs
# ------------------------------------------------------------------------------#

visualize_HCR <- function(HCRs, max_spwn=400000, int=1000) {
  
  out <- NULL
  
  for(j in HCRs){
    HCR = j
    
    for(i in seq(1,max_spwn, int)) {
      run.size = i
      
      # apply harvest control rules
      if(is.na(run.size)){run.size <- 0}
      if(run.size > 999000) {run.size <- 1000000} 
      if(HCR == "no.fishing"){HR.all <- 0}
      if(HCR == "status.quo"){ 
        catch <- ifelse(run.size<=42500, 0, run.size-42500)
        HR.all <- ifelse(run.size==0, 0, catch/run.size)
        if(HR.all > 0.8){       ## Add ER cap (80%) 
          catch <- run.size*0.8
          HR.all <- catch/run.size }}
      if(HCR == "status.quo.cap"){
        catch <- ifelse(run.size<=42500, 0, run.size-42500)
        HR.all <- ifelse(run.size==0, 0, catch/run.size)
        if(HR.all > 0.4){       ## Lower ER cap (40%) 
          catch <- run.size*0.4
          HR.all <- catch/run.size }}
      if(HCR == "moratorium"){
        catch <- ifelse(run.size<=71000, 0, run.size-71000)
        HR.all <- ifelse(run.size==0, 0, catch/run.size)
        if(HR.all > 0.8){       ## ER cap (80%) 
          catch <- run.size*0.8
          HR.all <- catch/run.size }}
      if(HCR == "moratorium.cap"){
        catch <- ifelse(run.size<=71000, 0, run.size-71000)
        HR.all <- ifelse(run.size==0, 0, catch/run.size)
        if(HR.all > 0.4){       ## lower ER cap (40%) 
          catch <- run.size*0.4
          HR.all <- catch/run.size }}
      if(grepl("fixed.ER", HCR)){
        if(run.size==0){ER <- 0}
        catch <- run.size*ER
        HR.all <- catch/run.size}
      if(HCR == "alternative"){
        if(run.size <= 19000) catch <- 0
        if(run.size >= 158333) catch <- run.size*0.4
        if(run.size > 19000 & run.size < 158333){
          dat <- data.frame(R=c(19000,158333), ER=c(0,0.4))
          lin <- lm(ER ~ R, data=dat)
          ER <- coef(lin)[1] + coef(lin)[2]*run.size
          catch <- run.size*ER
        }
        HR.all <- catch/run.size
      }
      itn <- data.frame(HCR = HCR,
                        run_size=run.size,
                        HR = HR.all,
                        catch=catch)
      out <- bind_rows(out, itn)
    }
  }
  return(out)
}





