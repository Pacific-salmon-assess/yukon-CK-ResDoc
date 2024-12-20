#functions

make_stan_data <- function(cu) {
  
  fn_data <- sp_har[sp_har$cu == cu, ]
  
  S_obs <- sp_har$spwn
  H_obs <- sp_har$N-sp_har$spwn
  H_obs[H_obs<=0] <- 0.01 # Replace any 0's with 0.01, otherwise you get ln(0) which breaks
  
  a_min <- 4
  a_max <- 7 
  nyrs <- nrow(fn_data) #number of spawning years
  A <- a_max - a_min + 1 #total age classes
  nRyrs <- nyrs + A - 1 #number of recruitment years, i.e. add in unobserved age classes at start to predict 1st year of spawners
  
  stan.data <- list("nyrs" = nyrs,
                    "a_min" = a_min,
                    "a_max" = a_max,
                    "A" = A,
                    "nRyrs" = nyrs + A - 1,
                    "A_obs" = A_obs,
                    "S_obs" = S_obs,
                    "H_obs" = H_obs,
                    "S_cv" = fn_data$spwn_cv,
                    "H_cv" = fn_data$harv_cv)
  
  return(stan.data)
}

my.ggsave <- function(filename = default_name(plot), plot = last_plot(), 
                      width= 9, height = 5.562, dpi= 180){
  ggsave(filename=filename, plot = last_plot(), width=width, height=height, dpi=dpi, bg="white")
}

# benchmark functions --------------------------------------------------------------------
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
