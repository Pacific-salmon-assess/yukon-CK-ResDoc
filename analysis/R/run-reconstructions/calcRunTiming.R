# calcRunTiming.R

calcRunTiming <- function( n_dtg=NULL )
{
  if( is.null(n_dtg) )
  {
    load("data/chinookYkData.Rdata")
    n_dtg <- apply( chinookYkData$n_sdtg[ , ,-(1:3), ], 2:4, sum )
  }
  
  days <- as.numeric(names(n_dtg[ ,1,1]))
  nD <- length(days)
  nT <- dim(n_dtg)[2]
  nG <- dim(n_dtg)[3]
  
  mu_tg <- matrix( data=NA, nrow=nT, ncol=nG )
  sd_tg <- matrix( data=NA, nrow=nT, ncol=nG )
  
  for( g in 1:nG )
  {
    for( t in 1:nT )
    {
      z <- NULL
      if( sum(n_dtg[ ,t,g],na.rm=1)>0 )
      {
        for( d in 1:nD ) 
        {
          if( !is.na(n_dtg[d,t,g]) )
            z <- append(z,rep(days[d],n_dtg[d,t,g]))
        }
        mu_tg[t,g] <- mean(z)
        sd_tg[t,g] <- sd(z)
      }
    }
  }
  return( list( mu_tg=mu_tg, sd_tg=sd_tg ) )
}