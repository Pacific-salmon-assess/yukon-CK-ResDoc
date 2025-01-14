#include <TMB.hpp>
#include <iostream>

// isNA() - replicates is.na() from R
template<class Type>
bool isNA(Type x){
  return R_IsNA(asDouble(x));
}

// square()
template <class Type> 
Type square(Type x){
  return pow(x,2);
}

// Bounded inverse logit
template <class Type> 
Type binvlogit(Type x)
{
  Type lb = -1;
  Type ub = 1;
  return invlogit(x)*(ub-lb) + lb;
}

// posfun() - replicates Dave Fournier's ADMB function
// Compares a quantity x to a threshold epsilon and if below the threshold, 
// increases x above the threshold and increments a penalty
template<class Type>
Type posfun(Type x, Type eps, Type &pen){
  pen += CppAD::CondExpLt(x, eps, Type(0.01) * pow(x-eps,2), Type(0.));
  return CppAD::CondExpGe(x, eps, x, eps/(Type(2.0)-x/eps));
}

// Objective function ------------------------------------------------------ //
template<class Type>
Type objective_function<Type>::operator() ()
{
  // Access multivariate distributions
  using namespace density;

  // Data and controls
  DATA_ARRAY(n_sdtg);     // Observed abundance by stock
  DATA_ARRAY(E_dtg);      // Daily counts
  DATA_VECTOR(I_t);       // Run size indices from mark-recapture
  DATA_VECTOR(day_d);     // Model days (Julian)
  DATA_VECTOR(CV_t);      // Model days (Julian)
  DATA_INTEGER(runRW);
  DATA_SCALAR(runSD);

  // Model dimensions
  int nS = n_sdtg.dim(0);    // Number of stocks
  int nD = n_sdtg.dim(1);    // Number of days
  int nT = n_sdtg.dim(2);    // Number of time steps
  int nG = n_sdtg.dim(3);    // Number of fleets

  // Parameters
  PARAMETER_ARRAY(lnRunSize_st);  // Run size (log scale)
  PARAMETER_VECTOR(lnArrivMu_s);  // Mean date of arrival in initial year (log scale)
  //PARAMETER(lnArrivSD);  // SD around mean Julian date of arrival (log scale)
  PARAMETER_VECTOR(lnArrivSD_s);  // SD around mean Julian date of arrival (log scale)
  PARAMETER_ARRAY(arrivErr_st);   // Process error in mean date of arrival
  //PARAMETER(lnErrSD);    // Process error standard deviation (log scale)
  PARAMETER_VECTOR(lnErrSD_s);    // Process error standard deviation (log scale)
  PARAMETER_ARRAY(logitCor_ss);   // Correlation matrix for process errors (log scale)
  PARAMETER_ARRAY(lnqE_sg);       // Count catchability (log scale)
  PARAMETER_VECTOR(lnqI_s);       // Index catchability
  PARAMETER_ARRAY(lnDisp_tg);

  // Transformed parameters
  array<Type> runSize_st(nS,nT);          // Run size
  vector<Type> arrivMu_s = exp(lnArrivMu_s);  // SD around mean Julian date of arrival
  vector<Type> arrivSD_s = exp(lnArrivSD_s);  // SD around mean Julian date of arrival
  vector<Type> errSD_s = exp(lnErrSD_s);

  // Latent variables
  vector<Type> runSize_t(nT);
  array<Type> mu_st(nS,nT);        // Mean Julian date of arrival
  array<Type> N_dt(nD,nT);         // Daily aggregate abundance
  array<Type> N_dst(nD,nS,nT);     // Daily abundance by stock
  array<Type> rho_dst(nD,nS,nT);   // Proportion of stock arriving each day
  matrix<Type> cor_ss(nS,nS);
  vector<Type> mrSD_t(nT);

  // Model predictions
  array<Type> Ihat_dtg(nD,nT,nG);     // Abundance indices
  vector<Type> mrIhat_t(nT);     // Abundance indices
  array<Type> Phat_sdtg(nS,nD,nT,nG); // Proportions by stock

  // Likelihood objects
  array<Type> nllI_tg(nT,nG);  // Abundance index likelihood
  vector<Type> nllP_g(nG);  // Stock composition likelihood
  Type nllMR = 0;
  Type nlp = 0;             // Prior penalty
  Type varPen = 0;

  // Initialize accumulating variables
  N_dt.fill(0);
  nllI_tg.fill(0);
  nllP_g.fill(0);
  
  // Intialize random walks
  mu_st.col(0) = exp(lnArrivMu_s);

  // In-season abundance dynamics
  for( int t=0; t<nT; t++ )
  {
    // Transform run size
    runSize_st.col(t) = exp( lnRunSize_st.col(t) );

    runSize_t(t) = runSize_st.col(t).sum();

    for( int s=0; s<nS; s++ )
    {

      //nlp -= dnorm( lnRunSize_st(s,t), Type(9), Type(15), TRUE );

      for( int s2=0; s2<nS; s2++ )
        cor_ss(s,s2) = binvlogit(logitCor_ss(s,s2));

      // Random walk in arrival timing
      if( t>0 )
        mu_st(s,t) = mu_st(s,t-1)*exp(arrivErr_st(s,t-1));

      // Relative daily arrivals
      for( int d=0; d<nD; d++ )
        rho_dst(d,s,t) = exp( -square(day_d(d)-mu_st(s,t))/(2*square(arrivSD_s(s))) );
      
      // Convert relative daily arrivals to proportions
      rho_dst.col(t).col(s) /= rho_dst.col(t).col(s).sum();

      // Daily abundance
      N_dst.col(t).col(s) = runSize_st(s,t)*rho_dst.col(t).col(s);

      // Add to aggregate abundance
      N_dt.col(t) += N_dst.col(t).col(s);

      // Store predicted abundance
      for( int g=0; g<nG; g++ )
        Ihat_dtg.col(g).col(t) += exp(lnqE_sg(s,g))*N_dst.col(t).col(s);

    } // next s
    
    // Mark-recapture predicted index
    mrIhat_t(t) = (exp(lnqI_s)*runSize_st.col(t)).sum();

  } // next t


  // Observation model ----------------------------------------------------- //

  for( int t=0; t<nT; t++ )
  {
    for( int g=0; g<nG; g++ )
    {
      for( int d=0; d<nD; d++ )
      {
        // Likelihood for abundance observations
        if( !isNA(E_dtg(d,t,g)) )
        {
          // Neg-binomial will fail if prediction is too close to 0
          // Keep it above 1e-8
          Type mu = posfun( Ihat_dtg(d,t,g), Type(1e-8), varPen );
          // Neg-binomial variance
          Type var = mu+mu*mu*exp(lnDisp_tg(t,g));
//          if( t==16 & g==1 & d<60 )
//           Type var = mu + Type(1e-4);
          // Neg-binomial likelihood
          nllI_tg(t,g) -= dnbinom_robust(E_dtg(d,t,g),log(mu),log(var-mu),TRUE);
        }

        // Predicted stock composition
        for( int s=0; s<nS; s++ )
          Phat_sdtg(s,d,t,g) = exp(lnqE_sg(s,g))*N_dst(d,s,t)/Ihat_dtg(d,t,g);

        // Multinomial likelihood for composition data
        if( !isNA(n_sdtg(0,d,t,g)) )
        {
          vector<Type> N_s = n_sdtg.col(g).col(t).col(d);
          vector<Type> P_s = Phat_sdtg.col(g).col(t).col(d);
          nllP_g(g) -= dmultinom( N_s, P_s, TRUE );
        }
      }
    }

    // Poisson likelihood for mark-recapture indices
    if( !isNA(I_t(t)) )
    {
      mrSD_t(t) = sqrt( log(square(CV_t(t))+1) );
      nllMR -= dnorm( log(I_t(t)),
                      log(mrIhat_t(t)),
                      mrSD_t(t),
                      TRUE );
    }

  }
  
  // Priors ---------------------------------------------------------------- //

  // Diagonal matrix with process error SDs on diagonal
  matrix<Type> D(nS,nS);
  D.fill(0);
  D.diagonal() = errSD_s;

  // Covariance matrix
  matrix<Type> cov_ss = D*cor_ss*D;

  // Multivariate density for process errors
  MVNORM_t<Type> errDens(cov_ss);

  // Negative-log penalty for process errors
  for( int t=0; t<(nT-1); t++ )
  {
    nlp += errDens(arrivErr_st.col(t));
    
    if( runRW==1 )
    {
      if( t < 20 ) // Comment this line out to apply run size prior in every year
      {
          vector<Type> runErr = lnRunSize_st.col(t+1)-lnRunSize_st.col(t);
          nlp -= dnorm( runErr, Type(0), runSD, TRUE ).sum();
      }
    }
  }

  // Total objective function
  Type objFun = nllI_tg.sum() +
                nllP_g.sum() +
                nllMR +
                nlp +
                1e3*varPen;

  // REPORT SECTION -------------------------------------------------------- //
  ADREPORT(runSize_st);
  ADREPORT(runSize_t);
  ADREPORT(errSD_s);
  ADREPORT(mu_st);

  // Data and controls
  REPORT(n_sdtg);
  REPORT(E_dtg);
  REPORT(I_t);
  REPORT(day_d);
  REPORT(CV_t);
  REPORT(mrSD_t);

  // Dimensions
  REPORT(nS);
  REPORT(nG);
  REPORT(nT);
  REPORT(nD);

  // Parameters
  REPORT(lnRunSize_st);
  REPORT(lnArrivMu_s);
  REPORT(lnArrivSD_s);
  REPORT(arrivErr_st);
  REPORT(logitCor_ss);
  REPORT(cor_ss);
  REPORT(lnErrSD_s);
  REPORT(lnqI_s);
  REPORT(lnDisp_tg);

  // Latent parameters
  REPORT(runSize_t);
  REPORT(runSize_st);
  REPORT(N_dt);
  REPORT(N_dst);
  REPORT(mu_st);
  REPORT(rho_dst);
  REPORT(lnqE_sg);
  REPORT(arrivSD_s);
  REPORT(errSD_s);

  // Model predictions
  REPORT(Ihat_dtg);
  REPORT(mrIhat_t);
  REPORT(Phat_sdtg);

  // Objective function
  REPORT(objFun);
  REPORT(nllI_tg);
  REPORT(nllP_g);
  REPORT(nllMR);
  REPORT(nlp);
  REPORT(cov_ss);
  REPORT(varPen);

  return objFun;
}




