// Stan version of state-space spawner-recruitment model with AR-1 process variation

data{
  int T;           // number of calender years
  vector[T] S_obs; // observed spawners
  vector[T] H_obs; // observed harvest
  vector[T] S_cv;  // spawner observation error CV
  vector[T] H_cv;  // harvest observation error CV
  real pSmax_mean;    // mean on prior for Smax
  real pSmax_sig;     // SD on prior for Smax
}

transformed data{
real ln_beta_pr;
real ln_beta_pr_sig;

ln_beta_pr_sig=sqrt(log(1+((1/pSmax_sig)*(1/pSmax_sig))/((1/pSmax_mean)*(1/pSmax_mean)))); //this converts sigma on the untransformed scale to a log scale
ln_beta_pr=log(1/pSmax_mean)-0.5*ln_beta_pr_sig*ln_beta_pr_sig; //convert smax prior to per capita slope - transform to log scale with bias correction

}

parameters{
  vector<lower=0>[T] lnR;             // log recruitment states
  real<lower=0> ln_alpha;                 // log Ricker a
  real ln_beta;                           // Ricker b - let it go negative
  real<lower=0> sigma_R;                  // process error
  real<lower=0> sigma_R0;                 // process error for first a.max years with no spawner link
  real<lower=-1,upper=1> phi;             // lag-1 correlation in process error
  real lnresid_0;                         // first residual for lag-1 correlation in process error
  real<lower=0> mean_ln_R0;               // "true" mean log recruitment in first a.max years with no spawner link
  vector<lower=0.01,upper=0.99>[T] U;  // harvest rate
}

transformed parameters{
  vector<lower=0>[T] N;              // run size states (vestigial from age structured model)
  vector<lower=0>[T] S;              // spawner states
  vector[T] lnS;                     // log spawner states
  vector<lower=0>[T] C;              // catch states
  vector[T] lnC;                     // log catch states
  vector<lower=0>[T] R;             // recruitment states
  vector[T] lnresid;                // log recruitment residuals
  vector[T] lnRm_1;                 // log recruitment states in absence of lag-one correlation
  vector[T] lnRm_2;                 // log recruitment states after accounting for lag-one correlation
  real<lower=0> sigma_R_corr;           // log-normal bias-corrected process error
  real beta;

  beta = exp(ln_beta);
  R = exp(lnR);
  sigma_R_corr = (sigma_R*sigma_R)/2;

  // Calculate returns, spawners and catch by return year
  for(t in 1:T) {
    N[t] = R[t]; //sort of redundant leftover from when there's age structure
    S[t] = N[t] * (1 - U[t]);
    lnS[t] = log(S[t]);
    C[t] = N[t] * U[t];
    lnC[t] = log(C[t]);
  }

  // Ricker SR with AR1 process on log recruitment residuals for years with brood year spawners
  for (t in 1:T) {
    lnresid[t] = 0.0;
    lnRm_1[t] = 0.0;
    lnRm_2[t] = 0.0;
  }

  for (t in 2:T) {
    lnRm_1[t] = lnS[t-1] + ln_alpha - beta * S[t-1];
    lnresid[t] = lnR[t] - lnRm_1[t];
  }

  lnRm_2[2] =  lnRm_1[2] + phi * lnresid_0;

  for (t in 3:T) {
    lnRm_2[t] =  lnRm_1[t] + phi * lnresid[t-1];
  }
}

model{
  // Priors
  ln_alpha ~ normal(1,2);
  ln_beta ~ normal(ln_beta_pr,ln_beta_pr_sig); //per capita capacity parameter - wide prior
  sigma_R ~ normal(1,2);
  lnresid_0 ~ normal(0,20);
  mean_ln_R0 ~ normal(0,20);
  sigma_R0 ~ inv_gamma(2,1);

  // Likelihoods
  // Process model
  lnR[1] ~ normal(mean_ln_R0, sigma_R0);
  lnR[2:T] ~ normal(lnRm_2[2:T], sigma_R_corr);

  // Observation model
  for(t in 1:T){
    U[t] ~ beta(1,1);
    H_obs[t] ~ lognormal(lnC[t], sqrt(log((H_cv[t]^2) + 1)));
    S_obs[t] ~ lognormal(lnS[t], sqrt(log((S_cv[t]^2) + 1)));
  }
}

//posterior predicitive check ala <https://mc-stan.org/docs/stan-users-guide/posterior-predictive-checks.html>
generated quantities{
  array[T] real H_rep;
  array[T] real S_rep;

  H_rep = lognormal_rng(lnC, sqrt(log((H_cv^2) + 1)));
  S_rep = lognormal_rng(lnS, sqrt(log((S_cv^2) + 1)));
}
