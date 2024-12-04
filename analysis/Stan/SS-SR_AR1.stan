// Stan version of age-structured state-space spawner-recruitment model with AR-1 process variation (adapted from Fleischman et al. CJFAS. 2013)

data{
  int nyrs;           // number of calender years
  int a_min;          // minimum age class
  int a_max;          // maximum age class
  int A;              // number of age classes
  int nRyrs;          // number of recruitment years
  int A_obs[nyrs, A]; // observed age composition in counts by age class
  vector[nyrs] S_obs; // observed spawners
  vector[nyrs] H_obs; // observed harvest
  vector[nyrs] S_cv;  // spawner observation error CV
  vector[nyrs] H_cv;  // harvest observation error CV
}

parameters{
  vector<lower=0>[nRyrs] lnR;             // log recruitment states
  real<lower=0> lnalpha;                  // log Ricker a
  real<lower=0> beta;                     // Ricker b
  real<lower=0> sigma_R;                  // process error
  real<lower=0> sigma_R0;                 // process error for first a.max years with no spawner link
  real<lower=-1,upper=1> phi;             // lag-1 correlation in process error
  real lnresid_0;                         // first residual for lag-1 correlation in process error
  real<lower=0> mean_ln_R0;               // "true" mean log recruitment in first a.max years with no spawner link
  vector<lower=0.01,upper=0.99>[nyrs] U;  // harvest rate
  vector<lower=0,upper=1>[3] prob;        // maturity schedule probs
  real<lower=0,upper=1> D_scale;          // governs variability of age proportion vectors across cohorts
  matrix<lower=0.01>[nRyrs, A] g;         // individual year/age class gamma variates for generating age at maturity proportions
}

transformed parameters{
  vector<lower=0>[nyrs] N;              // run size states
  vector<lower=0>[nyrs] S;              // spawner states
  vector[nyrs] lnS;                     // log spawner states
  vector<lower=0>[nyrs] C;              // catch states
  vector[nyrs] lnC;                     // log catch states
  vector<lower=0>[nRyrs] R;             // recruitment states
  vector[nRyrs] lnresid;                // log recruitment residuals
  vector[nRyrs] lnRm_1;                 // log recruitment states in absence of lag-one correlation
  vector[nRyrs] lnRm_2;                 // log recruitment states after accounting for lag-one correlation
  matrix<lower=0>[nyrs, A] N_ta;        // returns by age matrix
  matrix<lower=0, upper=1>[nRyrs, A] p; // age at maturity proportions
  vector<lower=0,upper=1>[4] pi;        // maturity schedule probs
  real<lower=0> D_sum;                  // inverse of D_scale which governs variability of age proportion vectors across cohorts
  vector<lower=0>[A] Dir_alpha;         // Dirichlet shape parameter for gamma distribution used to generate vector of age-at-maturity proportions
  matrix<lower=0, upper=1>[nyrs, A] q;  // age composition by year/age classr matrix

  // Maturity schedule: use a common maturation schedule to draw the brood year specific schedules
  pi[1] = prob[1];
  pi[2] = prob[2] * (1 - pi[1]);
  pi[3] = prob[3] * (1 - pi[1] - pi[2]);
  pi[4] = 1 - pi[1] - pi[2] - pi[3];
  D_sum = 1/D_scale^2;

  for (a in 1:A) {
    Dir_alpha[a] = D_sum * pi[a];
    for (y in 1:nRyrs) {
      p[y,a] = g[y,a]/sum(g[y,]);
    }
  }

  R = exp(lnR);

  // Calculate the numbers at age matrix as brood year recruits at age (proportion that matured that year)
  for (t in 1:nyrs) {
    for(a in 1:A){
      N_ta[t,a] = R[t+A-a] * p[t+A-a,a];
    }
  }

  // Calculate returns, spawners and catch by return year
  for(t in 1:nyrs) {
    N[t] = sum(N_ta[t,1:A]);
    S[t] = N[t] * (1 - U[t]);
    lnS[t] = log(S[t]);
    C[t] = N[t] * U[t];
    lnC[t] = log(C[t]);
  }

  // Calculate age proportions by return year
  for (t in 1:nyrs) {
    for(a in 1:A){
      q[t,a] = N_ta[t,a]/N[t];
    }
  }

  // Ricker SR with AR1 process on log recruitment residuals for years with brood year spawners
  for (i in 1:nRyrs) {
    lnresid[i] = 0.0;
    lnRm_1[i] = 0.0;
    lnRm_2[i] = 0.0;
  }

  for (y in (A+a_min):nRyrs) {
    lnRm_1[y] = lnS[y-a_max] + lnalpha - beta * S[y-a_max];
    lnresid[y] = lnR[y] - lnRm_1[y];
  }

  lnRm_2[A+a_min] =  lnRm_1[A+a_min] + phi * lnresid_0;

  for (y in (A+a_min+1):nRyrs) {
    lnRm_2[y] =  lnRm_1[y] + phi * lnresid[y-1];
  }
}

model{
  // Priors
  lnalpha ~ normal(0,3);
  beta ~ normal(0,1);
  sigma_R ~ normal(0,2);
  lnresid_0 ~ normal(0,20);
  mean_ln_R0 ~ normal(0,20);
  sigma_R0 ~ inv_gamma(2,1); 
  prob[1] ~ beta(1,1);
  prob[2] ~ beta(1,1);
  prob[3] ~ beta(1,1);
  D_scale ~ beta(1,1);

  // Likelihoods
  // Gamma variates for each year and age class which are used to determine age at maturity proportions
  for (y in 1:nRyrs) {
    for (a in 1:A) {
      //g[y,a] ~ gamma(Dir_alpha[a],1);
      target += gamma_lpdf(g[y,a]|Dir_alpha[a], 1);
    }
  }

  // First `a.max` years of recruits, for which there is no spawner link
  lnR[1:a_max] ~ normal(mean_ln_R0, sigma_R0);

  // State model
  lnR[(A+a_min):nRyrs] ~ normal(lnRm_2[(A+a_min):nRyrs], sigma_R);

  // Observation model
  for(t in 1:nyrs){
  //A_obs[t,1:A]) ~ multinomial(q[t,1:A]);
    target += multinomial_lpmf(A_obs[t,1:A]|to_vector(q[t,1:A]));
    U[t] ~ beta(1,1);
    H_obs[t] ~ lognormal(lnC[t], sqrt(log((H_cv[t]^2) + 1)));
    S_obs[t] ~ lognormal(lnS[t], sqrt(log((S_cv[t]^2) + 1)));
  }
}

generated quantities {

}
