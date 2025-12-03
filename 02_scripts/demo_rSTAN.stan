// Simple direct model - let the data determine the relationship
data {
  int<lower=1> N;
  int<lower=1> N_obs;
  int<lower=0> N_miss;
  int<lower=1> N_months;
  int<lower=1, upper=N> obs_idx[N_obs];
  int<lower=1, upper=N> miss_idx[N_miss];
  vector[N_obs] pds_obs;
  vector[N_obs] pds_se_obs;
  vector[N] rcn_obs;
  vector<lower=0>[N] rcn_se_obs;
  int<lower=1, upper=N_months> month[N];
}

parameters {
  // Separate beta and intercept for each month
  vector[N_months] alpha_month;             // intercept per month
  vector[N_months] beta_month;              // slope per month
  
  // Latent PDS for missing observations
  vector<lower=0>[N_miss] pds_latent_miss;
  
  // Residual variation
  real<lower=0> sigma;
}
transformed parameters {
  vector[N] pds_expected;
  
  // Linear model: PDS = alpha + beta * RCN (separate params per month)
  for (t in 1:N) {
    pds_expected[t] = alpha_month[month[t]] + beta_month[month[t]] * rcn_obs[t];
  }
}

model {
  // Very weak priors
  alpha_month ~ normal(0, 200);
  beta_month ~ normal(0, 2);
  sigma ~ normal(0, 100);
  
  // Likelihood for observed PDS
  for (i in 1:N_obs) {
    int t = obs_idx[i];
    pds_obs[i] ~ normal(pds_expected[t], sqrt(square(pds_se_obs[i]) + square(sigma)));
  }
  
  // Prior for missing PDS
  pds_latent_miss ~ normal(pds_expected[miss_idx], sigma);
}

generated quantities {
  vector[N] pds_full;
  vector[N] pds_pred;
  vector[N] pds_mean;
  real mean_beta = mean(beta_month);
  real mean_alpha = mean(alpha_month);
  
  pds_mean = pds_expected;
  
  // Combine observed and missing
  for (i in 1:N_obs) {
    pds_full[obs_idx[i]] = pds_obs[i];
  }
  for (j in 1:N_miss) {
    pds_full[miss_idx[j]] = pds_latent_miss[j];
  }
  
  // Posterior predictive
  for (t in 1:N) {
    pds_pred[t] = normal_rng(pds_expected[t], sigma);
  }
}
