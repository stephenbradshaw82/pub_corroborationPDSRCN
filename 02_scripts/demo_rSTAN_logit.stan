// Hierarchical logit-proportion calibration model

// Purpose:
//   Calibrate sparse PDS fishing-boating survey observations against continuous
//   RCN remote-camera boat counts, while enforcing the contextual constraint
//   that expected PDS is a fraction of RCN.
//
// Core model:
//   PDS_t = q_t * RCN_t
//   logit(q_t) = alpha_0 + month_effect[month_t] + year_effect[year_t]
//
// This model intentionally has no area effect as each site is modelled separately.
// Month and year effects are partially pooled through hierarchical priors.
//
// Reporting:
// pds_est_constrained or pds_full_constrained as the primary imputed time series by the construction (0 <= PDS <= RCN).
// pds_pred_raw used for posterior predictive checks only, because realised observation-error draws can exceed RCN.

data {
  int<lower=1> N;                         // all time records to predict over
  int<lower=1> N_obs;                     // records with observed PDS
  int<lower=0> N_miss;                    // records without observed PDS
  int<lower=1> N_months;
  int<lower=1> N_years;

  int<lower=1, upper=N> obs_idx[N_obs];
  int<lower=1, upper=N> miss_idx[N_miss]; // retained for data compatibility

  vector<lower=0>[N_obs] pds_obs;         // observed PDS estimates/counts
  vector<lower=0>[N_obs] pds_se_obs;      // standard error of PDS observation

  vector<lower=0>[N] rcn_obs;             // matched-window RCN counts for calibration,
                                          // or target-window RCN counts for prediction
  vector<lower=0>[N] rcn_se_obs;          // retained for data compatibility; not used

  int<lower=1, upper=N_months> month[N];
  int<lower=1, upper=N_years> year_idx[N];
}

parameters {
  // Overall mean fishing-boating proportion on logit scale.
  real alpha_0;

  // Hierarchical month effects, non-centred.
  real<lower=0> sigma_month;
  vector[N_months] z_month;

  // Hierarchical year effects, non-centred.
  real<lower=0> sigma_year;
  vector[N_years] z_year;

  // Extra residual variation on log1p(PDS) scale.
  real<lower=0> sigma_log;
}

transformed parameters {
  vector[N_months] month_effect_raw;
  vector[N_years] year_effect_raw;
  vector[N_months] month_effect;
  vector[N_years] year_effect;

  vector[N] logit_q;
  vector<lower=0, upper=1>[N] q;
  vector<lower=0>[N] pds_expected;
  vector[N] mu_log1p_pds;

  month_effect_raw = sigma_month * z_month;
  year_effect_raw = sigma_year * z_year;

  // Sum-to-zero centring improves identifiability of alpha_0.
  month_effect = month_effect_raw - mean(month_effect_raw);
  year_effect = year_effect_raw - mean(year_effect_raw);

  for (t in 1:N) {
    logit_q[t] = alpha_0 + month_effect[month[t]] + year_effect[year_idx[t]];
    q[t] = inv_logit(logit_q[t]);

    // Expected PDS is constrained to be between 0 and RCN.
    pds_expected[t] = q[t] * rcn_obs[t];

    // Observation model is applied on log1p scale to keep the likelihood positive and allow errors to scale with abundance.
    mu_log1p_pds[t] = log1p(pds_expected[t]);
  }
}

model {
  // Priors.
  // This prior centres the average fishing-boating proportion at 0.50 before seeing the data, while still allowing wide movement.
  // Readers could modify to site specific priors if justification exists.
  alpha_0 ~ normal(logit(0.50), 1.5);

  // Hierarchical variation in month and year effects on the logit-proportion scale.
  sigma_month ~ normal(0, 1);
  z_month ~ normal(0, 1);

  sigma_year ~ normal(0, 0.75);
  z_year ~ normal(0, 1);

  // Residual process variation on log1p(PDS) scale.
  sigma_log ~ normal(0, 0.5);

  // Likelihood for observed PDS.
  // PDS observation error is approximately propagated to the log1p scale via
  // the delta method: se[log1p(PDS)] ~= se[PDS] / (PDS + 1).
  for (i in 1:N_obs) {
    int t = obs_idx[i];
    real log_pds_i = log1p(pds_obs[i]);
    real log_pds_se_i = fmax(pds_se_obs[i] / (pds_obs[i] + 1.0), 1e-6);

    log_pds_i ~ normal(
      mu_log1p_pds[t],
      sqrt(square(log_pds_se_i) + square(sigma_log))
    );
  }
}

generated quantities {
  vector[N] q_est;                        // estimated fishing-boating proportion
  vector[N] pds_est_constrained;          // constrained expected PDS = q * RCN
  vector[N] pds_full_constrained;         // observed PDS where observed; constrained means where missing

  vector[N] pds_pred_raw;                 // positive posterior predictive draw; may exceed RCN
  vector[N] pds_pred_capped;              // capped posterior predictive draw for operational use

  real mean_q;
  real mean_month_effect;
  real mean_year_effect;

  q_est = q;
  pds_est_constrained = pds_expected;
  mean_q = mean(q_est);
  mean_month_effect = mean(month_effect);
  mean_year_effect = mean(year_effect);

  // Primary imputed series: keep observed PDS where available and use the
  // constrained model expectation where PDS is missing.
  for (t in 1:N) {
    pds_full_constrained[t] = pds_est_constrained[t];
  }

  for (i in 1:N_obs) {
    pds_full_constrained[obs_idx[i]] = pds_obs[i];
  }

  // Posterior predictive draws for model checking. The expected value is
  // constrained below RCN, but noisy realised draws on the log scale can exceed
  // RCN. Use pds_pred_raw for posterior predictive checks and pds_pred_capped
  // only where an operational product must obey PDS <= RCN at every draw.
  for (t in 1:N) {
    real draw = exp(normal_rng(mu_log1p_pds[t], sigma_log)) - 1.0;
    pds_pred_raw[t] = fmax(draw, 0.0);
    pds_pred_capped[t] = fmin(pds_pred_raw[t], rcn_obs[t]);
  }
}
