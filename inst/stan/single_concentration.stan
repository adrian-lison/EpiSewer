functions {
  #include functions/helper_functions.stan
  #include functions/link.stan
  #include functions/time_series.stan
  #include functions/approx_count_dist.stan
  #include functions/renewal.stan
  #include functions/dist_normal.stan
  #include functions/dist_lognormal.stan
  #include functions/dist_gamma.stan
  #include functions/dist_beta.stan
  #include functions/hurdle.stan
  #include functions/pcr_noise.stan
}
data {
  //# Measurements ----
  int<lower=0> n_measured; // number of replicates
  vector<lower=0>[n_measured] measured_concentrations; // measured concentrations
  vector<lower=0>[n_measured] n_averaged; // number of averaged technical replicates per measurement (is vector for vectorization)
  vector<lower=0>[n_measured] dPCR_total_partitions; // total number of partitions in dPCR
  int<lower=0, upper=3> obs_dist; // Parametric distribution for observation likelihood: 0 (default) for gamma, 1 for log-normal, 2 for truncated normal, 3 for normal

  //# Concentration ----
  array[2] real concentration_prior; // prior for true concentration

  //# Coefficient of variation (CV) of measurements ----
  int<lower=0, upper=2> cv_type; // 0 for constant, 1 for dPCR, 2 for constant_var
  array[2] real nu_upsilon_a_prior; // prior for pre-PCR CV
  int<lower=0, upper=1> total_partitions_observe; // 0 for not observed, 1 for observed
  array[cv_type == 1 ? 2 : 0] real nu_upsilon_b_mu_prior; // prior for parameter 2 of CV formula (number of partitions). Scaled by 1e-4 for numerical efficiency.
  array[cv_type == 1 && total_partitions_observe!=1 ? 2 : 0] real nu_upsilon_b_cv_prior; // prior on partition number variation (coefficient of variation)
  array[cv_type == 1 ? 2 : 0] real nu_upsilon_c_prior; // prior for parameter 3 of CV formula (partition size*(scaling factor, i.e. exp_conc_assay/exp_conc_ww)). Scaled by 1e+5 for numerical efficiency.
  array[cv_type == 1 ? 1 : 0] int <lower=0, upper=1> cv_pre_type; // 0 for gamma, 1 for log-normal
  array[cv_type == 1 ? 1 : 0] int <lower=0, upper=1> cv_pre_approx_taylor; // 0 for no Taylor expansion approximation, 1 for Taylor expansion approximation

  //# Limit of detection ----
  // LOD_model = 0: no LOD
  // LOD_model = 1: assumed LOD, LOD_scale provided
  // LOD_model = 2: estimated LOD based on dPCR model, needs dPCR parameters
  int<lower=0, upper=2> LOD_model;
  array[(LOD_model == 1) ? 1 : 0] real<lower=0> LOD_scale;
  real<lower=0, upper=1> LOD_drop_prob; // probability threshold for non-detection below which log likelihood contributions of observed concentrations are dropped from LOD model
}
transformed data {

  // number of averaged technical replicates per date
  real n_averaged_median = quantile(n_averaged, 0.5);

  // number of total partitions per measurement per date
  real total_partitions_median = quantile(dPCR_total_partitions, 0.5);

  // Upper relevant bound for LOD model
  real conc_drop_prob;
  if (LOD_model == 0) {
    conc_drop_prob = positive_infinity();
  } else {
    real LOD_expected_scale;
    if (LOD_model == 1) {
      LOD_expected_scale = LOD_scale[1];
    } else if (LOD_model == 2) {
      LOD_expected_scale = (
      (total_partitions_observe ? total_partitions_median : (nu_upsilon_b_mu_prior[1] * 1e4)) *
      nu_upsilon_c_prior[1] * 1e-5 *
      n_averaged_median
      );
    }
    conc_drop_prob = -log(LOD_drop_prob)/LOD_expected_scale; // concentrations above this value are irrelevant for LOD model (probability of non-detection is virtually zero)
  }

  int n_zero = num_zeros(measured_concentrations);
  int n_dropLOD = num_zeros(fmax(0, conc_drop_prob - measured_concentrations));
  array[n_zero] int<lower=0> i_zero;
  array[n_measured - n_zero] int<lower=0> i_nonzero;
  array[n_measured - n_dropLOD] int<lower=0> i_LOD;
  array[n_measured - n_zero - n_dropLOD] int<lower=0> i_nonzero_LOD;
  int i_z = 0;
  int i_nz = 0;
  int i_lod = 0;
  int i_nzs = 0;
  for (n in 1:n_measured) {
    if (measured_concentrations[n] == 0) {
      i_z += 1;
      i_zero[i_z] = n;
    } else {
      i_nz += 1;
      i_nonzero[i_nz] = n;
      if (measured_concentrations[n] < conc_drop_prob) {
        i_nzs += 1;
        i_nonzero_LOD[i_nzs] = n;
      }
    }
    if (measured_concentrations[n] < conc_drop_prob) {
        i_lod += 1;
        i_LOD[i_lod] = n;
    }
  }
}
parameters {
  real<lower=0> true_concentration;
  // Coefficient of variation of likelihood for measurements
  real<lower=0> nu_upsilon_a;
  array[(cv_type == 1) && (nu_upsilon_b_mu_prior[2] > 0) ? 1 : 0] real<lower=0> nu_upsilon_b_mu;
  array[(cv_type == 1) && total_partitions_observe!=1 && (nu_upsilon_b_cv_prior[2] > 0) ? 1 : 0] real<lower=0> nu_upsilon_b_cv;
  vector[(cv_type == 1) && total_partitions_observe!=1 ? n_measured : 0] nu_upsilon_b_noise_raw;
  array[(cv_type == 1) && nu_upsilon_c_prior[2] > 0 ? 1 : 0] real<lower=0> nu_upsilon_c;
}
transformed parameters {
  vector<lower=0>[(cv_type == 1) && total_partitions_observe!=1 ? n_measured : 0] nu_upsilon_b; // total partitions per measurement
  array[LOD_model > 0 ? 1 : 0] vector<lower=0>[n_measured] LOD_hurdle_scale;
  vector[n_measured] concentration = rep_vector(true_concentration, n_measured);
  vector[n_measured] p_zero_log = rep_vector(negative_infinity(), n_measured);
  vector[n_measured] p_zero = rep_vector(0, n_measured);

  if ((cv_type == 1) && total_partitions_observe!=1) {
    nu_upsilon_b = total_partitions_noncentered(
      param_or_fixed(nu_upsilon_b_mu, nu_upsilon_b_mu_prior),
      param_or_fixed(nu_upsilon_b_cv, nu_upsilon_b_cv_prior),
      nu_upsilon_b_noise_raw
    );
  }

 if (LOD_model > 0) {
    if (LOD_model == 1) {
      LOD_hurdle_scale[1] = rep_vector(LOD_scale[1], n_measured);
    } else if (LOD_model == 2) {
      LOD_hurdle_scale[1] = (
      n_averaged .*
      (total_partitions_observe ? dPCR_total_partitions : nu_upsilon_b * 1e4) *
       param_or_fixed(nu_upsilon_c, nu_upsilon_c_prior) * 1e-5
      );
    }
    p_zero_log[i_LOD] = log_hurdle_exponential(
        concentration[i_LOD],
        LOD_hurdle_scale[1][i_LOD], // LOD scale (c * m * n)
        cv_type == 1 ? nu_upsilon_a : 0, // nu_pre (pre-PCR CV)
        cv_type == 1 ? cv_pre_type[1] : 0 // Type of pre-PCR CV
        );
    p_zero = exp(p_zero_log);
  }

    // CV of each observation as a function of concentration
  vector[n_measured] cv;
  if (cv_type == 0) { // constant cv
    cv = rep_vector(nu_upsilon_a, n_measured);
  } else if (cv_type == 1) { // dPCR
    cv = cv_dPCR_pre(
      concentration, // lambda (concentration)
      nu_upsilon_a, // nu_pre (pre-PCR CV)
      (total_partitions_observe ? dPCR_total_partitions : nu_upsilon_b * 1e4), // m (number of partitions)
      param_or_fixed(nu_upsilon_c, nu_upsilon_c_prior) * 1e-5, // c (conversion factor)
      n_averaged, // n (number of averaged replicates)
      cv_pre_type[1], // Type of pre-PCR CV
      cv_pre_approx_taylor[1] // Should taylor approximation be used?
      );
  } else if (cv_type == 2) { // constant variance
    cv = (
      nu_upsilon_a * mean(measured_concentrations) /
      concentration
      );
  }

  vector[n_measured - n_zero] mean_conditional = concentration[i_nonzero] ./ (1-p_zero[i_nonzero]);
  vector[n_measured - n_zero] cv_conditional = sqrt(cv[i_nonzero]^2 .* (1-p_zero[i_nonzero]) - p_zero[i_nonzero]);
}
model {
  // Priors

  // True concentration
  true_concentration ~ normal(concentration_prior[1], concentration_prior[2]) T[0, ]; // truncated normal

  // Prior on cv of likelihood for measurements
  nu_upsilon_a ~ normal(nu_upsilon_a_prior[1], nu_upsilon_a_prior[2]) T[0, ]; // truncated normal
  if (cv_type == 1) {
    target += normal_prior_lpdf(nu_upsilon_b_mu | nu_upsilon_b_mu_prior, 0); // truncated normal
    target += normal_prior_lpdf(nu_upsilon_c | nu_upsilon_c_prior, 0); // truncated normal
    if (total_partitions_observe != 1) {
      target += normal_prior_lpdf(nu_upsilon_b_cv | nu_upsilon_b_cv_prior, 0); // truncated normal
      nu_upsilon_b_noise_raw ~ std_normal();
    }
  }

  // Likelihood
  {
    // limit of detection
    if (LOD_model > 0) {
      target += sum(p_zero_log[i_zero]); // below-LOD probabilities
      target += sum(log1m_exp(p_zero_log[i_nonzero_LOD])); // above-LOD probabilities
    }

    // measurements
    if (obs_dist == 0) {
      target += gamma3_lpdf(
        measured_concentrations[i_nonzero] |
        mean_conditional, // expectation
        cv_conditional // coefficient of variation
      );
    } else if (obs_dist == 1) {
      target += lognormal5_lpdf(
        measured_concentrations[i_nonzero] |
        mean_conditional, // log expectation
        cv_conditional // coefficient of variation
      );
    } else if (obs_dist == 2) {
      target += normal2_lpdf(
        measured_concentrations[i_nonzero] |
        mean_conditional, // expectation
        cv_conditional, // coefficient of variation,
        0 // truncate at zero
      );
    } else if (obs_dist == 3) {
      target += normal2_lpdf(
        measured_concentrations |
        concentration ./ (1-p_zero), // expectation
        sqrt(cv^2 .* (1-p_zero) - p_zero) // coefficient of variation
      );
    } else {
      reject("Distribution not supported.");
    }
  }
}
generated quantities {
  real predicted_concentration;
  real meas_conc;
  real above_LOD;

  above_LOD = bernoulli_rng(1-p_zero[1]);

  if (obs_dist == 0) {
    meas_conc = gamma3_rng(mean_conditional, cv_conditional)[1];
  } else if (obs_dist == 1) {
    meas_conc = lognormal5_rng(mean_conditional, cv_conditional)[1];
  } else if (obs_dist == 2) {
    meas_conc = normal2_rng(mean_conditional, cv_conditional, 0)[1]; // truncated at zero
  } else if (obs_dist == 3) {
    meas_conc = normal2_rng(mean_conditional, cv_conditional)[1];
  } else {
    reject("Distribution not supported.");
  }
  predicted_concentration = above_LOD * meas_conc;
}
