functions {
  #include functions/helper_functions.stan
  #include functions/link.stan
  #include functions/time_series.stan
  #include functions/approx_count_dist.stan
  #include functions/renewal.stan
  #include functions/normal2.stan
  #include functions/lognormal2.stan
  #include functions/gamma2.stan
  #include functions/hurdle.stan
  #include functions/pcr_noise.stan
}
data {
  int<lower=0> T; // number of total days in measured time span
  int<lower=0> n_samples; // number of samples from different dates
  int<lower=0> n_measured; // number of different measurements (including replicates)
  array[n_samples] int<lower=1, upper=T> sample_to_date; // index mapping samples to dates
  array[n_measured] int<lower=1, upper=n_samples> measure_to_sample; // index mapping measurements to samples
  vector<lower=0>[n_measured] measured_concentrations; // measured concentrations
  vector<lower=0>[n_measured] n_averaged; // number of averaged technical replicates per measurement (is vector for vectorization)
  int<lower=1> w; // composite window: how many past days the samples cover,
  // e.g. 1 for individual day samples, 7 for weekly composite samples, ...

  real<lower=0> load_mean; // mean load shed per person
  int<lower=0, upper=1> load_vari; // model individual-level variation in shedding loads
  array[load_vari ? 2 : 0] real nu_zeta_prior; // Prior on individual-level load coefficient of variation

  vector<lower=0>[T] flow; // flow rate for normalization of measurements

  // Sample date effects model
  int<lower=0> K; // number of sample date predictors
  matrix[T, K] X; // sample date predictor design matrix
  array[K > 0 ? 2 : 0] real eta_prior; // prior for sample date effects

  // Pre-replicate noise
  int<lower=0, upper=1> pr_noise; // Pre-replicate noise: Model variation before replication stage?
  array[pr_noise ? 2 : 0] real nu_psi_prior; // Prior on coefficient of variation fpr pre-replicate noise

  // Coefficient of variation (CV) of lognormal likelihood of measurements
  int<lower=0, upper=2> cv_type; // 0 for constant, 1 for ddPCR, 2 for constant_var
  array[2] real nu_upsilon_a_prior; // prior for pre-PCR CV
  real nu_upsilon_b_fixed;
  array[cv_type == 1 && nu_upsilon_b_fixed < 0 ? 2 : 0] real nu_upsilon_b_prior; // prior for parameter 2 of CV formula (number of droplets). Scaled by 1e-4 for numerical efficiency.
  real nu_upsilon_c_fixed;
  array[cv_type == 1 && nu_upsilon_c_fixed < 0 ? 2 : 0] real nu_upsilon_c_prior; // prior for parameter 3 of CV formula (droplet size*(scaling factor, i.e. exp_conc_assay/exp_conc_ww)). Scaled by 1e+5 for numerical efficiency.
  array[cv_type == 1 ? 1 : 0] int <lower=0, upper=1> cv_pre_type; // 0 for gamma, 1 for log-normal
  array[cv_type == 1 ? 1 : 0] int <lower=0, upper=1> cv_pre_approx_taylor; // 0 for no Taylor expansion approximation, 1 for Taylor expansion approximation

  // Limit of detection
  // LOD_model = 0: no LOD
  // LOD_model = 1: assumed LOD, LOD_scale provided
  // LOD_model = 2: estimated LOD based on ddPCR model, needs ddPCR parameters
  int<lower=0, upper=2> LOD_model;
  array[(LOD_model == 1) ? 1 : 0] real<lower=0> LOD_scale;

  // Residence time distribution
  int<lower=0> D; // last day of shedding
  vector[D + 1] residence_dist; // residence time distribution
  // --> probability for residence of zero days (same day arrival at sampling site) comes first

  // Shedding load distribution
  int<lower=0> S; // last day of shedding
  vector[S + 1] shedding_dist; // shedding load distribution
  // --> probability for shedding today comes first

  // Incubation period distribution
  int<lower=0> L; // maximum delay
  vector[L + 1] incubation_dist; // incubation period distribution
  // --> probability for a delay of zero comes first

  // Generation interval distribution
  int<lower=1> G; // maximum generation interval
  vector<lower=0>[G] generation_dist; // generation interval distribution
  // --> probability for a delay of one comes first (zero excluded)

  // Hyperpriors for seeding infections
  int<lower=0, upper=1> seeding_model; // 0 for fixed, 1 for random walk
  array[2] real iota_log_seed_intercept_prior;
  array[seeding_model == 1 ? 2 : 0] real iota_log_seed_sd_prior;

  // Stochastic (=1) or deterministic (=0) renewal process?
  int<lower=0, upper=1> I_sample;
  int<lower=0, upper=I_sample> I_overdispersion; // whether to model overdispersion via a negative binomial
  real I_xi_fixed; // fixed overdispersion parameter
  array[I_overdispersion && I_xi_fixed < 0 ? 2 : 0] real I_xi_prior; // prior on the overdispersion parameter

  // Hyperpriors for smoothing
  array[2] real R_level_start_prior;
  array[2] real R_trend_start_prior;
  array[2] real R_sd_prior;

  // Exponential smoothing priors / configuration
  int<lower=0> ets_diff; // order of differencing
  int<lower=0, upper=1> ets_noncentered; // use non-centered parameterization?
  real ets_alpha_fixed; // fixed value used if non-negative
  array[ets_alpha_fixed < 0 ? 2 : 0] real<lower=0> ets_alpha_prior;
  real ets_beta_fixed; // fixed value used if non-negative
  array[ets_beta_fixed < 0 ? 2 : 0] real<lower=0> ets_beta_prior;
  real ets_phi_fixed; // fixed value used if non-negative
  array[ets_phi_fixed < 0 ? 2 : 0] real<lower=0> ets_phi_prior;

  // Link function and corresponding hyperparameters
  // first element: 0 = inv_softplus, 1 = scaled_logit
  // other elements: hyperparameters for the respective link function
  array[4] real R_link;

  // Parametric distribution for observation likelihood
  // 1 (default) = truncated normal
  // 2 = lognormal
  int<lower=1, upper=2> obs_dist;
}
transformed data {
  vector[G] gi_rev = reverse(generation_dist);
  vector[L + 1] inc_rev = reverse(incubation_dist);
  vector[S + 1] shed_rev_log = log(reverse(shedding_dist));
  vector[D + 1] residence_rev_log = log(reverse(residence_dist));
  vector[T] flow_log = log(flow);

  // Upper relevant bound for LOD model
  real LOD_expected_scale;
  if (LOD_model == 1) {
    LOD_expected_scale = LOD_scale[1];
  } else if (LOD_model == 2) {
    LOD_expected_scale = (
    (nu_upsilon_b_fixed < 0 ? nu_upsilon_b_prior[1] : nu_upsilon_b_fixed) * 1e4 *
    (nu_upsilon_c_fixed < 0 ? nu_upsilon_c_prior[1] : nu_upsilon_c_fixed) * 1e-5
    );
  }
  real LOD_irrelevant = -log(1e-4)/LOD_expected_scale; // concentrations above this value are irrelevant for LOD model (probability of non-detection is virtually zero)

  int n_zero = num_zeros(measured_concentrations);
  int notsmall = num_zeros(fmax(0, LOD_irrelevant - measured_concentrations));
  array[n_zero] int<lower=0> i_zero;
  array[n_measured - n_zero] int<lower=0> i_nonzero;
  array[n_measured - n_zero - notsmall] int<lower=0> i_nonzero_small;
  int i_z = 0;
  int i_nz = 0;
  int i_nzs = 0;
  for (n in 1:n_measured) {
    if (measured_concentrations[n] == 0) {
      i_z += 1;
      i_zero[i_z] = n;
    } else {
      i_nz += 1;
      i_nonzero[i_nz] = n;
      if (measured_concentrations[n] < LOD_irrelevant) {
        i_nzs += 1;
        i_nonzero_small[i_nzs] = n;
      }
    }
  }

  // number of averaged technical replicates per date
  real n_averaged_median = quantile(n_averaged, 0.5);
  vector[T] n_averaged_all = rep_vector(n_averaged_median, T);
  for (i in 1:n_measured) {
    // note that if several measurements per sample exist,
    // the number of replicates of the last one will be used for that date
    n_averaged_all[sample_to_date[measure_to_sample[i]]] = n_averaged[i];
  }

  if (R_link[1] < 0 || R_link[1] > 1) {
    reject("Link function must be one of inv_softplus (0) or scaled_logit (1)");
  }
}
parameters {
  // log(R) time series prior
  real R_level_start; // starting value of the level
  real R_trend_start; // starting value of the trend
  real<lower=0> R_sd; // standard deviation of additive errors
  vector<multiplier=(ets_noncentered ? R_sd : 1)>[L + S + D + T - G - 1] R_noise; // additive errors

  // exponential smoothing / innovations state space process for log(R)
  array[ets_alpha_fixed < 0 ? 1 : 0] real<lower=0, upper=1> ets_alpha; // smoothing parameter for the level
  array[ets_beta_fixed < 0 ? 1 : 0] real<lower=0, upper=1> ets_beta; // smoothing parameter for the trend
  array[ets_phi_fixed < 0 ? 1 : 0] real<lower=0, upper=1> ets_phi; // dampening parameter of the trend

  // seeding
  real iota_log_seed_intercept;
  array[seeding_model == 1 ? 1 : 0] real<lower=0> iota_log_seed_sd;
  vector<multiplier=(seeding_model == 1 ? iota_log_seed_sd[1] : 1)>[seeding_model == 1 ? G - 1 : 0] iota_log_ar_noise;

  // realized infections
  array[I_overdispersion && I_xi_fixed < 0 ? 1 : 0] real<lower=0> I_xi; // positive to ensure identifiability
  vector<lower=0>[I_sample ? L + S + D + T : 0] I; // realized number of infections

  // individual-level shedding load variation
  array[load_vari ? 1 : 0] real<lower=0> nu_zeta; // coefficient of variation of individual-level load
  vector[load_vari ? S + D + T : 0] zeta_raw; // realized shedding load (non-centered)

  // sample date effects
  vector[K] eta;

  // Coefficient of variation of lognormal likelihood for measurements
  array[pr_noise ? 1 : 0] real<lower=0> nu_psi; // pre-replicaton coefficient of variation
  vector<multiplier = (pr_noise ? nu_psi[1] : 1)>[pr_noise ? n_samples : 0] psi; // realized concentration before replication stage
  real<lower=0> nu_upsilon_a;
  array[(cv_type == 1) && (nu_upsilon_b_fixed < 0) ? 1 : 0] real<lower=0> nu_upsilon_b;
  array[(cv_type == 1) && (nu_upsilon_c_fixed < 0) ? 1 : 0] real<lower=0> nu_upsilon_c;
}
transformed parameters {
  vector[L + S + D + T - G] R; // effective reproduction number
  vector[L + S + D + T] iota; // expected number of infections
  vector[S + D + T] lambda; // expected number of shedding onsets
  vector<lower = 0>[load_vari ? S + D + T : 0] zeta; // realized shedding load
  vector[T] pi_log; // log expected daily loads
  vector[T] kappa_log; // log expected daily concentrations
  vector[n_samples] rho_log; // log expected concentrations in (composite) samples
  array[LOD_model > 0 ? 1 : 0] vector<lower=0>[n_measured] LOD_hurdle_scale;

  // Innovations state space process implementing exponential smoothing
  R = apply_link(holt_damped_process(
    [R_level_start, R_trend_start]',
    ets_alpha_fixed < 0 ? ets_alpha[1] : ets_alpha_fixed,
    ets_beta_fixed < 0 ? ets_beta[1] : ets_beta_fixed,
    ets_phi_fixed < 0 ? ets_phi[1] : ets_phi_fixed,
    R_noise,
    ets_diff
  ), R_link);

  // seeding
  if (seeding_model == 0) {
    iota[1 : G] = exp(rep_vector(iota_log_seed_intercept, G));
  } else if (seeding_model == 1) {
    iota[1 : G] = exp(random_walk([iota_log_seed_intercept]', iota_log_ar_noise, 0));
  }
  // renewal process
  if (I_sample) {
    iota[(G + 1) : (L + S + D + T)] = renewal_process_stochastic(
      (L + S + D + T - G), R, G, gi_rev, I);
  } else {
    iota[(G + 1) : (L + S + D + T)] = renewal_process_deterministic(
      (L + S + D + T - G), R, G, gi_rev, iota);
  }

  // convolution from infections to shedding onsets (expected)
  lambda = convolve(inc_rev, I_sample ? I : iota)[(L + 1) : (L + S + D + T)];

  // calculation of total loads shed each day (expected)
  vector[D + T] omega_log;
  if (load_vari) {
    zeta = softplus(gamma_sum_approx(nu_zeta[1], lambda, zeta_raw), 10); // softplus as soft >0 constraint
    omega_log = log_convolve(
        shed_rev_log, // shedding load distribution
        log(load_mean) + log(zeta) // total load shed
        )[(S + 1) : (S + D + T)];
  } else {
    omega_log = log_convolve(
        shed_rev_log, // shedding load distribution
        log(load_mean) + log(lambda) // total load shed
        )[(S + 1) : (S + D + T)];
  }

  // calculation of total loads at sampling site (expected)
  if (D>0) {
    pi_log = log_convolve(
      residence_rev_log, // residence time distribution
      omega_log
      )[(D + 1) : (D + T)];
  } else {
    pi_log = omega_log;
  }

  // calculation of concentrations at measurement site by day (expected)
  // --> adjusted for flow and for date of sample effects
  if (K > 0) {
    kappa_log = pi_log - flow_log + X * eta;
  } else {
    kappa_log = pi_log - flow_log;
  }

  // concentrations in (composite) samples
  if (w > 1) { // multi-day composite samples
    for (i in 1 : n_samples) {
        rho_log[i] = log_sum_exp(
          kappa_log[(sample_to_date[i] - w + 1) : sample_to_date[i]]
          ) - log(w);
    }
  } else { // individual day samples
     for (i in 1 : n_samples) {
        rho_log[i] = kappa_log[sample_to_date[i]];
     }
  }

  // scale for LOD hurdle model
  if (LOD_model == 1) {
    LOD_hurdle_scale[1] = rep_vector(LOD_scale[1], n_measured);
  } else if (LOD_model == 2) {
    LOD_hurdle_scale[1] = (
    (nu_upsilon_b_fixed < 0 ? nu_upsilon_b[1] : nu_upsilon_b_fixed) * 1e4 *
    (nu_upsilon_c_fixed < 0 ? nu_upsilon_c[1] : nu_upsilon_c_fixed) * 1e-5 *
    n_averaged
    );
  }
}
model {
  // Priors

  // Innovations state space model
  ets_coefficient_priors_lp(ets_alpha, ets_alpha_fixed, ets_alpha_prior,
                            ets_beta, ets_beta_fixed, ets_beta_prior,
                            ets_phi, ets_phi_fixed, ets_phi_prior);
  R_level_start ~ normal(R_level_start_prior[1], R_level_start_prior[2]); // starting prior for level
  R_trend_start ~ normal(R_trend_start_prior[1], R_trend_start_prior[2]); // starting prior for trend
  R_sd ~ normal(R_sd_prior[1], R_sd_prior[2]) T[0, ]; // truncated normal
  R_noise ~ normal(0, R_sd); // Gaussian noise

  // Seeding
  iota_log_seed_intercept ~ normal(iota_log_seed_intercept_prior[1], iota_log_seed_intercept_prior[2]);
  if (seeding_model == 1) {
    iota_log_seed_sd[1] ~ normal(iota_log_seed_sd_prior[1], iota_log_seed_sd_prior[2]) T[0, ]; // truncated normal
    iota_log_ar_noise ~ normal(0, iota_log_seed_sd[1]); // Gaussian noise
  }

  // Sampling of infections
  if (I_sample) {
    if (I_overdispersion) {
      if (I_xi_fixed < 0) {
        I_xi[1] ~ normal(I_xi_prior[1], I_xi_prior[2]) T[0, ]; // truncated normal
        I[1 : (L + S + D + T)] ~ normal(iota, sqrt(iota .* (1 + iota * (I_xi[1] ^ 2)))); // approximates negative binomial
      } else {
        I[1 : (L + S + D + T)] ~ normal(iota, sqrt(iota .* (1 + iota * (I_xi_fixed ^ 2)))); // approximates negative binomial
      }
    } else {
      I[1 : (L + S + D + T)] ~ normal(iota, sqrt(iota)); // approximates Poisson
    }
  }

  // Prior on individual-level shedding load variation
  if (load_vari) {
    nu_zeta[1] ~ normal(nu_zeta_prior[1], nu_zeta_prior[2]) T[0, ]; // truncated normal
    zeta_raw ~ std_normal();
  }

  // Prior on sample date effects
  if (K > 0) {
    eta ~ normal(eta_prior[1], eta_prior[2]);
  }

  // Prior on cv of lognormal likelihood for measurements
  nu_upsilon_a ~ normal(nu_upsilon_a_prior[1], nu_upsilon_a_prior[2]) T[0, ]; // truncated normal
  if (cv_type == 1) {
    if (nu_upsilon_b_fixed < 0) {
      nu_upsilon_b ~ normal(nu_upsilon_b_prior[1], nu_upsilon_b_prior[2]) T[0, ]; // truncated normal
    }
    if (nu_upsilon_c_fixed < 0) {
      nu_upsilon_c ~ normal(nu_upsilon_c_prior[1], nu_upsilon_c_prior[2]) T[0, ]; // truncated normal
    }
  }

  // Likelihood
  {
    // accounting for pre-replicate noise
    vector[n_measured] concentrations;
    if (pr_noise) {
      // Prior on cv of pre-replication concentrations
      nu_psi[1] ~ normal(nu_psi_prior[1], nu_psi_prior[2]) T[0, ]; // truncated normal
      // lognormal distribution modeled as normal on the log scale
      // with mean = rho (mean_log = rho_log) and cv = nu_psi
      target += lognormal_log_lpdf(psi | rho_log, nu_psi[1]);
      concentrations = psi[measure_to_sample];
    } else {
      concentrations = rho_log[measure_to_sample];
    }

    // concentration on unit scale
    vector[n_measured - n_zero] concentrations_unit = exp(concentrations);

    // limit of detection
    if (LOD_model > 0) {
      // below-LOD probabilities for zero measurements
      target += sum(log_hurdle_exponential(
        concentrations_unit[i_zero],
        LOD_hurdle_scale[1][i_zero], // LOD scale (c * m * n)
        cv_type == 1 ? nu_upsilon_a : 0, // nu_pre (pre-PCR CV)
        cv_type == 1 ? cv_pre_type[1] : 0 // Type of pre-PCR CV
        ));
      // above-LOD probabilities for non-zero measurements
      target += sum(log1m_exp(log_hurdle_exponential(
        concentrations_unit[i_nonzero_small],
        LOD_hurdle_scale[1][i_nonzero_small], // LOD scale (c * m * n)
        cv_type == 1 ? nu_upsilon_a : 0, // nu_pre (pre-PCR CV)
        cv_type == 1 ? cv_pre_type[1] : 0 // Type of pre-PCR CV
      )));
    }

    // CV of each observation as a function of concentration
    vector[n_measured - n_zero] cv;
    if (cv_type == 0) { // constant cv
      cv = rep_vector(nu_upsilon_a, n_measured - n_zero);
    } else if (cv_type == 1) { // ddPCR
      cv = cv_ddPCR_pre(
        concentrations_unit[i_nonzero],
        nu_upsilon_a,
        (nu_upsilon_b_fixed < 0 ? nu_upsilon_b[1] : nu_upsilon_b_fixed) * 1e4,
        (nu_upsilon_c_fixed < 0 ? nu_upsilon_c[1] : nu_upsilon_c_fixed) * 1e-5,
        n_averaged[i_nonzero],
        cv_pre_type[1],
        cv_pre_approx_taylor[1]
        );
    } else if (cv_type == 2) { // constant variance
      cv = (
        nu_upsilon_a * mean(measured_concentrations[i_nonzero]) /
        concentrations_unit[i_nonzero]
        );
    }

    // measurements
    if (obs_dist == 1) {
      target += normal2_lpdf(
        measured_concentrations[i_nonzero] |
        concentrations_unit[i_nonzero], // expectation
        cv, // coefficient of variation
        0 // truncate at zero
      );
    }
    else if (obs_dist == 2) {
      target += lognormal4_lpdf(
        measured_concentrations[i_nonzero] |
        concentrations[i_nonzero], // log expectation
        cv // coefficient of variation
      );
    } else {
      reject("Distribution not supported.");
    }
  }
}
generated quantities {
  // predicted measurements
  // note that we here assume the same measurement variance as from composite samples,
  // which may be smaller than that of hypothetical daily measurements
  vector[T] predicted_concentration;
  {
    vector[T] pre_repl;
    vector[T] exp_pre_repl;
    vector[T] cv;
    vector[T] above_LOD; // will be a vector of 0s and 1s
    if (pr_noise) {
      pre_repl = to_vector(lognormal_log_rng(kappa_log, nu_psi[1]));
    } else {
      pre_repl = kappa_log;
    }

    exp_pre_repl = exp(pre_repl);

    if (LOD_model > 0) {
     vector[T] LOD_hurdle_scale_all;
      // scale for LOD hurdle model
      if (LOD_model == 1) {
        LOD_hurdle_scale_all = rep_vector(LOD_scale[1], T);
      } else if (LOD_model == 2) {
        LOD_hurdle_scale_all = (
        (nu_upsilon_b_fixed < 0 ? nu_upsilon_b[1] : nu_upsilon_b_fixed) * 1e4 *
        (nu_upsilon_c_fixed < 0 ? nu_upsilon_c[1] : nu_upsilon_c_fixed) * 1e-5 *
        n_averaged_all
        );
      }
      above_LOD = to_vector(bernoulli_rng(
        1-exp(log_hurdle_exponential(
          exp_pre_repl, // lambda (concentration)
          LOD_hurdle_scale_all,
          cv_type == 1 ? nu_upsilon_a : 0, // nu_pre (pre-PCR CV)
          cv_type == 1 ? cv_pre_type[1] : 0 // Type of pre-PCR CV
          ))
          ));
    } else {
      above_LOD = rep_vector(1, T);
    }

    if (cv_type == 0) {
      cv = rep_vector(nu_upsilon_a, T);
    } else if (cv_type == 1) {
      cv = cv_ddPCR_pre(
        exp_pre_repl,
        nu_upsilon_a,
        (nu_upsilon_b_fixed < 0 ? nu_upsilon_b[1] : nu_upsilon_b_fixed) * 1e4,
        (nu_upsilon_c_fixed < 0 ? nu_upsilon_c[1] : nu_upsilon_c_fixed) * 1e-5,
        n_averaged_all,
        cv_pre_type[1],
        cv_pre_approx_taylor[1]
        );
    } else if (cv_type == 2) {
      cv = (
        nu_upsilon_a * mean(measured_concentrations[i_nonzero]) /
        exp_pre_repl
        );
    }

    vector[T] meas_conc;
    if (obs_dist == 1) {
      meas_conc = normal2_lb_rng(exp_pre_repl, cv, 0);
    }
    else if (obs_dist == 2) {
      meas_conc = lognormal4_rng(pre_repl, cv);
    } else {
      reject("Distribution not supported.");
    }
    predicted_concentration = above_LOD .* meas_conc;
  }
}
