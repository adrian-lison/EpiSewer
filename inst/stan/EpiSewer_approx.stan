functions {
  #include functions/helper_functions.stan
  #include functions/link.stan
  #include functions/time_series.stan
  #include functions/approx_count_dist.stan
  #include functions/renewal.stan
  #include functions/dist_normal.stan
  #include functions/dist_lognormal.stan
  #include functions/dist_gamma.stan
  #include functions/dist_loggamma.stan
  #include functions/dist_gamma_sum.stan
  #include functions/dist_beta.stan
  #include functions/dist_gev.stan
  #include functions/hurdle.stan
  #include functions/pcr_noise.stan
}
data {
  // Measurements ----
  int<lower=0> T; // number of total days in measured time span
  int<lower=0> h; // number of days to forecast
  int<lower=0> n_samples; // number of samples from different dates
  int<lower=0> n_measured; // number of different measurements (including replicates)
  array[n_samples] int<lower=1, upper=T> sample_to_date; // index mapping samples to dates
  array[n_measured] int<lower=1, upper=n_samples> measure_to_sample; // index mapping measurements to samples
  vector<lower=0>[n_measured] measured_concentrations; // measured concentrations
  vector<lower=0>[n_measured] n_averaged; // number of averaged technical replicates per measurement (is vector for vectorization)
  vector<lower=0>[n_measured] dPCR_total_partitions; // total number of partitions in dPCR
  int<lower=1> w; // composite window: how many past days the samples cover, e.g. 1 for individual day samples, 7 for weekly composite samples, ...
  int<lower=0, upper=3> obs_dist; // Parametric distribution for observation likelihood: 0 (default) for gamma, 1 for log-normal, 2 for truncated normal, 3 for normal

  // Flow ----
  vector<lower=0>[T+h] flow; // flow rate for normalization of measurements

  // Sample date effects model ----
  int<lower=0> K; // number of sample date predictors
  matrix[T+h, K] X; // sample date predictor design matrix
  array[K > 0 ? 2 : 0] real eta_prior; // prior for sample date effects

  // Pre-replicate noise ----
  int<lower=0, upper=1> pr_noise; // Pre-replicate noise: Model variation before replication stage?
  array[pr_noise ? 2 : 0] real nu_psi_prior; // Prior on coefficient of variation for pre-replicate noise

  // Coefficient of variation (CV) of measurements ----
  int<lower=0, upper=2> cv_type; // 0 for constant, 1 for dPCR, 2 for constant_var
  array[2] real nu_upsilon_a_prior; // prior for pre-PCR CV
  int<lower=0, upper=1> total_partitions_observe; // 0 for not observed, 1 for observed
  array[cv_type == 1 ? 2 : 0] real nu_upsilon_b_mu_prior; // prior for parameter 2 of CV formula (number of partitions). Scaled by 1e-4 for numerical efficiency.
  array[cv_type == 1 && total_partitions_observe!=1 ? 2 : 0] real nu_upsilon_b_cv_prior; // prior on partition number variation (coefficient of variation)
  array[cv_type == 1 ? 2 : 0] real nu_upsilon_c_prior; // prior for parameter 3 of CV formula (partition size*(scaling factor, i.e. exp_conc_assay/exp_conc_ww)). Scaled by 1e+5 for numerical efficiency.
  array[cv_type == 1 ? 1 : 0] int <lower=0, upper=1> cv_pre_type; // 0 for gamma, 1 for log-normal
  array[cv_type == 1 ? 1 : 0] int <lower=0, upper=1> cv_pre_approx_taylor; // 0 for no Taylor expansion approximation, 1 for Taylor expansion approximation

  // Limit of detection ----
  // LOD_model = 0: no LOD
  // LOD_model = 1: assumed LOD, LOD_scale provided
  // LOD_model = 2: estimated LOD based on dPCR model, needs dPCR parameters
  int<lower=0, upper=2> LOD_model;
  array[(LOD_model == 1) ? 1 : 0] real<lower=0> LOD_scale;
  real<lower=0, upper=1> LOD_drop_prob; // probability threshold for non-detection below which log likelihood contributions of observed concentrations are dropped from LOD model

  // Residence time ----
  int<lower=0> D; // maximum residence time
  vector[D + 1] residence_dist; // residence time distribution
  // --> probability for residence of zero days (same day arrival at sampling site) comes first

  // Shedding ----
  real<lower=0> load_mean; // mean load shed per person
  int<lower=1> S; // maximum number of days with shedding
  vector[S + 1] shedding_dist; // shedding load distribution
  // --> probability for shedding today comes first
  int<lower=0, upper=1> load_vari; // model individual-level variation in shedding loads?
  array[load_vari ? 2 : 0] real nu_zeta_prior; // prior on coefficient of variation of individual-level load
  int<lower = 0> n_zeta_normal_approx;
  int<lower = 0> n_zeta_exact;
  array[load_vari ? n_zeta_normal_approx : 0] int<lower = 1, upper = S + D + T> zeta_normal_approx; // dates on which zeta can be approximated
  array[load_vari ? n_zeta_exact : 0] int<lower = 1, upper = S + D + T> zeta_exact; // dates on which zeta should not be approximated

  // Outliers ----
  int<lower=0, upper=1> outliers;
  array[outliers ? 3 : 0] real<lower=0> epsilon_prior;

  // Incubation period ----
  int<lower=1> L; // maximum incubation period
  vector[L + 1] incubation_dist; // incubation period distribution
  // --> probability for a delay of zero comes first

  // Generation interval ----
  int<lower=1> G; // maximum generation interval
  vector<lower=0>[G] generation_dist; // generation interval distribution
  // --> probability for a delay of one comes first (zero excluded)

  // Seeding of infections ----
  int<lower=0> se; // seeding extension (used when time series starts with many non-detects)
  int<lower=0, upper=1> seeding_model; // 0 for fixed, 1 for random walk seeding
  array[2] real iota_log_seed_intercept_prior;
  array[seeding_model == 1 ? 2 : 0] real iota_log_seed_sd_prior;
  row_vector[(G+se)] iota_log_seed_trend_reg;

  // Infection noise ----
  int<lower=0, upper=1> I_sample; // Stochastic (=1) or deterministic (=0) renewal process?
  int<lower=0, upper=I_sample> I_overdispersion; // whether to model overdispersion via a negative binomial
  array[I_overdispersion ? 2 : 0] real I_xi_prior; // prior on the overdispersion parameter

  // Basis spline (bs) configuration for smoothing infections
  // Sparse bs matrix: columns = bases (bs_n_basis), rows = time points (L+S+T-(G+se))
  int<lower=1> bs_n_basis; // number of B-splines
  vector[bs_n_basis - 1] bs_dists; // distances between knots
  int<lower=0> bs_n_w; // number of nonzero entries in bs matrix
  vector[bs_n_w] bs_w; // nonzero entries in bs matrix
  array[bs_n_w] int bs_v; // column indices of bs_w
  array[L + S + D + T + 1 - (G+se) + h] int bs_u; // row starting indices for bs_w plus padding
  array[2] real inf_ar_sd_prior; // sd hyperprior for random walk on log bs coeffs
  real<lower=0, upper=1> inf_smooth;
  real<lower=0, upper=1> inf_trend_smooth;
  real<lower=0, upper=1> inf_trend_dampen;

  // Reproduction number ----
  // Note that in this model, Rt is not estimated as a parameter but derived
  // directly from the estimated infection trajectory as a computed quantity
  int<lower=1> R_w; // R smoothing window (for compatibility with Cori et al.)
}
transformed data {
  vector[G] gi_rev = reverse(generation_dist);
  vector[L + 1] inc_rev = reverse(incubation_dist);
  vector[S + 1] shed_rev_log = log(reverse(shedding_dist));
  vector[D + 1] residence_rev_log = log(reverse(residence_dist));
  vector[T+h] flow_log = log(flow);

  // number of averaged technical replicates per date
  real n_averaged_median = quantile(n_averaged, 0.5);
  vector[T+h] n_averaged_all = rep_vector(n_averaged_median, T+h);
  for (i in 1:n_measured) {
    // note that if several measurements per sample exist,
    // the number of replicates of the last one will be used for that date
    n_averaged_all[sample_to_date[measure_to_sample[i]]] = n_averaged[i];
  }

  // number of total partitions per measurement per date
  real total_partitions_median = quantile(dPCR_total_partitions, 0.5);
  vector[T+h] total_partitions_all = rep_vector(total_partitions_median, T+h);
  for (i in 1:n_measured) {
    // note that if several measurements per sample exist,
    // the total partitions of the last one will be used for that date
    total_partitions_all[sample_to_date[measure_to_sample[i]]] = dPCR_total_partitions[i];
  }

  // Number of steps between spline coefficients
  int bs_n = to_int(sum(bs_dists));
  array[bs_n_basis] int bs_select = to_int(to_array_1d(
    cumulative_sum(append_row(1, bs_dists))
    ));

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
  array[n_measured] int<lower=0> i_all;
  array[n_zero] int<lower=0> i_zero;
  array[n_measured - n_zero] int<lower=0> i_nonzero;
  array[n_measured - n_dropLOD] int<lower=0> i_LOD;
  array[n_measured - n_zero - n_dropLOD] int<lower=0> i_nonzero_LOD;
  int i_z = 0;
  int i_nz = 0;
  int i_lod = 0;
  int i_nzs = 0;
  for (n in 1:n_measured) {
    i_all[n] = n;
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
  int n_dropzero = (obs_dist == 3 ? 0 : n_zero); // only include zeros for non-truncated normal
  array[n_measured - n_dropzero] int<lower=0> i_include;
  if (n_dropzero > 0) {
    i_include = i_nonzero;
  } else {
    i_include = i_all;
  }

  // half R smoothing window
  int R_w_half; // half of R smoothing window
  if ((R_w % 2)==0) {
    reject("Smoothing window size for R must be an odd number.");
  } else {
    if (R_w == 1) {
      R_w_half = 0;
    } else {
      R_w_half = (R_w - 1) %/% 2;
    }
  }
}
parameters {
  // log(iota) time series prior
  real<lower=0> inf_ar_sd; // sd for random walk on log bs coeffs
  vector[bs_n] bs_coeff_noise_raw; // additive errors (non-centered)

  // seeding
  real iota_log_seed_intercept;
  array[seeding_model == 1 ? 1 : 0] real<lower=0> iota_log_seed_sd;
  vector<multiplier=(seeding_model == 1 ? iota_log_seed_sd[1] : 1)>[seeding_model == 1 ? (G+se) - 1 : 0] iota_log_ar_noise;

  // realized infections
  array[I_overdispersion && (I_xi_prior[2] > 0) ? 1 : 0] real<lower=0> I_xi; // positive to ensure identifiability
  vector[I_sample ? L + S + D + T : 0] I_raw; // infection noise

  // individual-level shedding load variation
  array[load_vari && nu_zeta_prior[2] > 0 ? 1 : 0] real<lower=0> nu_zeta; // coefficient of variation of individual-level load
  vector[load_vari ? n_zeta_exact : 0] zeta_log_exact; // realized shedding load
  vector[load_vari ? n_zeta_normal_approx : 0] zeta_raw_approx; // realized shedding load (non-centered)

  // sample date effects
  vector[K] eta;

  // Coefficient of variation of likelihood for measurements
  array[pr_noise ? 1 : 0] real<lower=0> nu_psi; // pre-replicaton coefficient of variation
  vector<multiplier = (pr_noise ? nu_psi[1] : 1)>[pr_noise ? n_samples : 0] psi; // realized concentration before replication stage
  real<lower=0> nu_upsilon_a;
  array[(cv_type == 1) && (nu_upsilon_b_mu_prior[2] > 0) ? 1 : 0] real<lower=0> nu_upsilon_b_mu;
  array[(cv_type == 1) && total_partitions_observe!=1 && (nu_upsilon_b_cv_prior[2] > 0) ? 1 : 0] real<lower=0> nu_upsilon_b_cv;
  vector[(cv_type == 1) && total_partitions_observe!=1 ? n_measured : 0] nu_upsilon_b_noise_raw;
  array[(cv_type == 1) && nu_upsilon_c_prior[2] > 0 ? 1 : 0] real<lower=0> nu_upsilon_c;

  vector<lower=0>[outliers ? D + T : 0] epsilon; // external noise at low concentrations;
}
transformed parameters {
  vector[bs_n_basis] bs_coeff; // Basis spline coefficients
  vector[L + S + D + T] iota; // expected number of infections
  vector[h] iota_forecast; // spline-based forecast of expected infections
  vector<lower=0>[I_sample ? L + S + D + T : 0] I; // realized number of infections
  vector[L + S + D + T] I_noise_correction; // correction for approximate infection noise
  vector[load_vari ? S + D + T : 0] zeta_log; // realized shedding load
  vector[S + D + T] lambda; // expected number of shedding onsets
  vector[D + T] omega_log; // log expected daily loads in catchment
  vector[T] pi_log; // log expected daily loads at sampling site
  vector[T] kappa_log; // log expected daily concentrations
  vector[n_samples] rho_log; // log expected concentrations in (composite) samples
  vector<lower=0>[(cv_type == 1) && total_partitions_observe!=1 ? n_measured : 0] nu_upsilon_b; // total partitions per measurement
  array[LOD_model > 0 ? 1 : 0] vector<lower=0>[n_measured] LOD_hurdle_scale;

  // seeding
  {
    vector[(G+se)] iota_log_seed;
    if (seeding_model == 0) {
      iota_log_seed = rep_vector(iota_log_seed_intercept, (G+se));
    } else if (seeding_model == 1) {
      iota_log_seed = random_walk([iota_log_seed_intercept]', iota_log_ar_noise, 0);
    }
    iota[1 : (G+se)] = exp(iota_log_seed);

    // Spline smoothing of infections
    bs_coeff = holt_damped_process(
      [iota_log_seed[(G+se)], iota_log_seed_trend_reg * iota_log_seed]',
      inf_smooth,
      inf_trend_smooth,
      inf_trend_dampen,
      bs_coeff_noise_raw * inf_ar_sd, 0)[bs_select];

    vector[L + S + D + T - (G+se) + h] iota_all = apply_link(csr_matrix_times_vector(
      L + S + D + T - (G+se) + h, bs_n_basis, bs_w, bs_v, bs_u, bs_coeff
     ), rep_array(2, 1)); // log link
    iota[(G+se+1) : (L + S + D + T)] = iota_all[1 : (L + S + D + T - (G+se))];
    if (h > 0) {
      iota_forecast = iota_all[(L + S + D + T - (G+se) + 1) : (L + S + D + T - (G+se) + h)];
    }
  }

  if (I_sample) {
    vector[L + S + D + T] I_noise;
    if (I_overdispersion) {
        I_noise = I_raw .* sqrt(iota .* (1 + iota * (param_or_fixed(I_xi, I_xi_prior) ^ 2))); // approximates negative binomial
    } else {
      I_noise = I_raw .* sqrt(iota); // approximates Poisson
    }
    I_noise_correction[1:(G+se)] = rep_vector(0, (G+se));
    I_noise_correction[((G+se) + 1):(L + S + D + T)] = renewal_noise_correction(
      (L + S + D + T - (G+se)), G, se, gi_rev, I_noise
      );
    I[1 : (L + S + D + T)] = softplus(iota + I_noise_correction + I_noise, 10);
  } else {
    I_noise_correction = rep_vector(0, L + S + D + T);
  }

  // convolution from infections to shedding onsets (expected)
  lambda = convolve(inc_rev, I_sample ? I : iota)[(L + 1) : (L + S + D + T)];

  // calculation of total loads shed each day (expected)
  if (load_vari) {
    // Shedding load variation
    zeta_log[zeta_exact] = zeta_log_exact;
    zeta_log[zeta_normal_approx] = gamma_sum_log_approx(
      param_or_fixed(nu_zeta, nu_zeta_prior),
      lambda[zeta_normal_approx],
      zeta_raw_approx
    );

    omega_log = log_convolve(
        shed_rev_log, // shedding load distribution
        log(load_mean) + zeta_log // total load shed
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

  if (outliers) {
    pi_log = log_sum_exp(pi_log, log(load_mean) + log(epsilon)); // additive outlier component
  }

  // calculation of concentrations at measurement site by day (expected)
  // --> adjusted for flow and for date of sample effects
  if (K > 0) {
    kappa_log = pi_log - flow_log[1:T] + X[1:T] * eta;
  } else {
    kappa_log = pi_log - flow_log[1:T];
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

  if ((cv_type == 1) && total_partitions_observe!=1) {
    nu_upsilon_b = total_partitions_noncentered(
      param_or_fixed(nu_upsilon_b_mu, nu_upsilon_b_mu_prior),
      param_or_fixed(nu_upsilon_b_cv, nu_upsilon_b_cv_prior),
      nu_upsilon_b_noise_raw
    );
  }

  // scale for LOD hurdle model
  if (LOD_model == 1) {
    LOD_hurdle_scale[1] = rep_vector(LOD_scale[1], n_measured);
  } else if (LOD_model == 2) {
    LOD_hurdle_scale[1] = (
    n_averaged .*
    (total_partitions_observe ? dPCR_total_partitions : nu_upsilon_b * 1e4) *
    param_or_fixed(nu_upsilon_c, nu_upsilon_c_prior) * 1e-5
    );
  }
}
model {
  // Priors

  // Infections spline smoothing
  inf_ar_sd ~ normal(inf_ar_sd_prior[1], inf_ar_sd_prior[2]) T[0, ]; // truncated normal
  bs_coeff_noise_raw ~ std_normal(); // Gaussian noise

  // Seeding
  iota_log_seed_intercept ~ normal(iota_log_seed_intercept_prior[1], iota_log_seed_intercept_prior[2]);
  if (seeding_model == 1) {
    iota_log_seed_sd[1] ~ normal(iota_log_seed_sd_prior[1], iota_log_seed_sd_prior[2]) T[0, ]; // truncated normal
    iota_log_ar_noise ~ normal(0, iota_log_seed_sd[1]); // Gaussian noise
  }

  // Sampling of infections
  if (I_sample) {
    if (I_overdispersion) {
        target += normal_prior_lpdf(I_xi | I_xi_prior, 0); // truncated normal
    }
    I_raw[1 : (L + S + D + T)] ~ std_normal();
  }

  // Prior on individual-level shedding load variation
  if (load_vari) {
    target += normal_prior_lpdf(nu_zeta | nu_zeta_prior, 0); // truncated normal
    target += gamma3_sum_log_lpdf(zeta_log_exact | 1, param_or_fixed(nu_zeta, nu_zeta_prior), lambda[zeta_exact]); // gamma sum (exact)
    zeta_raw_approx ~ std_normal(); // non-centered noise for gamma sum (approx)
  }

  // Prior for additive outlier component
  if (outliers) {
    target += gev_lpdf(epsilon | epsilon_prior[1], epsilon_prior[2], epsilon_prior[3]); // generalized extreme value distribution
  }

  // Prior on sample date effects
  if (K > 0) {
    eta ~ normal(eta_prior[1], eta_prior[2]);
  }

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
    // accounting for pre-replicate noise
    vector[n_measured] concentration_log;
    if (pr_noise) {
      // Prior on cv of pre-replication concentrations
      nu_psi[1] ~ normal(nu_psi_prior[1], nu_psi_prior[2]) T[0, ]; // truncated normal
      // lognormal distribution modeled as normal on the log scale
      // with mean = rho (mean_log = rho_log) and cv = nu_psi
      target += lognormal_log_lpdf(psi | rho_log, nu_psi[1]);
      concentration_log = psi[measure_to_sample];
    } else {
      concentration_log = rho_log[measure_to_sample];
    }

    // concentration on unit scale
    vector[n_measured] concentration = exp(concentration_log);
    vector[n_measured] p_zero_log = rep_vector(negative_infinity(), n_measured);
    vector[n_measured] p_zero = rep_vector(0, n_measured);

    // limit of detection
    if (LOD_model > 0) {
      p_zero_log[i_LOD] = log_hurdle_exponential(
        concentration[i_LOD],
        LOD_hurdle_scale[1][i_LOD], // LOD scale (c * m * n)
        cv_type == 1 ? nu_upsilon_a : 0, // nu_pre (pre-PCR CV)
        cv_type == 1 ? cv_pre_type[1] : 0 // Type of pre-PCR CV
        );
      target += sum(p_zero_log[i_zero]); // below-LOD probabilities
      target += sum(log1m_exp(p_zero_log[i_nonzero_LOD])); // above-LOD probabilities
      p_zero = exp(p_zero_log);
    }

    // CV of each observation as a function of concentration
    vector[n_measured] cv;
    if (cv_type == 0) { // constant cv
      cv[i_include] = rep_vector(nu_upsilon_a, n_measured - n_zero);
    } else if (cv_type == 1) { // dPCR
      cv[i_include] = cv_dPCR_pre(
        concentration[i_include], // lambda (concentration)
        nu_upsilon_a, // nu_pre (pre-PCR CV)
        (total_partitions_observe ? dPCR_total_partitions[i_include] : nu_upsilon_b[i_include] * 1e4), // m (number of partitions)
        param_or_fixed(nu_upsilon_c, nu_upsilon_c_prior) * 1e-5, // c (conversion factor)
        n_averaged[i_include], // n (number of averaged replicates)
        cv_pre_type[1], // Type of pre-PCR CV
        cv_pre_approx_taylor[1] // Should taylor approximation be used?
        );
    } else if (cv_type == 2) { // constant variance
      cv[i_include] = (
        nu_upsilon_a * mean(measured_concentrations[i_include]) /
        concentration[i_include]
        );
    }

    vector[n_measured] mean_conditional;
    vector[n_measured] cv_conditional;
    mean_conditional[i_include] = concentration[i_include] ./ (1-p_zero[i_include]);
    cv_conditional[i_include] = sqrt(cv[i_include]^2 .* (1-p_zero[i_include]) - p_zero[i_include]);

    // measurements
    if (obs_dist == 0) {
      target += gamma3_lpdf(
        measured_concentrations[i_include] |
        mean_conditional[i_include], // expectation
        cv_conditional[i_include] // coefficient of variation
      );
    } else if (obs_dist == 1) {
      target += lognormal5_lpdf(
        measured_concentrations[i_include] |
        mean_conditional[i_include], // log expectation
        cv_conditional[i_include] // coefficient of variation
      );
    } else if (obs_dist == 2) {
      target += normal2_lpdf(
        measured_concentrations[i_include] |
        mean_conditional[i_include], // expectation
        cv_conditional[i_include], // coefficient of variation,
        0 // truncate at zero
      );
    } else if (obs_dist == 3) {
      target += normal2_lpdf(
        measured_concentrations[i_include] |
        mean_conditional[i_include], // expectation
        cv_conditional[i_include] // coefficient of variation
      );
    } else {
      reject("Distribution not supported.");
    }
  }
}
generated quantities {
  array[L + S + D + T - (G+se)] real R; // effective reproduction number

  // predicted measurements and forecasts
  // note that we here assume the same measurement variance as from composite samples,
  // which may be smaller than that of hypothetical daily measurements
  vector[T] predicted_concentration;
  vector[T] predicted_concentration_norm; // flow-normalized
  vector[h] predicted_concentration_forecast;
  vector[h] predicted_concentration_forecast_norm; // flow-normalized
  vector[h] R_forecast;
  vector[h] I_forecast;
  vector[h] lambda_forecast;
  vector[h] omega_log_forecast;
  vector[h] pi_log_forecast;
  vector[h] kappa_log_forecast;

  // reproduction number
  {
    vector[L + S + D + T] infs;
    vector[L + S + D + T - (G+se)] infness;
    if (I_sample) {
      infs = I;
    } else {
      infs = iota;
    }
    infness = infectiousness((L + S + D + T - (G+se)), G, se, gi_rev, infs);

    vector[L + S + D + T] new_infs_w = softplus(iota + I_noise_correction, 10);
    int max_t = L + S + D + T - (G+se);
    for (t in (R_w_half+1):(max_t-R_w_half)) {
      R[t] = sum(new_infs_w[((G+se)+t-R_w_half):((G+se)+t+R_w_half)]) ./
      sum(infness[(t-R_w_half):(t+R_w_half)]);
    }
  }

  // concentrations
  {
    vector[T+h] concentration_log;
    vector[T+h] concentration;
    vector[T+h] cv_all;
    vector[T+h] above_LOD; // will be a vector of 0s and 1s

    // Prediction for days until present
    if (pr_noise) {
      concentration_log[1:T] = to_vector(lognormal_log_rng(kappa_log, nu_psi[1]));
    } else {
      concentration_log[1:T] = kappa_log;
    }

    // Forecasting for days beyond present
    if (h>0) {

      // Forecasting of infections
      vector[h] I_noise_correction_forecast;
      if (I_sample) {
        vector[G+h] I_noise;
        vector[h] I_raw_forecast = std_normal_n_rng(h);
        I_noise[1:G] = I[(L + S + D + T - G + 1):(L + S + D + T)] - iota[(L + S + D + T - G + 1):(L + S + D + T)] - I_noise_correction[(L + S + D + T - G + 1):(L + S + D + T)];
        if (I_overdispersion) {
            I_noise[(G+1):(G+h)] = I_raw_forecast .* sqrt(iota_forecast .* (1 + iota_forecast * (param_or_fixed(I_xi, I_xi_prior) ^ 2))); // approximates negative binomial
        } else {
          I_noise[(G+1):(G+h)] = I_raw_forecast .* sqrt(iota_forecast); // approximates Poisson
        }
          I_noise_correction_forecast = renewal_noise_correction(
          h, G, 0, gi_rev, I_noise
          );
        I_forecast = softplus(iota_forecast + I_noise_correction_forecast + I_noise[(G+1):(G+h)], 10);
      } else {
        I_noise_correction_forecast = rep_vector(0, h);
        I_forecast = iota_forecast;
      }

      // Forecasting of R (back-calculation)
      vector[G + R_w - 1 + h] infs = append_row(I[(L + S + D + T - G - (R_w - 1) + 1):(L + S + D + T)], I_forecast);
      vector[R_w - 1 + h] infness = infectiousness(R_w - 1 + h, G, 0, gi_rev, infs);
      vector[G + R_w - 1 + h] new_infs_w = softplus(
          append_row(iota[(L + S + D + T - G - (R_w - 1) + 1):(L + S + D + T)], iota_forecast) +
          append_row(I_noise_correction[(L + S + D + T - G - (R_w - 1) + 1):(L + S + D + T)], I_noise_correction_forecast),
          10);

      for (t in 1:R_w_half) {
        int t_index = L + S + D + T - G - R_w_half + t;
        R[t_index] = sum(new_infs_w[(G+t):(G+t+R_w-1)]) ./
        sum(infness[t:(t+R_w-1)]);
      }
      for (t in 1:(h-R_w_half)) {
        R_forecast[t] = sum(new_infs_w[(G+t+R_w_half):(G+t+R_w_half+R_w-1)]) ./
        sum(infness[(R_w_half+t):(R_w_half+t+R_w-1)]);
      }

      // Forecasting of symptom onsets
      lambda_forecast = convolve(
        inc_rev, append_row((I_sample ? I : iota)[((L + S + D + T)+1-L):(L + S + D + T)], I_forecast)
        )[(L+1):(L+h)];

      // Forecasting of total loads
      if (load_vari) {
        vector[h] zeta_log_forecast = log(to_vector(gamma_rng(
          lambda_forecast/param_or_fixed(nu_zeta, nu_zeta_prior)^2,
          1/param_or_fixed(nu_zeta, nu_zeta_prior)^2
          )));
        omega_log_forecast = log_convolve(
            shed_rev_log, // shedding load distribution
            log(load_mean) + append_row(zeta_log[((S + D + T)+1-S):(S + D + T)], zeta_log_forecast) // total load shed
            )[(S+1):(S+h)];
      } else {
        omega_log_forecast = log_convolve(
            shed_rev_log, // shedding load distribution
            log(load_mean) + log(append_row(lambda[((S + D + T)+1-S):(S + D + T)], lambda_forecast)) // total load shed
            )[(S+1):(S+h)];
      }

      // Forecasting of total loads at sampling site (expected)
      if (D>0) {
        pi_log_forecast = log_convolve(
          residence_rev_log, // residence time distribution
          append_row(omega_log[((D + T)+1-D):(D + T)], omega_log_forecast)
          )[(D+1):(D+h)];
      } else {
        pi_log_forecast = omega_log_forecast;
      }

      // Forecasting of concentrations at measurement site by day (expected)
      // --> adjusted for flow and for date of sample effects
      if (K > 0) {
        kappa_log_forecast = pi_log_forecast - flow_log[(T+1):(T+h)] + X[(T+1):(T+h)] * eta;
      } else {
        kappa_log_forecast = pi_log_forecast - flow_log[(T+1):(T+h)];
      }

      // Forecasting of pre-replication concentrations
      if (pr_noise) {
        concentration_log[(T+1):(T+h)] = to_vector(lognormal_log_rng(kappa_log_forecast, nu_psi[1]));
      } else {
        concentration_log[(T+1):(T+h)] = kappa_log_forecast;
      }
    }

    concentration = exp(concentration_log);

    vector[T+h] nu_upsilon_b_all;
    if (cv_type == 1) {
      if (total_partitions_observe == 1) {
        nu_upsilon_b_all = total_partitions_all / 1e4;
      } else {
        nu_upsilon_b_all = total_partitions_noncentered(
          param_or_fixed(nu_upsilon_b_mu, nu_upsilon_b_mu_prior),
          param_or_fixed(nu_upsilon_b_cv, nu_upsilon_b_cv_prior),
          std_normal_n_rng(T+h)
          );
        for (i in 1:n_measured) {
          nu_upsilon_b_all[sample_to_date[measure_to_sample[i]]] = nu_upsilon_b[i];
        }
      }
    }

    vector[T+h] p_zero_all;
    if (LOD_model > 0) {
      vector[T+h] LOD_hurdle_scale_all;
      // scale for LOD hurdle model
      if (LOD_model == 1) {
        LOD_hurdle_scale_all = rep_vector(LOD_scale[1], T+h);
      } else if (LOD_model == 2) {
        LOD_hurdle_scale_all = (
        n_averaged_all .*
        nu_upsilon_b_all * 1e4 *
        param_or_fixed(nu_upsilon_c, nu_upsilon_c_prior) * 1e-5
        );
      }
      p_zero_all = exp(log_hurdle_exponential(
          concentration, // lambda (concentration)
          LOD_hurdle_scale_all,
          cv_type == 1 ? nu_upsilon_a : 0, // nu_pre (pre-PCR CV)
          cv_type == 1 ? cv_pre_type[1] : 0 // Type of pre-PCR CV
          ));
      p_zero_all = trim_or_reject_ub(
        p_zero_all,
        1-1e-5, // trim to almost 1
        1.01 // throw error when significantly above 1
      );
      above_LOD = to_vector(bernoulli_rng(1-p_zero_all));
    } else {
      p_zero_all = rep_vector(0, T+h);
      above_LOD = rep_vector(1, T+h);
    }

    if (cv_type == 0) {
      cv_all = rep_vector(nu_upsilon_a, T+h);
    } else if (cv_type == 1) {
      cv_all = cv_dPCR_pre(
        concentration, // lambda (concentration)
        nu_upsilon_a, // nu_pre (pre-PCR CV)
        nu_upsilon_b_all * 1e4, // m (number of partitions)
        param_or_fixed(nu_upsilon_c, nu_upsilon_c_prior) * 1e-5, // c (conversion factor)
        n_averaged_all, // n (number of averaged replicates) for all dates
        cv_pre_type[1], // Type of pre-PCR CV
        cv_pre_approx_taylor[1] // Should taylor approximation be used?
        );
    } else if (cv_type == 2) {
      cv_all = (
        nu_upsilon_a * mean(measured_concentrations[i_nonzero]) /
        concentration
        );
    }

    vector[T+h] mean_conditional_all = concentration ./ (1-p_zero_all);
    vector[T+h] cv_conditional_all = sqrt(cv_all^2 .* (1-p_zero_all) - p_zero_all);

    // correct potentially slightly negative approximations
    cv_conditional_all = trim_or_reject_lb(
      cv_conditional_all,
      1e-5, // trim to almost zero
      -0.01 // throw error when significantly below zero
    );

    vector[T+h] meas_conc;
    if (obs_dist == 0) {
      meas_conc = gamma3_rng(mean_conditional_all, cv_conditional_all);
    } else if (obs_dist == 1) {
      meas_conc = lognormal5_rng(mean_conditional_all, cv_conditional_all);
    } else if (obs_dist == 2) {
      meas_conc = normal2_rng(mean_conditional_all, cv_conditional_all, 0); // truncated at zero
    } else if (obs_dist == 3) {
      meas_conc = normal2_rng(mean_conditional_all, cv_conditional_all);
    } else {
      reject("Distribution not supported.");
    }
    predicted_concentration = above_LOD[1:T] .* meas_conc[1:T];
    predicted_concentration_norm = predicted_concentration .* flow[1:T];
    if (h>0) {
      predicted_concentration_forecast = above_LOD[(T+1):(T+h)] .* meas_conc[(T+1):(T+h)];
      predicted_concentration_forecast_norm = predicted_concentration_forecast .* flow[(T+1):(T+h)];
    }
  }
}
