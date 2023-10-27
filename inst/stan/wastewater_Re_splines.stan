functions {
  #include functions/helper_functions.stan
  #include functions/link.stan
  #include functions/time_series.stan
  #include functions/approx_count_dist.stan
  #include functions/renewal.stan
  #include functions/lognormal2.stan
  #include functions/gamma2.stan
  #include functions/hurdle.stan
}
data {
  int<lower=0> T; // number of total days in measured time span
  int<lower=0> n_samples; // number of samples from different dates
  int<lower=0> n_measured; // number of different measurements (including replicates)
  array[n_samples] int<lower=1, upper=T> sample_to_date; // index mapping samples to dates
  array[n_measured] int<lower=1, upper=n_samples> measure_to_sample; // index mapping measurements to samples
  vector<lower=0>[n_measured] measured_concentrations; // measured concentrations
  int<lower=1> w; // composite window: how many past days the samples cover,
  // e.g. 1 for individual day samples, 7 for weekly composite samples, ...

  real<lower=0> load_mean; // mean load shed per person
    int<lower=0, upper=1> load_vari; // model individual-level variation in shedding loads
  array[load_vari ? 2 : 0] real nu_prior; // Prior on individual-level load coefficient of variation

  vector<lower=0>[T] flow; // flow rate for normalization of measurements

  // Sample date effects model
  int<lower=0> K; // number of sample date predictors
  matrix[T, K] X; // sample date predictor design matrix
  array[K > 0 ? 2 : 0] real eta_prior; // prior for sample date effects

  // Noise
  int<lower=0, upper=1> pr_noise; // Pre-replicate noise: Model variation before replication step?
  array[pr_noise ? 2 : 0] real tau_prior; // Prior on coefficient of variation
  array[2] real cv_prior; // Prior for coefficient of variation of lognormal likelihood for measurements

  // Limit of detection
  real<lower=0> LOD;
  real<lower=0> LOD_sharpness;

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

  // Hyperpriors of random walk for seeding infections
  array[2] real iota_log_ar_start_prior;
  array[2] real iota_log_ar_sd_prior;

  // Stochastic (=1) or deterministic (=0) renewal process?
  int<lower=0, upper=1> I_sample;
  int<lower=0, upper=I_sample> I_overdispersion; // whether to model overdispersion via a negative binomial
  array[I_overdispersion ? 2 : 0] real I_xi_prior; // prior on the overdispersion parameter

  // Basis spline (bs) configuration for smoothing R
  // Sparse bs matrix: columns = bases (bs_n_basis), rows = time points (L+S+T-G)
  int<lower=1> bs_n_basis; // number of B-splines
  vector[bs_n_basis - 1] bs_dists; // distances between knots
  int<lower=0> bs_n_w; // number of nonzero entries in bs matrix
  vector[bs_n_w] bs_w; // nonzero entries in bs matrix
  array[bs_n_w] int bs_v; // column indices of bs_w
  array[L + S + D + T - G + 1] int bs_u; // row starting indices for bs_w plus padding
  array[2] real bs_coeff_ar_start_prior; // start hyperprior for random walk on log bs coeffs
  array[2] real bs_coeff_ar_sd_prior; // sd hyperprior for random walk on log bs coeffs

  // Link function and corresponding hyperparameters
  // first element: 0 = inv_softplus, 1 = scaled_logit
  // other elements: hyperparameters for the respective link function
  array[4] real R_link;
}
transformed data {
  vector[G] gi_rev = reverse(generation_dist);
  vector[L + 1] inc_rev = reverse(incubation_dist);
  vector[S + 1] shed_rev_log = log(reverse(shedding_dist));
  vector[D + 1] residence_rev_log = log(reverse(residence_dist));
  vector[T] log_flow = log(flow);

  int n_zero = num_zeros(measured_concentrations);
  array[n_zero] int<lower=0> i_zero;
  array[n_measured - n_zero] int<lower=0> i_nonzero;
  int i_z = 0;
  int i_nz = 0;
  for (n in 1:n_measured) {
    if (measured_concentrations[n] == 0) {
      i_z += 1;
      i_zero[i_z] = n;
    } else {
      i_nz += 1;
      i_nonzero[i_nz] = n;
    }
  }

  if (R_link[1] < 0 || R_link[1] > 1) {
    reject("Link function must be one of inv_softplus (0) or scaled_logit (1)");
  }
}
parameters {
  // log(R) time series prior
  real bs_coeff_ar_start; // intercept for random walk on log bs coeffs
  real<lower=0> bs_coeff_ar_sd; // sd for random walk on log bs coeffs
  vector[bs_n_basis - 1] bs_coeff_noise_raw; // additive errors (non-centered)

  // realized infections
  real iota_log_ar_start;
  real<lower=0> iota_log_ar_sd;
  vector<multiplier=iota_log_ar_sd>[G - 1] iota_log_ar_noise;
  array[I_overdispersion ? 1 : 0] real<lower=0> I_xi; // positive to ensure identifiability
  vector<lower=0>[I_sample ? L + S + D + T : 0] I; // realized number of infections

  // individual-level shedding load variation
  array[load_vari ? 1 : 0] real<lower=0> nu; // coefficient of variation of individual-level load
  vector[load_vari ? S + D + T : 0] zeta; // realized shedding load

  // sample date effects
  vector[K] eta;

  // Scale of lognormal likelihood for measurements
  array[pr_noise ? 1 : 0] real<lower=0> tau; // pre-replicaton coefficient of variation
  vector<multiplier = (pr_noise ? tau[1] : 1)>[pr_noise ? n_samples : 0] psi; // realized noise before replication step
  real<lower=0> cv;
}
transformed parameters {
  vector[bs_n_basis - 1] bs_coeff_noise; // additive errors
  vector[bs_n_basis] bs_coeff; // Basis spline coefficients
  vector[L + S + D + T - G] R; // effective reproduction number
  vector[L + S + D + T] iota; // expected number of infections
  vector[S + D + T] lambda; // expected number of shedding onsets
  vector[T] kappa_log; // log expected daily loads
  vector[T] pi_log; // log expected daily concentrations
  vector[n_samples] rho_log; // log expected concentrations in (composite) samplesles

  // Spline smoothing of R
  bs_coeff_noise = bs_coeff_noise_raw .* (bs_coeff_ar_sd * sqrt(bs_dists));
  bs_coeff = random_walk([bs_coeff_ar_start]', bs_coeff_noise, 0);

  R = apply_link(csr_matrix_times_vector(
    L + S + D + T - G, bs_n_basis, bs_w, bs_v, bs_u, bs_coeff
    ), R_link);

  // renewal process
  iota[1 : G] = exp(random_walk([iota_log_ar_start]', iota_log_ar_noise, 0)); // seeding
  if (I_sample) {
    iota[(G + 1) : (L + S + D + T)] = renewal_process_stochastic(
      (L + S + D + T - G), R, G, gi_rev, I);
  } else {
    iota[(G + 1) : (L + S + D + T)] = renewal_process_deterministic(
      (L + S + D + T - G), R, G, gi_rev, iota);
  }

  // convolution from infections to shedding onsets (expected)
  lambda = convolve(inc_rev, I_sample ? I : iota)[(L + 1) : (L + S + D + T)];

  // calculation of loads at measurement site by day (expected)
  // total load shed
  vector[S + D + T] load_realised_log;
  if (load_vari) {
    load_realised_log = log(zeta);
  } else {
    load_realised_log = log(load_mean) + log(lambda);
  }
  // distribution of load over time
  // this first convolves with the shedding load distribution and then
  // with the residence time distribution
  if (D>0) {
    kappa_log = log_convolve(
      residence_rev_log, // residence time distribution
      log_convolve(
        shed_rev_log, // shedding load distribution
        load_realised_log // total load shed
        )[(S + 1) : (S + D + T)]
      )[(D + 1) : (D + T)];
  } else {
    kappa_log = log_convolve(
      shed_rev_log, // shedding load distribution
      load_realised_log // total load shed
      )[(S + 1) : (S + T)];
  }

  // calculation of concentrations at measurement site by day (expected)
  // --> adjusted for flow and for date of sample effects
  if (K > 0) {
    pi_log = kappa_log - log_flow + X * eta;
  } else {
    pi_log = kappa_log - log_flow;
  }

  // concentrations in (composite) samples
  if (w > 1) { // multi-day composite samples
    for (i in 1 : n_samples) {
        rho_log[i] = log_sum_exp(pi_log[(sample_to_date[i] - w + 1) : sample_to_date[i]]) - log(w);
    }
  } else { // individual day samples
     for (i in 1 : n_samples) {
        rho_log[i] = pi_log[sample_to_date[i]];
     }
  }
}
model {
  // Priors

  // R spline smoothing
  bs_coeff_ar_start ~ normal(bs_coeff_ar_start_prior[1], bs_coeff_ar_start_prior[2]); // starting prior
  bs_coeff_ar_sd ~ normal(bs_coeff_ar_sd_prior[1], bs_coeff_ar_sd_prior[2]) T[0, ]; // truncated normal
  bs_coeff_noise_raw ~ std_normal(); // Gaussian noise

  // Sampling of infections
  iota_log_ar_start ~ normal(iota_log_ar_start_prior[1], iota_log_ar_start_prior[2]);
  iota_log_ar_sd ~ normal(iota_log_ar_sd_prior[1], iota_log_ar_sd_prior[2]) T[0, ]; // truncated normal
  iota_log_ar_noise ~ normal(0, iota_log_ar_sd); // Gaussian noise
  if (I_sample) {
    if (I_overdispersion) {
      I_xi[1] ~ normal(I_xi_prior[1], I_xi_prior[2]) T[0, ]; // truncated normal
      I[1 : (L + S + D + T)] ~ normal(iota, sqrt(iota .* (1 + iota * (I_xi[1] ^ 2)))); // approximates negative binomial
    } else {
      I[1 : (L + S + D + T)] ~ normal(iota, sqrt(iota)); // approximates Poisson
    }
  }

  // Prior on individual-level shedding load variation
  if (load_vari) {
    nu[1] ~ normal(nu_prior[1], nu_prior[2]); // truncated normal
    target += gamma3_sum_lpdf(zeta | load_mean, nu[1], lambda);
  }

  // Prior on sample date effects
  if (K > 0) {
    eta ~ normal(eta_prior[1], eta_prior[2]);
  }

  // Prior on scale of lognormal likelihood for measurements
  if (pr_noise) {
    tau[1] ~ normal(tau_prior[1], tau_prior[2]); // truncated normal
    psi ~ normal(0, sqrt(log1p(tau[1]^2))); // from tau(=cv) to sigma (lognormal)
  }
  cv ~ normal(cv_prior[1], cv_prior[2]); // truncated normal

  // Likelihood
  {
    // accounting for pre-replicate noise
    vector[n_measured] concentrations;
    if (pr_noise) {
      concentrations = (rho_log + psi)[measure_to_sample];
    } else {
      concentrations = rho_log[measure_to_sample];
    }

    // limit of detection
    if (LOD>0) {
     // below-LOD probabilities for zero measurements
    target += sum(log_inv_logit(hurdle_smooth(
      concentrations[i_zero], LOD, LOD_sharpness
      )));
      // above-LOD probabilities for non-zero measurements
    target += sum(log1m_inv_logit(hurdle_smooth(
      concentrations[i_nonzero], LOD, LOD_sharpness
      )));
    }

    // measurements
    target += lognormal4_lpdf(
      measured_concentrations[i_nonzero] |
        concentrations[i_nonzero],
      cv
      );
  }
}
generated quantities {
  // predicted measurements
  // note that we here assume the same measurement variance as from composite samples,
  // which may be smaller than that of hypothetical daily measurements
  vector[T] predicted_concentration;
  {
    vector[T] pre_repl;
    vector[T] above_LOD; // will be a vector of 0s and 1s
    if (pr_noise) {
      pre_repl = to_vector(normal_rng(pi_log, sqrt(log1p(tau[1]^2))));
    } else {
      pre_repl = pi_log;
    }
    if (LOD>0) {
     above_LOD = to_vector(bernoulli_rng(inv_logit(-hurdle_smooth(
      pre_repl, LOD, LOD_sharpness
      ))));
    } else {
      above_LOD = rep_vector(1, T);
    }
    predicted_concentration = above_LOD .* to_vector(lognormal4_rng(pre_repl, cv));
  }
}
