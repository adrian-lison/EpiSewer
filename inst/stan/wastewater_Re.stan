functions {
  #include functions/helper_functions.stan
  #include functions/time_series.stan
  #include functions/approx_count_dist.stan
  #include functions/renewal.stan
  #include functions/lognormal2.stan
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

  vector<lower=0>[T] flow; // flow rate for normalization of measurements

  // Sample date effects model
  int<lower=0> K; // number of sample date predictors
  matrix[T, K] X; // sample date predictor design matrix
  array[K > 0 ? 2 : 0] real eta_prior; // prior for sample date effects

  // Noise
  int<lower=0, upper=1> pre_replicate_noise; // Model variation before replication step?
  array[pre_replicate_noise ? 2 : 0] real tau_prior; // Prior on variation
  array[2] real sigma_prior; // Prior for scale of lognormal likelihood for measurements

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
}
transformed data {
  vector[G] gi_rev = reverse(generation_dist);
  vector[L + 1] inc_rev = reverse(incubation_dist);
  vector[S + 1] shed_rev_log = log(reverse(shedding_dist));
  vector[D + 1] residence_rev_log = log(reverse(residence_dist));
  vector[T] log_flow = log(flow);

  array[n_measured - num_zeros(measured_concentrations)] int<lower=0> i_nonzero;
  int i_nz = 0;
  for (n in 1:n_measured) {
    if (measured_concentrations[n] == 0) continue;
    i_nz += 1;
    i_nonzero[i_nz] = n;
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

  // realized infections
  real iota_log_ar_start;
  real<lower=0> iota_log_ar_sd;
  vector<multiplier=iota_log_ar_sd>[G - 1] iota_log_ar_noise;
  array[I_overdispersion ? 1 : 0] real<lower=0> I_xi; // positive to ensure identifiability
  vector<lower=0>[I_sample ? L + S + D + T : 0] I; // realized number of infections

  // sample date effects
  vector[K] eta;

  // Scale of lognormal likelihood for measurements
  array[pre_replicate_noise ? 1 : 0] real<lower=0> tau; // pre-replicaton variation
  vector<multiplier = (pre_replicate_noise ? tau[1] : 1)>[pre_replicate_noise ? n_samples : 0] psi; // realized noise before replication step
  real<lower=0> sigma;
}
transformed parameters {
  vector[L + S + D + T - G] R; // effective reproduction number
  vector[L + S + D + T] iota; // expected number of infections
  vector[S + D + T] lambda_log; // log expected number of symptom onsets
  vector[T] kappa_log; // log expected daily loads
  vector[T] pi_log; // log expected daily concentrations
  vector[n_samples] rho_log; // log expected concentrations in (composite) samples

  // ETS/Innovations state space process on log scale, starting value 1 on unit scale
  R = softplus(
    holt_damped_process(
      [R_level_start, R_trend_start]',
      ets_alpha_fixed < 0 ? ets_alpha[1] : ets_alpha_fixed,
      ets_beta_fixed < 0 ? ets_beta[1] : ets_beta_fixed,
      ets_phi_fixed < 0 ? ets_phi[1] : ets_phi_fixed,
      R_noise,
      0
    ), 4);

  // infections and renewal process
  iota[1 : G] = exp(random_walk([iota_log_ar_start]', iota_log_ar_noise, 0)); // seeding
  if (I_sample) {
    iota[(G + 1) : (L + S + D + T)] = renewal_process_stochastic(
      (L + S + D + T - G), R, G, gi_rev, I);
  } else {
    iota[(G + 1) : (L + S + D + T)] = renewal_process_deterministic(
      (L + S + D + T - G), R, G, gi_rev, iota);
  }

  // convolution from infections to log symptom onsets (expected)
  lambda_log = log(convolve(inc_rev, I_sample ? I : iota)[(L + 1) : (L + S + D + T)]);

  // calculation of loads at measurement site by day (expected)
  // this first convolves with the shedding load distribution and then
  // with the residence time distribution
  kappa_log = log_convolve(
    residence_rev_log,
    log_convolve(shed_rev_log, lambda_log)[(S + 1) : (S + D + T)]
    )[(D + 1) : (D + T)];

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

  // Innovations state space model
  ets_coefficient_priors_lp(ets_alpha, ets_alpha_fixed, ets_alpha_prior,
                            ets_beta, ets_beta_fixed, ets_beta_prior,
                            ets_phi, ets_phi_fixed, ets_phi_prior);
  R_level_start ~ normal(R_level_start_prior[1], R_level_start_prior[2]); // starting prior for level
  R_trend_start ~ normal(R_trend_start_prior[1], R_trend_start_prior[2]); // starting prior for trend
  R_sd ~ normal(R_sd_prior[1], R_sd_prior[2]) T[0, ]; // truncated normal
  R_noise ~ normal(0, R_sd); // Gaussian noise

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

  // Prior on sample date effects
  if (K > 0) {
    eta ~ normal(eta_prior[1], eta_prior[2]);
  }

  // Prior on scale of lognormal likelihood for measurements
  if (pre_replicate_noise) {
    tau[1] ~ normal(tau_prior[1], tau_prior[2]); // truncated normal
    psi ~ normal(0, tau[1]);
  }
  sigma ~ normal(sigma_prior[1], sigma_prior[2]); // truncated normal


  // Likelihood
  if (pre_replicate_noise) {
    target += lognormal3_lpdf(
      measured_concentrations[i_nonzero] | (rho_log+psi)[measure_to_sample][i_nonzero],
      sigma
      );
  } else {
    target += lognormal3_lpdf(
      measured_concentrations[i_nonzero] | rho_log[measure_to_sample][i_nonzero],
      sigma
      );
  }
}
generated quantities {
  // predicted measurements
  // note that we here assume the same measurement variance as from composite samples,
  // which may be smaller than that of hypothetical daily measurements
  array[T] real predicted_concentration;
  if (pre_replicate_noise) {
    vector[T] pre_repl = to_vector(normal_rng(pi_log, tau[1])); // pre-replication
    predicted_concentration = lognormal3_rng(pre_repl, sigma);
  } else {
    predicted_concentration = lognormal3_rng(pi_log, sigma);
  }
}
