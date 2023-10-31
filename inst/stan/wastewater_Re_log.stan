functions {
  #include functions/helper_functions.stan
  #include functions/time_series.stan
  #include functions/lognormal2.stan
  #include functions/approx_count_dist.stan
  #include functions/renewal.stan
}

data {
  int T; // number of total days in measured time span
  
  int n_measured; // number of different measurements
  vector<lower=0>[n_measured] measured_concentrations; // measured concentrations
  array[n_measured] int<lower=1,upper=T> measured_dates; // when these concentrations where measured
  int<lower=1> w; // composite window: how many past days the samples cover,
  // e.g. 1 for individual day samples, 7 for a weekly composite sample, ...
  vector<lower=0>[T] flow; // flow rate for normalization of measurements
  
  // Prior for scale of lognormal likelihood for measurements
  array[2] real sigma_prior;
  
    // Shedding load distribution
  int S; // last day of shedding
  vector[S+1] shedding_dist;// shedding load distribution
  // --> probability for a shedding today comes first
  
  // Incubation period distribution
  int L; // maximum delay
  vector[L+1] incubation_dist;// incubation period distribution
  // --> probability for a delay of zero comes first
  
  // Generation interval distribution
  int G;  // maximum generation interval
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
  vector[G] gi_rev_log;
  vector[L+1] inc_rev_log;
  vector[S+1] shed_rev_log;
  vector[T] log_flow;

  gi_rev_log = log(reverse(generation_dist));
  inc_rev_log = log(reverse(incubation_dist));
  shed_rev_log = log(reverse(shedding_dist));
  
  log_flow = log(flow);
}

parameters {
  // log(R) time series prior
  real R_level_start; // starting value of the level
  real R_trend_start; // starting value of the trend
  real<lower=0> R_sd; // standard deviation of additive errors
  vector<multiplier=(ets_noncentered ? R_sd : 1)>[L+S+T-G-1] R_noise; // additive errors
  
  // exponential smoothing / innovations state space process for log(R)
  array[ets_alpha_fixed < 0 ? 1 : 0] real<lower=0,upper=1> ets_alpha; // smoothing parameter for the level
  array[ets_beta_fixed < 0 ? 1 : 0] real<lower=0,upper=1> ets_beta; // smoothing parameter for the trend
  array[ets_phi_fixed < 0 ? 1 : 0] real<lower=0,upper=1> ets_phi; // dampening parameter of the trend
  
  // realized infections
  real iota_log_ar_start;
  real<lower=0> iota_log_ar_sd;
  vector<multiplier=iota_log_ar_sd>[G-1] iota_log_ar_noise;
  array[I_overdispersion ? 1 : 0] real<lower=0> I_xi; // positive to ensure identifiability
  vector[I_sample ? L+S+T : 0] I_log;  // realized number of infections
  
  // Scale of lognormal likelihood for measurements
  real<lower=0> sigma;
}

transformed parameters {
  vector[L+S+T-G] R_log; // effective reproduction number
  vector[L+S+T] iota_log; // expected number of infections
  vector[S+T] lambda_log; // log expected number of symptom onsets
  vector[T] pi_log; // log expected daily concentrations
  vector[n_measured] rho_log; // log expected concentrations in composite samples

  // ETS/Innovations state space process on log scale, starting value 1 on unit scale
  R_log = holt_damped_process(
      [R_level_start, R_trend_start]',
      ets_alpha_fixed < 0 ? ets_alpha[1] : ets_alpha_fixed,
      ets_beta_fixed < 0 ? ets_beta[1] : ets_beta_fixed,
      ets_phi_fixed < 0 ? ets_phi[1] : ets_phi_fixed,
      R_noise, 0);
      
  // renewal process
  iota_log[1:G] = random_walk([iota_log_ar_start]', iota_log_ar_noise, 0); // seeding
  if (I_sample) {
    iota_log[(G+1):(L+S+T)] = log_renewal_process_stochastic(L+S+T-G, R_log, G, gi_rev_log, I_log);
  } else {
    iota_log[(G+1):(L+S+T)] = log_renewal_process_deterministic(L+S+T-G, R_log, G, gi_rev_log, iota_log);
  }
  
  // convolution from infections to log symptom onsets (expected)
  lambda_log = log_convolve(inc_rev_log, I_sample ? I_log : iota_log)[(L+1):(L+S+T)];
  
  // calculation of concentrations at measurement site by day (expected)
  pi_log = log_convolve(shed_rev_log, lambda_log)[(S+1):(S+T)] - log_flow;
  
  // calculation of concentrations in composite samples
  for (i in 1:n_measured) {
    rho_log[i] = log_sum_exp(pi_log[(measured_dates[i]-w+1):measured_dates[i]]) - log(w);
  }
}

model {
  // Priors
  
  // Innovations state space model
  ets_coefficient_priors_lp(
    ets_alpha, ets_alpha_fixed, ets_alpha_prior,
    ets_beta, ets_beta_fixed, ets_beta_prior,
    ets_phi, ets_phi_fixed, ets_phi_prior
  );
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
      target += approx_negative_binomial_log_lpdf(I_log[1:(L+S+T)] | iota_log, I_xi[1]); // approximates negative binomial
    } else {
      target += approx_poisson_log_lpdf(I_log[1:(L+S+T)] | iota_log); // approximates Poisson
    }
  }
  
  // Prior on scale of lognormal likelihood for measurements
  sigma ~ normal(sigma_prior[1], sigma_prior[2]);  // truncated normal

  // Likelihood
  target += lognormal3_lpdf(measured_concentrations + 0.01 | rho_log, sigma);
}

generated quantities {
  // predicted measurements
  // note that we here assume the same measurement variance as from composite samples,
  // which may be smaller than that of hypothetical daily measurements
  array[T] real predicted_concentration = lognormal3_rng(pi_log, sigma);
  vector[L+S+T-G] R = exp(R_log); // effective reproduction number
  vector[L+S+T] iota = exp(iota_log); // expected number of infections
  vector[I_sample ? L+S+T : 0] I = exp(I_log);
}
