/**
* Stochastic renewal process over T time steps
*
* @param T Number of time steps
*
* @param R Effective reproduction number
*
* @param G Maximum generation interval
*
* @param gi_rev The generation interval distribution in reversed format
*
* @param I The sampled infections (see model block). These should include the
* seeded infections from time (1-G):0, the vector thus has length G+T
*
* @return The expected infections from time 1:T
*/
vector renewal_process_stochastic(int T, vector R, int G, vector gi_rev, vector I) {
  vector[T] iota;
  for (t in 1:T) {
    iota[t] = R[t] * dot_product(gi_rev, I[t:(G+t-1)]);
  }
  return(iota);
}

vector log_renewal_process_stochastic(int T, vector R_log, int G, vector gi_rev_log, vector I_log) {
  vector[T] iota_log;
  for (t in 1:T) {
    iota_log[t] = R_log[t] + log_dot_product(gi_rev_log, I_log[t:(G+t-1)]);
  }
  return(iota_log);
}

array[] vector renewal_process_stochastic_noncentered(int T, vector R, int G, vector gi_rev, vector iota, vector I_raw, real I_xi) {
  vector[G+T] iota_tmp = iota;
  vector[G+T] I;
  real I_xi_squared = square(I_xi); // (I_xi = 0) => no overdispersion
  I[1:G] = iota_tmp[1:G] + I_raw[1:G] .* sqrt(iota_tmp[1:G] .* (1 + iota_tmp[1:G]*I_xi_squared));
  for (t in 1:T) {
    iota_tmp[G+t] = R[t] * dot_product(gi_rev, I[t:(G+t-1)]);
    I[G+t] = iota_tmp[G+t] + I_raw[G+t] * sqrt(iota_tmp[G+t] * (1 + iota_tmp[G+t]*I_xi_squared));
  }
  return {iota_tmp, I};
}

array[] vector log_renewal_process_stochastic_noncentered(int T, vector R_log, int G, vector gi_rev_log, vector iota_log, real I_xi, vector I_log_raw) {
   vector[G+T] iota_log_tmp = iota_log;
   vector[G+T] I_log;
   I_log[1:G] = approx_negative_binomial_log_noncentered(iota_log_tmp[1:G], I_xi, I_log_raw[1:G]);
   for (t in 1:T) {
     iota_log_tmp[G+t] = R_log[t] + log_dot_product(gi_rev_log, I_log[t:(G+t-1)]);
     I_log[G+t] = approx_negative_binomial_log_noncentered(iota_log_tmp[G+t], I_xi, I_log_raw[G+t]);
   }
   return {iota_log_tmp, I_log};
}

vector renewal_process_deterministic(int T, vector R, int G, vector gi_rev, vector iota) {
  vector[G+T] iota_tmp = iota;
  for (t in 1:T) {
    iota_tmp[G+t] = R[t] * dot_product(gi_rev, iota_tmp[t:(G+t-1)]);
  }
  return(iota_tmp[(G+1):(G+T)]);
}

vector log_renewal_process_deterministic(int T, vector R_log, int G, vector gi_rev_log, vector iota_log) {
  vector[G+T] iota_log_tmp = iota_log;
  for (t in 1:T) {
    iota_log_tmp[G+t] = R_log[t] + log_dot_product(gi_rev_log, iota_log_tmp[t:(G+t-1)]);
  }
  return(iota_log_tmp[(G+1):(G+T)]);
}

// Approximation for autocorrelated infection noise
// added to expected infection time series
vector renewal_noise_correction(int T, int G, vector gi_rev, vector I_noise) {
  vector[G+T] noise_tmp = I_noise;
  for (t in 1:T) {
    noise_tmp[G+t] = I_noise[G+t] + dot_product(gi_rev, noise_tmp[t:(G+t-1)]);
  }
  return((noise_tmp-I_noise)[(G+1):(G+T)]);
}

vector infectiousness(int T, int G, vector gi_rev, vector I) {
  vector[T] infectiousness;
  for (t in 1:T) {
    infectiousness[t] = dot_product(gi_rev, I[t:(G+t-1)]);
  }
  return(infectiousness);
}

vector log_infectiousness(int T, int G, vector gi_rev_log, vector I_log) {
  vector[T] infectiousness;
  for (t in 1:T) {
    infectiousness[t] = log_dot_product(gi_rev_log, I_log[t:(G+t-1)]);
  }
  return(infectiousness);
}

/**
* Convolution of a time series
*
* @param f The weight function, e.g. the incubation period distribution
*
* @param g The time series to be convolved, e.g. infections
*
* @return The convolved time series. The first length(f)-1 elements are NA
* because the convolved values can only be computed starting from length(f).
*/
vector convolve(vector f, vector g) {
  int f_length = num_elements(f);
  int g_length = num_elements(g);
  vector[g_length] fg;
  for (t in f_length:g_length) {
    fg[t] = dot_product(f, g[(t-f_length+1):t]);
  }
  return(fg);
}

  /**
* Convolution of a time series for T time steps on log scale
**
* @param f The weight function, e.g. the incubation period distribution
*
* @param g The time series to be convolved, e.g. infections
*
* @return The convolved log time series. The first length(f)-1 elements are NA
* because the convolved values can only be computed starting from length(f).
*/
vector log_convolve(vector f, vector g) {
  int f_length = num_elements(f);
  int g_length = num_elements(g);
  vector[g_length] fg;
  for (t in f_length:g_length) {
    fg[t] = log_dot_product(f, g[(t-f_length+1):t]);
  }
  return(fg);
}
