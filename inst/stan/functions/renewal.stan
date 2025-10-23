/**
* Stochastic renewal process over T time steps
*
* @param T Number of time steps
*
* @param R Effective reproduction number
*
* @param G Maximum generation interval
*
* @param se Length of shedding phase extension (beyond G)
*
* @param gi_rev The generation interval distribution in reversed format
*
* @param I The sampled infections (see model block). These should include the
* seeded infections from time (1-G):0, the vector thus has length G+T
*
* @return The expected infections from time 1:T
*/
vector renewal_process_stochastic(int T, vector R, int G, int se, vector gi_rev, vector I) {
  vector[T] iota;
  for (t in 1:T) {
    iota[t] = R[t] * dot_product(gi_rev, I[(t+se):(G+se+t-1)]);
  }
  return(iota);
}

vector log_renewal_process_stochastic(int T, vector R_log, int G, int se, vector gi_rev_log, vector I_log) {
  vector[T] iota_log;
  for (t in 1:T) {
    iota_log[t] = R_log[t] + log_dot_product(gi_rev_log, I_log[(t+se):(G+se+t-1)]);
  }
  return(iota_log);
}

/*
* // Example of how  non-centered renewal process can be modeled in EpiSewer:
*  array[2] vector[L + S + D + T] renewal_res = renewal_process_stochastic_noncentered(
*    (L + S + D + T - G), R, G, gi_rev, iota, I_raw, I_overdispersion ? I_xi[1] : 0
*  ); // this returns an array with 2 elements
*  iota = renewal_res[1]; // expected
*  I = renewal_res[2]; // realized
*
* // Then in model block:
*  I_raw[1 : (L + S + D + T)] ~ std_normal();
*
* // Note however that the non-centered parameterization does not work well for
* // the renewal model (longer runtime and higher chance of divergent transitions).
*/
array[] vector renewal_process_stochastic_noncentered(int T, vector R, int G, int se, vector gi_rev, vector iota, vector I_raw, real I_xi) {
  vector[(G+se)+T] iota_tmp = iota;
  vector[(G+se)+T] I;
  real I_xi_squared = square(I_xi); // (I_xi = 0) => no overdispersion
  I[1:(G+se)] = iota_tmp[1:(G+se)] + I_raw[1:(G+se)] .* sqrt(iota_tmp[1:(G+se)] .* (1 + iota_tmp[1:(G+se)]*I_xi_squared));
  for (t in 1:T) {
    iota_tmp[(G+se)+t] = R[t] * dot_product(gi_rev, I[(t+se):(G+se+t-1)]);
    I[(G+se)+t] = iota_tmp[(G+se)+t] + I_raw[(G+se)+t] * sqrt(iota_tmp[(G+se)+t] * (1 + iota_tmp[(G+se)+t]*I_xi_squared));
  }
  return {iota_tmp, I};
}

array[] vector log_renewal_process_stochastic_noncentered(int T, vector R_log, int G, int se, vector gi_rev_log, vector iota_log, real I_xi, vector I_log_raw) {
   vector[(G+se)+T] iota_log_tmp = iota_log;
   vector[(G+se)+T] I_log;
   I_log[1:(G+se)] = approx_negative_binomial_log_noncentered(iota_log_tmp[1:(G+se)], I_xi, I_log_raw[1:(G+se)]);
   for (t in 1:T) {
     iota_log_tmp[(G+se)+t] = R_log[t] + log_dot_product(gi_rev_log, I_log[(t+se):(G+se+t-1)]);
     I_log[(G+se)+t] = approx_negative_binomial_log_noncentered(iota_log_tmp[(G+se)+t], I_xi, I_log_raw[(G+se)+t]);
   }
   return {iota_log_tmp, I_log};
}

vector renewal_process_deterministic(int T, vector R, int G, int se, vector gi_rev, vector iota) {
  vector[(G+se)+T] iota_tmp = iota;
  for (t in 1:T) {
    iota_tmp[(G+se)+t] = R[t] * dot_product(gi_rev, iota_tmp[(t+se):(G+se+t-1)]);
  }
  return(iota_tmp[((G+se)+1):((G+se)+T)]);
}

vector log_renewal_process_deterministic(int T, vector R_log, int G, int se, vector gi_rev_log, vector iota_log) {
  vector[(G+se)+T] iota_log_tmp = iota_log;
  for (t in 1:T) {
    iota_log_tmp[(G+se)+t] = R_log[t] + log_dot_product(gi_rev_log, iota_log_tmp[(t+se):(G+se+t-1)]);
  }
  return(iota_log_tmp[((G+se)+1):((G+se)+T)]);
}

// simulation of a renewal process
array[] vector renewal_process_stochastic_sim_rng(int T, vector R, int G, int se, vector gi_rev, vector iota, real I_xi) {
  vector[(G+se)+T] iota_tmp = iota;
  vector[(G+se)+T] I;
  real I_xi_squared = square(I_xi); // (I_xi = 0) => no overdispersion
  I[1:(G+se)] = normal_lb_rng(iota_tmp[1:(G+se)], sqrt(iota_tmp[1:(G+se)] .* (1 + iota_tmp[1:(G+se)]*I_xi_squared)), 0);
  for (t in 1:T) {
    iota_tmp[(G+se)+t] = R[t] * dot_product(gi_rev, I[(t+se):(G+se+t-1)]);
    I[(G+se)+t] = normal_lb_rng(iota_tmp[(G+se)+t], sqrt(iota_tmp[(G+se)+t] * (1 + iota_tmp[(G+se)+t]*I_xi_squared)), 0);
  }
  return {iota_tmp, I};
}

// Approximation for autocorrelated infection noise
// added to expected infection time series
vector renewal_noise_correction(int T, int G, int se, vector gi_rev, vector I_noise) {
  vector[(G+se)+T] noise_tmp = I_noise;
  for (t in 1:T) {
    noise_tmp[(G+se)+t] = I_noise[(G+se)+t] + dot_product(gi_rev, noise_tmp[(t+se):(G+se+t-1)]);
  }
  return((noise_tmp-I_noise)[((G+se)+1):((G+se)+T)]);
}

vector infectiousness(int T, int G, int se, vector gi_rev, vector I) {
  vector[T] infness;
  for (t in 1:T) {
    infness[t] = dot_product(gi_rev, I[(t+se):(G+se+t-1)]);
  }
  return(infness);
}

vector log_infectiousness(int T, int G, int se, vector gi_rev_log, vector I_log) {
  vector[T] infness;
  for (t in 1:T) {
    infness[t] = log_dot_product(gi_rev_log, I_log[(t+se):(G+se+t-1)]);
  }
  return(infness);
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

/**
* Find growth rate for a given reproduction number and discrete generation time
* distribution
*
* @param R The effective reproduction number
*
* @param g The generation time distribution
*
* @return The calculated growth rate
*/
real get_growth_rate(real R, vector g) {
    real r = 0.1;  // initial guess
    int G = num_elements(g);
    vector[G] ks = linspaced_vector(G, 1, G);  // time indices

    for (i in 1:5) {
      vector[G] ek = exp(-r * ks);
      real S = dot_product(ek, g);
      real S_prime = -dot_product(ks .* ek, g);
      real f = inv(S) - R;                // f(r) = 1/S - R
      real f_prime = (-1) / (S * S) * S_prime; // f'(r) = -1/S^2 * S'
      r -= f / f_prime;  // Newton-Raphson update
    }
    return r;
  }

real get_R_from_growth_rate(real r, vector g) {
    int G = num_elements(g);
    vector[G] ks = linspaced_vector(G, 1, G);
    vector[G] ek = exp(-r * ks);
    real S = dot_product(ek, g);
    return inv(S);
  }
