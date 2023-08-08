/** ---------------------------------------------------------------------------
Time series process functions
---------------------------------------------------------------------------- */
 
  /**
  * Vectorized AR(1) process with coefficient of 1
  * Note that if the first value of the process is already included in the
  * increments_standardized vector, then you can provide a start_value = 0
  */
  vector random_walk(vector start_values, vector increments, int diff_order);
  
  vector random_walk(vector start_values, vector increments, int diff_order) {
    if (diff_order == 0) {
       return cumulative_sum(append_row(start_values[1], increments));
     } else {
      vector[diff_order] next_start = start_values[2:(diff_order+1)];
      int next_n = num_elements(increments) + diff_order - 1;
      vector[next_n] diffs = random_walk(next_start, increments, diff_order - 1);
      return cumulative_sum(append_row(start_values[1], diffs));
    }
  }
  
  /**
  * Vectorized MA(q) process with all coefficients = 1
  * Note that if the first value of the process is already included in the
  * increments_standardized vector, then you can provide a start_value = 0
  */
  vector simple_ma(vector start_values, vector increments, int Q, int diff_order);
  
  vector simple_ma(vector start_values, vector increments, int Q, int diff_order) {
    if (diff_order == 0) {
      int n = num_elements(increments)+1;
      vector[n] normal_cumsum = cumulative_sum(append_row(0, increments));
      vector[n] lagged_cumsum = cumulative_sum(append_row(rep_vector(0, Q+2), increments))[1:n];
      vector[n] window_sum = normal_cumsum - lagged_cumsum;
      return start_values[1] + window_sum; // start_value represents mean here
     } else {
      vector[diff_order] next_start = start_values[2:(diff_order+1)];
      int next_n = num_elements(increments) + diff_order - 1;
      vector[next_n] diffs = simple_ma(next_start, increments, Q, diff_order - 1);
      return cumulative_sum(append_row(start_values[1], diffs));
    }
  }
  
  /**
  * Damped Holt's method as an innovation state space model, vectorized implementation
  * Linear method with additive noise and without seasonality
  */
  vector holt_damped_process(vector start_values, real alpha, real beta_star, real phi, vector noise, int diff_order);
  
  vector holt_damped_process(vector start_values, real alpha, real beta_star, real phi, vector noise, int diff_order) {
    if (diff_order == 0) {
      int n = 1 + num_elements(noise); // n = 1 (for start value) + length of noise
      vector[n] y; // observations
      vector[n] epsilons = append_row(0, noise); // innovations
      
      // Compute level values
      vector[n] level = start_values[1] + alpha * cumulative_sum(append_row(0, epsilons))[1:n];
      
      // Compute trend values
      vector[n] trend;
      real beta = alpha * beta_star; // below, we always need the product of alpha and beta*
      if (phi == 0) {
        // special case: no trend
        trend = rep_vector(0, n);
      } else if (phi == 1) {
        // special case: trend, no dampening
        trend = start_values[2] + cumulative_sum(beta * cumulative_sum(append_row(0, epsilons))[1:n]);
      } else {
        // general case: trend and dampening
        vector[n] b = rep_vector(0, n);
        for (t in 2:n) { b[t] = phi * b[t - 1] + beta * epsilons[t - 1]; }
        trend = start_values[2] + phi * cumulative_sum(b);
      }
      
      // Compute observations 
      y = level + trend + epsilons;
      return y;
    } else {
      vector[diff_order+1] next_start = start_values[2:(diff_order+2)];
      int next_n = num_elements(noise) + diff_order;
      // Compute holt_damped_process at one diff_order lower 
      vector[next_n] diffs = holt_damped_process(next_start, alpha, beta_star, phi, noise, diff_order - 1);
      // Add up differences
      return cumulative_sum(append_row(start_values[1], diffs));
    }
  }
  
  void ets_coefficient_priors_lp(
    array[] real ets_alpha, real ets_alpha_fixed, array[] real ets_alpha_prior,
    array[] real ets_beta, real ets_beta_fixed, array[] real ets_beta_prior,
    array[] real ets_phi, real ets_phi_fixed, array[] real ets_phi_prior
  ) {
    if(ets_alpha_fixed < 0) {
      ets_alpha[1] ~ beta(ets_alpha_prior[1], ets_alpha_prior[2]); 
    }
    if(ets_beta_fixed < 0) {
      ets_beta[1] ~ beta(ets_beta_prior[1], ets_beta_prior[2]); 
    }
    if(ets_phi_fixed < 0) {
       // dampening needs a tight prior, roughly between 0.8 and 0.98
      ets_phi[1] ~ beta(ets_phi_prior[1], ets_phi_prior[2]);
    }
  }